

#const XYIterGMMorLevel=Union{XYIterLevel, XYIterGMM}

iterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl, RegressionType::Val;
  ) where {TM, TV, T} = sum(
    (Θ(Xv, RegressionType) .- Xv.ṽ).^2)


function currentloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl,
   absμₖ::TM) where {TM<:CuMatrix, TV<:CuVector, T}

  someones::CuVector{T} = CUDA.ones(T, size(Θ.A,1))
  return sum((absμₖ*someones .- Xv.ṽ).^2)
end



function testiterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl;
  RegressionType::Val=Val(:choleskyzygote)) where {TM, TV, T}

  vhattest = iterloss(deepcopy(Θ), Xv, RegressionType)

  vhatver = sum((Θ(Xv, RegressionType) .- Xv.ṽ) .^ 2)

  @assert vhatver ≈ vhattest
end

function itermodel(panel::AbstractDataFrame,  zs::ZSpec,::Val{:ols},
    ::Type{T}=PARAM[:itertype],
    ::Type{TV}=PARAM[:itergpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:itergpu] ? CUDA.CuMatrix{T} : Matrix{T},
    PredictionType::Val = Val(PARAM[:predictiontype]),
    RegressionType::Val = Val(PARAM[:iterregressiontype]),
    SolveType::Val = Val(PARAM[:itersolvetype]),
    ; kwargs... #in case we want to use kwargs in other methods
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  (T == Float32) && (PARAM[:minweightmagnitude]<10^(-7)) && @warn(
    "T=$T while PARAM[:minweightmagnitude]=$(PARAM[:minweightmagnitude]).
    Therefore the minweight is effectivelly 0- is this really what you want?")
  CUDA.allowscalar(false)

  Zms::MeasureSpec = zs.Zms
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRtw::Vector{Symbol} = Fweights(Zms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  FW::Vector{Symbol} = Zms.FW
  FWgroup = Fgroupedcontrols(Zms)

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, T,TV,TM,
    weightdata = (;Fw, FRtw, FRLw), Fvol=Zms.Fvol; zs, FW, FWgroup)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  Rtws::TM = p.dat[:FRtw]
  RLws::TM = p.dat[:FRLw]
  W::TM = p.dat[:FW]
  Wgroup = p.dat[:FWgroup]
  ts::Vector{Int} = p.ts


  v::TV = p.dat[:Fvol]

  #println("ts: $ts")
  Xv = XYIterControl(ws, Rtws, RLws, v, ts, W=W, Wgroup=Wgroup)

  #construct the parameter
  Θ = VolumePartsIter(dims.T, dims.K, ts, T)

  local λ⁺::T = T(-1)

  #for debugging purposes- will trhow an error at some point if true
  PARAM[:runstacktraceoniter] && solveA!(iterloss,Θ,Xv,SolveType)

  #main optimization routine
  try
    λ⁺ = solveA!(iterloss,Θ,Xv,SolveType)
  catch err
    errmessage = "$err"
    if length(errmessage) > 20_000 #sometimes we get really long error messages
      errmessage = errmessage[1:20_000]
    end
    #@warn "Solver failed with error $errmessage"
    @error "Solver failed with error $errmessage\n" exception=(err, catch_backtrace())

  finally
    results::DataFrame = formresults(panel, Zms, Θ, Xv, Val{:iter}())

    #recompute the loss if needed
    if λ⁺<0
      throw("λ⁺ < 0!! This probably means the gpu failed. If thats not what happened
      it shouldn't be too hard to fix up the code and allow for this scenario. If it is,
      Just restart Julia and rerun.")
      λ⁺ = updateA₁!(Θ.g|>TV|>deepcopy, iterloss, Θ, Xv, RegressionType, Val{:identity}())
      updateAfromGAₜ!(Θ,Xv,1) #update to the final values of A
      #currentloss(Θ, Xv, projectvolumefromA(Θ,Xv))
      λ⁺ = lossfromA(Θ, Xv, RegressionType)
    end



    λ⁺ = isfinite(λ⁺) ? Int(round(min(λ⁺, typemax(Int)÷10))) : T(-99.0)

    return λ⁺, results
  end

  @assert false #should never reach here due to the finalizer
end

function itermodel(panel::AbstractDataFrame,  zs::ZSpec,::Val{:bootstrap},
    ::Type{T}=PARAM[:itertype],
    ::Type{TV}=PARAM[:iterbootstrapgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:iterbootstrapgpu] ? CUDA.CuMatrix{T} : Matrix{T},
    PredictionType::Val = Val(PARAM[:predictiontype]),
    RegressionType::Val = Val(PARAM[:iterregressiontype]);
    funds, kwargs...
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  (T == Float32) && throw(
    "If you really want to use Float32 using the bootstrap approach,
    you can disable this.")
  CUDA.allowscalar(false)

  Zms::MeasureSpec = zs.Zms
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRtw::Vector{Symbol} = Fweights(Zms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  FW::Vector{Symbol} = Zms.FW
  FWgroup = Fgroupedcontrols(Zms)

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, T,TV,TM,
    weightdata = (;Fw, FRtw, FRLw), Fvol=Zms.Fvol; zs, FW, FWgroup)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  Rtws::TM = p.dat[:FRtw]
  RLws::TM = p.dat[:FRLw]
  W::TM = p.dat[:FW]
  Wgroup = p.dat[:FWgroup]
  ts::Vector{Int} = p.ts


  v::TV = p.dat[:Fvol]

  #println("ts: $ts")
  Xv = XYIterControl(ws, Rtws, RLws, v, ts, W=W, Wgroup=Wgroup)

  #construct the parameter
  Θ = VolumePartsIter(dims.T, dims.K, ts, T, TM, TV)

  G = extractGfromfunds(TM; funds, zs, tidx, ts)
  println("dims G: $(size(G))
    dims Θ.G: $(size(Θ.G))")
  Θ.G .= G

  local λ⁺::T = T(-1)

  #for debugging purposes- will trhow an error at some point if true
  PARAM[:runstacktraceoniter] && solveA!(iterloss,Θ,Xv,
    Val{:bootstrap}())

  #main optimization routine
  #try
    λ⁺ = solveA!(iterloss,Θ,Xv,Val{:bootstrap}())
  #=catch err
    errmessage = "$err"
    if length(errmessage) > 20_000 #sometimes we get really long error messages
      errmessage = errmessage[1:20_000]
    end
    #@warn "Solver failed with error $errmessage"
    @error "Solver failed with error $errmessage\n" exception=(err, catch_backtrace())

  finally
    results::DataFrame = formresults(panel, Zms, Θ, Xv, Val{:iter}())

    #recompute the loss if needed
    if λ⁺<0
      throw("λ⁺ < 0!! This probably means the gpu failed. If thats not what happened
      it shouldn't be too hard to fix up the code and allow for this scenario. If it is,
      Just restart Julia and rerun.")
      λ⁺ = updateA₁!(Θ.g|>TV|>deepcopy, iterloss, Θ, Xv, RegressionType, Val{:identity}())
      updateAfromGAₜ!(Θ,Xv,1) #update to the final values of A
      #currentloss(Θ, Xv, projectvolumefromA(Θ,Xv))
      λ⁺ = lossfromA(Θ, Xv, RegressionType)
    end=#


    results::DataFrame = formresults(panel, Zms, Θ, Xv, Val{:iter}())
    results = addfundcolstoresults(results; funds, zs, tidx, ts)

    λ⁺ = isfinite(λ⁺) ? Int(round(min(λ⁺, typemax(Int)÷10))) : T(-99.0)

    return λ⁺, results
  #end

  @assert false #should never reach here due to the finalizer
end

#extracts the G from the fund file
function extractGfromfunds(::Type{TM}; funds, zs, tidx, ts, fundsource=PARAM[:iterfundsource],
  fundassetsprefix = PARAM[:iterfundassetprefix]) where TM <: AbstractMatrix
  ms = zs.originalms
  f = funds[fundsource] |> deepcopy
  alldates = sort(keys(tidx) |> collect)

  @assert issorted(f, :date)
  funddates = f.date

  #verify we have coverage of the growth rates over the estimation period
  @assert (setdiff(alldates, funddates) |> isempty) "
    setdiff(alldates,funddates): $(setdiff(alldates,funddates))"

  cols = ms.Fξs .|> (s->Symbol(fundassetsprefix, s))

  cleanfunds = completecases(f, cols)

  A = f[f.date .|> (d->insorted(d, alldates)), cols] |> Matrix{Float64}
  G = (A[2:end, :] ./ A[1:(end-1),:])' |> TM

  #NOTE: this approach won't work with smoothing- but wouldn't smoothing
  #break the accounting relationship?
  return G
end

function addfundcolstoresults(results; funds, zs, tidx, ts, fundsource=PARAM[:iterfundsource],
  fundassetsprefix = PARAM[:iterfundassetprefix],
  additionalfundscols = PARAM[:iteradditionalfundscols]) where TM <: AbstractMatrix
  ms = zs.originalms
  f = funds[fundsource] |> deepcopy
  alldates = sort(keys(tidx) |> collect)

  @assert issorted(f, :date)
  funddates = f.date

  #verify we have coverage of the growth rates over the estimation period
  @assert (setdiff(alldates, funddates) |> isempty) "
    setdiff(alldates,funddates): $(setdiff(alldates,funddates))"

  cols = [ms.Fξs .|> (s->Symbol(fundassetsprefix, s)); additionalfundscols]

  cleanfunds = completecases(f, cols)

  colstoadd = f[f.date .|> (d->insorted(d, alldates)), [:date; cols;]] |> deepcopy

  newcolnames = cols .|> (s)->Symbol(fundsource, s)
  rename!(colstoadd, Dict(cols .=> newcolnames))

  results = leftjoin(results, colstoadd, on=:date, validate=(true,true))



  #NOTE: this approach won't work with smoothing- but wouldn't smoothing
  #break the accounting relationship?
  return results
end
