
function optkeywords()
  #this will hold the arguments
  optkwargs::Dict{Symbol,Real} = Dict()
  kws = [:f_tol, :g_tol, :x_tol, :iterations, :Δ]
  sizehint!(optkwargs, length(kws))

  paramkw(s) = Symbol(PARAM[:optprefix], :_ , s)

  for kw ∈ kws
    arg = PARAM[paramkw(kw)]
    if !(arg===nothing)
      push!(optkwargs, kw=>arg)
    end
  end

  return optkwargs
end

function getoptbounds(zs::ZSpec,
  Θ::VolumePartsDLSq, dims::NamedTuple,
  PredictionType::Val = Val(PARAM[:predictiontype]);
  minvol::Float64=10^-6)

  Abounds = PARAM[:abounds]
  V₀bounds = PARAM[:v0bounds]

  #this willl hold the bounds
  ZAbounds::Matrix{Tuple{Float64,Float64}} = Matrix{Tuple{Float64,Float64}}(undef, dims.K, dims.T)
  ZV₀bounds::Vector{Tuple{Float64,Float64}}= Vector{Tuple{Float64,Float64}}(undef, dims.T-1)

  function no0log(x)
    (x<0) && error("no0log cannot have a negative input (did you use a negative V₀ bound?)")
    return log(max(minvol,x))
  end

  #need to scale the bounds at the same level as the coefficients
  σᵥ::Float64 = zs.scale[zs.Zms.Fvol]
  scaledV₀bounds = (t-> no0log(t / σᵥ)).(V₀bounds)

  #scale the intercept term
  for t ∈ 1:(dims.T-1)
    ZV₀bounds[t] = scaledV₀bounds
  end

  #scale the A matrix coefficeint bounts
  for (k, FZcoef) ∈ enumerate(zs.Zms.Fξs)
    FZw::Symbol = zs.Zms[FZcoef].X[:Fw]
    σ_kw = zs.scale[FZw]
    scaledAbounds = (t->t*σ_kw/σᵥ).(Abounds)
    for t ∈ 1:dims.T
      ZAbounds[k,t] =  scaledAbounds
    end
  end

  #combine the bounds
  bounds=[vec(ZAbounds); ZV₀bounds]

  return bounds
end



#generates an in-place prediction function
function leastsquaresobj(Θ::AbstractVolumeParts, X::NamedTuple, v::AbstractVector,
  PredictionType::Val)

  function obj!(out::AbstractVector, Θ::AbstractVolumeParts)
    out .= Θ(X,v,PredictionType) .- v
    return nothing
  end

  return obj!
end

function optmodel(panel::AbstractDataFrame,  zs::ZSpec, ::Val{:levellsq},
    ::Type{T}=PARAM[:opttype],
    ::Type{TV}=Vector{T},
    ::Type{TM}=Matrix{T},
    PredictionType::Val = Val(PARAM[:predictiontype]),
    optkwargs = optkeywords(),
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "invalid predictiontype $PredictionType for model type $(PARAM[:fluxmodeltype])")
  (T == Float64) || @warn("T=$T, which may lead to numerical instability.
    Is this really what you want?")

  Zms::MeasureSpec = zs.Zms
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  essentialfields = [Fw; FRLw; Zms.Fvol; ]

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, Zms,T,TV,TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=Zms.Fvol)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  ts::Vector{Int} = p.ts
  X = (ws=ws, RLws=RLws)
  v::TV = p.dat[:Fvol]

  #construct the parameter
  Θ = VolumePartsDLSq(dims.T, dims.K, ts, T)


  J, J! = JVolumePartsDLSq(Θ, X, v, T, PredictionType)
  testjacobian(Θ, X, v, T, J = J, J! = J!, dims=dims)
  obj! = leastsquaresobj(Θ, X, v, PredictionType)

  out = similar(v)
  obj!(out,Θ)
  vhatver = Θ(X,v,PredictionType)
  @assert out .+ v ≈ vhatver

  #this is not used in the optimization, but keeps track of our objective
  function loss(Θ::AbstractVolumeParts, X::NamedTuple,y::AbstractVector) where TX<:AbstractMatrix

    #keep count of each iteration
    absμ = Θ(X, y, PredictionType)

    sumsq = sum((absμ .- y) .^ 2)

    return sumsq
  end

  @info "Beginning optimization with $(optkwargs)"
  opt = optimize!(LeastSquaresOptim.LeastSquaresProblem(x=Θ, f! = obj!, g! = J!, J=J, output_length=Ncleanpanel),
    LevenbergMarquardt(LeastSquaresOptim.LSMR()); optkwargs...)

  @info "opt: \n$opt"
  λ⁺ = loss(Θ, X,v)

  #compute the full jacobian with AD
  Θver = VolumePartsLSq(dims.T, dims.K, ts, T)
  Θver .= Θ
  ∇Θver = gradient(()->loss(Θver,X,v),Flux.params(Θver))
  rmsver = sqrt(mean(vec∇(Θver, ∇Θver).^2))
  @info "lst squares gradient rms: $rmsver"


  #=J!(J,Θ)

  Jcolsums = vec(sum(J, dims=1))
  @info "opt:\n$opt"
  @info "Jcolsums: $Jcolsums"
  @info "L₁(Jcolsums)=$(sum(abs.(Jcolsums)))"=#

  results::DataFrame = formresults(panel, Zms, Θ, PredictionType)

  return Int(round(λ⁺)), results
end


function optmodel(panel::AbstractDataFrame,  zs::ZSpec,::Val{:levellsqbb},
    ::Type{T}=PARAM[:opttype],
    ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{T} : Matrix{T},
    PredictionType::Val = Val(PARAM[:predictiontype]);
    evalmodels::Bool = PARAM[:optbbcomparemodels],
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "invalid predictiontype $PredictionType for model type $(PARAM[:fluxmodeltype])")
  (T == Float64) || @warn("T=$T, which may lead to numerical instability.
    Is this really what you want?")
  CUDA.allowscalar(false)

  Zms::MeasureSpec = zs.Zms
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  essentialfields = [Fw; FRLw; Zms.Fvol; ]

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, Zms,T,TV,TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=Zms.Fvol)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  ts::Vector{Int} = p.ts
  X = (ws=ws, RLws=RLws)
  v::TV = p.dat[:Fvol]

  #construct the parameter
  Θ = VolumePartsDLSq(dims.T, dims.K, ts, T)

  testprediction(Θ, X, v, PredictionType)

  #least squares objective funciton
  function loss(Θ::AbstractVolumeParts, X::NamedTuple,y::AbstractVector)

    #keep count of each iteration
    absμ = Θ(X, y, PredictionType)

    sumsq = sum((absμ .- y) .^ 2)

    return sumsq
  end

  function loss(Θ::AbstractVolumeParts, X::NamedTuple,y::CuVector{T})

    #keep count of each iteration
    absμ = Θ(X, y, PredictionType)

    krn(absμᵢ,yᵢ) = (absμᵢ - yᵢ)*(absμᵢ - yᵢ)
    sumsq = CUDA.sum(krn.(absμ,y))

    return sumsq
  end

  #this just puts the problem in the form used by bboptimize
  function loss(Π::AbstractVector)
    Θ.Π .= TV(Π)
    return loss(Θ, X, v)
  end

  bounds = getoptbounds(zs, Θ, dims, PredictionType)

  if evalmodels
    @warn "Testing models! This could take a while...."
    optbbcomparemodels(
      bounds=bounds,
      loss=loss,
      Θ=Θ,
      X=X,
      y=v,)
  end

  opt, Ψ = optbbestimate(
    bounds=bounds,
    loss=loss,
    Θ=Θ,
    X=X,
    y=v,
  )

  λ⁺ = loss(Θ, X,v)

  #const a differentiable version with nromal data types of parameters
  Θver = VolumePartsLSq(dims.T, dims.K, Θ.ts, T)
  Θver .= Vector(Θ.Π)
  @assert loss(Θver,X,v) ≈ λ⁺

  optbbresults(loss=loss, opt=opt, Ψ=Ψ, Θ=Θ, Θver=Θver, X=X, y=v)
  results::DataFrame = formresults(panel, Zms, Θver, PredictionType)

  return Int(round(min(λ⁺, typemax(Int)÷10))), results
end
