function mcmcmodel(panel::AbstractDataFrame,  zs::ZSpec,::Val{:randomeffects},
    ::Type{T}=PARAM[:mcmctype],
    ::Type{TV}=PARAM[:mcmcgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:mcmcgpu] ? CUDA.CuMatrix{T} : Matrix{T},
    prior::Symbol=PARAM[:mcmcprior],
    Nburn::Int=PARAM[:mcmcnburn],
    Ndraws::Int=PARAM[:mcmcndraws]
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  @assert PARAM[:predictiontype] ∈ [:leveldemean, :level]
  (T == Float64) || @warn("T=$T, which may lead to numerical instability.
    Is this really what you want?")
  CUDA.allowscalar(false)

  Zms::MeasureSpec = zs.Zms
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  essentialfields = [Fw; FRLw; Zms.Fvol; ]

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, zs,T,TV,TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=Zms.Fvol)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  ts::Vector{Int} = p.ts
  Ascale::Vector{T} = p.Ascale

  v::TV = p.dat[:Fvol]

  #println("ts: $ts")
  Xv = XYIter(ws, RLws, v, ts, Ascale)


  #construct the main MCMC object for book keeping and parameters
  hyper=PARAM[prior]
  initial=(
    A₀ = ones(T, dims.K),
    G = ones(T, dims.K, dims.T-1),
    σ² = ones(T, dims.T-1),
    ω = zeros(T, dims.T-1)
    )
  cc = CapacityChains(T, Nburn=Nburn, Ndraws=Ndraws, paramstemplate=initial, hyper=hyper)
  #construct the parameter
  Θ = VolumePartsMCMC(TM, TV, T; initial, cc, Xv, ts)
  #printmln(Θ.Ψ)
  createmcmc(Θ)

  #for debugging purposes- will trhow an error at some point if true
  PARAM[:runstacktraceonmcmc] #=&& solveA!(iterloss,Θ,Xv,_Xv, SolveType)

  #main optimization routine
  try
    λ⁺ = solveA!(iterloss,Θ,Xv,_Xv, SolveType)
  catch err
    errmessage = "$err"
    if length(errmessage) > 20_000 #sometimes we get really long error messages
      errmessage = errmessage[1:20_000]
    end
    @warn "Solver failed with error $errmessage"
  finally
    results::DataFrame = formresults(panel, Zms, Θ, Val{:iter}())

    #recompute the loss if needed
    if λ⁺<0
      updateA₁!(Θ.g|>deepcopy, iterloss, Θ, Xv, _Xv, RegressionType, Val{:identity}())
      updateAfromGAₜ!(Θ,Xv, _Xv,1) #update to the final values of A
      λ⁺ = currentloss(Θ, Xv, _Xv, projectvolumefromA(Θ,Xv,_Xv))
    end



    if Ncleanpanel < 500_000 && PARAM[:testiterresults]
      Aexpanded = view(Θ.A', Θ.ts, :)
      LAexpanded = view(Θ.A', Θ.tsM1, :)
      need2demean = genabsμₖ.(Aexpanded, LAexpanded, Xv.ws, Xv.RLws)
      factored = Xv.M * need2demean

      Fcoefs::Vector{Symbol} = Fcoefficients0(zs.Zms)
      Fwcoefs::Vector{Symbol} = (s->Symbol(:w,s)).(Fcoefs)
      FRLwcoefs::Vector{Symbol} = (s->Symbol(:RLw,s)).(Fcoefs)
      FAcoefs::Vector{Symbol} = (s->Symbol(:A,s)).(Fcoefs)
      FLAcoefs::Vector{Symbol} = (s->Symbol(:LA,s)).(Fcoefs)
      Ffcoefs::Vector{Symbol} = (s->Symbol(:f,s)).(Fcoefs)

      dfdebug = hcat(
        DataFrame(ts=Θ.ts, tsM1=Θ.tsM1, v=Xv.v, vtilde=Xv.ṽ),
        DataFrame(Xv.ws, Fwcoefs),
        DataFrame(Xv.RLws, FRLwcoefs),
        DataFrame(Aexpanded, FAcoefs),
        DataFrame(LAexpanded, FLAcoefs),
        DataFrame(factored, Ffcoefs))

      dfdebug |> CSV.write(
        "$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_opt$(iterresultsname(λ⁺|>round|>Int)).csv")
    end

    return Int(round(min(λ⁺, typemax(Int)÷10))), results
  end=#

  @assert false #should never reach here due to the finalizer
end

function updateΨ!(Θ::AbstractVolumePartsMCMC{TM, TV, T}, ::Val{:randomeffects}) where {TM,TV,T}
  #@info "TM:$TM, TV:$TV, T:$T"
  updateA₀!(Θ)
end

#this will create the mcmc chain
function createmcmc(Θ::AbstractVolumePartsMCMC; stateupdater::Val = Val{PARAM[:mcmcstateupdater]}())
  @unpack cc, Ψ = Θ
  @unpack Nburn = cc

  Niter = Nburn + size(cc.value,1)

  for i ∈ 1:Niter
    (i==Nburn+1) && @assert(cc.idx==1)
    updateΨ!(Θ, stateupdater)
    (i<Niter) && updateandincrement!(cc, Ψ, verify=true)
  end
end
