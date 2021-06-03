

#creates a model using flux that uses levels (as opposed to differences)
#WARNING- DO NOT DELETE BELOW
function fluxlevelmlemodel(panel::AbstractDataFrame,  zs::ZSpec,
    ::Type{T}=PARAM[:fluxtype],
    ::Type{TV}=PARAM[:fluxgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:fluxgpu] ? CUDA.CuMatrix{T} : Matrix{T},
    PredictionType::Val = Val(PARAM[:predictiontype])) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "invalid predictiontype $PredictionType for model type $(PARAM[:fluxmodel])")

  #(PARAM[:fluxgpu] && (T≠Float32)) && error("use Float32 for fluxtype when using gpu")

  CUDA.allowscalar(false)
  Zms::MeasureSpec = zs.Zms

  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, Zms, T, TV, TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=Zms.Fvol)
  dims = p.dims
  tidx::Dict = p.tidx
  cleanpanel::SubDataFrame = p.cleanpanel
  Ncleanpanel::Int = nrow(cleanpanel)
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  ts::Vector{Int} = p.ts
  v::TV = p.dat[:Fvol]
  X = (ws=ws, RLws=RLws)
  #construct the parameter
  Θ = VolumePartsMLE(dims.T, dims.K, ts, T)
  Π=Flux.params(Θ)

  #define the loss function (in this case, log-likelihood)
  #Npanelones = ones(Npanel) |> fluxgpu
  iter::Int = 0
  function loss(X::NamedTuple, y::AbstractVector)
    θ = view(exp.(Θ.τ),Θ.expand.τidx)
    #keep count of each iteration
    iter += 1
    absμ = Θ(X, y, PredictionType)

    #θ = exp.(τlong) #here τ=logθ
    k = absμ ./ θ  #absμ ./ θ # μ = kθ,

    if (sum(θ .> 0) ≠ Npanel) || (sum(k .> 0) ≠ Npanel)
      return Inf #this will cause the program to stop
    end

    loglike = sum(vecgammalogpdf(k, θ, y))

    return -loglike
  end

  #NOTE:" cuda optimized version
  function loss(X::NamedTuple, y::CuVector)
    θ = CuVector(exp.(Θ.τ)[Θ.expand.τidx])
    #keep count of each iteration
    iter += 1
    absμ = Θ(X, y, PredictionType)

    #θ = CUDA.exp_fast.(τlong) #here τ=logθ
    k = absμ ./ θ # CUDA.div_fast  crashes for some reason

    if (CUDA.sum(θ .> 0) ≠ Ncleanpanel) || (CUDA.sum(k .> 0) ≠ Ncleanpanel)
      return Inf #this will cause the program to stop
    end

    loglike = CUDA.sum(vecgammalogpdf(k, θ, y))

    #println("minval: $minval, loglike: $loglike")
    return -loglike
  end

  #use this function to throw data from the main block over to
  #a seperate function to report error information
  debugcb() = debuginf(cleanpanel, tidx,  Θ, Matrix(ws), Matrix(RLws), Vector(v))

  @info "computing test gradient. Testing function: $(sum(Θ(X, v, Val{:level}())))"
  gs = gradient(()->loss(X, v),Π)
  #@info "gs[A]: $(gs[Θ.A])"

  test∇VolumeParts!(loss, X, v, Θ)

  @info "beginning training with $(length(Π)) parameter vectors"

  λ⁺=solvevol!(OPT_INDEX[PARAM[:fluxopt]], T;
    loss=loss,
    Θ=Θ,
    X=X, y=v,
    debugcb = debugcb)


  results::DataFrame = formresults(panel, Zms, Θ, PredictionType)

  return (Int(round(λ⁺)), results)
end


function cudadifsq(a,b)
  Δ = a - b
  return Δ*Δ
end



function fluxlevellsqmodel(panel::AbstractDataFrame,  zs::ZSpec,
    ::Type{T}=PARAM[:fluxtype],
    ::Type{TV}=PARAM[:fluxgpu] ? CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:fluxgpu] ? CuMatrix{T} : Matrix{T};
    #::Type{PredictionType} = Val(PARAM[:predictiontype])
    PredictionType::Val = Val(PARAM[:predictiontype])#Capacity.eval(:(Val(PARAM[:predictiontype])))
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}
  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "invalid predictiontype $PredictionType for model type $(PARAM[:fluxmodel])")
  #(PARAM[:fluxgpu] && (T≠Float32)) && error("use Float32 for fluxtype when using gpu")

  CUDA.allowscalar(false)
  Zms::MeasureSpec = zs.Zms

  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, Zms,T, TV, TM,
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
  Θ = VolumePartsLSq(dims.T, dims.K, ts, T)
  Π=Flux.params(Θ)
  length(Π) == 0 && error("No parameters noticed by Flux!")

  #define the loss function (in this case, log-likelihood)
  #Npanelones = ones(Npanel) |> fluxgpu
  #ptype2 = Val{:levellog}()

  function loss(X::NamedTuple,y::AbstractVector) where TX<:AbstractMatrix

    #keep count of each iteration
    absμ = Θ(X, y, PredictionType)

    sumsq = sum((absμ .- y) .^ 2)

    return sumsq
  end

  #NOTE:" cuda optimized version
  function loss(X::NamedTuple, y::CuVector)

    #keep count of each iteration
    absμ = Θ(X, y, PredictionType)

    sumsq = CUDA.sum(cudadifsq.(absμ,y))

    #println("minval: $minval, loglike: $loglike")
    return sumsq
  end

  #use this function to throw data from the main block over to
  #a seperate function to report error information
  debugcb() = debuginf(cleanpanel, tidx,  Θ, Matrix(ws), Matrix(RLws), Vector(v))

  @info "computing test gradient. Testing function: $(sum(Θ(X, v, PredictionType)))"


  test∇VolumeParts!(loss, X, v, Θ)

  @info "beginning training with $(length(Π)) parameter vectors"

  λ⁺=solvevol!(OPT_INDEX[PARAM[:fluxopt]], T;
    loss=loss,
    Θ=Θ,
    X=X,
    y=v,
    debugcb = debugcb)


  results::DataFrame = formresults(panel, Zms, Θ, PredictionType)

  return Int(round(λ⁺)), results
end
