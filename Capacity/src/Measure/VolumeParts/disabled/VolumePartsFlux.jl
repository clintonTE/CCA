


###########∇VolumeParts
#this struct/function computes the Jacobian
#vec∇(Θ::VolumeParts, ∇Θ) = reduce(vcat, (p->vec(∇Θ[getproperty(Θ,p)])).(Θ.Πelements))

function ∇VolumeParts!(loss, X::NamedTuple, y, Θ::AbstractVolumeParts)
  function ∇Θ!(out::AbstractVector, Θ)
    ∇Θ = gradient(()->loss(X,y), Flux.params(Θ))
    out .= vec∇(Θ, ∇Θ)
    return Θ
  end

  return ∇Θ!
end

Flux.@functor VolumePartsLSq
#Flux.trainable(a::VolumePartsLSq) = (a.A, a.V₀,)
Flux.params(Θ::VolumePartsLSq) = Flux.params(Θ.A, Θ.V₀,)

Flux.@functor VolumePartsIter
Flux.params(Θ::VolumePartsIter) = Flux.params(Θ.G,)

#WARNING- OPT_INDEX was originally in Capacity.jl- best to leave it here if possible, but it may
#need to be defined earlier
OPT_INDEX = Dict{Symbol,Any}(
  :adam=>Flux.ADAM, #47s (9.7->3.7)
  :radam=>Flux.RADAM, #fail
  :rmsprop=>Flux.RMSProp, #44s (8.4->2.7)
  :nadam=>Flux.NADAM, #70s (9.7->3.6)
  :adamw=>Flux.ADAMW,# 43s, (9.1->3.5)
  :momentum=>Flux.Momentum, #fail
  :nesterov=>Flux.Nesterov,#fail
  :descent=>Flux.Descent, #fail
  #:adadelta=>Flux.ADADelta, #fail #WARNING: uses non-standard learning parameterization
  :adamax=>Flux.AdaMax, #51s (8.8->3.1)
  :adagrad=>Flux.ADAGrad, #44s (8.1->-1.9(!))
  :amsgrad=>Flux.AMSGrad) #42s->(9.6->-0.3)



function test∇VolumeParts!(loss, X::NamedTuple, y, Θ::VolumePartsLSq)
  @info "test ∇VolumeParts!"
  Π = Flux.params(Θ) #formatted parameters
  @info "length: $(length(Π))"
  gs = gradient(()->loss(X,y),Π)

  rms(v) = sqrt(mean(v.^2))
  vrms(m) = rms(vec(m))

  @info gs[Θ.A]

  @info "test gradrms: Θ.A-$(vrms(gs[Θ.A])); Θ.V₀-$(vrms(gs[Θ.V₀]))"

  ∇Θ! = ∇VolumeParts!(loss, X,y, Θ)

  #out .= reduce(vcat, [∇Θ[Θ.A], vec(∇Θ[Θ.V₀])])#vec∇(Θ, ∇Θ)
  out = similar(Θ.Π)
  ∇Θ!(out, Θ)
  NvecA::Int = length(vec(Θ.A))
  NvecV₀::Int = length(vec(Θ.V₀))
  @info "test ∇Θ!: Θ.A-$(rms(out[1:NvecA]));" *
    " Θ.V₀-$(vrms(out[(NvecA+1):(NvecA + NvecV₀)]))"
end
