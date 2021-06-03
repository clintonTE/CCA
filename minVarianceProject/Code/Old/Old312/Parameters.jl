
abstract type AbstractParameters end

#holds the main parameter object
mutable struct Parameters <: AbstractParameters
  σ²G::Float64
  ζ²G::Float64
  ζ²P::Float64
  SG::Float64
  SGP::Float64
end

#named constructor for testing purposes
function Parameters(;
  σ²G::Float64 = D_σ²G,
  ζ²G::Float64 = D_ζ²G,
  ζ²P::Float64 = D_ζ²P,
  SG::Float64 = D_SG,
  SGP::Float64 = D_SGP)::Parameters

  return Parameters(σ²G, ζ²G, ζ²P, SG, SGP)
end

function Base.print(Θ::Parameters)
  return  "σ²G = $(Θ.σ²G); ζ²G = $(Θ.ζ²G); ζ²P = $(Θ.ζ²P); SG = $(Θ.SG); SGP = $(Θ.SGP)"
end

Base.println(Θ::Parameters)::String = Base.print(Θ) * "\n"

copy2vector(Θ::Parameters)::Vector{Float64} = [Θ.σ²G, Θ.ζ²G, Θ.ζ²P, Θ.SG, Θ.SGP]

Base.names(Θ::Parameters) = ["sigma2G", "zeta2G", "zeta2P", "SG", "SGP"]

#holds the priors
struct Priors <: AbstractParameters
  θG::Float64
  δ²G::Float64
  αG::Float64
  βG::Float64
  αP::Float64
  βP::Float64
  θSG::Float64
  δ²SG::Float64
  θSGP::Float64
  δ²SGP::Float64
end

#named constructor for testing purposes
function Priors(;  θG::Float64 = D_θG,
  δ²G::Float64 = D_δ²G,
  αG::Float64 = D_αG,
  βG::Float64 = D_βG,
  αP::Float64 = D_αP,
  βP::Float64 = D_βP,
  θSG::Float64 = D_θSG,
  δ²SG::Float64 = D_δ²SG,
  θSGP::Float64 = D_θSGP,
  δ²SGP::Float64 = D_δ²SGP)::Priors

  return Priors(θG,  δ²G,  αG,  βG,  αP,  βP,  θSG,  δ²SG,  θSGP,  δ²SGP)
end
