
abstract type AbstractParameters end

#holds the main parameter object
mutable struct Parameters <: AbstractParameters
  σ²G::Float64
  ζ²G::Float64
  ζ²P::Float64
  ζGP::Float64
  SG::Float64
  SGP::Float64 #not an actual
end

#named constructor for testing purposes
function Parameters(;
  σ²G::Float64 = D_σ²G,
  ζ²G::Float64 = D_ζ²G,
  ζ²P::Float64 = D_ζ²P,
  ζGP::Float64 = D_ζGP,
  SG::Float64 = D_SG,
  SGP::Float64 = D_SGP)::Parameters

  return Parameters(σ²G, ζ²G, ζ²P, ζGP, SG, SGP)
end

function Base.print(Θ::Parameters)
  return  "σ²G = $(Θ.σ²G); ζ²G = $(Θ.ζ²G); ζ²P = $(Θ.ζ²P);
    ζGP = $(Θ.ζGP); SG = $(Θ.SG); SGP = $(Θ.SGP)"
end

Base.println(Θ::Parameters)::String = Base.print(Θ) * "\n"

copy2vector(Θ::Parameters)::Vector{Float64} = [Θ.σ²G, Θ.ζ²G, Θ.ζ²P, Θ.ζGP, Θ.SG, Θ.SGP]

Base.names(Θ::Parameters) = ["sigma2G", "zeta2G", "zeta2P", "zetaGP", "SG", "SGP"]

#holds the priors
struct Priors <: AbstractParameters
  θG::Float64
  δ²G::Float64
  Ψ::Matrix{Float64}
  ν::Float64
  θSG::Float64
  δ²SG::Float64
  θSGP::Float64
  δ²SGP::Float64
end

#named constructor for testing purposes
function Priors(;priortype::Symbol = :default)

  local θG::Float64
  local δ²G::Float64
  local Ψ::Matrix{Float64}
  local ν::Float64
  local θSG::Float64
  local δ²SG::Float64
  local θSGP::Float64
  local δ²SGP::Float64

  if priortype == :default || priortype == :informative
    θG = D_θG
    δ²G = D_δ²G
    Ψ = D_Ψ
    ν = D_ν
    θSG = D_θSG
    δ²SG = D_δ²SG
    θSGP = D_θSGP
    δ²SGP = D_δ²SGP
  elseif priortype == :uninformative
    θG = U_θG
    δ²G = U_δ²G
    Ψ = U_Ψ
    ν = U_ν
    θSG = U_θSG
    δ²SG = U_δ²SG
    θSGP = U_θSGP
    δ²SGP = U_δ²SGP
  else
    error("Prior type $priortype not recognized.")
  end

  return Priors(θG, δ²G, Ψ, ν, θSG, δ²SG, θSGP, δ²SGP)
end
