
abstract type AbstractPosteriors end

#holds the posterior distributions
mutable struct Posteriors <: AbstractPosteriors
  σ²G!::Function
  ζ²G!::Function
  ζ²P!::Function
  SG!::Function
  SGP!::Function
end

#NOTE: Each posterior distribution as a corresponding test version

function Posteriors(;
  σ²G!::Function = posteriorσ²G!,
  ζ²G!::Function = posteriorζ²G!,
  ζ²P!::Function = posteriorζ²P!,
  SG!::Function = posteriorSG!,
  SGP!::Function = posteriorSGP!)::Posteriors

  return Posteriors(σ²G!, ζ²G!, ζ²P!, SG!, SGP!)
end

function PosteriorsTest(;
  σ²G::Function = posteriorσ²Gtest!,
  ζ²G::Function = posteriorζ²Gtest!,
  ζ²P::Function = posteriorζ²Ptest!,
  SG::Function = posteriorSGtest!,
  SGP::Function = posteriorSGPtest!)::Posteriors

  return Posteriors(σ²G, ζ²G, ζ²P, SG, SGP)
end

#############
function posteriorσ²G!(Θ::Parameters, Π::Priors)::Nothing

  ζ²Gstar::Float64 = 1.0 / (1.0/Θ.ζ²G + 1.0/Θ.ζ²P + 1.0/Π.δ²G)

  μ::Float64 = (Θ.SG/Θ.ζ²G + Θ.SGP/Θ.ζ²P + Π.θG/Π.δ²G) * ζ²Gstar

  Θ.σ²G = rand(TruncatedNormal(μ, ζ²Gstar, 0.0, ∞))

  return nothing
end


#################
function posteriorSG!(Θ::Parameters, Π::Priors)::Nothing
  ζ²SGstar::Float64 = 1.0 / (1.0/Θ.ζ²G + 1.0/Π.δ²SG)

  μ::Float64 = (Θ.σ²G/Θ.ζ²G + Π.θSG/Π.δ²SG) * ζ²SGstar

  Θ.SG = rand(TruncatedNormal(μ, ζ²SGstar, 0.0, ∞))

  return nothing
end

###########################
function posteriorSGP!(Θ::Parameters, Π::Priors)::Nothing
  ζ²SGPstar::Float64 = 1.0 / (1.0/Θ.ζ²G + 1.0/Π.δ²SGP)

  μ::Float64 = (Θ.σ²G/Θ.ζ²G + Π.θSGP/Π.δ²SGP) * ζ²SGPstar

  #Θ.SGP = rand(TruncatedNormal(μ, ζ²SGPstar, 0.0, ∞))
  Θ.SGP = rand(Normal(μ, ζ²SGPstar))

  return nothing
end


#########################
function posteriorζ²G!(Θ::Parameters, Π::Priors)::Nothing
  α::Float64 = Π.αG + 0.5
  β::Float64 = Π.βG + (Θ.SG-Θ.σ²G)^2/2

  Θ.ζ²G = rand(InverseGamma(α, β))

  return nothing
end


#########################
function posteriorζ²P!(Θ::Parameters, Π::Priors)::Nothing
  α::Float64 = Π.αP + 0.5
  β::Float64 = Π.βP + (Θ.SGP-Θ.σ²G)^2/2

  Θ.ζ²P = rand(InverseGamma(α, β))

  return nothing
end

##############****************************#################################
##############****************************#################################
##############****************************#################################
##############****************************#################################

function posteriorσ²Gtest!(Θ::Parameters, Π::Priors)::Nothing

  Θ.σ²G = rand(Normal((Θ.SG/Θ.ζ²G+Θ.SGP/Θ.ζ²P+Π.θG/Π.δ²G)/(1/Θ.ζ²G+1/Θ.ζ²P+1/Π.δ²G),
    1.0/(1/Θ.ζ²G+1/Θ.ζ²P+1/Π.δ²G)))

  return nothing
end

function posteriorSGtest!(Θ::Parameters, Π::Priors)::Nothing

  Θ.SG = rand(Normal((Θ.σ²G/Θ.ζ²G+Π.θSG/Π.δ²SG)/(1/Θ.ζ²G+1/Π.δ²SG),
    1.0/(1/Θ.ζ²G+1/Π.δ²SG)))

  return nothing
end


function posteriorSGPtest!(Θ::Parameters, Π::Priors)::Nothing

  Θ.SGP = rand(Normal((Θ.σ²G/Θ.ζ²G+Π.θSGP/Π.δ²SGP)/(1/Θ.ζ²G+1/Π.δ²SGP),
    1.0/(1/Θ.ζ²G+1/Π.δ²SGP)))

  return nothing
end

function posteriorζ²Gtest!(Θ::Parameters, Π::Priors)::Nothing

  Θ.ζ²G = rand(InverseGamma(Π.αG+0.5, Π.βG + (Θ.SG - Θ.σ²G)^2/2))

  return nothing
end

function posteriorζ²Ptest!(Θ::Parameters, Π::Priors)::Nothing

  Θ.ζ²P = rand(InverseGamma(Π.αP + 0.5, Π.βP + (Θ.SGP - Θ.σ²G)^2/2))

  return nothing
end
