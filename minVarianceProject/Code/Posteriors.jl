
abstract type AbstractPosteriors end

#holds the posterior distributions
mutable struct Posteriors <: AbstractPosteriors
  σ²G!::Function
  Z!::Function
  SG!::Function
  SGP!::Function
end

#NOTE: Each posterior distribution as a corresponding test version

function Posteriors(;
  σ²G!::Function = posteriorσ²G!,
  Z!::Function = posteriorZ!,
  SG!::Function = posteriorSG!,
  SGP!::Function = posteriorSGP!)::Posteriors

  return Posteriors(σ²G!, Z!, SG!, SGP!)
end

function PosteriorsTest(;
  σ²G!::Function = posteriorσ²Gtest!,
  Z!::Function = posteriorZtest!,
  SG!::Function = posteriorSGtest!,
  SGP!::Function = posteriorSGPtest!)::Posteriors

  return Posteriors(σ²G!, Z!, SG!, SGP!)
end

#############
function posteriorσ²G!(Θ::Parameters, Π::Priors)::Nothing

  μG::Float64 = ((-Θ.SGP)*(-Θ.ζ²G + Θ.ζGP) - Θ.SG*(Θ.ζGP - Θ.ζ²P))/(Θ.ζ²G - 2*Θ.ζGP + Θ.ζ²P)
  γG2::Float64 = (-Θ.ζGP^2 + Θ.ζ²G*Θ.ζ²P)/(Θ.ζ²G - 2*Θ.ζGP + Θ.ζ²P)
  ζ²Gstar::Float64 = 1/(1/γG2+1/Π.δ²G)
  μ::Float64 = ζ²Gstar * (μG/γG2+Π.θG/Π.δ²G)

  ζ²Gstar = max(ζ²Gstar, SMALL)

  try
    Θ.σ²G = rand(TruncatedNormal(μ,ζ²Gstar, 0.0, ∞))
  catch err
    println("ζ²Gstar: $ζ²Gstar, μ=$μ, Θ: $Θ")
    error(err)
  end

  return nothing
end

###########PosteriorsZ
function posteriorZ!(Θ::Parameters, Π::Priors)::Nothing

  s11::Float64 = (Θ.SG - Θ.σ²G)^2
  s12::Float64 = (Θ.SGP - Θ.σ²G)*(Θ.SG - Θ.σ²G)
  s22::Float64 = (Θ.SGP - Θ.σ²G)^2
  Σ::Matrix{Float64} = [s11 s12; s12 s22]

  Z::Matrix{Float64} = rand(InverseWishart(Π.ν + 1, Π.Ψ .+ Σ))
  Θ.ζ²G = Z[1,1]
  Θ.ζGP = Z[1,2]
  Θ.ζ²P = Z[2,2]
  return nothing
end


#################
function posteriorSG!(Θ::Parameters, Π::Priors)::Nothing
  μSG = Θ.σ²G+(Θ.ζGP/Θ.ζ²P)*(Θ.SGP-Θ.σ²G)
  γ²SG = Θ.ζ²G-(Θ.ζGP^2/Θ.ζ²P)
  ζ²SGstar::Float64 = 1.0 / (1.0/γ²SG+ 1.0/Π.δ²SG)

  μ::Float64 = (μSG/γ²SG + Π.θSG/Π.δ²SG) * ζ²SGstar

  ζ²SGstar = max(ζ²SGstar, SMALL)


  Θ.SG = rand(TruncatedNormal(μ,ζ²SGstar, 0.0, ∞))

  return nothing
end

##################
function posteriorSGP!(Θ::Parameters, Π::Priors)::Nothing
  μSGP = Θ.σ²G+(Θ.ζGP/Θ.ζ²G)*(Θ.SG-Θ.σ²G)
  γ²SGP = Θ.ζ²P-(Θ.ζGP^2/Θ.ζ²G)
  ζ²SGPstar::Float64 = 1.0 / (1.0/γ²SGP+ 1.0/Π.δ²SGP)

  μ::Float64 = (μSGP/γ²SGP + Π.θSGP/Π.δ²SGP) * ζ²SGPstar

  #println("SGTest## μSG: $μSGP  γ²SG: $γ²SGP ζ²SGstar: $ζ²SGPstar")
  ζ²SGPstar = max(ζ²SGstar, SMALL)
  Θ.SGP = rand(TruncatedNormal(μ, ζ²SGPstar, 0.0, ∞))

  return nothing
end


##############****************************#################################
##############****************************#################################
##############****************************#################################
##############****************************#################################

function posteriorσ²Gtest!(Θ::Parameters, Π::Priors)::Nothing

  μG::Float64 = ((-Θ.SGP)*(-Θ.ζ²G + Θ.ζGP) - Θ.SG*(Θ.ζGP - Θ.ζ²P))/(Θ.ζ²G - 2*Θ.ζGP + Θ.ζ²P)
  γG2::Float64 = (-Θ.ζGP^2 + Θ.ζ²G*Θ.ζ²P)/(Θ.ζ²G - 2*Θ.ζGP + Θ.ζ²P)



  ζ²Gstar::Float64 = 1/(1/γG2+1/Π.δ²G)
  μ::Float64 = ζ²Gstar * (μG/γG2+Π.θG/Π.δ²G)

  Θ.σ²G = rand(TruncatedNormal(μ, ζ²Gstar, 0.0, ∞))
  return nothing
end

function posteriorZtest!(Θ::Parameters, Π::Priors)::Nothing

  Σ::Matrix{Float64} =
    [(Θ.SG-Θ.σ²G)^2 (Θ.SG-Θ.σ²G)*(Θ.SGP-Θ.σ²G); (Θ.SG-Θ.σ²G)*(Θ.SGP-Θ.σ²G) (Θ.SGP-Θ.σ²G)^2]

  Z::Matrix{Float64} = rand(InverseWishart(Π.ν + 1, Π.Ψ + Σ))
  Θ.ζ²G = Z[1,1]
  Θ.ζGP = Z[1,2]
  Θ.ζ²P = Z[2,2]
  return nothing
end

function posteriorSGtest!(Θ::Parameters, Π::Priors)::Nothing
  ρSG::Float64 = Θ.ζGP/(Θ.ζ²G*Θ.ζ²P)^0.5
  μSG::Float64 = Θ.σ²G+(Θ.ζ²G/Θ.ζ²P)^0.5*ρSG*(Θ.SGP-Θ.σ²G)
  γ²SG::Float64 = (1-ρSG^2)*Θ.ζ²G

  ζ²SGstar::Float64 = 1/(1/γ²SG+1/Π.δ²SG)
  μ::Float64 = ζ²SGstar*(μSG/γ²SG+Π.θSG/Π.δ²SG)



  Θ.SG = rand(TruncatedNormal(μ, ζ²SGstar, 0.0, ∞))

  return nothing
end

function posteriorSGPtest!(Θ::Parameters, Π::Priors)::Nothing
  ρSGP::Float64 = Θ.ζGP/(Θ.ζ²G*Θ.ζ²P)^0.5
  μSGP::Float64 = Θ.σ²G+(Θ.ζ²P/Θ.ζ²G)^0.5*ρSGP*(Θ.SG-Θ.σ²G)
  γ²SGP::Float64 = (1-ρSGP^2)*Θ.ζ²P

  ζ²SGPstar::Float64 = 1/(1/γ²SGP+1/Π.δ²SGP)
  μ::Float64 = ζ²SGPstar*(μSGP/γ²SGP+Π.θSGP/Π.δ²SGP)

  #println("SGTest## μSG: $μSGP  γ²SG: $γ²SGP ζ²SGstar: $ζ²SGPstar")

  Θ.SGP = rand(TruncatedNormal(μ, ζ²SGPstar, 0.0, ∞))

  return nothing
end
