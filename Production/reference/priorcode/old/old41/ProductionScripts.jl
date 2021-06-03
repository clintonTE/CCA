using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end
using Revise, Finometrics, Random

try #seems to be a bug with running projects locally
  using Production
catch err
  using Production
end

Random.seed!(11)

const REDEFINE = true

#this is to avoid annoying redefinition warnings when possible
if REDEFINE || !(@isdefined DEFINED)
  const STATES_PER_VAR = 20
  const TOL = 1.0e-10
  const D_TOL = 1e-7
  const MAX_ITER = 10000

  #***************competitive parameters**************************
  const COM_τ = [0.0, 0.0]
  const COM_δ = 0.2
  const COM_Θ = (0.7, 0.3)
  const COM_α = 1.05
  const COM_β = 0.98
  const COM_γ = Float64(1/3)
  const COM_ζ = (1.0, 0.0)

  const COM_p = [0.81 0.68; 0.19 0.32]
  const COM_M = 2

  #some heuristics for the maximum values
  const COM_KMAX = (COM_Θ[1]/(COM_p[1,1]^(1/COM_α)*COM_δ))^(1/(1-COM_γ))
  const COM_YMAX = COM_KMAX^COM_γ*COM_Θ[1]/COM_p[1,2]^(1/COM_α)

  const COM_YHSTATES = collect(range(TOL, stop=COM_YMAX, length=STATES_PER_VAR))
  const COM_KSTATES = collect(range(TOL, stop=COM_KMAX, length=STATES_PER_VAR))


  #***************monopolistic parameters**************************
  const MON_τ = COM_τ
  const MON_δ = COM_δ
  const MON_Θ = COM_Θ
  const MON_α = COM_α
  const MON_β = COM_β
  const MON_γ = COM_γ
  const MON_ζ = (3.0, 1.0)


  const MON_p = COM_p
  const MON_M = COM_M

  #technically this sould be based on a sepearte
  const MON_YHSTATES = COM_YHSTATES
  const MON_KSTATES = COM_KSTATES

  const DEFINED = true
end

function productionscript(firmtype::Symbol)::Nothing
  local F::Firm
  local Φ::Problem

  #the problem parameters
  local τ::Vector{Float64}
  local δ::Float64
  local Θ::NTuple{2,Float64}
  local α::Float64
  local β::Float64
  local γ::Float64
  local kstates::Vector{Float64}
  local yhstates::Vector{Float64}
  local ζ::NTuple{2,Float64}
  local p::Matrix{Float64}

  #assign values based on the scenario
  if firmtype ∈ (:competitive, :COM, :com)
    τ= COM_τ
    δ= COM_δ
    Θ= COM_Θ
    α= COM_α
    β= COM_β
    γ= COM_γ
    kstates= COM_KSTATES
    yhstates= COM_YHSTATES
    ζ= COM_ζ
    p= COM_p
  elseif firmtype ∈ (:monpolistic, :MON, :mon)
    τ = MON_τ
    δ = MON_δ
    Θ = MON_Θ
    α = MON_α
    β= MON_β
    γ= MON_γ
    kstates = MON_KSTATES
    yhstates = MON_YHSTATES
    ζ= MON_ζ
    p = MON_p
  else
    error("firmtype $firmtype not recognized")
  end

  println("\n#**********RUNNING $firmtype using $STATES_PER_VAR states/var***********#")
  F = Firm(τ=τ, δ=δ, Θ=Θ, α=α, β=β, γ=γ, kstates=kstates, yhstates=yhstates, ζ=ζ, p=p)
  Φ = Problem(F, p, dtol = D_TOL, maxiter = MAX_ITER)

  return nothing
end

@time productionscript(:COM)
