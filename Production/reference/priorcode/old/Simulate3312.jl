
struct Simulation
  Φ::Problem
  nsims::Int #number of sims
  burnin::Int
  D::Vector{Float64} #dividend realizations
  V::Vector{Float64} #value realizations
  Y::Vector{Float64} #production
  Yₕ::Vector{Float64} #high state production
  Yₗ::Vector{Float64} #low state production
  K::Vector{Float64} #capital
end


function Simulation(Φ; nsims::Int = DEF_N_SIMS, burnin::Int)

  local s::Int = 1 #current state (index)
  local ls::Int #lagged state
  local iψ::Int = length(iψ) ÷ 2 #current policy state
  local liψ::Int #lagged policy state
  local ph::Float64


  local weights::Vector{Vector{Float64}} = (s::Int->vec(p[:,s])).(1:Φ.S)

  local V::Vector{Float64} = Vector{Float64}(undef, nsims) #value realizations
  local Y::Vector{Float64} = Vector{Float64}(undef, nsims) #current production
  local Yh::Vector{Float64} = Vector{Float64}(undef, nsims) #high state production
  local Yl::Vector{Float64} = Vector{Float64}(undef, nsims) #low state production
  local K::Vector{Float64} = Vector{Float64}(undef, nsims) #capital

  local samplefrom = (1,2,3)

  for i ∈ (-burnin):nsims
    #select the current state
    s = sample(samplefrom, weights[s])

    if i ≥ 1
      ph::Float64 =
      V[i] = Φ.F.V[iψ, s] #get the current value
      K[i] = ψ[1][iψ] #get the current capital level
      Yh[i] = ψ[2][iψ] #get the high state production
      Yl[i] = getyl(Yh[i],K[i], weights[s][1], 1 -weights[s][1], Φ.F) #derive low state production
      Y[i] = i==1 ? Yh[i] : Yl[i]
    end

    iψ = Φ.Π[iψ, s]
  end



end
