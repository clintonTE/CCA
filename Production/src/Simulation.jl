
struct Simulation
  Φ::Problem
  nsims::Int #number of sims
  burnin::Int
  S::Vector{Int}
  D::Vector{Float64} #dividend realizations
  V::Vector{Float64} #value realizations
  Y::Vector{Float64} #production
  Yh::Vector{Float64} #high state production
  Yl::Vector{Float64} #low state production
  K::Vector{Float64} #capital
  εh::Vector{Float64}
  εl::Vector{Float64}
  I::Vector{Float64} #investment realizations
end


function Simulate!(Ψ::Simulation)
  local s::Int = 1 #current state (index)
  local sₗ::Int = 1#lagged state
  local iψ::Int = Int(round(sqrt(Ψ.Φ.N)/2)) #heuristic to ensure a feasable starting value
  local iψₜ₁::Int #the next policy
  local ph::Float64
  local pl::Float64
  local weights::Vector{Weights} = (s::Int->Weights(vec(Ψ.Φ.p[:,s]))).(1:Ψ.Φ.S)
  local samplefrom::Vector{Int} = [1,2]

  #guess randomly until we get a guaranteed feasable allocation
  local ctr::Int = 1
  while ((maximum(Ψ.Φ.F.R[:,iψ,s,1]) ≤ 0.0) || (maximum(Ψ.Φ.F.R[:,iψ,s,2]) ≤ 0.0)) && ctr < 10_000
    iψ = rand(1:Ψ.Φ.N)
    ctr+=1
  end

  (ctr ≥ 10_000) && error("Could not find feasable starting point")

  for i ∈ (-Ψ.burnin):Ψ.nsims
    #select the current state
    sₗ = s
    s = sample(samplefrom, weights[s])
    iψₜ₁ = Ψ.Φ.F.Π[iψ, sₗ, s] #get the policy as an integer

    if i ≥ 1
      Ψ.S[i] = -s + 2 #in a high state S[i] will be 1, in a low state it will be 0
      ph = weights[sₗ][1]
      Ψ.V[i] = Ψ.Φ.F.V[iψ, sₗ, s] #get the current value
      Ψ.K[i] = Ψ.Φ.F.ψ[1][iψ] #get the current capital level
      Ψ.Yh[i] = Ψ.Φ.F.ψ[2][iψ] #get the high state production

      Ψ.Yl[i] = getyl(Ψ.Yh[i],Ψ.K[i], ph, Ψ.Φ.F) #derive low state production
      Ψ.Y[i] = s==1 ? Ψ.Yh[i] : Ψ.Yl[i]

      Ψ.εh[i] = Ψ.Yh[i] / (Ψ.K[i])^Ψ.Φ.F.γ
      Ψ.εl[i]  = Ψ.Yl[i] / (Ψ.K[i])^Ψ.Φ.F.γ

      Ψ.D[i] = Ψ.Φ.F.R[iψ, iψₜ₁, sₗ, s]

      Ψ.I[i] = Ψ.Φ.F.ψ[1][iψₜ₁] - Ψ.K[i]*(1-Ψ.Φ.F.δ)

    end

    (maximum(Ψ.Φ.F.R[:,iψ,sₗ,s]) ≤ 0.0) && @warn("Non-feasable policy implmented at iteration $i")
    iψ = iψₜ₁
  end

  return Ψ
end

function Simulation(Φ;
  nsims::Int = error("nsims is required"),
  burnin::Int = nsims,
  simulate::Bool = true)::Simulation

  local S::Vector{Int} = Vector{Int}(undef, nsims) #dividend realizations
  local D::Vector{Float64} = Vector{Float64}(undef, nsims) #dividend realizations
  local V::Vector{Float64} = Vector{Float64}(undef, nsims) #value realizations
  local Y::Vector{Float64} = Vector{Float64}(undef, nsims) #current production
  local Yh::Vector{Float64} = Vector{Float64}(undef, nsims) #high state production
  local Yl::Vector{Float64} = Vector{Float64}(undef, nsims) #low state production
  local K::Vector{Float64} = Vector{Float64}(undef, nsims) #capital
  local εh::Vector{Float64} = Vector{Float64}(undef, nsims) #productivity allocations
  local εl::Vector{Float64} = Vector{Float64}(undef, nsims) #productivity allocations
  local I::Vector{Float64} = Vector{Float64}(undef, nsims) #investment

  Ψ::Simulation = Simulation(Φ, nsims, burnin, S, D, V, Y, Yh, Yl, K, εh, εl, I)

  #simulates if true, otherwise just returns the allocated shell
  simulate && Simulate!(Ψ)

  return Ψ
end
