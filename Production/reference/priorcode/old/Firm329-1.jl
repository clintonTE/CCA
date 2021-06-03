
struct Firm{T <: NTuple}
  P::Function #inverse demand function
  ϕ::Function #production as a function of capital
  τ::T #tax
  α::Float64
  γ::Float64
  δ::Float64 #depreciation
  Θ::T #technology per state

  ψ::Vector{NTuple} #kroniker product of all possible state values
  π::Matrix{CartesianIndex} #policy states
  v::Matrix{Float64} #value states

  R::Array{Float64,3}
  M::Array{T,3}
end


function getlinearprice(ζ::NTuple{2,Float64})::Function
  local P::Function

  ζ[2] > 0.0 && @warn "Upward sloping demand curve! If not intentional, reparameterize ζ[1]"

  #set the function, allowing for a linear demand curve
  P(Y::Float64) = ζ[1] + ζ[2] * Y
  return P
end

#gets the cobb-douglas produciton function
function getCD(γ::Float64)::Function
  ϕ(k::Float64)::Float64 = k^γ

  return ϕ
end

function getyl(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  kern::Float64 = (k^(F.γ*F.α)/pl) - (ph/pl)* (yh/F.Θ[1])^F.α
  yl::Float64 = (kern > 0.0) ? (kern)^(1/F.α)*F.Θ[2] : 0.0

  return yl
end

function getYlcheck(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  Yl::Float64 = ((1/pl)*(k^(F.γ*F.α)-ph*(yh/F.Θ[1])^F.α))^(1/F.α)*F.Θ[2]
end

function getyh(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  kern::Float64 = (k^(F.γ*F.α)/pl) - (ph/pl)* (yh/F.Θ[1])^F.α
  yl::Float64 = (kern > 0.0) ? (kern)^(1/F.α)*F.Θ[2] : 0.0

  return yl
end

function getyhcheck(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  Yl::Float64 = ((1/pl)*(k^(F.γ*F.α)-ph*(yh/F.Θ[1])^F.α))^(1/F.α)*F.Θ[2]
end


function getkernel(pₜ::Float64, F::T where T<:Firm)::T
end

function Firm(;τ::T = DEF_τ,
  δ::Float64 = DEF_δ,
  Θ::T = DEF_Θ,
  α::Float64=DEF_α,
  γ::Float64 = DEF_γ,
  kstates::Vector{Float64} = DEF_KSTATES,
  yhstates::Vector{Float64} = DEF_YSTATES,
  ζ::T= DEF_ζ,
  p::Matrix{Float64} = DEF_p)::Firm{T} where T<:NTuple #for performance reasons, don't use vectors

  ψ::Vector{NTuple} = combine(kstates, yhstates)
  N::Int = length(ψ)
  S::Int = length(Θ)

  P::Function = getlinearprice(ζ)
  ϕ::Function = getCD(γ)

  π::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef, N, S)
  v::Matrix{Float64} = zeros(N,S)
  R::Array{Float64,3} = Array{Float64,3}(undef, N, N, S)
  M::Array{T,3} = fill(T(zeros(S)), N, N, S) #preallocate

  F::Firm{T} = Firm{T}(P, ϕ, τ, α, γ, δ, Θ, ψ, π, v, R, M)

  #println("getYl: ",getyl(1.4, 2.3, .3, .9, F)) #NOTE: Checks this formula
  #println("getYlcheck: ",getylcheck(1.4, 2.3, .3, .9, F)) #NOTE: Checks this formula

  #now build up R
  for s::Int ∈ 1:S
    for r ∈ 1:N, c ∈ 1:N
      (kₜ₁::Float64, yₜ₁s::Float64) =  ψ[r]
      (kₜ::Float64, yₜ::Float64) =  ψ[c]
      τₜ::Float64 = τ[s]

      #get the low state for the subsequent period
      yₜ₁l::Float64 = getyl(yₜ₁, kₜ₁, p[1,s], p[2,s], F::Firm)

      #need to get the current low state production if we are in the low state
      pₜ::Float64 = P(yₜ)

      dₜ::Float64 = pₜ*yₜ*(1-τₜ)-(kₜ₁-(1-δ)*kₜ)+τₜ*δ*kₜ


      F.R[r,c,s] = dₜ


    end
  end

  return F
end
