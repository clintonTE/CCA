
struct Firm
  F::Function #production as a function of capital
  τ::Vector{Float64} #tax
  δ::Float64 #depreciation

  Θ::Vector{Float64} #technology per state
  π::Matrix{CartesianIndex} #policy states
  v::Matrix{Float64} #value states

  R::Array{Float64,3}
end

#gets the cobb-douglas produciton function
function getCD(γ::Float64)::Function
  F(k::Float64)::Float64 = k^γ

  return F
end

function Firm(τ::Vector{Float64} = DEF_τ,
  δ::Float64 = DEF_δ,
  Θ::Vector{Float64} = DEF_Θ,
  γ::Float64 = DEF_γ,
  kstates::Vector{Float64} = DEF_KSTATES,
  ystates::Vector{Float64} = DEF_YSTATES)

  kystates = combine(kstates, ystates)
  N::Int = length(kystates)
  M::Int = length(Θ)
  F = getCD(γ)

  π::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef, N, M)
  v::Matrix{Float64} = Matrix{Float64}(undef, N, M)
  R::Array{Float64,3} = Array{Float64,3}(undef, N, N, M)

  return Firm(F, τ, δ, Θ, π, v, R)
end
