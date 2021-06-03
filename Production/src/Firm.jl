

struct Firm{T <: NTuple}
  ζ::NTuple #inverse demand function
  τ::Vector{Float64} #tax
  α::Float64
  β::Float64 #time preference parameter
  γ::Float64 #production share
  δ::Float64 #depreciation
  Θ::T #technology per state

  ψ::Vector{Vector{U} where U <: Union{Float64, Int}} #kroniker product of all possible k and yh

  #indexing: [current k & yh, past s, current s]
  Π::Array{Int, 3} #policy states
  V::Array{Float64,3} #value states

  #holds the return in each state
  #indexing: [future k & yh, current k & yh, past s, current s]
  R::Array{Float64,4}

  #holds the pricing kernel in each state
  #indexing: [future k & yh, current k & yh, past s, current s, future s]
  M::Array{Float64,5}
end

function getyl(yh::Float64, k::Float64, ph::Float64, F::Firm)
  local pl::Float64 = 1.0 - ph
  local kern::Float64 = (k^(F.γ*F.α)/pl) - (ph/pl) * (yh/F.Θ[1])^F.α
  local yl::Float64 = (kern ≥ 0.0) ? (kern)^(1/F.α)*F.Θ[2] : NaN

  return yl
end


#pass the following parameters to get the pricing kernel
function getkernel(s::Int, yₜ::Float64, yₜ₁h::Float64, yₜ₁l::Float64, p::Matrix{Float64},
  F::T where T<:Firm)::NTuple{2,Float64}

  #derive the remaining parameters
  #hopefully this will aid in debugging
  ζ₀::Float64 = F.ζ[1]
  ζ₁::Float64 = F.ζ[2]
  Pₜ::Float64 = ζ₀ - ζ₁*yₜ

  pₜ₁h::Float64 = p[1,s]
  pₜ₁l::Float64 = 1 - pₜ₁h

  Θₜ::Float64 = F.Θ[s]
  Θh::Float64 = F.Θ[1]
  Θl::Float64 = F.Θ[2]

  τh::Float64 = F.τ[1]
  τl::Float64 = F.τ[2]

  #1/(F.ζ[1]-F.ζ[2]*yₜ₁h) * (yₜ₁h/yₜ)^(F.α-1) * (F.Θ[1]/F.Θ[s])^(-F.α)
  #now compute the intermediate quantities
  ξʰ::Float64 = (1/(ζ₀ - ζ₁*yₜ₁h)) * (yₜ₁h/yₜ)^(F.α-1) * (Θh/Θₜ)^(-F.α)
  ξˡ::Float64 = (1/(ζ₀ - ζ₁*yₜ₁l)) * (yₜ₁l/yₜ)^(F.α-1) * (Θl/Θₜ)^(-F.α)
  ξ::Float64 = ξʰ*pₜ₁h + ξˡ*pₜ₁l

  κʰ::Float64 = (ζ₀ - 2*ζ₁*yₜ₁h)*τh
  κˡ::Float64 = (ζ₀ - 2*ζ₁*yₜ₁l)*τl

  Mdenom::Float64 = ξ*(1-τh)*(1-τl) + pₜ₁l*κˡ*ξˡ*(1-τh) + pₜ₁h*κʰ*ξʰ*(1-τl)

  Mₜ₁h::Float64 = max(Pₜ*ξʰ*(1-τl)/(Mdenom),0.0) #no arbitrage implies no negative kernal
  Mₜ₁l::Float64 = max(Pₜ*ξˡ*(1-τh)/(Mdenom),0.0)

  rendog::Float64 = 1.0/(pₜ₁h*Mₜ₁h + pₜ₁l*Mₜ₁l)
  Mₜ₁h *= F.β * rendog
  Mₜ₁l *= F.β * rendog

  return (Mₜ₁h, Mₜ₁l)

end


function fillRCol!(F::T, p::Matrix{Float64}, sₗ::Int, s::Int, c::Int)::T where T<:Firm
  N::Int = length(F.ψ[1])
  @assert s ≤ size(p,1)
  @assert c ≤ N

  for r ∈ 1:N #2x performance gain or more
    #make some type assertions to help with debugging and optimization
    local τₜ::Float64
    local yₜ::Float64 #will hold the true produciton, not just the high state production
    local yₜ₁l::Float64 #future low stae productivity (fully specified by high state prod and state)
    local yₜl::Float64 #current low state productivitiy

    local pₜ₁h::Float64 #probability of transitioning to a high state
    local pₜ₁l::Float64 #probability of transitioning to a low state
    local pₜh::Float64 #probability of having transitioned from a high state

    local dₜ::Float64 #the dividend payment
    local Mₜ₁::NTuple{2, Float64} #pricing kernel given the current state

    local Pₜ::Float64 #current price
    local Pₜ₁h::Float64 #future price in high state
    local Pₜ₁l::Float64 #future price in low state

    local allocationₜ::Float64 #used to check the feasability of the current allocaiton
    local allocationₜ₁::Float64 #used to check the feasability of the future allocaiton

    local kₜ₁::Float64 #kt+1
    local yₜ₁h::Float64 #yt+1 high
    local kₜ::Float64 #kt
    local yₜh::Float64 #yt high

    local feasable::Bool

    kₜ₁ = F.ψ[1][c]
    yₜ₁h = F.ψ[2][c]
    kₜ = F.ψ[1][r]
    yₜh = F.ψ[2][r]

    τₜ = F.τ[s]

    # compute the future low state production from the high state production
    pₜ₁h = p[1,s]
    pₜ₁l = 1.0 - pₜ₁h
    yₜ₁l = getyl(yₜ₁h, kₜ₁, pₜ₁h, F::Firm)

    #compute the current low income state
    pₜh = (sₗ == 1) ? p[1,1] : p[1,2]
    yₜl = getyl(yₜh, kₜ, pₜh, F::Firm)

    if s == 1 #need to recover the previously selected low productivity
      yₜ = yₜh
    else
      yₜ = yₜl
    end

    yₜ = max(yₜ, N_TOL)
    yₜ₁l = max(yₜ₁l, N_TOL)
    yₜ₁h = max(yₜ₁h, N_TOL)


    #prices
    Pₜ = F.ζ[1] - F.ζ[2] * yₜ
    Pₜ₁h = F.ζ[1] - F.ζ[2] * yₜ₁h
    Pₜ₁l = F.ζ[1] - F.ζ[2] * yₜ₁l

    if (isnan(yₜl) || isnan(yₜ₁l) || (Pₜ<0) || (Pₜ₁h<0) || (Pₜ₁l<0)) #make sure the choice is feasable
      feasable = false
    else
      dₜ = Pₜ*yₜ*(1-τₜ)-(kₜ₁-(1-F.δ)*kₜ)+τₜ*F.δ*kₜ
      feasable = (dₜ ≥ 0) #no negative dividends
    end

    if feasable
      Mₜ₁ = getkernel(s, yₜ, yₜ₁h, yₜ₁l, p, F)
    else
      Mₜ₁ = (0.0,0.0) #non-feasable states have no future value
      dₜ = PENALTY #large penalty for non-feasable values
    end

    F.R[r, c, sₗ, s] = dₜ
    F.M[r, c, sₗ, s, :] .= Mₜ₁

  end

  return F
end

function refreshR!(F::T, p::Matrix{Float64})::T where T<:Firm
  N::Int = length(F.ψ[1])
  S::Int = length(F.Θ)

  for sₗ ∈ 1:S, s ∈ 1:S, c ∈ 1:N
    fillRCol!(F, p, sₗ, s, c)
  end

  #drange = (f->Int(round(f))).(range(5, stop=N, length=10))
  #display(F.R[(max(end-10,1)):end ,(max(end-10,1)):end ,1,1] .- F.R[(max(end-10,1)):end ,(max(end-10,1)):end ,1,2])
  #display(F.R[(max(end-10,1)):end ,(max(end-10,1)):end ,1,1] .- F.R[(max(end-10,1)):end ,(max(end-10,1)):end ,1,2])
  return F
end


function Firm(;
    ζ::NTuple= error("ζ is required"),
    τ::Vector{Float64} = error("τ is required"),
    δ::Float64 = error("δ is required"),
    Θ::T = error("Θ is required"),
    α::Float64=error("α is required"),
    β::Float64 = error("β is required"),
    γ::Float64 = error("γ is required"),
    kstates::Vector{Float64} = error("kstates is required"),
    yhstates::Vector{Float64} = error("yhstates is required"),
    p::Matrix{Float64} = error("p is required"))::Firm{T} where T<:NTuple #for performance reasons, don't use vectors

  #rules that need to be followed for this approach to work
  @assert length(τ) == length(Θ) #tax is identical between the two states
  @assert sum(p) == 2.0 #probabilities add to 1
  @assert ζ[2] ≥ 0.0 #demand curve is weakly downward sloping

  S::Int = length(Θ)
  ψ::Vector{Vector{U} where U <: Union{Float64, Int}} = combine2vectors(kstates, yhstates)

  #display(ψ)

  N::Int = length(ψ[1])


  Π::Array{Int, 3} = Array{Int, 3}(undef, N, S, S)
  V::Array{Float64,3} = zeros(N,S,S)
  R::Array{Float64,4} = Array{Float64,4}(undef, N, N, S, S)
  M::Array{Float64,5} = Array{Float64,5}(undef, N, N, S, S, S) #preallocate

  F::Firm{T} = Firm{T}(ζ, τ, α, β, γ, δ, Θ, ψ, Π, V, R, M)

  #build up R
  refreshR!(F, p)

  #println("sums of columns $(sum(F.R, dims=1))")
  #println("sums of rows $(sum(F.R, dims=1))")

  #println("getYl: ",getyl(1.4, 2.3, .3, .9, F)) #NOTE: Checks this formula
  #println("getYlcheck: ",getylcheck(1.4, 2.3, .3, .9, F))

  #println("getkernel: ",getkernel(1, 2.5, 2.7, 3.4, p, F)) #NOTE: Checks this formula
  #println("getkernel: ",getkernelcheck(1, 2.5, 2.7, 3.4, p, F), "\n")
  #display(F.M[:,:,2,2])
  return F
end


#=
function getkernelcheck(s::Int, yₜ::Float64, yₜ₁h::Float64, yₜ₁l::Float64, p::Matrix{Float64},
  F::T where T<:Firm)

  #define derived parameters
  pₜ₁h = p[1,s]
  pₜ₁l = 1 - pₜ₁h
  Pₜ::Float64 = F.ζ[1] - F.ζ[2]*yₜ


  ξʰ::Float64 = 1/(F.ζ[1]-F.ζ[2]*yₜ₁h) * (yₜ₁h/yₜ)^(F.α-1) * (F.Θ[1]/F.Θ[s])^(-F.α)
  ξˡ::Float64 = 1/(F.ζ[1]-F.ζ[2]*yₜ₁l) * (yₜ₁l/yₜ)^(F.α-1) * (F.Θ[2]/F.Θ[s])^(-F.α)
  ξ::Float64 = pₜ₁h*ξʰ + pₜ₁l*ξˡ
  κʰ::Float64 = (F.ζ[1]-2*F.ζ[2]*yₜ₁h)*F.τ[1]
  κˡ::Float64 = (F.ζ[1]-2*F.ζ[2]*yₜ₁l)*F.τ[2]

  Mₜ₁l::Float64 = (Pₜ*ξˡ*(1-F.τ[1])/
    (ξ*(1-F.τ[1])*(1-F.τ[2]) + pₜ₁l*ξˡ*κˡ*(1-F.τ[1]) + pₜ₁h*ξʰ*κʰ*(1-F.τ[2])))
  Mₜ₁h::Float64 = ((Pₜ-pₜ₁l*κˡ*Mₜ₁l)*ξʰ /
    (ξ*(1-F.τ[1]) + pₜ₁h*ξʰ*κʰ))

  rendog::Float64 = 1.0 / (Mₜ₁h*pₜ₁h + pₜ₁l*Mₜ₁l)
  Mₜ₁h *= F.β * rendog
  Mₜ₁l *= F.β * rendog
  #println("Mdenom: $((ξ*(1-F.τ[1])*(1-F.τ[2]) + pₜ₁l*ξˡ*κˡ*(1-F.τ[1]) + pₜ₁h*ξʰ*κʰ*(1-F.τ[2])))")
  return (Mₜ₁h, Mₜ₁l)
end=#



#=
function getylcheck(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  Yl::Float64 = ((1/pl)*(k^(F.γ*F.α)-ph*(yh/F.Θ[1])^F.α))^(1/F.α)*F.Θ[2]
end=#
