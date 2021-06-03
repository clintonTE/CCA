

struct Firm{T <: NTuple}
  ζ::NTuple #inverse demand function
  τ::T #tax
  α::Float64
  γ::Float64 #production share
  δ::Float64 #depreciation
  Θ::T #technology per state

  ψ::Vector{NTuple} #kroniker product of all possible state values
  π::Matrix{CartesianIndex} #policy states
  v::Matrix{Float64} #value states

  R::Array{Float64,3}
  M::Array{T,3}
end

function getyl(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  kern::Float64 = (k^(F.γ*F.α)/pl) - (ph/pl)* (yh/F.Θ[1])^F.α
  yl::Float64 = (kern ≥ 0.0) ? (kern)^(1/F.α)*F.Θ[2] : NaN

  return yl
end




#pass the following parameters to get the pricing kernel
function getkernel(s::Int, yₜ::Float64, yₜ₁h::Float64, yₜ₁l::Float64, p::Matrix{Float64},
  F::T where T<:Firm)::NTuple{3,Float64}

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

  Mₜ₁h::Float64 = Pₜ*ξʰ*(1-τl)/(Mdenom)
  Mₜ₁l::Float64 = Pₜ*ξˡ*(1-τh)/(Mdenom)

  return (Mₜ₁h, Mₜ₁l, Mₜ₁l)

end


function Firm(;  ζ::NTuple= DEF_ζ,
  τ::T = DEF_τ,
  δ::Float64 = DEF_δ,
  Θ::T = DEF_Θ,
  α::Float64=DEF_α,
  γ::Float64 = DEF_γ,
  kstates::Vector{Float64} = DEF_KSTATES,
  yhstates::Vector{Float64} = DEF_YSTATES,
  p::Matrix{Float64} = DEF_p)::Firm{T} where T<:NTuple #for performance reasons, don't use vectors

  #rules that need to be followed for this approach to work
  @assert Θ[2] == Θ[3] #tech level has to be the same between states
  @assert τ[2] == τ[3] #tax is identical between the two states
  @assert sum(p[:,2] .≠ p[:,3]) == 0 #transition probabilities are the same if originating from low
  @assert sum(p) == 3.0 #probabilities add to 1
  @assert ζ[2] ≥ 0.0 #demand curve is weakly downward sloping

  #make some type assertions to help with debugging and optimization
  local τₜ::Float64
  local yₜ::Float64 #will hold the true produciton, not just the high state production
  local yₜ₁l::Float64 #future low stae productivity (fully specified by high state prod and state)
  local yₜl::Float64 #current low state productivitiy

  local pₜ₁h::Float64 #probability of transitioning to a high state
  local pₜ₁l::Float64 #probability of transitioning to a low state
  local pₜh::Float64 #probability of having transitioned from a high state
  local pₜl::Float64 #probability of having transitioned from a low state

  local dₜ::Float64 #the dividend payment
  local Mₜ₁::NTuple{3, Float64} #pricing kernel given the current state

  local Pₜ::Float64 #current price
  local Pₜ₁h::Float64 #future price in high state
  local Pₜ₁l::Float64 #future price in low state

  local allocationₜ::Float64 #used to check the feasability of the current allocaiton
  local allocationₜ₁::Float64 #used to check the feasability of the future allocaiton

  ψ::Vector{NTuple} = combine(kstates, yhstates)
  #println(ψ)
  N::Int = length(ψ)
  S::Int = length(Θ)

  π::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef, N, S)
  v::Matrix{Float64} = zeros(N,S)
  R::Array{Float64,3} = Array{Float64,3}(undef, N, N, S)
  M::Array{T,3} = fill(T(zeros(S)), N, N, S) #preallocate

  F::Firm{T} = Firm{T}(ζ, τ, α, γ, δ, Θ, ψ, π, v, R, M)

  #println("getYl: ",getyl(1.4, 2.3, .3, .9, F)) #NOTE: Checks this formula
  #println("getYlcheck: ",getylcheck(1.4, 2.3, .3, .9, F))

  #println("getkernel: ",getkernel(1, 2.5, 2.7, 3.4, p, F)) #NOTE: Checks this formula
  #println("getkernel: ",getkernelcheck(1, 2.5, 2.7, 3.4, p, F), "\n")

  #println("getkernel: ",getkernel(2, 2.5, 2.7, 3.4, p, F)) #NOTE: Checks this formula
  #println("getkernel: ",getkernelcheck(2, 2.5, 2.7, 3.4, p, F), "\n")

  #println("getkernel: ",getkernel(3, 2.5, 2.7, 3.4, p, F)) #NOTE: Checks this formula
  #println("getkernel: ",getkernelcheck(3, 2.5, 2.7, 3.4, p, F), "\n")

  #now build up R
  for s::Int ∈ 1:S
    for r ∈ 1:N, c ∈ 1:N
      #println("$(ψ[r])")
      (kₜ₁::Float64, yₜ₁h::Float64) =  ψ[r]
      (kₜ::Float64, yₜh::Float64) =  ψ[c]
      τₜ = τ[s]

      # compute the future low state production from the high state production
      pₜ₁h = p[1,s]
      pₜ₁l = 1.0 - p[1,s]
      yₜ₁l = getyl(yₜ₁h, kₜ₁, pₜ₁h, pₜ₁l, F::Firm)

      #compute the current low income state
      pₜh = (s == 2) ? p[1,1] : p[1,2]
      pₜl = 1.0 - pₜh
      yₜl = getyl(yₜh, kₜ, pₜh, 1 - pₜh, F::Firm)

      if s == 1 #need to recover the previously selected low productivity
        yₜ = yₜh
      else
        yₜ = yₜl
      end

      yₜ = max(yₜ, N_TOL)
      yₜ₁l = max(yₜ₁l, N_TOL)
      yₜ₁h = max(yₜ₁h, N_TOL)


      #prices
      Pₜ = ζ[1] - ζ[2] * yₜ
      Pₜ₁h = ζ[1] - ζ[2] * yₜ₁h
      Pₜ₁l = ζ[1] - ζ[2] * yₜ₁l

      if isnan(yₜl) || isnan(yₜ₁l) || (Pₜ<0) || (Pₜ₁h<0) || (Pₜ₁l<0)
        Mₜ₁ = (0.0,0.0,0.0)
        dₜ = 0.0
      else


      #allocationₜ = (pₜh*(yₜh/Θ[1])^α + pₜl*(yₜl/Θ[2])^α)^(1/α)

      #check for negative prices or a non-feasable current allocation
      #if (Pₜ < 0) || (kₜ^γ + TOL_CONSTRAINT < allocationₜ)
      #  yₜ = 0.0
      #end
      #compute the dividend and record as the return
        dₜ = Pₜ*yₜ*(1-τₜ)-(kₜ₁-(1-δ)*kₜ)+τₜ*δ*kₜ


      #check for future non-feasability and/or negative prices

      #allocationₜ₁ = (pₜ₁h*(yₜ₁h/Θ[1])^α + pₜ₁l*(yₜ₁l/Θ[2])^α)^(1/α)

      #if (yₜ ≤ 0) || (Pₜ₁h < 0) || (Pₜ₁l < 0) || (kₜ₁^γ + TOL_CONSTRAINT < allocationₜ₁)
      #  Mₜ₁ = (0.0,0.0,0.0)
      #else
        #println("(s, yₜ, yₜ₁h, yₜ₁l): $((s, yₜ, yₜ₁h, yₜ₁l))")
        Mₜ₁ = getkernel(s, yₜ, yₜ₁h, yₜ₁l, p, F)
      #end
      end

      #if sum((isnan).([Mₜ₁...;dₜ])) ≠ 0
        #println("ERROR: Mₜ₁: $Mₜ₁, dₜ: $dₜ")
      #  println("(s, yₜ, yₜ₁h, yₜ₁l, p):", "$((s, yₜ, yₜ₁h, yₜ₁l, p))")
      #end
      @assert sum((isnan).([Mₜ₁...;dₜ])) == 0 #no NaNs allowed

      F.R[r,c,s] = dₜ
      F.M[r,c,s] = Mₜ₁


    end
  end

  #display(F.M[:,:,1])

  return F
end



#=
function getkernelcheck(s::Int, yₜ::Float64, yₜ₁h::Float64, yₜ₁l::Float64, p::Matrix{Float64},
  F::T where T<:Firm)::NTuple{3,Float64}

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

  #println("Mdenom: $((ξ*(1-F.τ[1])*(1-F.τ[2]) + pₜ₁l*ξˡ*κˡ*(1-F.τ[1]) + pₜ₁h*ξʰ*κʰ*(1-F.τ[2])))")
  return (Mₜ₁h, Mₜ₁l, Mₜ₁l)
end=#

#=
function getylcheck(yh::Float64, k::Float64, ph::Float64, pl::Float64, F::Firm)
  Yl::Float64 = ((1/pl)*(k^(F.γ*F.α)-ph*(yh/F.Θ[1])^F.α))^(1/F.α)*F.Θ[2]
end=#
