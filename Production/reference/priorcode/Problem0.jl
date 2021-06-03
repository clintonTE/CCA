

mutable struct Problem
  F::Firm

  N::Int
  S::Int
  p::Matrix{Float64}

  dtol::Float64
  maxiter::Int
end

function updatefirmτ!(Φ::Problem, τ::Vector{Float64})::Problem
  Φ.F.τ .= τ
  refreshR!(Φ.F, Φ.p)
  findpolicy!(Φ::Problem)
  return Φ
end

function findpolicy!(Φ::Problem)::Problem
  local rslices::Matrix{SubArray} #view of r in each state
  local vslices::Matrix{SubArray} #view of v in each state
  local πslices::Matrix{SubArray} #view of the policy in each state
  local mslices::Array{SubArray,3} #view of the pricing kernel in each state and each future state
  local iter::Int = 0
  local derr::Float64 = 1.0
  local vrow::Matrix{Float64} = Matrix{Float64}(undef,1,Φ.N) #captures updated values
  local πrow::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef,1,Φ.N) #captures updated policy
  local vold::Array{Float64,3} = deepcopy(Φ.F.V) #holds the previous value for comparison purposes

  #each view here maps all current states to a single future state
  rslices = Matrix{SubArray}(undef, Φ.S, Φ.S)
  mslices = Array{SubArray,3}(undef, Φ.S, Φ.S, Φ.S)
  vslices = Matrix{SubArray}(undef, Φ.S, Φ.S) #will a view for each lagged and current state
  πslices = Matrix{SubArray}(undef, Φ.S, Φ.S)

  #make our convenience views
  for sₗ ∈ 1:Φ.S, s ∈ 1:Φ.S
    rslices[sₗ, s] = view(Φ.F.R, :, :, sₗ, s)
    vslices[sₗ, s] = view(Φ.F.V, :, sₗ, s)
    πslices[sₗ, s] = view(Φ.F.Π, :, sₗ, s)
    mslices[sₗ, s,:] .= ((sₜ₁::Int)->view(Φ.F.M, :, :, sₗ, s, sₜ₁)).(1:Φ.S)
  end

  #these are useful for the BLAS libraries
  C::Matrix{Float64} = zeros(Φ.N, Φ.N)
  β::Float64 = 0.0
  N1::Vector{Float64} = ones(Φ.N)
  N²::Int = Φ.N*Φ.N

  while (derr > Φ.dtol) && (iter < Φ.maxiter)
    for sₗ ∈ 1:Φ.S, s ∈ 1:Φ.S #need to use vec since findmax returns column vectors

      vrow, πrow = findmax(rslices[sₗ, s] .+
        (Φ.p[1, s] .* mslices[sₗ, s,1] .* vslices[s, 1]) .+
        (Φ.p[2, s] .* mslices[sₗ, s,2] .* vslices[s, 2]), dims=1)

      #need to do this to fix the dimensionality of the returned values from findmax
      vslices[sₗ, s] .= vec(vrow)
      πslices[sₗ, s] .= vec(πrow)

    end

    derr = maximum((abs).(vold .- Φ.F.V))
    vold .= Φ.F.V

    iter += 1

    (iter % 500 == 0) && println("iter: $iter, derr: $derr")

  end

  (iter == Φ.maxiter) && @warn("Convergence failure iter: $iter, derr: $derr")


  return Φ
end


function Problem(F::Firm{NTuple{2,Float64}}, p::Matrix{Float64};
    dtol::Float64=DEF_D_TOL, maxiter::Int = DEF_MAX_ITER,
    verbose::Bool = length(F.ψ[1]) < 1000)::Problem

  N::Int = length(F.ψ[1])
  S::Int = size(p,1)

  Φ::Problem = Problem(F, N, S, p, dtol, maxiter)
  findpolicy!(Φ)

  #display(Φ.F.V)
  if verbose
    highpolicy::Matrix{Float64} = Matrix{Int}(undef, Φ.N, Φ.S)
    lowpolicy::Matrix{Float64} = Matrix{Int}(undef, Φ.N, Φ.S)

    highvalue::Matrix{Float64} = Matrix{Float64}(undef, Φ.N, Φ.S)
    lowvalue::Matrix{Float64} = Matrix{Float64}(undef, Φ.N, Φ.S)
    for s ∈ 1:Φ.S
      highpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,1,s]) > 0.0) ? (Φ.F.Π[n,1,s])[1] : -1.0).(1:Φ.N)
      lowpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,2,s]) > 0.0) ? (Φ.F.Π[n,2,s])[1] : -1.0).(1:Φ.N)

      highvalue[:,s] = ((n)->(maximum(Φ.F.R[:,n,1,s]) > 0.0) ? (Φ.F.V[n,1,s]) : -1.0).(1:Φ.N)
      lowvalue[:,s]  = ((n)->(maximum(Φ.F.R[:,n,2,s]) > 0.0) ? (Φ.F.V[n,2,s]) : -1.0).(1:Φ.N)

      highvalue[:,s] .= (v::Float64->round(v,digits=1)).(highvalue[:,s])
      lowvalue[:,s] .= (v::Float64->round(v,digits=1)).(lowvalue[:,s])
    end

    for n ∈ 1:Φ.N
      k::Float64 = round(Φ.F.ψ[1][n], digits=2)
      yh::Float64 = round(Φ.F.ψ[2][n], digits=2)

      ylgh::Float64 = getyl(Φ.F.ψ[2][n], Φ.F.ψ[1][n], Φ.p[1,1], Φ.F)
      ylgl::Float64 = getyl(Φ.F.ψ[2][n], Φ.F.ψ[1][n], Φ.p[1,2], Φ.F)

      ylgh = round(ylgh, digits=2)
      ylgl = round(ylgl, digits=2)

      #println("n: $n; k,yh: ($k, $yh); yl, v,p|sₗ=h: $ylgh, $(highvalue[n,:]), $(highpolicy[n,:])")
      println("n: $n; k,yh: ($k, $yh); yl, p|sₗ=h: $ylgh, $(highpolicy[n,:]); p|sₗ=l: $ylgl, $(lowpolicy[n,:] )")
    end
  end

  return Φ
end
