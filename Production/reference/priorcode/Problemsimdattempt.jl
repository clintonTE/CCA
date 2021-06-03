

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

#creates pointers for valid widths
#finding valid widths:
function validsimd(i::Int)
  try
    maximum(ones(i+10)[VecRange{i}(1)])
    return true
  catch
    return false
  end
end


function getsimdptrs(N::Int; start::Int = 0, validsimdwidths::Vector{Int} = deepcopy(VALID_SIMD_WIDTHS))

  #this will hold the simd pointers
  ptrs = Vector{VecRange}()

  #need these to be ordered
  sort!(validsimdwidths)
  targetwidth::Int = pop!(validsimdwidths)
  sizehint!(ptrs, targetwidth + (N ÷ targetwidth))

  local ctr::Int = 1
  while ctr::Int ≤ N
    if N - (ctr-1) ≥ targetwidth
      push!(ptrs, VecRange{targetwidth}(0) + ctr)
      ctr += targetwidth
      #println(ctr)
    else
      targetwidth = pop!(validsimdwidths)
      #println(targetwidth)
    end
  end

  @assert N== (ctr - 1)

  return ptrs
end

#gets the maximum by column of a square matrix. Storage space must be provided
function colmaxsquare2!(A::AbstractMatrix{Float64},
  V::AbstractVector{Float64}, ptrs::Vector{VecRange})::Nothing

  N::Int = size(A,2)
  @assert length(V) == N
  @assert size(A,1) == N

  #local maxv::Vector{Float64}
  #local maxi::Vector{Int}

  #@fastmath @inbounds for c ∈ 1:N
  itemp::Vector{Int} = collect(1:N)

  for c ∈ 1:N
    V[c] = maximum((p->maximum(A[(c-1)*N + p])).(ptrs))
  end

  return nothing
end

#gets the maximum by column of a square matrix. Storage space must be provided
function colmaxsquare!(A::AbstractMatrix{Float64},
  V::AbstractVector{Float64},
  Π::AbstractVector{Int})::Nothing

  N::Int = size(A,2)
  @assert length(V) == N
  @assert length(Π) == N
  @assert size(A,1) == N

  #local maxv::Vector{Float64}
  #local maxi::Vector{Int}

  #@fastmath @inbounds for c ∈ 1:N
  for c ∈ 1:N

    local maxv::Float64 = A[1,c]
    local maxi::Int = 1

    @fastmath @inbounds @simd for r ∈ 2:N
      (maxv < A[r,c]) && (maxv = A[r,c]; maxi=r)
    end

    V[c] = maxv
    Π[c] = maxi
  end

  return nothing
end


function findpolicy!(Φ::Problem)::Problem
  local iter::Int = 0
  local derr::Float64 = 1.0
  #local vrow::Matrix{Float64} = Matrix{Float64}(undef,1,Φ.N) #captures updated values
  #local πrow::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef,1,Φ.N) #captures updated policy
  local vold::Array{Float64,3} = deepcopy(Φ.F.V) #holds the previous value for comparison purposes

  #each view here maps all current states to a single future state
  local rslices = Matrix{SubArray{Float64}}(undef, Φ.S, Φ.S)
  local mslices = Array{SubArray{Float64},3}(undef, Φ.S, Φ.S, Φ.S)
  local vslices = Matrix{SubArray{Float64}}(undef, Φ.S, Φ.S) #will a view for each lagged and current state
  local voldslices = Matrix{SubArray{Float64}}(undef, Φ.S, Φ.S) #will a view for each lagged and current state
  local πslices = Matrix{SubArray{Int}}(undef, Φ.S, Φ.S)
  local ptrs::Vector{VecRange} =  getsimdptrs(Φ.N)

  #make our convenience views
  for sₗ ∈ 1:Φ.S, s ∈ 1:Φ.S
    rslices[sₗ, s] = view(Φ.F.R, :, :, sₗ, s)
    vslices[sₗ, s] = view(Φ.F.V, :, sₗ, s)
    voldslices[sₗ, s] = view(vold, :, sₗ, s)
    πslices[sₗ, s] = view(Φ.F.Π, :, sₗ, s)
    mslices[sₗ, s,:] .= ((sₜ₁::Int)->view(Φ.F.M, :, :, sₗ, s, sₜ₁)).(1:Φ.S)
  end

  while (derr > Φ.dtol) && (iter < Φ.maxiter)
    vold .= Φ.F.V
    for sₗ ∈ 1:Φ.S, s ∈ 1:Φ.S #need to use vec since findmax returns column vectors


      colmaxsquare2!(rslices[sₗ, s] .+
        (Φ.p[1, s] .* mslices[sₗ, s,1] .* vslices[s, 1]) .+
        (Φ.p[2, s] .* mslices[sₗ, s,2] .* vslices[s, 2]),
        vslices[sₗ, s], ptrs)
      #vrow, πrow = findmax(rslices[sₗ, s] .+
      #  (Φ.p[1, s] .* mslices[sₗ, s,1] .* vslices[s, 1]) .+
      #  (Φ.p[2, s] .* mslices[sₗ, s,2] .* vslices[s, 2]), dims=1)

      #need to do this to fix the dimensionality of the returned values from findmax
      #vslices[sₗ, s] .= vec(vrow)
      #πslices[sₗ, s] .= vec(πrow)

    end
    derr = maximum((abs).(vold .- Φ.F.V))


    iter += 1

    (iter % 500 == 0) && println("iter: $iter, derr: $derr")

  end

  M::Matrix{Float64} = Matrix{Float64}(undef, Φ.N, Φ.N)

  for sₗ ∈ 1:Φ.S, s ∈ 1:Φ.S
    println("sl: $sₗ, $s")
    M .= rslices[sₗ, s] .+
      (Φ.p[1, s] .* mslices[sₗ, s,1] .* voldslices[s, 1]) .+
      (Φ.p[2, s] .* mslices[sₗ, s,2] .* voldslices[s, 2])
    display(M[:,2])
    display(Φ.F.V[:,sₗ,s])
    println((c->findfirst(isequal(Φ.F.V[c,sₗ,s]),M[:,c])).(1:Φ.N))
    πslices[sₗ, s] .= (c->findfirst(isequal(Φ.F.V[c,sₗ,s]),M[:,c])[1]).(1:Φ.N)
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
      highpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,1,s]) > 0.0) ? (Φ.F.Π[n,1,s]) : -1.0).(1:Φ.N)
      lowpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,2,s]) > 0.0) ? (Φ.F.Π[n,2,s]) : -1.0).(1:Φ.N)

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

      println("n: $n; k,yh: ($k, $yh); yl, v,p|sₗ=h: $ylgh, $(highvalue[n,:]), $(highpolicy[n,:])")
      #println("n: $n; k,yh: ($k, $yh); yl, p|sₗ=h: $ylgh, $(highpolicy[n,:]); p|sₗ=l: $ylgl, $(lowpolicy[n,:] )")
      #println(Φ.F.Π[n,1,1])
    end
  end

  return Φ
end
