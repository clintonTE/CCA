

mutable struct Problem
  F::Firm

  Nr::Int
  Nc::Int
  S::Int
  p::Matrix{Float64}

  dtol::Float64
  maxiter::Int
end

function updatefirmτ!(Φ::Problem, τ::Vector{Float64})::Problem
  Φ.F.τ .= τ
  refreshR!(Φ.F, Φ.p)

  return Φ
end

function findpolicy!(Φ::Problem)::Problem
  local rslices::Vector{SubArray} #view of r in each state
  local vslices::Vector{SubArray} #view of v in each state
  local πslices::Vector{SubArray} #view of the policy in each state
  local mslices::Matrix{SubArray} #view of the pricing kernel in each state and each future state
  local iter::Int = 0
  local derr::Float64 = 1.0
  local vrow::Matrix{Float64} = Matrix{Float64}(undef,1,Φ.Nc) #captures updated values
  local πrow::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef,1,Φ.Nc) #captures updated policy
  local vold::Matrix{Float64} = deepcopy(Φ.F.V) #holds the previous value for comparison purposes
  local RVt::Matrix{Float64} = Matrix{Float64}(undef,Φ.Nr,Φ.Nc) #used to build the choice set for each state

  #each view here maps all current states to a single future state
  rslices = Vector{SubArray}(undef, Φ.S)
  mslices = Matrix{SubArray}(undef, Φ.S, Φ.S)
  vslices = Matrix{SubArray}(undef, Φ.S, Φ.S) #will a view for each lagged and current state
  πslices = Vector{SubArray}(undef, Φ.S)

  #make our convenience views
  for s ∈ 1:Φ.S
    rslices[s] = view(Φ.F.R, :, :,  s)
    vslices[:, s] = ((sₗ::Int)->view(Φ.F.V,  s)).(1:Φ.S)
    πslices[s] = view(Φ.F.Π, :, s)
    mslices[s,:] .= ((sₜ₁::Int)->view(Φ.F.M, :, :, s, sₜ₁)).(1:Φ.S)
  end

  #display(rslices[1] .+
  #  Φ.p[1, 1] .* mslices[1,1] .* vslices[1] .+
  #  Φ.p[2, 1] .* mslices[1,2] .* vslices[2] .+
  #  Φ.p[3, 1] .* mslices[1,3] .* vslices[3])


  while (derr > Φ.dtol) && (iter < Φ.maxiter)
    for s ∈ 1:Φ.S #need to use vec since findmax returns column vectors
      RVt .= rslices[s]
      vrow, πrow = findmax(rslices[s] .+
        (Φ.p[1, s] .* mslices[s,1] .* vslices[1]') .+
        (Φ.p[2, s] .* mslices[s,2] .* vslices[2]'), dims=1)

      #need to do this to fix the dimensionality of the returned values from findmax
      vslices[s] .= vec(vrow)
      πslices[s] .= vec(πrow)
      #isplay(vrow)
      #=println(maximum(Φ.p[1, s] .* mslices[s,1] .+
        Φ.p[2, s] .* mslices[s,2] .+
        Φ.p[3, s] .* mslices[s,3]))=#
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
    dtol::Float64=DEF_D_TOL, maxiter::Int = DEF_MAX_ITER)::Problem
  Nr::Int = length(F.ψr[1])
  Nc::Int = length(F.ψc[1])
  S::Int = size(p,1)

  Φ::Problem = Problem(F, Nr, Nc, S, p, dtol, maxiter)
  findpolicy!(Φ)

  #display(Φ.F.V)
  if Φ.Nr < 1000 #checks if the policy function is feasable to display
    highpolicy::Matrix{Float64} = Matrix{Int}(undef, Φ.Nr, Φ.S)
    lowpolicy::Matrix{Float64} = Matrix{Int}(undef, Φ.Nr, Φ.S)

    highvalue::Matrix{Float64} = Matrix{Float64}(undef, Φ.Nr, Φ.S)
    lowvalue::Matrix{Float64} = Matrix{Float64}(undef, Φ.Nr, Φ.S)
    for s ∈ 1:Φ.S
      highpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,s]) > 0.0) ? (Φ.F.Π[n,s])[1] : -1.0).(1:Φ.Nr)
      lowpolicy[:,s] .= ((n)->(maximum(Φ.F.R[:,n,s]) > 0.0) ? (Φ.F.Π[n,s])[1] : -1.0).((Φ.Nr+1):Φ.Nc)

      highvalue[:,s] = ((n)->(maximum(Φ.F.R[:,n,s]) > 0.0) ? (Φ.F.V[n,s]) : -1.0).(1:Φ.Nr)
      lowvalue[:,s]  = ((n)->(maximum(Φ.F.R[:,n,s]) > 0.0) ? (Φ.F.V[n,s]) : -1.0).((Φ.Nr+1):Φ.Nc)

      highvalue[:,s] .= (v::Float64->round(v,digits=1)).(highvalue[:,s])
      lowvalue[:,s] .= (v::Float64->round(v,digits=1)).(lowvalue[:,s])
    end

    for n ∈ 1:Φ.Nr
      k::Float64 = round(Φ.F.ψc[1][n], digits=2)
      yh::Float64 = round(Φ.F.ψc[2][n], digits=2)

      ylgh::Float64 = getyl(Φ.F.ψc[2][n], Φ.F.ψc[1][n], Φ.p[1,1], 1- Φ.p[1,1], Φ.F)
      ylgl::Float64 = getyl(Φ.F.ψc[2][n], Φ.F.ψc[1][n], Φ.p[1,2], 1- Φ.p[1,2], Φ.F)

      ylgh = round(ylgh, digits=2)
      ylgl = round(ylgl, digits=2)

      println("n: $n; k,yh: ($k, $yh); yl, v,p|sₗ=h: $ylgh, $(highvalue[n,:]), $(highpolicy[n,:])")
      #println("n: $n; k,yh: ($k, $yh); yl, p|sₗ=h: $ylgh, $(highpolicy[n,:]); p|sₗ=l: $ylgl, $(lowpolicy[n,:] )")
    end
  end

  return Φ
end
