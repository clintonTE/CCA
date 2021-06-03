using Revise

using DataFrames, Statistics, StatsBase, LinearAlgebra, BenchmarkTools, Random, CUDA
CUDA.allowscalar(false)
Random.seed!(11)



function StatsBase.partialcor(X::AbstractArray{T,2}, Z::AbstractArray{T,2};
    zincludesintercept = false) where T
  N,K = size(X)
  @assert size(X,1) == size(Z,1)

  local Xe::Matrix
  if zincludesintercept
    Xₑ = X - Z * (cholesky!(Symmetric(Z'*Z))\(Z'*X))
  else
    Z1 = hcat(ones(T,N), Z)
    Xₑ = X - Z1 * (cholesky!(Symmetric(Z1'*Z1))\(Z1'*X))
  end
  return cor(Xₑ)
end

function StatsBase.pairwisepartialcor(X::AbstractArray{T,2}, Z::AbstractArray{T,2};
    zincludesintercept = false) where T
  N,K = size(X)
  @assert size(X,1) == size(Z,1)

  local Xe::Matrix
  if zincludesintercept
    Xₑ = X - Z * (cholesky!(Symmetric(Z'*Z))\(Z'*X))
  else
    Z1 = hcat(ones(T,N), Z)
    Xₑ = X - Z1 * (cholesky!(Symmetric(Z1'*Z1))\(Z1'*X))
  end
  return cor(Xₑ)
end

function testpartialcor(;N=252, K=500, R = 4)
  X = rand(N, K)
  Z = rand(N, R)
  #slow, manual version
  function partialcorver(X::AbstractArray{T,2}, Z::AbstractArray{T,2}) where T
    N,K = size(X)
    @assert size(X,1) == size(Z,1)

    ρ = ones(T, K,K)
    for c ∈ 1:K
      for r ∈ 1:(c-1)
        ρ[r,c] = partialcor(X[:,c], X[:,r],Z)
      end
    end

    ρ .= Symmetric(ρ)
    return ρ
  end

  #create an orthogonalized version for testing purposes
  Z1 = hcat(ones(N), Z)
  Xr = X .- Z1*(cholesky!(Z1'*Z1)\(Z1'*X))

  @info "testing partial correlations"
  @assert (partialcor(Xr, Z) .≈ cor(Xr)) |> all
  @assert (partialcor(X,Z) .≈ partialcor(X,Z1,zincludesintercept=true)) |> all
  #@btime partialcor($X, $Z1,zincludesintercept=true)
  @btime partialcor($X, $Z)

  if K ≤ 1000
    @assert (partialcorver(X,Z) .≈ cor(Xr)) |> all
    @assert (partialcor(X,Z) .≈ partialcorver(X,Z)) |> all
  end

end



@time testpartialcor(N=252, K=5000, R=4)
