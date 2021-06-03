using Revise
using CUDA, LinearAlgebra, Zygote, BenchmarkTools



#=
NOTE: below is verbatim adjoint
https://github.com/FluxML/Zygote.jl/blob/4262dbaaa22ddc55a5389ddc5fd696059797c6a7/src/lib/array.jl
@adjoint function cholesky(Σ::Union{StridedMatrix, Symmetric{<:Real, <:StridedMatrix}}; check = true)
  C = cholesky(Σ, check = check)
  return C, function(Δ::NamedTuple)
    issuccess(C) || throw(PosDefException(C.info))
    U, Ū = C.U, Δ.factors
    Σ̄ = similar(U.data)
    Σ̄ = mul!(Σ̄, Ū, U')
    Σ̄ = copytri!(Σ̄, 'U')
    Σ̄ = ldiv!(U, Σ̄)
    Σ̄ = BLAS.trsm!('R', 'U', 'T', 'N', one(eltype(Σ)), U.data, Σ̄)
    Σ̄[diagind(Σ̄)] ./= 2
    return (UpperTriangular(Σ̄),)
  end
end=#

Zygote.@adjoint function \(A::Cholesky, B::AbstractMatrix)
  Y, back = Zygote.pullback((U, B)->U \ (U' \ B), A.U, B)
  return Y, function(Ȳ)
    Ā_factors, B̄ = back(Ȳ)
    return ((uplo=nothing, info=nothing, factors=Ā_factors), B̄)
  end
end

#NOTE: new adjoint
Zygote.@adjoint function cholesky(Σ::CuMatrix; check = true)
  C = cholesky(Σ, check = check)
  return C, function(Δ::NamedTuple)
    issuccess(C) || throw(PosDefException(C.info))
    U, Ū = C.U, Δ.factors
    Σ̄ = similar(U.data)
    Σ̄ = mul!(Σ̄, Ū, U')
    Σ̄ = Zygote.copytri!(Σ̄, 'U')
    Σ̄ = ldiv!(U, Σ̄)
    Σ̄ = CUBLAS.trsm!('R', 'U', 'T', 'N', one(eltype(Σ)), U.data, Σ̄)
    Σ̄[diagind(Σ̄)] ./= 2
    return (UpperTriangular(Σ̄),)
  end
end


function ∇XtXinvXty(X,y)
  tf(X,y) = sum(cholesky(X'*X)\ (X' * y))
  tf(X,y)
  ∇tf = gradient((X)->tf(X,y),X)[1]

end

function choltest(N=50, K=5)
  X = rand(N,K)
  y = rand(N)

  dX = X |> CuMatrix{Float64}
  dy = y |> CuVector{Float64}

  ∇chol = ∇XtXinvXty(X,y)
  d∇chol = ∇XtXinvXty(dX,dy)
  @assert d∇chol |> Matrix ≈ ∇chol
  @info "cholesky test passed"
end

function cholbench(N=10^7, K=3)
  X = rand(N,K)
  y = rand(N)

  dX = X |> CuMatrix{Float64}
  dy = y |> CuVector{Float64}

  @info "benching CPU"
  ∇chol = @btime ∇XtXinvXty($X,$y)

  @info "benching gpu"
  d∇chol = CUDA.@sync @btime ∇XtXinvXty($dX,$dy)
  @assert d∇chol |> Matrix ≈ ∇chol
  @info "cholesky test passed"
end


choltest()
