#using Revise
#using CUDA, LinearAlgebra, BenchmarkTools, Zygote, SparseArrays, UnPack
#import Base.*
#abstract type AbstractAnnhilator{T} end

#creates a pre-calculated projection object for the dummy matrix in a fixed effect regression
#this approach should balance speed and performance
struct FixedAnnihilator{TF,T} <: AbstractAnnhilator{T}
  F::TF
  FCFt::TF

  FixedAnnihilator(F::TF,FCFt::TF) where {TF} = new{TF, eltype(TF)}(F,FCFt)
end

FixedAnnihilator(F::TF) where {TF<:SparseMatrixCSC} = FixedAnnihilator(F,cholesky!(Matrix(F'*F))\I*F' |> TF)
function FixedAnnihilator(F::TF, ::Type{<:Any}, ::Type{T} = eltype(TF)) where {TF,TA,T}
  (maximum(F) - T(1) == minimum(F) == zero(T)) || error("Fixed annihilation matrix only
    makes sense with a dummy matrix")
  return FixedAnnihilator(SparseMatrixCSC{T}(F))
end

#version for cuarrays
function FixedAnnihilator(F::TF, ::Type{<:CuArray}, ::Type{T} = eltype(TF)) where {TF,T}
  (maximum(F) - T(1) == minimum(F) == zero(T)) || error("Fixed annihilation matrix only
    makes sense with a dummy matrix")
    M = FixedAnnihilator(F |> Matrix{T}, Matrix{T})
    Fsp = M.F |> CUSPARSE.CuSparseMatrixCSR
    FCFtsp = M.FCFt |> CUSPARSE.CuSparseMatrixCSR

  return FixedAnnihilator(Fsp,FCFtsp)
end



#function *(M::FixedAnnihilator,X::TX, ::Type{T}=eltype(TX))::TX where {T, TX<:AbstractArray}
#  return X .- M.F*(M.FCFt*X)
  #CUSPARSE.mm!('N',T(-1.0),M.F, M.FCFt* X, T(1.0),deepcopy(X),'O')
#end

function *(M::FixedAnnihilator,X::TX, ::Type{T}=eltype(TX))::TX where {T, TX<:AbstractArray}
  #return deepcopy(X)
  return X .- M.F*(M.FCFt*X)
end

#the below is really fast but fragile, hopefully it will get less fragile over time
function *(M::FixedAnnihilator, X::TX, ::Type{T}=eltype(TX))::TX where {T, TX<:CuMatrix}
  #return deepcopy(X)
  return CUSPARSE.mm!('N','N',T(-1.0),M.F, M.FCFt* X, T(1.0),deepcopy(X),'O')
end

#special inplace version for a slight performance boost
function multiply!(M::FixedAnnihilator, X::TX, ::Type{T}=eltype(TX))::TX where {T, TX<:CuMatrix}
  #return X
  return CUSPARSE.mm!('N','N',T(-1.0),M.F, M.FCFt* X, T(1.0),X,'O')
end

#the above high speed operation works much better on a matrix than a vector
function *(M::FixedAnnihilator, y::Ty, ::Type{T}=eltype(Ty))::Ty where {T, Ty<:CuVector}
  #return deepcopy(y)
  return y .- M.F*(M.FCFt*y)#CUSPARSE.mv!('N',T(-1.0),M.F, M.FCFt * y, T(1.0),deepcopy(y),'O')
end

function multiply!(M::FixedAnnihilator, y::Ty, ::Type{T}=eltype(Ty))::Ty where {T, Ty<:CuVector}
  #return y
  return CUSPARSE.mv!('N',T(-1.0),M.F, M.FCFt * y, T(1.0),y,'O')
end

#fallback for the in-place multiply
function multiply!(M::FixedAnnihilator, X::TX, ::Type{T}=eltype(TX)) where {TX,T}
  #return X
  X .= M*X
  return X
end

function testFixedAnnihilator(F::TF, ::Type{T}=eltype(TF)) where {TF,T}
  M = FixedAnnihilator(F)
  v = rand(T, size(F,1))
  X= rand(T, size(F,1), 5)
  @assert Matrix(M.F) ≈ F
  @assert M*v ≈ v .- F*((cholesky!(Matrix(F'*F))\I) * F' * v)
  @assert M*X ≈ X .- F*((cholesky!(Matrix(F'*F))\I) * F' * X)

  return nothing
end

###############################
####Dense Annihilator
###############################


#creates a pre-calculated dense projection object
struct DenseAnnihilator{TW,T} <: AbstractAnnhilator{T}
  W::TW
  WtWinvWt::TW

  DenseAnnihilator(W::TW,WtWinvWt::TW) where {TW} = new{TW, eltype(TW)}(W,WtWinvWt)
#  DenseAnnihilator(::Nothing) where {TW} = new{Nothing, Nothing}(nothing, nothing)
end

function DenseAnnihilator(W::Matrix)

  try
    #return DenseAnnihilator(W, cholesky!(Hermitian(W'*W))\(W'))
    return DenseAnnihilator(W, svd(Symmetric(W'*W))\(W'))
  catch err
    @warn "ERROR! printing printmln(W'*W):"
    printmln(W'*W)
    error("Error: $err")
  end
end

function DenseAnnihilator(W::CuMatrix{T}) where T
  Wcpu = W |> Matrix{T}
  try
    #return DenseAnnihilator(W, cholesky!(Hermitian(W'*W))\(W'))
    return DenseAnnihilator(W, svd(Symmetric(Wcpu'*Wcpu))\(Wcpu') |> CuMatrix{T})
  catch err
    @warn "ERROR! printing printmln(W'*W):"
    printmln(W'*W)
    println("Printing diag of SVD")
    printmln(svd(W'*W).S |> Diagonal)
    println("Printing inverse:")
    printmln(inv(svd(W'*W |> Matrix)))
    if size(W,1) < 10^6
      W |> Matrix |> DataFrame |> CSV.write("$(PARAM[:testpath])\\error_W_$(size(W)).csv")
    end
    error("Error: $err")
  end
end

function *(M::DenseAnnihilator{TW,T}, X::TX)::TX where {TW, TX<:AbstractArray, T}

  return X - M.W*(M.WtWinvWt*X)
end

#special case for if we have no
#*(M::DenseAnnihilator{Nothing,Nothing}, X::TX)::TX where {TX<:AbstractArray} = X

#expand the annihilator matrix
Base.Matrix{T}(M::DenseAnnihilator{TW,T}) where {TW<:AbstractArray,T} = I - M.W * M.WtWinvWt
function CUDA.CuMatrix{T}(M::DenseAnnihilator{TW,T}) where {TW<:CuArray, T}
  @unpack W, WtWinvWt = M

  (N,K) = size(W)
  id = CuMatrix{T}(I, N, N)
  return id - W * WtWinvWt
end



function testdenseannihilator(::Type{TW}, ::Type{T}=eltype(TW); N=2500, K=10, bench=true) where {TW,T}
  W = rand(N, 2K) |> TW
  X = rand(N, K) |> TW
  M = DenseAnnihilator(W)
  A = M |> TW

  ver = X - W*(cholesky(W'*W)\(W'*X)) |> Matrix
  @assert ver ≈ ( (M * X) |> Matrix)
  @assert ver ≈ (A * X) |> Matrix
  verf(X) = A * X |> sum
  testf(X) = M * X |> sum
  @assert (gradient(verf, X)[1] |> Array) ≈ (gradient(testf, X)[1] |> Array)

  if bench
    print("Testing M*X for type $TW:")
    @btime CUDA.@sync $M*$X
    print("Testing A*X for type $TW:")
    @btime CUDA.@sync $A*$X
    print("Testing ∇M*X for type $TW:")
    @btime CUDA.@sync gradient($testf, $X)[1]
    print("Testing ∇A*X for type $TW:")
    @btime CUDA.@sync gradient($verf, $X)[1]
  end

  @info "Passed testdenseannihilator for type $TW"
end


  #test vectorized application of projections
function testdenseannihilator(
    blocks::Int, ::Type{TW}, ::Type{T}=eltype(TW); N=2500, K=10, bench=true)  where {TW,T}

  vecW = [rand(T, N + rand((-K^2):K^2), 2K) |> TW for b ∈ 1:blocks]
  vecX = [rand(T, size(vecW[b],1), K) |> TW  for b ∈ 1:blocks]
  vecM = (DenseAnnihilator).(vecW)

  for b in 1:blocks
    W = vecW[b]
    X = vecX[b]

    @assert (vecM[b] * vecX[b]) |> Matrix ≈ (X - W*(cholesky(W'*W)\(W'*X)) |> Matrix)
  end


  #  (W,X,zb)->Threads.@spawn((W*X) .* zb)).() .|> fetch .|> sum |> sum
  vecmult(z::Tz, ::Type{TW}) where {Tz, TW <: CuMatrix} = (
    (W,X,zb)->(W*X) .* zb |> sum).(vecM, vecX,z) |> sum

  vecmult(z::Tz, ::Type{TW}) where {Tz, TW <: Matrix} = ((W,X,zb)->(W*X) .* zb).(vecM, vecX,z) .|> sum |> sum
  ∇vecmult(z) = gradient(vecmult,z)[1]

  z = [rand(T, 1, K) |> TW  for b ∈ 1:blocks] #testing placeholder

  #warmup- make sure everything works
  CUDA.reclaim()
  @assert !(∇vecmult(z) === nothing)
  @assert reduce(vcat, ∇vecmult(z) .|> Array) .|> abs |> sum > 0.0

  if bench
    print("Testing vecmult for type $TW:")
    @btime CUDA.@sync $vecmult($z)
    print("Testing ∇vecmult for type $TW:")
    @btime CUDA.@sync $∇vecmult($z)
  end
end


function testdenseannihilator(blocks::Int = 1100; N=2500, K=10, bench=false, benchblocks=false)
  testdenseannihilator(Matrix{Float64}; N, K, bench)
  testdenseannihilator(CuMatrix{Float64}; N, K, bench)
  testdenseannihilator(blocks, Matrix{Float64}; N, K, bench=benchblocks)
  testdenseannihilator(blocks, CuMatrix{Float64}; N, K,bench= benchblocks)
end

#sleep(0.2)
#CUDA.allowscalar(false)
#testdenseannihilator(1100, bench=false, benchblocks=true, N=2500, K=10)
