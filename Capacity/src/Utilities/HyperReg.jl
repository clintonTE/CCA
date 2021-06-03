
#code is borrowed liberally from CUDA.jl/lib/cusolver/wrappers.jl
#                            and CUDA.jl/lib/cublas/wrappers.jl
for (rfname, rsname, T) in
    ((:cusolverDnSpotrfBatched,:cusolverDnSpotrsBatched,:Float32),
    (:cusolverDnDpotrfBatched,:cusolverDnDpotrsBatched,:Float64))
    @eval begin
function hyperreg(X::Vector{<:CuMatrix{$T}},
                      y::Vector{<:CuMatrix{$T}})
    cuuplo = CUDA.CUBLAS.cublasfill('U')

    if length(X) != length(y)
        throw(DimensionMismatch(""))
    end

    mx,nx = size(X[1])
    my,ny = size(y[1])

    if mx != my
        throw(ArgumentError("X Matrix must have same num of rows as y matrix"))
    end

    for (Xs, ys) in zip(X,y)
        if (size(Xs) != (mx, nx)) || (size(ys) != (my, ny))
            throw(DimensionMismatch("Dimensions of batched array entries must be invariant"))
        end
    end

    A = CUDA.CUBLAS.gemm_batched('T','N', X, X) #X'X
    B = CUDA.CUBLAS.gemm_batched('T','N', X, y) #X'y

    #Here we compute cholesky!(X'X)
    n = size(A[1],1)
    lda = max(1,stride(A[1],2))
    Aptrs = CUDA.CUBLAS.unsafe_batch(A)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUSOLVER.$rfname(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, Aptrs, lda, infoarray, length(A))
    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    # now solve C\(X'y) where C is the cholesky decomposition computed above
    nrhs = size(B[1])[2]
    ldb = max(1,stride(B[1],2))
    Bptrs = CUDA.CUBLAS.unsafe_batch(B)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUSOLVER.$rsname(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, nrhs, Aptrs, lda, Bptrs, ldb, infoarray, length(A))
    CUDA.CUSOLVER.unsafe_free!(Aptrs)
    CUDA.CUSOLVER.unsafe_free!(Bptrs)

    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    B
end end end

for (rfbuffer, rfname, rsname, T) in
    ((:cusolverDnSpotrf_bufferSize, :cusolverDnSpotrf, :cusolverDnSpotrs, :Float32),
    (:cusolverDnDpotrf_bufferSize, :cusolverDnDpotrf, :cusolverDnDpotrs, :Float64))
    @eval begin


function hyperreg(X::CuMatrix{$T},  y::CuMatrix{$T})
  cuuplo::Char = 'U'#CUDA.CUBLAS.cublasfill('U')

  if size(X,1) != length(y)
    throw(DimensionMismatch(""))
  end

  mx,nx = size(X)
  my,ny = size(y)

  if mx != my
      throw(ArgumentError("X Matrix must have same num of rows as y matrix"))
  end

  #@info "point 1"
  A = X' * X
  B = X' * y#CUDA.CUBLAS.gemm('T','N', X, y) #X'y

  size(A,1) == size(A,2) == size(B,1) || throw(ArgumentError("uh oh, something is wrong"))

  #@info "point 2"
  #Here we compute cholesky!(X'X)
  #first setup
  n = size(A,1)
  lda = max(1,stride(A,2))
  devinfo = CuArray{Cint}(undef, 1)
  CUDA.CUSOLVER.@workspace eltyp=$T size=CUDA.CUSOLVER.@argout(
    CUDA.CUSOLVER.$rfbuffer(CUDA.CUSOLVER.dense_handle(),
      cuuplo, n, A, lda, out(Ref{Cint}(0)))
    )[] buffer->begin
      CUDA.CUSOLVER.$rfname(CUDA.CUSOLVER.dense_handle(),
        cuuplo, n, A, lda, buffer, length(buffer), devinfo)
      end
  info = CUDA.CUSOLVER.@allowscalar devinfo[1]
  #CUDA.CUSOLVER.unsafe_free!(devinfo)

  if info != 0
    throw(ArgumentError(string("Invalid value during cholesky decomp at ",-info)))
  end

  #@info "point 3"
  # now solve C\(X'y) where C is the cholesky decomposition computed above
  nrhs = size(B,2)
  ldb = max(1,stride(B,2))

  devinfo = CuArray{Cint}(undef, 1)
  CUDA.CUSOLVER.$rsname(CUDA.CUSOLVER.dense_handle(), cuuplo, n, nrhs, A, lda, B, ldb, devinfo)

  info = CUDA.CUSOLVER.@allowscalar devinfo[1]


  if info != 0
    throw(ArgumentError,string("Invalid value solving triangular matrix at ",-info))
  end

  CUDA.CUSOLVER.unsafe_free!(devinfo)
  CUDA.CUSOLVER.unsafe_free!(A)

  B
end

#handles the case where we are given a vector for y instead of a matrix
hyperreg(X::CuMatrix{$T},  y::CuVector{$T}) = hyperreg(X,
  CUDA.unsafe_wrap(CuArray,CUDA.pointer(y,1),(length(y),1)))

end end


function regprojectedyonx(M::FixedAnnihilator, X::TX, projectedy) where TX
    projectedX = M * X
    return cholesky!(projectedX' * projectedX)\(projectedX' * projectedy)
end

function regprojectedyonx(M::FixedAnnihilator, X::TX, projectedy::CuMatrix) where {TX<:CuMatrix}
    projectedX = M * X
    return hyperreg(projectedX, projectedy)
end

regprojectedyonx(M::FixedAnnihilator, X::TX, projectedy::CuVector) where {TX<:CuMatrix} = (
    regprojectedyonx(M,X, reshape(projectedy,:,1) |> deepcopy))


######this fixes
#when this issue is resolved, can delete the below
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

#the rest of the below is testing code
#https://github.com/FluxML/Zygote.jl/issues/933
function ∇XtXinvXty(X,y)
  tf(X,y) = sum(cholesky(X'*X)\ (X' * y))
  tf(X,y) #works fine
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
