using Revise, CUDA, LinearAlgebra, Random

using CUDA.CUBLAS
using CUDA.CUBLAS: band, bandex

Random.seed!(11)
function mwe()

  local m = 20
  local n = 35
  local k = 13

  elty=Float32

  local alpha = rand(elty)
  local beta = rand(elty)

  local A = triu(rand(elty, m, m))
  local B = rand(elty,m,n)
  local C = zeros(elty,m,n)
  local dA = CuArray(A)
  local dB = CuArray(B)
  local dC = CuArray(C)
  local failed=false


  try
    C = alpha*(A\B)
    dC = copy(dB)
    CUBLAS.xt_trsm!('L','U','N','N',alpha,dA,dC)
    CUDA.synchronize()
    # move to host and compare
    h_C = Array(dC)
    @assert C ≈ h_C
  catch err
    @warn "xt_trsm! gpu failed!! error: $err"
    failed=true
  end

  try
    C  = alpha*(A\B)
    h_C = CUBLAS.xt_trsm('L','U','N','N',alpha,Array(dA),Array(dB))
    CUDA.synchronize()
    @assert C ≈ h_C
  catch err
    @warn "xt_trsm cpu failed!! error: $err"
    failed=true
  end

  return failed
end

for i ∈ 1:10^3
  if mwe()
    @info "Failed on iteration $i"
    break
  end
end
#=
try
    # generate matrices
    syrkx_A = rand(elty, n, k)
    syrkx_B = rand(elty, n, k)
    d_syrkx_A = CuArray(syrkx_A)
    d_syrkx_B = CuArray(syrkx_B)
    d_syrkx_C = CUBLAS.xt_syrkx('U','N',d_syrkx_A,d_syrkx_B)
    final_C = syrkx_A*transpose(syrkx_B)
    # move to host and compare
    h_C = Array(d_syrkx_C)
    @assert triu(final_C) ≈ triu(h_C)
catch err
    @warn "xt_syrkx gpu failed!! error: $err"
end=#
