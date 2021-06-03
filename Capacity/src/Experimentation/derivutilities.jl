using Revise,Zygote,LinearAlgebra,Random, ChainRules,BenchmarkTools
Random.seed!(1111)

#value of the function-used for testing purposes




#D(sum(A)) = [matrix of ones]
f_sA(A) = sum(A)
D_sA(A::TM, ::Type{T}=eltype(TM)) where {TM<:AbstractMatrix,T} = ones(T, size(A))

#D(At) = (noop- doesn't affect the input)
f_At(A) = A'
D_sAt(A::TM, ::Type{T}=eltype(TM)) where {TM,T} = D_A(A,T)

#D(AtA): let b=ones(K). Then D(AtA)=2A*bb'
f_AtA(A) = A'*A
function D_AtA(A::TM, ::Type{T}=eltype(TM)) where {TM,T}
  N,K = size(A)
  postmult = fill(T(2.0), K,K)
  #DA = D_sA(A)
  return A*postmult
end

#why doesn't this work?
function ChainRulesCore.rrule(::typeof(f_AtA), A::TM, ::Type{T}) where{TM,T}
  K::Int = size(A,2)
  f_AtA(A), (A)->A*fill(T(2.0), K,K)
end

function C_sAtA(A::TM, ::Type{T}=eltype(TM)) where {TM,T}
  N,K = size(A)
  AtA,pAtA  = rrule(*, A',A)
  sAtA,psAtA = rrule(sum, AtA)

  _, _psAtA = psAtA(T(1))
  _,_,_pAtA2 = pAtA(_psAtA)
  #_,_pAt = pAt(extern(_pAtA1))


  return T(2) .* _pAtA2
end


f_Aty(A,y) = A'*y
D_Aty(A::TM, y::AbstractVector{T}, ::Type{T}=eltype(A)) where {TM,T} =(
  reduce(hcat,[y for i in 1:size(A,2)]))



#D(AtAinv) = -AtAinv * D(AtA) * AtAinv
AtAinv(A) = cholesky(A'*A)\I
function D_AtAinv(A::TM, AtAinv::TM=cholesky(A'*A)\I, ::Type{T}=eltype(TM)) where {TM,T}
  N,K = size(A)
  return T(-1) .* A*AtAinv*fill(T(2.0),K,K)*AtAinv
  #return T(-1) .* DA*AtAinv*DAtA*AtAinv
end

#function rrule(::typeof(AtAinv), A::TM where TM, ::Type{T})
#  AtAinv(A),


f_AtAinvAty(A,y) = cholesky(X'*X)\y
#=D((X'X)^-1*X'y) wrt X:
  =D((X'X)^-1)*X'y+(X'X)^-1*D(X'y)
  =#

#=function D_AtAinvAty(A::TM, y::AbstractVector{T}, ::Type{T}=eltype(A)) where {TM,T}

  #build up least squares
  At, pAt = rrule(transpose, A)
  AtA,pAtA = rrule(*, At,A)
  cholAtA, pcholAtA = rrule(cholesky, AtA)
  Aty, pAty

  AtAinv = cholesky(A'*A)\I
  DAtAinv = D_AtAinv(A)
  Aty =A'*y
  DAty = D_Aty(A,y)

  term1 = DAtAinv*Aty
  term2 = AtAinv*DAty
  return  term1 .+ term2
end=#




function testD(f,D,args...; testname="", erroronfail::Bool = true)
  testval = D(args...)

  ∇(args) = gradient((args)->sum(f(args...)), args)
  verval = ∇(args)[1][1]

  local ε::Float64
  try
    ε =sum(abs.(testval .- verval))
  catch err
    msg = "Test $testname: failed probably due to dimensional mismatch: $err
      #args: $args
      testval: $testval
      verval: $verval"
    error(msg)
  end
  if args==verval
    @assert εrms==0
    @info "Test $testname: passed with a difference of 0"
  elseif testval≈verval
    @info "Test $testname: passed with absolute error of $(ε)"
  else
    msg = "Test $testname: failed with absolute error of $(ε)
      #args: $args
      testval: $testval
      verval: $verval"

    erroronfail && error(msg)
    @warn msg
  end


  return nothing
end

function testderivs(N=3,K=2)
  A= rand(N,K)
  y=rand(N)
  testD(f_sA, D_sA, A, testname="noop")
  testD(f_At, D_At, A, testname="transpose")
  testD(f_AtA, D_AtA, A, testname="AtA")
  testD(f_AtA, C_sAtA, A, testname="CAtA")
  #testD(f_AtA, C2_sAtA, A, testname="CAtA")

  testD(f_Aty, D_Aty, A, y, testname="Aty")

  g(A) = gradient((A)->sum(f_AtA(A)), A)
  @btime $g($A)[1][1]
  @btime C_sAtA($A)
  #@btime D_AtA($A)
  #testD(f_AtAinv, D_AtAinv, A, testname="AtAinv")
  #testD(f_AtAinvAty, D_AtAinvAty, A, y, testname="AtAinvAty")

end

testderivs(3,2)
