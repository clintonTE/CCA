using Revise
using LinearAlgebra,Random, BenchmarkTools, Zygote
Random.seed!(1111)
#TODO- adapt to the current problem under considration

#value of the function-used for testing purposes


f_Xa2(X) = ((a) -> (x->(x-a)^2).(X))
D_Xa2(X) = ((a)->2a .- 2X)

#projection matrix that is a function of a parameter
function f_PXay(f_Xa, X, y)
  return (
  function PXay(a)
    Xa = f_Xa(a)
    return Xa*(cholesky(Xa'*Xa)\(Xa'*y))
  end)
end

function D_PXayold(f_Xa, D_Xa, X, y)
  return (
  function DPXay(a)
    Xa = f_Xa(a)
    DXa = D_Xa(a)
    XatXa = Xa' * Xa
    XatXainv = cholesky(XatXa)\I
    return (DXa*XatXainv*(Xa'*y) -
      Xa*XatXainv*(DXa'*Xa+Xa'*DXa)*XatXainv * Xa'*y +
      Xa*XatXainv*(DXa'*y))
  end)
end

function D_PXay(f_Xa, D_Xa, X, y)
  return (
  function DPXay(a, X::TM = f_Xa(a), DaX::TM = D_Xa(a)) where TM

    K= size(X,2)
    cuI = TM(I, K, K) #need an identity matrix
    XtXinv = cholesky!(X'*X)\cuI
    XXtXinv = X*XtXinv
    XtXinvXtv = XtXinv * (X'*y)

    return DaX*XtXinvXtv -
      XXtXinv*(DaX'*X+X'*DaX)*XtXinvXtv +
      XXtXinv*(DaX'*y)
  end)
end



function testD(f,D,args...; testname="", erroronfail::Bool = true)
  testval = D(args...) |> sum

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
  X= rand(N,K)
  y=rand(N) + X*rand(K)
  h = 10^-4
  a = 10.0

  Xa2 = f_Xa2(X)
  DXa2 = D_Xa2(X)
  Pa2 = f_PXay(Xa2, X, y)
  DPa2 = D_PXay(Xa2, DXa2, X, y)

  sPa2(a) = sum(Pa2(a))
  fd = (sPa2(a+h) - sPa2(a))/h
  @info "fd: $fd"

  testD(Pa2, DPa2, a, testname="Pay")

end

@time testderivs(10^4,10)
