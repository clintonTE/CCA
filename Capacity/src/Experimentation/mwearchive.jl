

#=
IN_STREAM(p::String) = LZ4DecompressorStream(open(p))
IN_STREAM(F::Function, p::String) = open(F, LZ4DecompressorStream, p)

OUT_STREAM(p::String) = LZ4CompressorStream(open(p, "w"))
OUT_STREAM(F::Function, p::String) = open(F, LZ4CompressorStream, p, "w")

IN_FEATHER_STREAM(p::String) = LZ4DecompressorStream(open(p))
IN_FEATHER_STREAM(F::Function, p::String) = open(F, LZ4DecompressorStream, p)

OUT_FEATHER_STREAM(p::String) = LZ4CompressorStream(open(p, "w"))
OUT_FEATHER_STREAM(F::Function, p::String) = open(F, LZ4CompressorStream, p, "w")


function testout(N; instream = IN_STREAM, outstream = OUT_STREAM)
  tdf::DataFrame = DataFrame(mint = rand([missing, 1,2,3,4], N), mfloat = [rand(N-1); missing],
    #msymbol = rand([missing, :abc, :def, :dsaf], N),
    mstring = rand([missing, "abc", "def", "dsaf"], N))

  p = "tdf.feather.lz4"
  local tdfin::DataFrame

  #outstream("$p.jls.lz4") do s
  #  serialize(s, tdf)
  #end

  outstream(p) do io
    Feather.write(io, tdf)
  end
  @time instream("$p") do io
    tdfin = Feather.materialize(io)
  end
  #@time tdfin = Feather.materialize(p)


  #println("size after: $(size(tdfin))")

end

testout(1000)=#

#=function alluniquemwe(i::Int)
  df = DataFrame(x1=rand(i), x2=rand(i))

  print("allunique: ")
  @time println("$(allunique(eachrow(df)))")

  print("roughly equivalent to allunique: ")
  @time println("$(length(groupby(df, [:x1,:x2]))==size(df,1))")
end

alluniquemwe(10^6)=#

#=using StatsBase
function testq(v = [5,11,11,11,14,14,19,21,21,30])
  N::Int = length(v)
  println("ecdf: $((i->ecdf(v)(i)).(v))")
  println("rank: $(1 .- competerank(v, rev=true)./N)")
end
testq()=#

#=using StatsPlots, Distributions, Random
function testnormal(i)
  n = rand(Normal(), i)
  v1 = (abs).(rand(Normal(), i))
  v2 = (abs).(rand(Normal(), i))

  dif = v1 .- v2

  plot(
    qqplot(Normal, dif)
  )

end=#
#=using Distributions, StatsBase, Statistics

d = Laplace(2,1)
dn = Laplace(-2,1)

v = rand(d,10^7)
vn = rand(dn,10^7)

absv = (abs).(v)
absvn = (abs).(vn)

cdfv(x,v)=sum(v .< x)/length(v)

function tf(y,u, b=1.5)
  out = 1

  if -y < u
    out -= 0.5*exp((-y-u)/b)
  else
    out -= (1-0.5*exp(-(-y-u)/b))
  end

  if y < u
    out -= 1-0.5*exp((y-u)/b)
  else
    out -= 1-(1-0.5*exp(-(y-u)/b))
  end

  return out
end

function tf2(y,u, b=1.5)
  local out

  if y < u
    out = 0.5*exp((y-u)/b) - 0.5*exp((-y-u)/b)
  elseif -y < u
    out = 1 - 0.5*exp((-y-u)/b)-0.5*exp(-(y-u)/b)
  elseif u < -y
    out = 0.5*exp(-(-y-u)/b)-0.5*exp(-(y-u)/b)
  else
    @assert false
  end

  return out
end

function tf3(y,u, b=1.5)
  local out

  if y < abs(u)
    #out = 0.5*exp((y-u)/b) - 0.5*exp((-y-u)/b)
    out  = exp(-abs(u)/b)*sinh(y/b)
  elseif -y < u
    out = 1 - exp(-y/b)*cosh(abs(u)/b)
  else
    @assert false
  end

  return out
end
println()
#=println("\ntf(1,2):  $(tf(1,2))", ", tf(2.5,2): ", tf(2.5,2))
println("tf(1,-2): ", tf(1,-2), ", tf(2.5,-2): ", tf(2.5,-2))

println("tf2(1,2):  $(tf2(1,2))", ", tf2(2.5,2): ", tf2(2.5,2))
println("tf2(1,-2): ", tf2(1,-2), ", tf2(2.5,-2): ", tf2(2.5,-2))

println("*tf3(1,2):  $(tf3(1,2))", ", tf3(2.5,2): ", tf3(2.5,2))
println("*tf3(1,-2): ", tf3(1,-2), ", tf3(2.5,-2): ", tf3(2.5,-2))=#

function phi1(y,μ, b=1.5)
  local out
  (y < 0) && error("support only over positive numbers")
  if y ≤ abs(μ)
    out=1/b*exp(-abs(μ)/b)*cosh(y/b)
  elseif y>abs(μ)
    out=1/b*exp(-y/b)*cosh(abs(μ)/b)
  else
    @assert false
  end

  return out
end

phi2(y,μ, b=1.5) = 1/(2*b)*(exp(-abs(y-abs(μ))/b) + exp(-abs(-y-abs(μ))/b))

function phi3(y,μ, b=1.5)
  d = Laplace(abs(μ),b)

  return pdf(d,y) + pdf(d,-y)
end

function phi4(y,μ, b=1.5)
  dpos = Laplace(abs(μ),b)
  dneg = Laplace(-abs(μ),b)

  return pdf(dpos,y) + pdf(dneg,y)
end

function phi4cdf(y,μ, b=1.5)
  dpos = Laplace(abs(μ),b)
  dneg = Laplace(-abs(μ),b)

  return cdf(dpos,y) + cdf(dneg,y)
end

ty = collect(0.0:.0001:10)
vphi1 = (y->phi1(y,-2)).(ty)
vphi2 = (y->phi2(y,-2)).(ty)
vphi3 = (y->phi3(y,-2)).(ty)
vphi4 = (y->phi4(y,-2)).(ty)

maxdelta12 = maximum((vphi1 .- vphi2).^2)^0.5
maxdelta13 = maximum((vphi1 .- vphi3).^2)^0.5
maxdelta14 = maximum((vphi1 .- vphi4).^2)^0.5


println("Max dif12: $maxdelta12")
println("Max dif13: $maxdelta13")
println("Max dif14: $maxdelta14")
println("cdf: $(phi4cdf(Inf,2))")=#

#=using Revise
using Turing, Distributions, Tracker, Zygote, ReverseDiff

function testturing(n::Int, k::Int, backend)
  Turing.setadbackend(backend)

  turnprogress(true)

  ϵ = rand(n)
  X = rand(n,k)
  β::Vector{Float64} = collect(1:k)

  yₐ::Vector{Float64} = 1.0 .+ X*β .+ ϵ

  kzeros = zeros(k)

  @model testmod(y, ::Type{T}=Vector{Float32}) where {T} = begin
    σ² ~ truncated(Normal(0,k*10),0,Inf)
    b₀ ~ Normal(0,k)

    b = T(undef, k)
    #b ~ MvNormal(kzeros,k*10)
    b ~ filldist(Normal(0.0,k*10),k)

    μ = vec(X * b .+ b₀ )
    y ~ arraydist((Laplace).(μ,σ²))
  end

  m = testmod(yₐ)

  chain = sample(m, HMC(0.01, 100), 1000)
end

@time testturing(1000,5,:forwarddiff)=#




#####Turingdist
#=abstract type TuringMVDistribution <: ContinuousMultivariateDistribution end

#define a fully vectorized multi-variate gamma distribution
struct MvGamma{T<:Real,
  Gammas<:AbstractVector,
  Shape<:AbstractVector,
  Scale<:AbstractVector} <: TuringMVDistribution

  Γ::Gammas
  α::Shape
  θ::Scale
end

function MvGamma(α::AbstractVector{T}, θ::AbstractVector{T}) where {T<:Real}
  length(α) == length(θ) || throw(DimensionMismatch("The dimensions of α and θ are inconsistent."))

  Γ::Vector{Gamma{T}} = Gamma{T}.(α,θ)
  return MvGamma{T, typeof(Γ), typeof(α), typeof(θ)}(Γ, α, θ)
end

Base.eltype(::Type{<:MvLogNormal{T}}) where {T} = T

Base.length(d::MvGamma) = length(d.α)
shape(d::MvGamma) = d.α
scale(d::MvGamma) = d.θ
params(d::MvGamma) = (d.α, d.θ)

#Base.show(io::IO, d::MvGamma) =
#    show_multline(io, d, [(:dim, length(d)), (:α, shape(d)), (:θ, scale(d))])

mean(d::MvGamma) = d.α .* d.θ
var(d::MvGamma) = d.α .* d.θ .^ 2
skewness(d::MvGamma) = 2 ./ sqrt.(d.α)
kurtosis(d::MvGamma) = 6 ./ d.α

mgf(d::MvGamma, t::Real) = (1 .- t .* d.θ).^(-d.α)
cf(d::MvGamma, t::Real) = (1 .- im .* t .* d.θ).^(-d.α)

function _rand!(rng::AbstractRNG, d::MvGamma, x::AbstractVector{T}) where T<:Real
  x .= (γ->rand(γ)).(d.Γ)
end

function _rand!(d::MvGamma, A::DenseMatrix{T}) where T<:Real
  draws::Int = size(A,2)
  mvgammadraw!(x::AbstractVector{T}) = _rand!(d,x)

  mapslices(mvgammadraw!, A, dims=1)
end

function _logpdf!(r::AbstractArray, d::AbstractMvNormal, x::AbstractMatrix)
    sqmahal!(r, d, x)
    c0 = mvnormal_c0(d)
    for i = 1:size(x, 2)
        @inbounds r[i] = c0 - r[i]/2
    end
    r
end

function testmvgamma()
  d=MvGamma(ones(10),ones(10))

  rand(d)
end

testmvgamma()=#

###MWE for flux with MoreParam
#=using Revise, CUDAapi, Flux, BenchmarkTools
import CuArrays
struct MoreParam{Tc<:AbstractVector, Td<:AbstractVector}
  c::Tc
  d::Td
end

struct MWEStruct{Ta<:AbstractVector, Tb<:AbstractVector}
  a::Ta
  b::Tb
  moreparam::Tmoreparam
end

MWEStruct(N) = MWEStruct(rand(N), rand(N))

(m::MWEStruct)(x) = x

function mwe(N)
  Flux.@functor MWEStruct

  cpustruct = MWEStruct(N)
  println("type of cpu struct: ", typeof(cpustruct))
  gpustruct = cpustruct |> gpu
  println("type of gpu struct: ", typeof(gpustruct))
end

mwe(10)=#

###MWE for flux regarding first run
#=using Revise, Flux
import CuArrays

struct MWEStruct{Ta<:AbstractVector, Tb<:AbstractVector}
  a::Ta
  b::Tb
end

MWEStruct(N) = MWEStruct(rand(N), rand(N))

(m::MWEStruct)(x) = x

function mwe1(N)
  Flux.@functor MWEStruct

  cpustruct = MWEStruct(N)
  println("type of cpu struct: ", typeof(cpustruct))
  gpustruct = cpustruct |> gpu
  println("type of gpu struct: ", typeof(gpustruct))
end

function mwe2(N)
  Flux.@functor MWEStruct

  cpustruct = MWEStruct(N)
  println("type of cpu struct: ", typeof(cpustruct))
  gpustruct = cpustruct |> gpu
  println("type of gpu struct: ", typeof(gpustruct))
end

println("First run:")
mwe(10)
println("\nSecond run:")
mwe(10)=#

####MWE for CuArrays sum
#=using Flux,BenchmarkTools,LinearAlgebra
import CuArrays, CUDAnative, CUDAapi

CuArrays.allowscalar(false)
function mwesum(N)
  cpuv= rand(Float32, N)
  print("\nStandard cpu sum: ")
  @btime sum($cpuv)

  cuv=CuArrays.cu(cpuv)
  print("Standard cuda sum: ")
  @btime sum($cuv)

  print("Summing using allocation and dot: ")
  onesum(v) = dot(cuv, CuArrays.ones(N))
  @btime $onesum($cuv)

  cuvones = CuArrays.ones(N)
  print("Summing just using dot and a pre-allocated vector of 1s: ")
  onesum2(v) = dot(cuv, cuvones)
  @btime $onesum2($cuv)
end

sleep(0.5)
mwesum(10^8)=#



###MWE for trainable
#=using Revise, Flux
sleep(0.5)
struct Affine
  W
  b
  c
end

Affine(in::Integer, out::Integer) =
  Affine(randn(out, in), randn(out), rand())
(m::Affine)(x) = m.W * x .+ m.b.^m.c


Flux.@functor Affine
Flux.trainable(a::Affine) = (a.W,a.b,)

function testtrainable()
  W= rand(3,3)
  b=rand(3)
  a = Affine(W, b, rand())
  p=Flux.params(a)


  @info "num params: $(length(p))"
  x=rand(3)
  y=rand(3)
  loss(x,y) = sum((a(x) .- y).^2)

  gs = gradient(()->loss(x,y), p)

  println(gs[a.W][:, 1:2])
end

testtrainable()=#



##################Testing for cuda performance
#result: gpu use with views is fastest
#=using Revise, CUDAapi, Flux, BenchmarkTools
import CuArrays
CuArrays.allowscalar(false)
function testagg(v1, m1, v2, m2)
  uk = abs.(m1 .* v1' .- m2 .* v2')
  return sum(uk, dims=1)
end

function testcpustd(N, Nlong, K)
  v1 = rand(Float32, N)
  v2 = rand(Float32, N)
  ts = rand(1:N, Nlong)

  v1long = rand(Float32, Nlong)
  v2long = rand(Float32, Nlong)

  m1 = rand(Float32, K, N)
  m2 = rand(Float32, K, N)

  m1long = rand(Float32, K, Nlong)
  m2long = rand(Float32, K, Nlong)

  function stdagg()
    v1long .= v1[ts]
    v2long .= v2[ts]
    m1long .= m1[:, ts]
    m2long .= m2[:, ts]

    #testagg($v1long, $m1long, $v2long, $m2long)
    testagg(v1long, m1long, v2long, m2long)
  end

  print("testing cpustd: ")
  #stdagg()
  @btime $stdagg()
end

function testgpustd(N, Nlong, K)
  v1 = rand(N) |> gpu
  v2 = rand(N) |> gpu
  ts = rand(1:N, Nlong)

  v1long = rand(Nlong) |> gpu
  v2long = rand(Nlong) |> gpu

  m1 = rand(K, N) |> gpu
  m2 = rand(K, N) |> gpu

  m1long = rand(K, Nlong) |> gpu
  m2long = rand(K, Nlong) |> gpu

  function stdagg()
    v1long .= v1[ts]
    v2long .= v2[ts]
    m1long .= m1[:, ts]
    m2long .= m2[:, ts]

    #@btime $testagg($v1long, $m1long, $v2long, $m2long)
    testagg(v1long, m1long, v2long, m2long)
  end
  #InteractiveUtils.@code_warntype stdagg()

  print("testing gpustd: ")
  #stdagg()
  @btime $stdagg()
end

function testcpuview(N, Nlong, K)
  v1 = rand(Float32, N)
  v2 = rand(Float32, N)
  ts = rand(1:N, Nlong)

  v1long = view(v1, ts)
  v2long = view(v2, ts)

  m1 = rand(Float32, K, N)
  m2 = rand(Float32, K, N)

  m1long = view(m1, :, ts)
  m2long = view(m2, :, ts)

  function stdagg()

    #@btime $testagg($v1long, $m1long, $v2long, $m2long)
    testagg(v1long, m1long, v2long, m2long)
  end

  print("testing cpuview: ")
  #stdagg()
  @btime $stdagg()
end

function testgpuview(N, Nlong, K)
  v1 = rand(Float32, N) |> gpu
  v2 = rand(Float32, N) |> gpu
  ts = rand(1:N, Nlong)

  v1long = view(v1, ts)
  v2long = view(v2, ts)

  m1 = rand(Float32, K, N) |> gpu
  m2 = rand(Float32, K, N) |> gpu

  m1long = view(m1, :, ts)
  m2long = view(m2, :, ts)

  function stdagg()

    #@btime $testagg($v1long, $m1long, $v2long, $m2long)
    testagg(v1long, m1long, v2long, m2long)
  end

  print("testing gpuview: ")
  #stdagg()
  @btime $stdagg()
end

function testgpuhybrid(N, Nlong, K)
  v1 = rand(Float32, N) |> gpu
  v2 = rand(Float32, N) |> gpu
  ts = rand(1:N, Nlong)

  sv1long = view(v1, ts)
  sv2long = view(v2, ts)

  m1 = rand(Float32, K, N) |> gpu
  m2 = rand(Float32, K, N) |> gpu

  sm1long = view(m1, :, ts)
  sm2long = view(m2, :, ts)

  v1long = rand(Nlong) |> gpu
  v2long = rand(Nlong) |> gpu
  m1long = rand(K, Nlong) |> gpu
  m2long = rand(K, Nlong) |> gpu

  function stdagg()

    #@btime $testagg($v1long, $m1long, $v2long, $m2long)
    testagg(CuArray(sv1long), CuArray(sm1long), CuArray(sv2long), CuArray(sm2long))
  end

  print("testing gpuhybrid: ")
  #stdagg()
  @btime $stdagg()
end

N=3000
Nlong = 3000*N
K = 10
sleep(0.4)
println("")
#testcpustd(N, Nlong, K)
#testgpustd(N, Nlong, K)

#testcpuview(N, Nlong, K)
#testgpuview(N, Nlong, K)

#testgpuhybrid(N, Nlong, K)=#


###This is  vecloggammapdf
#=
using Revise, CUDAapi, Flux, BenchmarkTools, CUDAnative, Distributions, Zygote, Flux, DiffRules
using SpecialFunctions
import CuArrays
CuArrays.allowscalar(false)

#GPU vectorized log gamma pdf
function vecgammalogpdf(k::TV, θ::TV, x::TV)::TV where TV<:CuArrays.CuVector
  Γs_part = cudagammalogpdf_part.(k,θ,x)

  return Γs_part .- cudaveclgamma(k)
end

#see https://math.stackexchange.com/questions/1441753/approximating-the-digamma-function
#and https://en.wikipedia.org/wiki/Digamma_function
function cudaapproxdigamma(x::Real)
  adj = 0.0f0
  ψ = x
  #the polynomial is only accurate for large x
  for i ∈ 1:3
    adj -= 1f0/ψ
    ψ += 1f0
  end

  ψ = (CUDAnative.log(ψ) - 1f0 / 2f0 / ψ - 1f0 / 12f0 / (ψ * ψ) +
    1f0/120f0 * CUDAnative.pow(ψ, -4) - 1f0/252f0 * CUDAnative.pow(ψ, -6) +
    1f0/240f0 * CUDAnative.pow(ψ, -8) - 5f0/660f0 * CUDAnative.pow(ψ, -10) +
    691f0/32760f0 * CUDAnative.pow(ψ, -12) - 1f0/12f0 * CUDAnative.pow(ψ, -14)) + adj

  return ψ
end

#does everything except compute the log gamma function
#adapted from StatsFuns
#https://github.com/JuliaStats/StatsFuns.jl/blob/master/src/distrs/gamma.jl
cudagammalogpdf_part(k,θ,x) = - k * CUDAnative.log(θ) + (k - 1f0) * CUDAnative.log(x) - x / θ

#vectorized log gamma function
cudaveclgamma(k) = (CUDAnative.lgamma).(k)

#define the adjoint
Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
Zygote.refresh()



#y_cudaapproxdigamma(y::Real, k::Real) = y * cudaapproxdigamma(k)
#y_cudavecdigamma(y, k) = (y_cudaapproxdigamma).(y, k)

#WARNING: below doesn't work
#=cudagammalogpdf(k,θ,x) = (-CUDAnative.lgamma(k) - k * CUDAnative.log(θ) +
  (k - 1f0) * CUDAnative.log(x) - x / θ)
Zygote.@adjoint CUDAnative.lgamma(k) = CUDAnative.lgamma(k), y->y*cudaapproxdigamma(k)
Zygote.refresh()

Zygote.@adjoint function vecgammalogpdf(k::TV, θ::TV, x::TV) where TV<:CuArrays.CuVector
  y, back = Zygote.broadcast_forward(CuArrays.cufunc(cudagammalogpdf), k, θ, x)
  return y, Δy -> (nothing, nothing, back(Δy)...)
end
Zygote.refresh()=#

#DiffRules.@define_diffrule CUDAnative.lgamma(k) = :(cudaapproxdigamma($k))
#eval(Zygote.ForwardDiff.unary_dual_definition(:CUDAnative, :lgamma))
#Zygote.refresh()


#Zygote.@adjoint cudaveclgamma(k) = Zygote.broadcasted(
#  Zygote.CuArrayStyle, CuArrays.cufunc(cudaapproxdigamma), k)

#Zygote.@adjoint CUDAnative.lgamma(k) = CUDAnative.lgamma(k), y->(y * cudaapproxdigamma(k), )
#Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
#Zygote.@adjoint function vecgammalogpdf(k::TV, θ::TV, x::TV)::TV where TV<:CuArrays.CuVector
#  y, back = Zygote.pullback(cudagammalogpdf_part, k,θ,x)
#end
#Zygote.refresh()


#= NOTE: Freezing this block as it works
cudagammalogpdf_part(k,θ,x) = - k * CUDAnative.log(θ) + (k - 1f0) * CUDAnative.log(x) - x / θ
cudaveclgamma(k) = (CUDAnative.lgamma).(k)
cudavecdigamma(k) = (cudaapproxdigamma).(k)

Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
Zygote.refresh()
=#


#lightly optimized cpu version for testing
function vecgammalogpdf(k::TV, θ::TV, x::TV)::TV where TV<:AbstractVector

  Γs::TV = -(SpecialFunctions.loggamma).(k) .- k .* log.(θ) .+ (k .- 1f0) .* log.(x) .- x ./ θ

  return Γs
end

slowvecgammalogpdf(k, θ, x) = ((kᵢ, θᵢ, xᵢ)->logpdf(Gamma(kᵢ,θᵢ), xᵢ)).(k, θ, x)

#test function for vecgammalogpdf
function testvecgammalogpdf(N; tol = 10^-6)
  τ = abs.(randn(Float32, N)) .+ .001f0
  absμ = abs.(randn(Float32, N)) .+ .001f0

  θ = 1.0f0 ./ τ
  k = absμ .* τ
  v = rand(Float32, N)

  cuθ = θ |> gpu
  cuk = k |> gpu
  cuv = v |> gpu

  print("cpu time:")
  @btime sum(slowvecgammalogpdf($k,$θ,$v))
  s = sum(slowvecgammalogpdf(k,θ,v))/N

  print("cpu time:")
  @btime sum(vecgammalogpdf($k,$θ,$v))
  sopt = sum(vecgammalogpdf(k,θ,v))/N

  print("gpu time:")
  @btime sum(vecgammalogpdf($cuk,$cuθ,$cuv))
  cus = sum(vecgammalogpdf(cuk,cuθ,cuv))/N


  s=isfinite(s) ? s : 0.0f0
  cus = isfinite(cus) ? cus : 0.0f0

  Δ = abs(cus-s)+abs(cus-sopt)
  (Δ < tol) || error("Δ $Δ greater then tol. s: $s, cus: $cus, sopt: $sopt")

  @info "Δ: $Δ"
end

#test it
function testgrad(N::Int)
  τ = abs.(randn(Float32, N)) .+ .001f0
  absμ = abs.(randn(Float32, N)) .+ .001f0

  θ = 1.0f0 ./ τ
  k = absμ .* τ
  v = rand(Float32, N)

  cuθ = θ |> gpu
  cuk = k |> gpu
  cuv = v |> gpu

  @info "testing function values (N=$N) cpu: $(sum(vecgammalogpdf(k, θ, v))) " *
    "gpu: $(sum(vecgammalogpdf(cuk, cuθ, cuv)))"
  @info "Additional check of value: $(sum(slowvecgammalogpdf(k, θ, v)))"

  @info "cpu version gradient N=$N"
  if N ≤ 10
    gs = gradient((k,θ)->sum(vecgammalogpdf(k, θ, v)),k,θ)
    @info "gs: $gs"
  else
    gs = @btime gradient((k,θ)->sum(vecgammalogpdf(k, θ, $v)),$k,$θ)
  end

  @info "gpu version gradient N=$N"
  if N ≤ 10
    cugs = gradient((cuk,cuθ)->sum(vecgammalogpdf(cuk, cuθ, cuv)),cuk,cuθ)
    @info "cugs: $cugs"
  else
    cugs = @btime gradient((cuk,cuθ)->sum(vecgammalogpdf(cuk, cuθ, $cuv)),$cuk,$cuθ)
  end
end

testgrad(5)
testgrad(10^6)=#
###log normal dpf
#sortof a cuda kernal
#=using Revise, BenchmarkTools, Distributions, Zygote, Flux, CuArrays, CUDAnative
function cudanormlogprop(μᵢ::Real, logτᵢ::Real, xᵢ::Real)
  μMxᵢ = μᵢ-xᵢ
  return -1f0 * μMxᵢ*μMxᵢ*CUDAnative.exp(logτᵢ) + logτᵢ
end

function vecnormlogprop(μ::TV, logτ::TV, x::TV)::TV where TV<:AbstractVector
  uMx = μ .- x
  return -1f0 .* uMx .* uMx .* exp.(logτ) .+ logτ
end


#Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
#Zygote.refresh()
function vecnormlogprop(μ::TV, logτ::TV, x::TV)::TV where TV<:CuArrays.CuVector
  return cudanormlogprop.(μ,logτ,x)
end
cudanormlogpdffromprop(x::Real) = 0.5f0*x - 0.5f0*log(2f0*3.1415927f0)
slowvecnormlogpdf(μ, σ, x) = ((μᵢ, σᵢ, xᵢ)->logpdf(Normal(μᵢ,σᵢ), xᵢ)).(μ, σ, x)

function vecnormlogpdf(μ::TV, logτ::TV, x::TV)::TV where TV<:CuArrays.CuVector
  return cudanormlogpdffromprop.(vecnormlogprop(μ, logτ, x))
end

function vecnormlogpdf(μ::TV, logτ::TV, x::TV)::TV where TV<:AbstractVector
  return 0.5f0 .* vecnormlogprop(μ, logτ, x) .- 0.5f0 .* log.(2f0 .* 3.1415927f0)
end

function testnormgrad(N::Int)
  τ = log.(abs.(randn(Float32, N)) .+ .001f0)
  μ = abs.(randn(Float32, N)) .+ .001f0

  σ = (1.0f0 ./ exp.(τ)).^0.5f0
  v = rand(Float32, N)

  cuμ = μ |> gpu
  cuτ = τ |> gpu
  cuσ = σ |> gpu
  cuv = v |> gpu

  @info "testing prop values (N=$N) cpu: $(sum(vecnormlogprop(μ, τ, v))) " *
    "gpu: $(sum(vecnormlogprop(cuμ, cuτ, cuv)))"
  @info "corrected cpu values: $(sum(vecnormlogpdf(μ, τ, v)))"
  @info "corrected gpu values: $(sum(vecnormlogpdf(cuμ, cuτ, cuv)))"
  @info "Additional check of value: $(sum(slowvecnormlogpdf(μ, σ, v)))"

  @info "cpu version gradient N=$N"
  if N ≤ 10
    gs = gradient((μ, τ)->sum(vecnormlogprop(μ, τ, v)),μ, τ)
    @info "gs: $gs"
  else
    gs = @btime gradient((μ,τ)->sum(vecnormlogprop(μ, τ, $v)),$μ,$τ)
  end

  @info "gpu version gradient N=$N"
  if N ≤ 10
    cugs = gradient((cuμ,cuτ)->sum(vecnormlogprop(cuμ, cuτ, cuv)),cuμ,cuτ)
    @info "cugs: $cugs"
  else
    cugs = @btime gradient((cuμ,cuτ)->sum(vecnormlogprop(cuμ, cuτ, $cuv)),$cuμ,$cuτ)
  end
end

testnormgrad(3)
testnormgrad(10^5)=#

###custom adjoint for standard vp
#=using Revise, Zygote, Flux, BenchmarkTools, CuArrays

struct mwestructidx{TA,TiA}
  A::TA
  iA::TiA
end

function mwestructidx(K::Int)
  A = reshape(collect(1:K^2), (K,K))
  iA = to_indices(A, ((i->(i-1)%K+1).(1:K^2), :))
  return mwestructidx(A, iA)
end
expand(m::mwestructidx) = view(m.A,m.iA...)
#expand(m::mwestructidx) = SubArray(m.A,m.iA)

Zygote.@adjoint SubArray(A, iA) = SubArray(A,iA), Zygote.∇getindex(A, iA)
Zygote.refresh()

function (m::mwestructidx)(x)
  sA = expand(m)
  sum(sA * x)
end

Flux.@functor mwestructidx
Flux.trainable(a::mwestructidx) = (a.A,)

function zygotemweidx(N=5)

  #=s=mwestructidx(N)
  x = collect(1:N)
  @info "test f: $(s(x))"

  #gs = gradient((A)->sum(A[parentindices(s.sA)...]*x),s.A)
  p=Flux.params(s)
  gs = @btime gradient(()->$s($x),$p)
  if N < 10
    @info "gradient: $(gs[s.A])"
  end=#
  K=N
  A = reshape(collect(1:K^2), (K,K))
  iA = to_indices(A, ((i->(i-1)%K+1).(1:K^2), :))
  x = collect(1:N)

  @info "New view gradient: "
  fview(Ain) = sum(view(Ain, iA...)*x)
  gs = @btime gradient($fview,$A)

end=#



#####
#=
function cudasoftabs(μᵢ::Real)
  ϵ = 1f-5
  return CUDAnative.sqrt(μᵢ*μᵢ-ϵ*ϵ) + ϵ
end

#standard sub-array format
using Revise, Zygote, Flux, BenchmarkTools, CuArrays, CUDAnative
CuArrays.allowscalar(false)

#@adjoint hadamardprod(A, x::CuArray)

struct mwestruct{TA,TiA,TsA}
  A::TA
  B::TA
  i::TiA
  sA::TsA
  sB::TsA
end

zygote_pls_notice_me(A, sA) = sA
Zygote.@adjoint zygote_pls_notice_me(A, sA) = sA, Zygote.∇getindex(A, parentindices(sA))
Zygote.refresh()

function mwestruct(K::Int, A::AbstractMatrix, B::AbstractMatrix)

  i = (i->(i-1)%K+1).(1:K^2)
  sA = view(A, :, (i->(i-1)%K+1).(1:K^2))
  sB = view(B, :, (i->(i-1)%K+1).(1:K^2))
  return mwestruct(A, B, i, sA, sB)
end

expand(m::mwestruct, ::Any) = zygote_pls_notice_me(m.A, m.sA), zygote_pls_notice_me(m.B, m.sB)
expand(m::mwestruct, ::CuArray) =  CuMatrix(m.A[:,m.i]), CuMatrix(m.B[:,m.i])


function (m::mwestruct)(x)
  sA, sB = expand(m,x)

  sum(abs.((sA .-  sA.*sB .- sB)  ./ (x .+ 1f0)))
end

function (m::mwestruct)(x::CuArray)
  sA, sB = expand(m,x)

  sum(cudasoftabs.(sA .-  sA.*sB .- sB)  ./ (x .+ 1f0))
end


Flux.@functor mwestruct
Flux.trainable(a::mwestruct) = (a.A,a.B)


function zygotemwe(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  B = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  A .-= mean(A)
  B .-= mean(B)
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)

  s=mwestruct(N, A, B)
  @info "cpu: test f: $(s(x))"


  #gs = gradient(()->s(x),Flux.params(s))
  if N < 10
    gs = gradient(()->s(x),Flux.params(s))
    @info "cpu gradient: $(values(gs.grads))" #WARNING- uncomment soon
  else
    gs = @btime gradient(()->$s($x),Flux.params($s))
  end

  cuA = A |> gpu
  cuB = B |> gpu
  cux = x |> gpu

  cus=mwestruct(N, A, B)
  @info "gpu: test f: $(cus(cux))"
  if N < 10
    cugs = gradient(()->cus(cux),Flux.params(cus))
    #sumgrad = sum(vcat((v->vec(v)).(values(cugs.grads))))
    #vcat((v->vec(v)).(values(cugs.grads)))
    @info "gpu gradient: $(values(cugs.grads))"
    #@info "gpu gradient: $(vcat((v->vec(v)).(values(cugs.grads))))"
  else
    cugs = @btime gradient(()->$cus($cux),Flux.params($cus))
  end

end


function zygotemwecheck(N=5)

  A = reshape(collect(1:N^2), (N,N))
  #sA = view(A, (i->(i-1)%N+1).(1:N^2), :)
  xin = collect(1:N)
  As = view(A, :, (i->(i-1)%N+1).(1:N^2))
  #display(As .* xin)
  f(x)=sum(view(A, :, (i->(i-1)%N+1).(1:N^2)) .* x)
  fA(x)=sum(A .* x)
  #@info "dims A: $(size(A)), dims As: $(size(As))"

  gs = gradient(()->f(xin), Flux.params(A))
  if N < 10
    @info "check gradient: $(gs[A])"
  end
end




N = 4
#@time zygotemweidx(N)
zygotemwe(N)
#zygotemwecheck(N)=#


#TODO: get one of these work arounds or something like them into capacity
#=using Revise, Zygote, BenchmarkTools, CuArrays, CUDAnative, Flux
function badnewsbears(N)

  A = reshape(collect(1:N^2), (N,N))
  Acu = A |> gpu

  x = vcat((_->collect(1:N^2)').(1:N)...)
  xp = hcat((_->collect(1:N^2)).(1:N)...)
  xcu = x |>gpu
  xpcu = xp |>gpu

  colindices = (i->(i-1)%N+1).(1:N^2)

  function f(A,x)
    sA = view(A, :, (i->(i-1)%N+1).(1:N^2))
    sum(sA .* x)
  end

  function f(A,x, colindices)
    As = A[:, colindices]
    sum(As .* x)
  end

  #function fA2(A,x)
  #  Afull = hcat((i->A[:,i]).((i->(i-1)%N+1).(1:N^2))...)
  #  sum(Afull .* x)
  #end

  @info "cpu value: $(f(A,x)); gpu value: $(f(Acu,xcu))"

  if N≤5
    gs = gradient((A,x)->f(A,x), A,x)
    @info "check cpu gradient: $(gs)"
  else
    gs = @btime gradient((A,x)->$f(A,x), $A,$x)
  end

  if N≤5
    gscu = gradient((Acu, xcu)->f(Acu, xcu, colindices), A, x)
    @info "check gpu gradient: $(gscu)"
  else
    gs = @btime gradient((Acu, xcu)->$f(Acu, xcu, $colindices), $Acu, $xcu)
  end
end

sleep(0.5)
badnewsbears(3)=#
#badnewsbears(100)
#=
using Revise, CuArrays, Flux
function cusubmwe(N)
  v = ones(Float32, N)
  vgpu = v|> gpu

  indices = (i->(i-1)%N + 1).(1:N^2)
  subv = view(v, indices)
  subvgpu = view(vgpu, indices)

  subv .+= 1
  subvgpu .+= 1

  @info "indices: $indices"
  @info "cpu: $(v)"
  @info "gpu: $(vgpu)"

end

N=3
cusubmwe(N)=#



####Experimentation with Hessians
#=using Revise, Zygote, Flux
struct AB{Ta,Tb}
  a::Ta
  b::Tb
end

AB(K) = AB(rand(),rand(K))

(m::AB)(X) = m.a .+ X*m.b

Flux.@functor AB

function testab(N::Int, K::Int)
  atrue=1.0
  btrue=repeat([3.0], K)

  #gen x, y, and noise
  X = rand(N,K)
  ϵ = K*rand(N) #more noise for more vars
  y = atrue .+ X*btrue .+ ϵ

  ab = AB(K)
  Π = Flux.params(ab)

  loss(X,y) = sum((y .- X*ab.b .- ab.a).^2)

  function lossab(v)
    m = AB(v[1],v[2:(K+1)])
    return sum((y-m(X)).^2)
  end


  gs = gradient(()->loss(X,y), Π)

  θ,re = Flux.destructure(Π)
  h = Zygote.hessian((v)->lossab(v), [ab.a; ab.b])
  display(h)
  #h = Zygote.hessian((a,b)->lossab(a,b), ab.a, ab.b)

  #println(h)
end

N=10
K=5
testab(N,K)=#

###WTF MWE
#=using Revise, Zygote, Flux, BenchmarkTools, CuArrays, CUDAnative
CuArrays.allowscalar(false)

#@adjoint hadamardprod(A, x::CuArray)

struct mwestruct{TA,TiA,TsA}
  A::TA
  i::TiA
  sA::TsA
end

zygote_pls_notice_me(A, sA) = sA
Zygote.@adjoint zygote_pls_notice_me(A, sA) = sA, Zygote.∇getindex(A, parentindices(sA))
Zygote.refresh()

function mwestruct(K::Int, A::AbstractMatrix)

  i = (i->(i-1)%K+1).(1:K^2)
  sA = view(A, :, (i->(i-1)%K+1).(1:K^2))
  return mwestruct(A, i, sA)
end

expand(m::mwestruct, ::Any) = zygote_pls_notice_me(m.A, m.sA)
expand(m::mwestruct, ::CuArray) =  CuMatrix(m.A[:,m.i])


function (m::mwestruct)(x)
  sA = expand(m,x)

  #sum(sA .* (sA .* sA ) .+  sB  ./ (x .+ 1f0))
  sum(sA .* x)
end

#=function (m::mwestruct)(x::CuArray)

end=#

Flux.@functor mwestruct
Flux.trainable(a::mwestruct) = (a.A,)


function zygotemwe(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)

  s=mwestruct(N, A)
  @info "cpu: test f: $(s(x))"

  #gs = gradient(()->s(x),Flux.params(s))
  if N < 10
    gs = gradient(()->s(x),Flux.params(s))
    @info "cpu gradient: $(values(gs.grads))" #WARNING- uncomment soon
  else
    gs = @btime gradient(()->$s($x),Flux.params($s))
  end

  cux = x |> gpu

  cus=mwestruct(N, A)
  p = Flux.params(cus)
  @info "gpu: test f: $(cus(cux))"
  if N < 10
    cugs = gradient(()->cus(cux),p)
    #sumgrad = sum(vcat((v->vec(v)).(values(cugs.grads))))
    #vcat((v->vec(v)).(values(cugs.grads)))
    @info "gpu gradient: $(values(cugs.grads))"
    #@info "gpu gradient: $(vcat((v->vec(v)).(values(cugs.grads))))"
  else
    cugs = @btime gradient(()->$cus($cux),Flux.params($cus))
  end

end

#=function zygotemwedirect(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)
  i = (i->(i-1)%N+1).(1:N^2)

  function asum(m)
    sum(m .* A[:,i])
  end

  @info "valuecpu: $(asum(x))"
  gs = gradient(()->asum(x),Flux.params(A))
  @info "cpu gradient: $(values(gs.grads))"


  cuA = A |> gpu
  cux = x |> gpu
  function asum(m::CuArray)
    sum(m .* CuMatrix(A[:,i]))
  end

  @info "valuegpu: $(asum(cux))"
  cugs = gradient(()->asum(cux),Flux.params(A))
  #sumgrad = sum(vcat((v->vec(v)).(values(cugs.grads))))
  #vcat((v->vec(v)).(values(cugs.grads)))
  @info "gpu gradient: $(values(cugs.grads))"


end

zygotemwedirect(3)=#
zygotemwe(3)
=#


###standard mwe2
#standard sub-array format

#=using Revise, Zygote, Flux, BenchmarkTools, CUDA, LinearAlgebra
CUDA.allowscalar(false)

###below utility could be used in Capacity
printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end


#this utility re-interprets the object(!)
#the object MUST remain referenced, or bad things will happen
#tries to do some basic safety checks, but its still possible to screw things up
function less_unsafe_wrap(::Type{Tdest};
  source::Vector = error("src vector is required"),
  start::Int = error("start is required"),
  dims::Tuple)::Tdest where Tdest<:AbstractArray

  Ndat = prod(dims)

  T = eltype(Tdest)
  (eltype(source) == T) || error("eltype of src must be of same type as T")

  #this has the effect of making sure we are getting space from somewhere that exists
  protowrap::SubArray = view(source, start:(start + Ndat -1))
  @assert length(protowrap) == Ndat

  #a final check
  p = pointer(source, start)
  wrap::Tdest = unsafe_wrap(Tdest, p, dims)
  @assert sum(vec(protowrap) .== vec(wrap)) == Ndat

  return wrap
end

struct mwestruct{TΠ, TA,TiA,TsA}
  Π::TΠ
  A::TA
  B::TA
  i::TiA
  sA::TsA
  sB::TsA
end

zygote_pls_notice_me(A, sA) = sA
Zygote.@adjoint zygote_pls_notice_me(A, sA) = sA, Zygote.∇getindex(A, parentindices(sA))
Zygote.refresh()

function mwestruct(K::Int, A::TA, B::TA) where TA<:AbstractMatrix

  i = (i->(i-1)%K+1).(1:K^2)
  Π = vcat(vec(A), vec(B))

  linkedA = less_unsafe_wrap(TA, source = Π, start = 1, dims=size(A))
  linkedB = less_unsafe_wrap(TA, source = Π, start = length(vec(A))+1, dims=size(B))
  sA = view(A, :, (i->(i-1)%K+1).(1:K^2))
  sB = view(B, :, (i->(i-1)%K+1).(1:K^2))

  @assert linkedA == A
  @assert linkedB == B
  return mwestruct(Π, linkedA, linkedB, i, sA, sB)
end

#expand(m::mwestruct, ::Matrix{<:AbstractFloat}) = zygote_pls_notice_me(m.A, m.sA), zygote_pls_notice_me(m.B, m.sB)
expand(m::mwestruct, ::Matrix{<:Real}) = Matrix(m.A[:,m.i]), Matrix(m.B[:,m.i])
#expand(m::mwestruct, ::Matrix{<:Real}) = CuMatrix(m.A[:,m.i]), CuMatrix(m.B[:,m.i])
expand(m::mwestruct, ::CuArray) =  CuMatrix(m.A[:,m.i]), CuMatrix(m.B[:,m.i])


function (m::mwestruct)(x)
  #error("got here")
  sA, sB = expand(m,x)
  sA.*sB.*x .+ sA .- 0.5f0*sB
end


Flux.@functor mwestruct
Flux.trainable(a::mwestruct) = (a.A,a.B)


mwevecgrad(a, gs) = reduce(vcat, [vec(gs[a.A]),vec(gs[a.B])])

function zygotemwe(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  B = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  A .-= mean(A)
  B .-= mean(B)
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)

  s=mwestruct(N, A, B)
  @info "cpu: test f: $(sum(s(x)))"

  #gs = gradient(()->s(x),Flux.params(s))
  if N < 10
    gs = gradient(()->sum(s(x)),Flux.params(s))
    vecgrad = mwevecgrad(s,gs)
    @info "cpu gradient: $(vecgrad)" #WARNING- uncomment soon
  else
    gs = @btime gradient(()->sum($s($x)),Flux.params($s))
  end

  cuA = CuMatrix(A)
  cuB = CuMatrix(B)
  cux = CuMatrix(x)

  cus=mwestruct(N, A, B)
  @info "gpu: test f: $(sum(cus(cux)))"
  if N < 10
    cugs = gradient(()->sum(cus(cux)),Flux.params(cus))
    #sumgrad = sum(vcat((v->vec(v)).(values(cugs.grads))))
    #vcat((v->vec(v)).(values(cugs.grads)))
    cuvecgrad = mwevecgrad(cus,cugs)
    @info "gpu gradient: $(cuvecgrad)"
    #@info "gpu gradient: $(vcat((v->vec(v)).(values(cugs.grads))))"
  else
    cugs = @btime gradient(()->sum($cus($cux)),Flux.params($cus))
  end

end

zygotemwe(2)=#


###for test jacob- uses the standard code too
#=function testjacob(N, ::Type{T}) where T<:Real

  A = Matrix{T}(reshape(collect(1:N^2), (N,N)))
  B = Matrix{T}(reshape(collect(1:N^2), (N,N)))
  A .-= mean(A)
  B .-= mean(B)
  x = hcat((_->Vector{T}(collect(1:N))).(1:N^2)...)
  xreal = Matrix{Real}(x)
  cux = x |> gpu

  s=mwestruct(N, A, B)

  Areal=Matrix{Real}(A)
  Breal=Matrix{Real}(B)
  sreal = mwestruct(N, Areal, Breal)
  @info "zygote test f: $(s(x)); forwardtestf  $(sreal(x))"

  "Inputs:"
  @info "A/B"
  printmln(s.A)
  @info "sA/sB"
  printmln(s.sA)
  @info "x"
  printmln(x)

  #NOTE: The function computes as follows
  # sA = [A A]; sB = [B B]
  # out = sA.*sB.*x .+ sA .- 0.5sB
  #     = [A A] .* [B B] .* x .+ [A A] .- 0.5 .* [B B]
  @info "Test Jacobian with Symbolic answer"
  if s(x) == s.sA.*s.sB.*x .+ s.sA .- 0.5s.sB
    DA = [B B] .* x .+ 1
    DB = [A A] .* x .- 0.5
    ∇::Vector{Vector{T}} = Vector{Vector{T}}(undef, 8)
    ∇[1] = vec([1 0 1 0; 0 0 0 0] .* DA)
    ∇[2] = vec([0 0 0 0; 1 0 1 0] .* DA)
    ∇[3] = vec([0 1 0 1; 0 0 0 0] .* DA)
    ∇[4] = vec([0 0 0 0; 0 1 0 1] .* DA)
    ∇[5] = vec([1 0 1 0; 0 0 0 0] .* DB)
    ∇[6] = vec([0 0 0 0; 1 0 1 0] .* DB)
    ∇[7] = vec([0 1 0 1; 0 0 0 0] .* DB)
    ∇[8] = vec([0 0 0 0; 0 1 0 1] .* DB)
    Jtest = hcat(∇...)
    printmln(Jtest)
  else
    @warn "A.*sB.*x .+ sA .- 0.5sB ≠ output, so symbolic check is meaningless as written"
  end
  #  DB = [A A] .* x .- 0.5

  #=function freal(Π,i)
    #sreal.Π .= Π
    sreal.A .= reshape(Π[1:N^2], N, N)
    sreal.B .= reshape(Π[(N^2+1):(2N^2)], N, N)

    @info "got here"

    #sA = view(Ar, :, (i->(i-1)%N+1).(1:N^2))
    #sB = view(Br, :, (i->(i-1)%N+1).(1:N^2))

    #out = vec(abs.((sA .-  sA.*sB .- sB)  ./ (x .+ 1f0)))
    out = sreal(xreal)
    return out[i]
  end

  function makeJforward(Π)
    Nrow = length(vec(sreal(x)))
    out = reduce(vcat, (i->ForwardDiff.gradient((Π)->freal(Π,i), Π)').(1:Nrow))

    return out
  end

  @info "Start with forward diff: "
  if N ≤ 5
    printmln(makeJforward(sreal.Π))
  else
    @btime $makeJforward($sreal.Π)
  end

  =#

  f(i) = vec(s(x))[i]

  #the below vector of matrices
  ι::Vector{CuMatrix{T}} = begin
    ιcalls = eachcol(diagm(ones(T,length(x))))
    (v-> CuMatrix{T}(reshape(v, size(cux)))).(ιcalls)
  end
  #printmln(ι)
    #reshape((j->j==i + 0f0).(1:length(m)), size(x)) |>gpu
  function cuf(i)
    m = s(cux)

    return sum(m .*ι[i])
  end

  function makeJgpu(Π)
    Nrow = length(vec(s(x)))
    g(i) = mwevecgrad(s, Zygote.gradient(()->cuf(i),Π))
    out = reduce(vcat, (i->g(i)').(1:Nrow))
  end

  function makeJ(Π)
    Nrow = length(vec(s(x)))
    g(i) = mwevecgrad(s, Zygote.gradient(()->f(i),Π))
    out = reduce(vcat, (i->g(i)').(1:Nrow))
  end

  @info "Zygote CPU: "
  if N ≤ 5
    JZygote = makeJ(Flux.params(s))
    printmln(JZygote )
    @assert Jtest == JZygote
  else
    Πparam = Flux.params(s)
    @btime $makeJ($Πparam)
  end

  @info "Zygote GPU: "
  if N ≤ 5
    JZygote = makeJgpu(Flux.params(s))
    printmln(JZygote )
    @assert Jtest == JZygote
  else
    Πparam = Flux.params(s)
    @btime $makeJgpu($Πparam)
  end

  function fdiff(Π)
    sreal.Π .= Π

    return sreal(x)
  end

  @info "ForwardDiff cpu"
  if N ≤ 5
    Π = Vector{Real}(sreal.Π)
    JFdif = (ForwardDiff.jacobian(fdiff,Π))
    printmln(JFdif)
    @assert Jtest == JFdif
  else
    Πparam = Flux.params(s)
    @btime $makeJ($Πparam)
  end

  #=function fdiffgpu(Π)
    m = sreal(cux)

    return sum(m .*ι[i])
  end

  @info "ForwardDiff cpu"
  if N ≤ 5
    Π = Vector{Real}(sreal.Π)
    JZygote = (ForwardDiff.jacobian(fdiffgpu,Π))
    printmln(JZygote)
    @assert Jtest == JZygote
  else
    Πparam = Flux.params(s)
    @btime $makeJ($Πparam)
  end=#


end

testjacob(2, Float32)=#
#v=rand(1000)
#m=rand(10,10)
#wrapm = less_unsafe_wrap(typeof(m), source=v, start=101, dims=size(m))
#wrapm2 = less_unsafe_wrap(typeof(m), source=v, start=201, dims=size(m))

#


###
#true mwe of the above
#this utility re-interprets the object(!)
#the object MUST remain referenced, or bad things will happen
#tries to do some basic safety checks, but its still possible to screw things up
#=
using Revise, Flux, Zygote, LinearAlgebra, ForwardDiff, BenchmarkTools

struct mwestruct{TΠ, TA}
  Π::TΠ
  A::TA
  B::TA
end

#constructor
function mwestruct(K::Int, A::TA, B::TA) where TA<:AbstractMatrix

  i = (i->(i-1)%K+1).(1:K^2)
  Π = vcat(vec(A), vec(B)) #this is the parameter vector

  startA = 1
  p = pointer(Π, startA)
  linkedA = unsafe_wrap(TA, p, size(A))

  startB = length(vec(A))+1
  p = pointer(Π, startB)
  linkedB = unsafe_wrap(TA, p, size(B))

  return mwestruct(Π, linkedA, linkedB)
end

#lets call this the loss function (much simpler than the non-mwe version)
function (m::mwestruct)(x)
  abs.((m.A .-  m.A.*m.B .- diag(m.B))  ./ (x .+ 1f0))
end

Flux.@functor mwestruct
Flux.trainable(a::mwestruct) = (a.A,a.B)

#this computes the gradient as a vector corresponding to the parameter vectorW
mwevecgrad(a, gs) = reduce(vcat, [vec(gs[a.A]),vec(gs[a.B])])

printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end

function zygotemwe(N)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  B = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  x = Matrix{Float32}(reshape(collect(1:N^2), (N,N))) .* 2
  s=mwestruct(N, A, B)

  Areal = Matrix{Real}(reshape(collect(1:N^2), (N,N)))
  Breal = Matrix{Real}(reshape(collect(1:N^2), (N,N)))
  sreal=mwestruct(N, Areal, Breal)


  function freal(Π,i)
    sreal.Π .= Π
    out = vec(sreal(x))
    return out[i]
  end

  function makeJforward(Π)

    Nrow = length(vec(sreal(x)))

    out = reduce(vcat, (i->ForwardDiff.gradient((Π)->freal(Π,i), Π)').(1:Nrow))

    return out
  end

  f(i) = vec(s(x))[i]

  function makeJ(Π)
    Nrow = length(vec(s(x)))
    g(i) = mwevecgrad(s, Zygote.gradient(()->f(i),Π))
    out = reduce(vcat, (i->g(i)').(1:Nrow))
  end


  @info "Start with forward diff: "
  if N ≤ 5
    printmln(makeJforward(s.Π))
  else
    @btime $makeJforward($s.Π)
  end

  @info "Zygote CPU: "
  if N ≤ 5
    printmln(makeJ(Flux.params(s)))
  else
    Πparam = Flux.params(s)
    @btime $makeJ($Πparam)
  end


  #@info "f(Π): $(f(s.Π))"

  #Π = Flux.params(s)
  #fJ = fjacobian(f, s.Π)
  #fJgeneric =

  #@info jacobian(s,x,Flux.params(s))
  #=gs = gradient(()->sum(s(x)),Flux.params(s))
  vecgrad = mwevecgrad(s,gs)
  @info "gradient: $(vecgrad)"=#
end

zygotemwe(20)=#

###
#=
struct mwestructv{TΠ, TA,TiA,TsA}
  Π::TΠ
  A::TA
  B::TA
  i::TiA
  sA::TsA
  sB::TsA
end

function mwestructv(K::Int, A::TA, B::TA) where TA<:AbstractMatrix

  AB = hcat(A,B)
  Π = vec(AB)

  viewA = view(AB, :, 1:size(A,2))
  viewB = view(AB, :, (size(viewA',2)+1):size(AB,2))

  @assert sum(viewA .== A) == length(vec(A))
  @assert sum(viewB .== B) == length(vec(B))

  i = (i->(i-1)%K+1).(1:K^2)
  sA = view(viewA', :, (i->(i-1)%K+1).(1:K^2))
  sB = view(viewB', :, (i->(i-1)%K+1).(1:K^2))
  return mwestructv(Π, viewA, viewB, i, sA, sB)
end
#expand(m::mwestructv, ::Any) = zygote_pls_notice_me(m.A, m.sA), zygote_pls_notice_me(m.B, m.sB)
#expand(m::mwestructv, ::CuArray) =  CuMatrix(permutedims(m.A)[:, m.i]), CuMatrix(permutedims(m.B)[:, m.i])
expand(m::mwestructv, ::Any) = zygote_pls_notice_me(m.A, m.A'), zygote_pls_notice_me(m.B, m.B')
expand(m::mwestructv, ::CuArray) =  CuMatrix(permutedims(m.A)), CuMatrix(permutedims(m.B))

function (m::mwestructv)(x)
  sA, sB = expand(m,x)


  sum(abs.((sA .-  sA.*sB .- sB)  ./ (x .+ 1f0)))
end

Flux.@functor mwestructv
Flux.trainable(a::mwestructv) = (a.A,a.B)

function zygotemwev(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N))')
  B = Matrix{Float32}(reshape(collect(1:N^2), (N,N))')
  A .-= mean(A)
  B .-= mean(B)
  #x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N)...)



  s=mwestructv(N, A, B)
  @info "####view version: ####\ncpu: test f: $(s(x))"


  #gs = gradient(()->s(x),Flux.params(s))
  if N < 10
    gs = gradient(()->s(x),Flux.params(s))
    vecgrad = mwevecgrad(s,gs)
    @info "cpu gradient: $(vecgrad)" #WARNING- uncomment soon
  else
    gs = @btime gradient(()->$s($x),Flux.params($s))
  end

  cuA = A |> gpu
  cuB = B |> gpu
  cux = x |> gpu

  cus=mwestructv(N, A, B)
  @info "gpu: test f: $(cus(cux))"
  if N < 10
    cugs = gradient(()->cus(cux),Flux.params(cus))
    #sumgrad = sum(vcat((v->vec(v)).(values(cugs.grads))))
    #vcat((v->vec(v)).(values(cugs.grads)))
    cuvecgrad = mwevecgrad(cus,cugs)
    @info "gpu gradient: $(cuvecgrad)"
  else
    cugs = @btime gradient(()->$cus($cux),Flux.params($cus))
  end
end

zygotemwev(3)=#

#=
using Revise
using Flux, BenchmarkTools, CUDA, CUDA, ForwardDiff, LinearAlgebra, Random
using SparseDiffTools, SparseArrays, SparsityDetection
CUDA.allowscalar(true)
Random.seed!(1111)

function tcudiff(N, ::Type{T} = Float32) where T<:Real
  A = rand(T, N,N)
  Asp = sparse((f->f>T(0.8) ? f : 0f0).(A))
  cuA = A |> gpu
  cuAsp = Asp |> gpu

  function f!(out, A)
    out .= A .+ A .* A .+ T(1) .+ A.^0.5 .+ A .*A .*A .+ T(0.5) .+ exp.(A)
    #out .= A .+ A'*A .+ T(1)
  end

  krn(x) = x + x*x + 1f0 + CUDA.sqrt(x) + x*x*x + 0.5f0 + CUDA.exp(x)
  function f!(out, A::CuMatrix)
    out .= krn.(A)
    #out .= A .+ A'*A .+ T(1)
  end

  function f(A)

    return A .+ A .* A .+ A.^0.5 .+ A .*A .*A .+ T(0.5) .+ exp.(A)
    #out .= A .+ A'*A .+ T(1)
  end

  function f(A::CuMatrix)

    return krn.(A)
    #out .= A .+ A'*A .+ T(1)
  end

  function fsp!(out, A)
    out .= A .+ A .* A .+ A.^0.5 .+ A .*A .*A
    #out .= A .+ A'*A .+ T(1)
  end

  #J = rand(T, N,N)
  #jac = (T).(sparse(pattern))
  #colors = matrix_colors(A)
  #cucolors = matrix_colors(cuA)
  #forwarddiff_color_jacobian!(jac, f, A, colorvec=colors)

  @info "test forward cpu"
  (N<5) && @info "test ∇f forward cpu: $(ForwardDiff.jacobian(f, A))"
  (N>5) && @btime ForwardDiff.jacobian($f, $A)

  @info "test color cpu (inplace)"
  J = rand(T, N^2, N^2)
  cache = SparseDiffTools.ForwardColorJacCache(f!,A, dx = similar(A))
  SparseDiffTools.forwarddiff_color_jacobian!(J, f!, A, cache)
  (N<5) && @info "test ∇f color cpu: $(J)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian!($J, $f!, $A, $cache)

  @info "test color gpu (inplace)"
  cuJ = J |> gpu
  cucache = SparseDiffTools.ForwardColorJacCache(f!,cuA, dx = similar(cuA))
  SparseDiffTools.forwarddiff_color_jacobian!(cuJ, f!, cuA, cucache)
  (N<5) && @info "test ∇f color gpu: $(cuJ)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian!($cuJ, $f!, $cuA, $cucache)


  #sparsity_pattern = jacobian_sparsity(fsp!,similar(Asp), Asp)
  #@info "test sparse cpu (inplace): $(fsp!(similar(Asp), Asp))"
  Jsp = spzeros(Float32, N^2, N^2)
  cachesp = SparseDiffTools.ForwardColorJacCache(fsp!,Asp, dx = similar(Asp))
  #SparseDiffTools.forwarddiff_color_jacobian!(Jsp, fsp!, Asp, cachesp)
  (N<5) && @info "test ∇f sparse cpu: $(Jsp)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian!($Jsp, $f!, $Asp, $cachesp)

  #=@info "test sparse gpu"
  #cuJ = J |> gpu
  cucache = SparseDiffTools.ForwardColorJacCache(f,cuA, dx = similar(cuA))
  cudx = similar(cuA)
  cuJ = SparseDiffTools.forwarddiff_color_jacobian(f, cuA, cucache, dx = cudx, jac_prototype=cuJ)
  (N<5) && @info "test ∇f sparse gpu: $(cuJ)"
  @btime SparseDiffTools.forwarddiff_color_jacobian($f, $cuA, $cucache, dx = $cudx, jac_prototype=$cuJ)=#
end

#tcudiff(5)=#

#=using Revise
using Flux, BenchmarkTools, CUDA, CUDA, ForwardDiff, LinearAlgebra, Random

function mwe(N, ::Type{T}=Float32) where T<:Real
  A::Matrix{T} = rand(T, N,N)
  cuA = A |> gpu

  function f!(out, A)
    out .= A .+ A .* A .+ 1f0
  end

  krn(x) = x + x*x + 1f0
  function f!(out, A::CuMatrix{Float32})
    out .= krn.(A)
  end

  function f(A)
    return A .+ A .* A .+ 1f0
  end

  function f(A::CuMatrix{Float32})
    return krn.(A)
  end

  J = rand(T, N^2, N^2)
  @info "test cpu (inplace)"
  cache = SparseDiffTools.ForwardColorJacCache(f!,A, dx = similar(A))
  SparseDiffTools.forwarddiff_color_jacobian!(J, f!, A, cache)
  (N<5) && @info "test ∇f cpu inplace: $(J)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian!($J, $f!, $A, $cache)

  @info "test cpu (out of place)"
  cacheoos = SparseDiffTools.ForwardColorJacCache(f,A, dx = similar(A))
  J = SparseDiffTools.forwarddiff_color_jacobian(f, A, cacheoos)
  (N<5) && @info "test ∇f cpu oop: $(J)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian($f, $A, $cacheoos)


  @info "test gpu (inplace)"
  cuJ = J |> gpu
  cucache = SparseDiffTools.ForwardColorJacCache(f!,cuA, dx = similar(cuA))
  SparseDiffTools.forwarddiff_color_jacobian!(cuJ, f!, cuA, cucache)
  (N<5) && @info "test ∇f gpu inplace: $(cuJ)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian!($cuJ, $f!, $cuA, $cucache)

  @info "test gpu (outofplace)"
  cucacheoop = SparseDiffTools.ForwardColorJacCache(f,cuA, dx = similar(cuA))
  cuJ = SparseDiffTools.forwarddiff_color_jacobian(f, cuA, cucacheoop)
  (N<5) && @info "test ∇f  gpu oop: $(cuJ)"
  (N>5) && @btime SparseDiffTools.forwarddiff_color_jacobian($f, $cuA, $cucacheoop)

end

mwe(12)=#

#####CUDA arrays
#standard sub-array format

#=using Revise, Zygote, Flux, BenchmarkTools, CUDA, LinearAlgebra
CUDA.allowscalar(false)

###below utility could be used in Capacity
printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end


#this utility re-interprets the object(!)
#the object MUST remain referenced, or bad things will happen
#tries to do some basic safety checks, but its still possible to screw things up
function less_unsafe_wrap(::Type{Tdest};
  source::AbstractVector = error("src vector is required"),
  start::Int = error("start is required"),
  dims::Tuple)::Tdest where Tdest<:AbstractArray

  Ndat = prod(dims)

  T = eltype(Tdest)
  (eltype(source) == T) || error("eltype of src must be of same type as T")

  #this has the effect of making sure we are getting space from somewhere that exists
  protowrap = view(source, start:(start + Ndat -1))
  @assert length(protowrap) == Ndat

  #a final check
  p = pointer(source, start)
  wrap::Tdest = unsafe_wrap(Tdest, p, dims)
  @assert sum(vec(protowrap) .== vec(wrap)) == Ndat

  return wrap
end

function less_unsafe_wrap(::Type{Tdest};
  source::AbstractVector = error("src vector is required"),
  start::Int = error("start is required"),
  dims::Tuple)::Tdest where Tdest<:CUDA.CuArray

  Ndat = prod(dims)

  T = eltype(Tdest)
  (eltype(source) == T) || error("eltype of src must be of same type as T")

  #this has the effect of making sure we are getting space from somewhere that exists
  protowrap = view(source, start:(start + Ndat -1))
  @assert length(protowrap) == Ndat

  #a final check
  p = CUDA.pointer(source, start)
  @info dims
  wrap::Tdest = CUDA.unsafe_wrap(CuArray, p, dims)
  @assert sum(vec(protowrap) .== vec(wrap)) == Ndat

  return wrap
end

struct mwestruct{TΠ, TA,TiA,TsA}
  Π::TΠ
  A::TA
  B::TA
  i::TiA
  sA::TsA
  sB::TsA
end

function mwestruct(K::Int, A::TA, B::TA) where TA<:AbstractMatrix

  i = (i->(i-1)%K+1).(1:K^2)
  Π = vcat((vec(A)), vec(B))

  linkedA = less_unsafe_wrap(TA, source = Π, start = 1, dims=size(A))
  linkedB = less_unsafe_wrap(TA, source = Π, start = length(vec(A))+1, dims=size(B))
  sA = view(A, :, (i->(i-1)%K+1).(1:K^2))
  sB = view(B, :, (i->(i-1)%K+1).(1:K^2))

  @assert linkedA == A
  @assert linkedB == B
  return mwestruct(Π, linkedA, linkedB, i, sA, sB)
end

noopkrn(x) = x
vecc(A) = noopkrn.(vec(A))
function mwestruct(K::Int, A::TA, B::TA) where TA<:CuMatrix

  i = (i->(i-1)%K+1).(1:K^2)
  Π = vcat(vecc(A), vecc(B))

  #@info TA

  linkedA = less_unsafe_wrap(TA, source = Π, start = 1, dims=size(A))
  linkedB = less_unsafe_wrap(TA, source = Π, start = length(vec(A))+1, dims=size(B))
  sA = view(A, :, i)
  sB = view(B, :, i)

  @assert linkedA == A
  @assert linkedB == B
  return mwestruct(Π, linkedA, linkedB, i, sA, sB)
end

#expand(m::mwestruct, ::Matrix{<:AbstractFloat}) = zygote_pls_notice_me(m.A, m.sA), zygote_pls_notice_me(m.B, m.sB)
expand(m::mwestruct, ::Matrix{<:Real}) = Matrix(m.A[:,m.i]), Matrix(m.B[:,m.i])
#expand(m::mwestruct, ::Matrix{<:Real}) = CuMatrix(m.A[:,m.i]), CuMatrix(m.B[:,m.i])
expand(m::mwestruct, ::CuArray) =  CuMatrix(m.A[:,m.i]), CuMatrix(m.B[:,m.i])

krn(sAᵢ,sBᵢ,xᵢ) = sAᵢ*sBᵢ*xᵢ  + sAᵢ - 0.5f0*sBᵢ
function (m::mwestruct)(x)
  #error("got here")
  sA, sB = expand(m,x)
  sA.*sB.*x .+ sA .- 0.5f0*sB
end

function (m::mwestruct)(x::CuArray)
  krn.(m.sA,m.sB,x)
end


mwevecgrad(a, gs) = reduce(vcat, [vec(gs[a.A]),vec(gs[a.B])])

function zygotemwe(N=5)
  A = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  B = Matrix{Float32}(reshape(collect(1:N^2), (N,N)))
  A .-= mean(A)
  B .-= mean(B)
  x = hcat((_->Vector{Float32}(collect(1:N))).(1:N^2)...)

  s=mwestruct(N, A, B)
  @info "cpu: test f:"
  printmln(s(x))

  cuA = CUDA.CuMatrix(A)
  cuB = CUDA.CuMatrix(B)

#  @assert typeof(cuA) <: CUDA.CuArray

  cus=mwestruct(N, cuA, cuB)
  @info "gpu: test f:"
  printmln(cus(cux))


end

zygotemwe(2)=#

###MWE for FLUX
#=using CUDA, Flux
function mwecuda()
  A = rand(100,10)

  @info "A: $(typeof(A))"
  @info "gpu(A): $(typeof(gpu(A)))"
  @info "cu(A): $(typeof(CUDA.cu(A)))"
  @assert typeof(CUDA.cu(A)) <: CUDA.CuArray
  @assert typeof(gpu(A)) <: CUDA.CuArray
end
mwecuda()=#

#=using Revise, CUDA, Flux, Zygote

function cudamwe()
  m = rand(100,1000)
  @info "$(typeof(CUDA.cu(m)))"
end

cudamwe()=#

#=using Revise,Flux, CUDA, Zygote
struct ZygoteMWE{TA<:AbstractMatrix, TB<:AbstractMatrix}
  A::TA
  B::TB
end

function (Θ::ZygoteMWE{TA,TB})() where {TA, TB}

  return sum(abs.(Θ.A .- Θ.B))
end

absA_M_B(a,b) = CUDA.abs(a-b)
function (Θ::ZygoteMWE{TA,TB})() where {TA<:Union{CUDA.CuArray,Flux.CuArrays.CuArray}, TB}

  return CUDA.sum(absA_M_B.(Θ.A, Θ.B))
end

Flux.@functor ZygoteMWE
Flux.trainable(a::ZygoteMWE) = (a.A,a.B)
vecgrad(a, gs) = reduce(vcat, [vec(gs[a.A]),vec(gs[a.B])])

function zygotemwe2(N)
  A = rand(N,N)
  B = rand(N,N)
  s = ZygoteMWE(A,B)
  ∇s =gradient(()->s(), Flux.params(s))
  @info "∇s: $(vecgrad(s,∇s))"

  fluxcuA = Flux.CuArrays.cu(A)
  fluxcuB = Flux.CuArrays.cu(B)
  fluxcus = ZygoteMWE(fluxcuA,fluxcuB)
  ∇fluxcus = gradient(()->fluxcus(), Flux.params(fluxcus))
  @info "∇fluxcus: $(vecgrad(fluxcus,∇fluxcus))"

  cuA = CUDA.cu(A)
  cuB = CUDA.cu(B)
  cus = ZygoteMWE(cuA,cuB)
  ∇cus = gradient(()->cus(), Flux.params(cus))
  @info "∇cus: $(vecgrad(cus,∇cus))"
end

zygotemwe2(2)=#

###Simple example of many regressions
#=using Revise, LinearAlgebra, BenchmarkTools, CUDA, Random
Random.seed!(1111)
printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end


function gels_batched2!(trans::Char,
                      A::Vector{<:CuMatrix{Float32}},
                      C::Vector{<:CuMatrix{Float32}})
    cutrans = CUDA.CUBLAS.cublasop(trans)
    if length(A) != length(C)
        throw(DimensionMismatch(""))
    end
    for (As,Cs) in zip(A,C)
        m,n = size(As)
        mC,nC = size(Cs)
    end
    m,n = size(A[1])
    if m < n
        throw(ArgumentError("System must be overdetermined"))
    end

    nrhs = size(C[1])[2]
    lda = max(1,stride(A[1],2))
    ldc = max(1,stride(A[1],2))
    Aptrs = CUDA.CUBLAS.unsafe_batch(A)
    Cptrs = CUDA.CUBLAS.unsafe_batch(C)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUBLAS.cublasSgelsBatched(
      CUDA.CUBLAS.handle(), cutrans, m, n, nrhs, Aptrs, lda, Cptrs, ldc, [info], infoarray, length(A))
    CUDA.CUBLAS.unsafe_free!(Cptrs)
    CUDA.CUBLAS.unsafe_free!(Aptrs)

    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    A, C, infoarray
end

function gels_batched2(trans::Char,
                     A::Vector{<:CuMatrix{Float32}},
                     C::Vector{<:CuMatrix{Float32}})
    gels_batched2!(trans, deepcopy(A), deepcopy(C))
end

function potrf_batched2!(uplo::Char,
                      A::Vector{<:CuMatrix{Float32}})
    cuuplo = CUDA.CUBLAS.cublasfill(uplo)

    m,n = size(A[1])
    for As in A
        ms,ns = size(As)
        if (size(As) != (m, n))
            throw(DimensionMismatch("Dimensions of batched array entries must be invariant"))
        end
    end
    if m != n
        throw(ArgumentError("Matrix must be square"))
    end


    lda = max(1,stride(A[1],2))
    Aptrs = CUDA.CUBLAS.unsafe_batch(A)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUSOLVER.cusolverDnSpotrfBatched(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, Aptrs, lda, infoarray, length(A))
    CUDA.CUSOLVER.unsafe_free!(Aptrs)

    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    A, infoarray
end

function potrs_batched2!(uplo::Char,
                      A::Vector{<:CuMatrix{Float32}},
                      B::Vector{<:CuMatrix{Float32}})
    cuuplo = CUDA.CUBLAS.cublasfill(uplo)

    if length(A) != length(B)
        throw(DimensionMismatch(""))
    end
    m,n = size(A[1])
    mb,nb = size(B[1])
    for (As, Bs) in zip(A,B)
        if (size(As) != (m, n)) || (size(Bs) != (mb, nb))
            throw(DimensionMismatch("Dimensions of batched array entries must be invariant"))
        end
    end
    if m != n
        throw(ArgumentError("A Matrix must be square"))
    end

    if m != mb
        throw(ArgumentError("B Matrix must have same num of rows as A matrix"))
    end

    nrhs = size(B[1])[2]
    lda = max(1,stride(A[1],2))
    ldb = max(1,stride(B[1],2))
    Aptrs = CUDA.CUBLAS.unsafe_batch(A)
    Bptrs = CUDA.CUBLAS.unsafe_batch(B)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUSOLVER.cusolverDnSpotrsBatched(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, nrhs, Aptrs, lda, Bptrs, ldb, infoarray, length(A))
    CUDA.CUSOLVER.unsafe_free!(Aptrs)
    CUDA.CUSOLVER.unsafe_free!(Bptrs)

    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    A, B, infoarray
end=#

#=function hyperreg(X::Vector{<:CuMatrix{Float32}},
                      y::Vector{<:CuMatrix{Float32}})
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
    CUDA.CUSOLVER.cusolverDnSpotrfBatched(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, Aptrs, lda, infoarray, length(A))
    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    nrhs = size(B[1])[2]
    ldb = max(1,stride(B[1],2))
    Bptrs = CUDA.CUBLAS.unsafe_batch(B)
    info  = zero(Cint)
    infoarray = CUDA.zeros(Cint, length(A))
    CUDA.CUSOLVER.cusolverDnSpotrsBatched(
      CUDA.CUSOLVER.dense_handle(), cuuplo, n, nrhs, Aptrs, lda, Bptrs, ldb, infoarray, length(A))
    CUDA.CUSOLVER.unsafe_free!(Aptrs)
    CUDA.CUSOLVER.unsafe_free!(Bptrs)

    if info != 0
        throw(ArgumentError,string("Invalid value at ",-info))
    end

    B
end=#

#code is borrowed liberally from CUDA.jl/lib/cusolver/wrappers.jl
#                            and CUDA.jl/lib/cublas/wrappers.jl
#=
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


function pad0(M::TM, trows::Int, ::Type{T} = eltype(TM)) where {TM<:AbstractMatrix, T}
  padded = vcat(M, TM(zeros(T, trows-size(M,1), size(M,2))))
  sM = view(padded, 1:size(M,1), 1:size(M,2))
  return (sM, padded)
end

function pad0(V::TV, trows::Int, ::Type{T} = eltype(TV)) where {TV<:AbstractVector, T}
  padded = vcat(V, TV(zeros(T, trows-length(V))))
  sV = view(padded, 1:length(V))
  return (sV, padded)
end

function manyregmwe(T=1000,N=5000,K=10)

  #generate the data
  manysmallX = [rand(rand((N÷2):N),K) for i in 1:T]
  manysmally = [rand(size(x,1)) for x in manysmallX]

  linest(X,y) = cholesky!(X'*X) \ (X'*y)
  bs = [Vector{Float64}(undef, K) for i in 1:T]

  #cpu version
  function cpurunreg()
    Threads.@threads for i in 1:T
      @inbounds bs[i] .= linest(manysmallX[i], manysmally[i])
    end
  end

  manysmallXcu = cu.((X->pad0(X,N)[2]).(manysmallX))
  manysmallycu = cu.((y->Matrix(reshape(pad0(y,N)[2], N, 1))).(manysmally))#cu.(manysmally)

  function gpurunreg()
    bscu = hyperreg(manysmallXcu,manysmallycu)
    return bscu
  end
  cpurunreg()
  ys = gpurunreg()

  @info("cpu1: ")
  @btime $cpurunreg()

  @info("gpu: ")
  @btime $gpurunreg()

  #make sure the solutions are equivelent
  cpusol = reduce(hcat, bs)
  gpusol = reduce(hcat, vec.(Matrix.(ys)))
  cpusol ≈ gpusol || error("cpu ≠ gpu")

  return nothing
end


manyregmwe()=#

#=using Revise, LinearAlgebra, BenchmarkTools, CUDA
function gpuattack(T=1000,N=5000,K=10)
  X = rand(N,K)
  y=rand(1000)

  XtX = X'*X
  Xout = deepcopy(XtX)

  Xcu = cu(X)
  XtXcu = cu(XtX)
  Xoutcu = cu(Xout)


  function cholcu(Xout, X)
    @inbounds Xout .= CUDA.cholesky(X'*X).U

    return nothing
  end

  @cuda CUDA.cholesky!(X'*X)
end

gpuattack()=#
