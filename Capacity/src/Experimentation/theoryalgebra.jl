using Revise
using LinearAlgebra, StatsBase, Statistics, Distributions, Random
#Random.seed!(11)

#for displaying matrices
printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end


#generates a well-defined covariance matrix
function randomcovar(K, N=10^4; displaymatrix=false)
  ϵ = rand(Uniform(-1,1), N,K)
  ψ = rand(Uniform(-1,1), N)
  β = (k->k*(-1)^k).(1:K)

  X = β' .* ψ
  @assert size(X) == (N,K)
  X .+= ϵ
  Σ = cov(X)

  displaymatrix && printmln(Σ)
  return Σ
end

function checkmatrixalgebra(dims::NamedTuple; testiter=10^5)
  atol = dims.K/testiter^0.5
  #generate the characteristic matrix
  Ξ = rand(Laplace(), dims.N, dims.K+1)
  #Ξ .-= mean(Ξ, dims=1)
  Ξ[:,end] = rand(dims.N) #last column is the market

  #signal parameters
  δi = rand(dims.K+1)
  Σθ = randomcovar(dims.K+1)
  σ²ϵ = rand(dims.N) .+ 0.0001
  #Σϵ = Diagonal(σ²ϵ)
  Σϵ = randomcovar(dims.N)

  Rf = 1.2
  P = rand(TriangularDist(1.0,11.0), dims.N)
  N1 = ones(dims.N)

  #distributions
  distθ = MvNormal(δi, Σθ)
  distϵ = MvNormal(zeros(dims.N), Σϵ)

  #test covariance of portfolio
  wi = rand(dims.N) .+ 0.01
  wi ./= sum(wi)
  ωi = wi ./ P
  vervec = [wi'* (Ξ * rand(distθ) ./ P .+ rand(distϵ) ./ P .- Rf * N1) for i in 1:testiter]
  vervar = var(vervec)
  testvar = ωi' * (Ξ*Σθ*Ξ' .+ Σϵ) * ωi
  msg  = "covariance of portfolio test with Δ $(abs(vervar-testvar)) and max atol $atol"
  isapprox(testvar, vervar, atol=atol) ? @info("passed $msg") : @warn("failed $msg")

  #test expectation
  vermean = mean(vervec)
  testmean = ωi' * (Ξ * δi .- P * Rf)
  msg  = "expectation of portfolio test with Δ $(abs(vermean-testmean)) and max atol $atol"
  isapprox(testmean, vermean, atol=atol) ? @info("passed $msg") : @warn("failed $msg")


  #test woodbury
  Σθinv = svd(Σθ)\I
  Σϵinv = Σϵ\I
  IN = Diagonal(I,dims.N)

  verω = (Ξ *Σθ * Ξ' + Σϵ)\I * (Ξ * δi .- P .* Rf)
  testw1 = Σϵinv*(I - Ξ*((Σθinv .+ Ξ'*Σϵinv*Ξ)\I)*Ξ'*Σϵinv)*(Ξ*δi .- P .* Rf)
  @assert verω ≈ testw1
  @info "passed test prewoodbury"

  si = (I - ((Σθinv .+ Ξ'*Σϵinv*Ξ)\I) * (Ξ'*Σϵinv * Ξ)) * δi
  κ = (I - Ξ*((Σθinv .+ Ξ'*Σϵinv*Ξ)\I) * (Ξ'*Σϵinv))
  testw2 = Σϵinv*(Ξ*si .- κ * P .* Rf)
  @assert verω ≈ testw2

  testw3 = P .* Σϵinv*(Ξ*si  .- κ * P .* Rf)
  #testw3 = Σϵinv*(Ξ*si .- κ * P .* Rf)
  ver3 = P .* verω
  @assert ver3 ≈ testw3
  @info "passed test woodbury"

  ###Market amalgamation test
  Ξᴷ = Ξ[:,1:(end-1)]
  ξₘ = Ξ[:,end] .* si[end] .- κ * P .* Rf
  siᴷ = si[1:(end-1)]
  testamalg1 = ωi'*Σϵinv*Ξᴷ*siᴷ .+ ωi'*Σϵinv*ξₘ

  Ξᴺ=deepcopy(Ξ)
  Ξᴺ[:,end] .-= κ * P .* Rf/si[end]
  testamalg2 = ωi'*Σϵinv*Ξᴺ*si

  veramalg = ωi'*Σϵinv*(Ξ*si .- κ * P .* Rf)

  @assert testamalg1 ≈ veramalg
  @assert testamalg2 ≈ veramalg
  @info "passed market amalgamation test"

  #Signal test
  verxi = verω .* P#Σϵinv * (Ξ*si - κ*P*Rf) .* P
  testxi = [Σϵinv * Ξ[:,k] .* si[k] .* P for k ∈ 1:(dims.K+1)]
  testxi[end] .= Σϵinv * (Ξ[:,end] .* si[end] .- κ*P .* Rf) .* P
  @assert verxi ≈ sum(testxi)
  @info "passed portfolio weight test. Portfolio weights (net): $(sum.(testxi))"


end

function checkwork()
  dims = (N=1000, K=5)
  checkmatrixalgebra(dims,testiter=10^5)
end

@time checkwork()
