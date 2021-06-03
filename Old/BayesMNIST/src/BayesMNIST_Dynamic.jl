using Revise
module BayesMNIST

using Statistics, LinearAlgebra, StatsBase, Flux, Flux.Data.MNIST, DistributionsAD
using Flux: softmax, onehotbatch#, params
using DynamicHMC, LogDensityProblems, TransformVariables, DynamicHMC.Diagnostics, Parameters
using CuArrays, Random, Distributions, Zygote, BenchmarkTools, InteractiveUtils#, Plots
using CUDAapi
if has_cuda()
    @info "CUDA is on"
    import CuArrays
    CuArrays.allowscalar(false)
end
import ForwardDiff

using Revise

#load the data
#the purpose is to explore flux in a Bayesian context

Random.seed!(11)

const D=784
const K=10

#Lasso Logistic Categorical regression
struct LCRegression{TY, TX, Tατ, Tθτ}
  Y::TY
  X::TX

  ατβ::Tατ
  θτβ::Tθτ

end


function (problem::LCRegression)(θ)
    #unpack parameters and data
    @unpack Y, X, ατβ, θτβ = problem
    @unpack βs, τβs = θ
    Ncovariates, Npoints = size(X,1),size(X,2)
    Nclasses = size(βs,2)

    @assert size(βs,1) == Ncovariates
    @assert size(βs,2) == length(τβs) == Nclasses

    #p(τβᵢ)
    logτβpri = sum([logpdf(Gamma(ατβ, θτβ), τβs[i]) for i ∈ 1:Nclasses])

    #p(βᵢ|τβᵢ)
    logβpri = sum([logpdf(MvNormal(zeros(Ncovariates), 1/τβs[i]), βs[:,i]) for i ∈ 1:Nclasses])

    μs = softmax(X' * βs, dims=2)
    #μs = (xᵢ->softmax(βs'*xᵢ)).(X)
    #  display(μs[1])
    #p(Y|β,τβ)
    #loglik = sum([logpdf(Categorical(μs[n, :]), Y[n]) for n ∈ 1:Npoints])
    #println(size(μs[1]))
    #println(size(Y[1]))
    #display(logpdf(Multinomial(1, μs[1]), Y[1]))
    loglik = sum([logpdf(Multinomial(1,μs[n,:]), Y[:,n]) for n ∈ 1:Npoints])


    return logτβpri + logβpri + loglik
end

#hyperpriors
const Πmwe = Dict{Symbol, Float64}(
  :ατβ=>1.0,
  :θτβ=>1.0
)

function LCtransform(p::LCRegression)
    return as((βs=as(Array, D+1, K),
      τβs=as(Array, asℝ₊, K)))
end


function dynamicmnistmwe(Npoints::Int = 1000,
    ::Type{V} = AbstractVector, ::Type{T}=Float32) where {V,T}
  imgs = MNIST.images()
  labels = MNIST.labels()

  #make each image into a feature vector
  #concatenate each feature vector into a matrix
  preX =hcat(float.(reshape.(imgs,:))...)

  #add an intercept
  preX = vcat(preX, ones(1,size(preX,2)))[:,1:Npoints]
  #X = (Vector).(eachcol(preX)) #|>gpu
  X = preX
  preY = onehotbatch(labels[1:Npoints], 0:9)
  Y = Matrix(preY)

  #Y=labels[1:5_000]

  #make the problem
  p = LCRegression(Y, X, Πmwe[:ατβ], Πmwe[:θτβ])

  #transform the support
  t = LCtransform(p)
  P = TransformedLogDensity(t, p)

  #differentiate via AD
  ∇P =  ADgradient(:Zygote, P)

  @info "beginning MCMC"

  testβs = rand(D+1,K)# |>gpu
  testτβs = rand(K)

  #InteractiveUtils.@code_warntype p((βs=testβs, τβs=testτβs))
  #@btime $p((βs=$testβs, τβs=$testτβs))
  results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 1000, reporter = LogProgressReport())
  serialize("results.jls", results)

  posterior = transform.(t, results.chain)

  β_posterior = first.(transform.(t, results.chain))
  display(mean(β_posterior))

end

#dynamicmnistmwe()

end # module
