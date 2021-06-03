using Revise
module BayesMNIST

using Flux, Flux.Data.MNIST, Statistics, LinearAlgebra, StatsBase
using Flux: onehotbatch, onecold, crossentropy, throttle, @epochs
using Flux: softmax, gradient#, params
using DynamicHMC, LogDensityProblems, TransformVariables, DynamicHMC.Diagnostics, Parameters
using CUDAapi, Random, Distributions, DataFrames, Zygote#, Plots
import ForwardDiff
if has_cuda()
    @info "CUDA is on"
    import CuArrays
    CuArrays.allowscalar(false)
end
using Revise

#load the data
#the purpose is to explore flux in a Bayesian context

Random.seed!(11)

const D=784
const K=10

#Lasso Logistic Categorical regression
struct LLCRegression{TY, TX, Tαλ, Tθλ, Tατ, Tθτ}
  Y::TY
  X::TX

  αλ::Tαλ
  θλ::Tθλ

  ατβ::Tατ
  θτβ::Tθτ

end

###Dimensionality
#=
We have N datapoints, D features, K classes
Y ∈

=#


function (problem::LLCRegression)(θ)
    #unpack parameters and data
    @unpack Y, X, αλ, θλ, ατβ, θτβ = problem
    @unpack βs, τβs, λ²s = θ
    Ncovariates, Npoints = size(X)
    Nclasses = length(λ²s)

    @assert size(βs,1) == Ncovariates
    @assert size(βs,2) == length(τβs) == length(λ²s) == Nclasses


    #p(λ²ᵢ|αλ,θλ)
    logλ²pri = sum([logpdf(Gamma(αλ, θλ), λ²s[i]) for i ∈ 1:Nclasses])


    logτβpri = sum([logpdf(Gamma(ατβ, θτβ), τβs[i]) for i ∈ 1:Nclasses])

    #p(βᵢ|λ²ᵢ, τβᵢ, βᵢ)
    logβpri = sum([logpdf(Laplace(0, 1/τβs[i]), λ²s[i]^0.5*sum((abs).(βs[:,i]))) for i ∈ 1:Nclasses])

    μs = softmax(X' * βs, dims=2)
    #loglik = sum(logpdf(Multinomial(Npoints, μs[n, :])))
    loglik = sum([logpdf(Multinomial(1, μs[n, :]), Y[n]) for n ∈ 1:Npoints])


    return logλ²pri + logτβpri + logβpri + loglik
end

#hyperpriors
const Π = Dict{Symbol, Float64}(
  :αλ=>1.0f0,
  :θλ=>1.0f0,

  :ατβ=>1.0f0,
  :θτβ=>1.0f0
)

function LLCtransform(p::LLCRegression)
    return as((βs=as(Array, size(p.X,1), K),
      τβs=as(Array, asℝ₊, K),
      λ²s=as(Array, asℝ₊, K)))
end


function dynamicmnist(::Type{V} = AbstractVector, ::Type{T}=Float64) where {V,T}
  imgs = MNIST.images()
  labels = MNIST.labels()

  CuArrays.allowscalar(false)

  #make each image into a feature vector
  #concatenate each feature vector into a matrix
  preX =hcat(float.(reshape.(imgs,:))...)

  #add an intercept
  preX = vcat(preX, ones(1,size(preX,2)))[:,1:1000]
  X = preX# |> gpu

  #creates a one-hot matrix of 0 vectors

  Y = Matrix{Float64}(onehotbatch(labels[1:1000], 0:9))

  #error("stop")
  #make the problem
  p = LLCRegression(Y, X, Π[:αλ], Π[:θλ], Π[:ατβ], Π[:θτβ])

  #transform the support
  t = LLCtransform(p)
  P = TransformedLogDensity(t, p)

  #differentiate via AD
  ∇P =  ADgradient(:ForwardDiff, P)

  @info "beginning MCMC"
  testβs = rand(D+1,K)# |> gpu
  testτβs = rand(K)
  testλ²s = rand(K)
  #p((βs=testβs, τβs=testτβs, λ²s=rand(K)))
  results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, 10, reporter = LogProgressReport()#,
    #initialization = (ϵ = 0.1, ),
    #warmup_stages = fixed_stepsize_warmup_stages()
    )

  posterior = transform.(t, results.chain)

  β_posterior = first.(transform.(t, results.chain))
  display(mean(β_posterior))

  λ_posterior = last.(transform.(t, results.chain))
  display(mean(λ_posterior) .^ 0.5)

  println("sizeX: $(size(X))")
  println("sizeY: $(size(Y))")

end



end # module
