using Revise
module BayesMNIST

using Flux, Flux.Data.MNIST, Statistics, LinearAlgebra, StatsBase
using Flux: onehotbatch, onecold, crossentropy, throttle, @epochs
using Flux: kldivergence
using CUDAapi, CuArrays, Random, Distributions, DataFrames#, Plots
using Revise

#load the data
#the purpose is to explore flux in a Bayesian context

Random.seed!(11)


  include("MNISTModels.jl")
  #include("TuringBNN.jl")

end # module
