using Revise

using Juno, UnPack
using CSV, DataFrames, Serialization, Dates, Statistics, StatsBase, BenchmarkTools
using LinearAlgebra, Distributions, Distributed, VegaLite, Random
using CodecZlib, CodecZstd, GLM, XLSX, InteractiveUtils, Base.Iterators
using Optim, Pkg, SparseArrays#, BlackBoxOptim
using MCMCChains
using CUDA, GPUArrays
using Transducers, BangBang, LineSearches, LoopVectorization,Zygote

using Finometrics
using CUDA: CuVector, CuMatrix, CuArray



function oddball(N::Int = 10_000_000)
  df = DataFrame(g=rand(1:(N÷100) |> collect, N), x=rand(N))
  preconditioner = demeanξ ∘ quantiles
  df[!,:upperlower] = Vector{Union{Missing,Float64}}(undef, size(df,1))

  gdf = groupby(df, :g)
  Threads.@threads for i ∈ 1:length(gdf)
    oddball(gdf[i], preconditioner)
  end
end

@inline quantiles(ξraw::AbstractVector{<:Union{Real,Missing}}) = (
  length(ξraw) .- competerank(ξraw,rev=true) .+ 1)./length(ξraw)

@inline demeanξ(ξraw::AbstractVector{<:Union{Real, Missing}}) = ξraw .- mean(skipmissing(ξraw))

function oddball(df::AbstractDataFrame, preconditioner; thresholds = (lower=0.33, upper=0.67))

  ξ = preconditioner(df[:,:x])
  ranks::Vector{Float64} = ξ |> quantiles
  df[:,:upperlower] .= ((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower))
  (minimum(df.upperlower) ≈ -1.0 * maximum(df.upperlower)) || error("Ranking Failed!!
    $((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower))
    typeof:$(typeof((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower)))")

end

@time oddball()
