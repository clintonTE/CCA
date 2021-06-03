using Revise

using DataFrames, Statistics, Distributions, StatsBase,
  LinearAlgebra, BenchmarkTools, Random, Optim, Zygote, UnPack, Random

MFloat64 = Union{Missing, Float64}

Random.seed!(11)

# a simple AR(1) model
struct AR1{TV<:AbstractVector, TM<:AbstractMatrix}
  ϕ::TV
  μ::TV
  Σ::TM
  ΣcholL::TM
end

#primary constructor, performs an integrity check on Σ integrity checks
function AR1(;ϕ::TV, μ::TV,  Σ::TM, ΣcholL::TM=cholesky(Σ).L |> Matrix) where {TV,TM}

  @assert (Σ .≈ ΣcholL * ΣcholL') |> all

  return AR1(ϕ, μ, Σ, ΣcholL)
end

centerar1(ar1::AR1) = AR1(ar1.ϕ, 0.0, ar1.Σ, ar1.ΣcholL)


function AR1(Xy::AbstractMatrix{T}) where T

  X̃y = @view Xy[2:end, :]
  LX̃y = @view Xy[1:(end-1), :]

  #compute the persistence as a regression
  #the special case is for a constant value in all terms, which would give ϕ=0.0, σ²=0.0, μ=c
  ϕ = ((x̃,Lx̃)->ifelse(maximum(x̃)-minimum(x̃) + 1.0 ≈ 1.0, 0.0, cov(x̃,Lx̃)/var(Lx̃))
    ).(eachcol(X̃y), eachcol(LX̃y))

  μ = mean(Xy, dims=1) |> vec

  ε = vcat((Xy[1,:] .- μ)', X̃y .- LX̃y .* ϕ' .- μ')
  Σ = cov(ε)
  ΣcholL = semicovcholl(ε)


  #println("μ: $μ, ϕ: $ϕ, b: $b, sum(Xy): $(sum(Xy)), Xy[1:10, :]: $(Xy[1:10, :])")
  #=global sigma = Σ
  display(Σ)
  println("b=$b
    ϕ: $ϕ")=#

  @assert length(μ) == length(ϕ) == size(Xy,2) == size(Σ,1)
  @assert (Σ .≈ ΣcholL * ΣcholL') |> all
  @assert (ϕ .≈ ((x̃,Lx̃)-> ifelse(maximum(x̃) - minimum(x̃) + 1.0 ≈ 1.0, 0.0,
    cov(x̃ .- mean(x̃), Lx̃ .- mean(Lx̃))/var(Lx̃ .- mean(Lx̃)))
    ).(eachcol(X̃y), eachcol(LX̃y))) |> all

  return AR1(;ϕ,μ,Σ,ΣcholL)
end


function simulatear1!(ar1::AR1, Xyt::AbstractMatrix;
    ) where TV<:AbstractVector
  @unpack ϕ, Σ, μ = ar1
  K::Int = length(ϕ)

  #first simulate the centered distribution
  #we will do this all in place

  #first generate the shocks
  #can't plug the covar matrix into MVnormal directly as its
  Xyt .= ar1.ΣcholL * rand!(MvNormal(K,1.0), Xyt)


  Xyt[:,1] .*= 1 ./ (1 .- ϕ.^2) .^ 0.5 # the first entry is simulated from the unconditional distribution

  #use a profiled breakpoint- large matrices its better to vectorize, small and not
  if K ≤ 3
    for t ∈ 2:size(Xyt,2)
      for r ∈ 1:size(Xyt,1)
        Xyt[r,t] += ϕ[r] * Xyt[r,t-1]
      end
    end
  else
    for t ∈ 2:size(Xyt,2)
      Xyt[:,t] .+= ϕ .* Xyt[:,t-1]
    end
  end
  Xyt .+= μ
  return nothing
end
#non-inplace version for testing
function simulatear1(ar1, N)
  res = Matrix{Float64}(undef, length(ar1.ϕ), N, )
  simulatear1!(ar1, res)
  return res'
end

function compareoncear1xy!(ar1::AR1, Xyt::AbstractMatrix,
  comparefunc::Tcomparefunc) where Tcomparefunc

  simulatear1!(ar1, Xyt)
  comparefunc(Xyt')
end

#this runs the simulations- basically a parametric bootstrap under the null
#statfunc is computed for each pass, while aggfunc is computed on the statistics
#at the end of the pass
function comparear1undernull(X::Tx,y::Ty,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-5,
    ) where {Tx<:AbstractMatrix{Float64}, Ty <: AbstractVector{Float64}, Tstatfunc, Taggfunc}

  Xy = hcat(X,y)
  K = size(X,2)
  ar1 = AR1(Xy)

  results = Dict(s => missings(Float64, Nsimulations, K) for s ∈ keys(statfunc(Xy)))
  @assert (k->ismissing.(results[k]) |> all).(keys(results)) |> all


  #this assumes stationary results are meaningless
  if checkstationarity && ((ar1.ϕ .≥ 1.0-stationaritytol) |> any)
    #display(ar1.Σ)
    #throw("ar1.Σ is not stationary!: eigen(ar1.Σ).values: $(eigen(ar1.Σ).values)")
    return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
  end

  Xyts = [deepcopy(Xy' |> Matrix{Float64}) for t ∈ 1:Threads.nthreads()]

  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for n ∈ 1:Nsimulations
  #WARNING: threading crashes it for some unknown reason
  #for n ∈ 1:Nsimulations

    t::Int = Threads.threadid()
    #res = compareoncear1xy!(xar1, xs[1], yar1, ys[1], statfunc)
    res = compareoncear1xy!(ar1, Xyts[t], statfunc;)

    #this structure allows us to compute multiple statistics
    #println(res)
    for s ∈ keys(res)

      @assert ismissing.(results[s][n, :]) |> all

      results[s][n,:] .= res[s] |> vec

    end
  end

  return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
end

comparear1undernull(X::Tx,y::Ty,args...; kwargs...
    ) where {Tx <: AbstractMatrix, Ty<:AbstractVector} = comparear1undernull(
    X|>Matrix{Float64},y|>Vector{Float64},args...; kwargs...)


#this is just a wrapper for the 2-vector case
function compareundernull(x::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc, comparefunc::Tcomparefunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-6
    ) where {Tx<:AbstractVector, Tstatfunc, Taggfunc, Tcomparefunc}


  return comparefunc(reshape(x,length(x),1) |> Matrix{Float64}, y, statfunc, aggfunc;
    Nsimulations, checkstationarity, stationaritytol)
end

#wrapper for the wrapper
comparear1undernull(x::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-4
    ) where {Tx<:AbstractVector, Tstatfunc, Taggfunc} = compareundernull(
      x,y,statfunc,aggfunc, comparear1undernull;
      Nsimulations, checkstationarity, stationaritytol)

#=function comparearima11undernull(X::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-4,
    dXts = (size(X,2) > 1 ?
      [similar(X[2:end, :]' |> Matrix{Float64}) for t ∈ 1:Threads.nthreads()] :
      [similar(X[2:end, :] |> vec |> Vector{Float64}) for t ∈ 1:Threads.nthreads()]),
    dys = [similar(y[2:end]) for t ∈ 1:Threads.nthreads()],
    ) where {Tx<:AbstractMatrix, Tstatfunc, Taggfunc}


  #difference everything
  dX = reduce(hcat, (xi->xi[2:end] .- xi[1:(end-1)]).(eachcol(X)))

  #handle the intercept as a special case
  for (dxi, xi) ∈ zip(eachcol(dX), eachcol(X))
    if maximum(xi) - minimum(xi) + 1.0 ≈ 1.0
      dxi .= mean(xi)
    end
  end

  dx = dX |> eachcol |> collect
  dy = y[2:end] .- y[1:(end-1)]

  K = length(dx)

  dxar1 = (AR1.(dx))# .|> centerar1
  dyar1 = AR1(dy)# |> centerar1

  results = Dict(s => missings(Float64, Nsimulations, K) for s ∈ keys(statfunc(dX,dy)))
  @assert (k->ismissing.(results[k]) |> all).(keys(results)) |> all

  #this assumes stationary results are meaningless
  if checkstationarity && max((dxar1i -> abs(dxar1i.ϕ)).(dxar1)..., abs(dyar1.ϕ)) ≥ 1.0-stationaritytol
    return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
  end

  ϕdx::Vector{Float64} = dxar1 .|> dxar1i->dxar1i.ϕ
  σdx::Vector{Float64} = dxar1 .|> dxar1i->dxar1i.σ²^0.5
  μdx::Vector{Float64} = dxar1 .|> dxar1i->dxar1i.μ

  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for n ∈ 1:Nsimulations
  #for n ∈ 1:Nsimulations

    t::Int = Threads.threadid()
    #res = compareoncear1xy!(xar1, xs[1], yar1, ys[1], statfunc)
    res = compareoncear1xy!(dXts[t], dyar1, dys[t], statfunc; ϕx=ϕdx, σx=σdx, μx=μdx)


    #this structure allows us to compute multiple statistics
    for s ∈ keys(res)

      @assert ismissing.(results[s][n, :]) |> all
      results[s][n,:] .= res[s]
    end
  end

  return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
end=#


#this runs the simulations- basically a parametric bootstrap under the null
#statfunc is computed for each pass, while aggfunc is computed on the statistics
#at the end of the pass
function comparearima11undernull(X::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-5,
    ) where {Tx<:AbstractMatrix, Tstatfunc, Taggfunc}

  Xy = hcat(X,y)
  dXy = Xy[2:end,:] .- Xy[1:(end-1),:]
  K = size(X,2)

  #handle the intercept as a special case
  for (dxi, xi) ∈ zip(eachcol(dXy), eachcol(Xy))
    if maximum(xi) - minimum(xi) + 1.0 ≈ 1.0
      dxi .= mean(xi)
    end
  end

  ar1 = AR1(dXy)

  results = Dict(s => missings(Float64, Nsimulations, K) for s ∈ keys(statfunc(dXy)))
  @assert (k->ismissing.(results[k]) |> all).(keys(results)) |> all

  #this assumes stationary results are meaningless
  if checkstationarity &&  ((ar1.ϕ .≥ 1.0-stationaritytol) |> any)
    return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
  end

  dXyts = [deepcopy(dXy' |> Matrix{Float64}) for t ∈ 1:Threads.nthreads()]

  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for n ∈ 1:Nsimulations
  #WARNING: threading crashes it for some unknown reason
  #for n ∈ 1:Nsimulations

    t::Int = Threads.threadid()
    #res = compareoncear1xy!(xar1, xs[1], yar1, ys[1], statfunc)
    res = compareoncear1xy!(ar1, dXyts[t], statfunc;)

    #this structure allows us to compute multiple statistics
    #println(res)
    for s ∈ keys(res)

      @assert ismissing.(results[s][n, :]) |> all
      results[s][n,:] .= res[s] |> vec
    end
  end

  return (; Dict(s=>aggfunc(results[s], s) for s ∈ keys(results))...)
end

comparearima11undernull(X::Tx,y::Ty,args...; kwargs...
    ) where {Tx <: AbstractMatrix, Ty<:AbstractVector} = comparearima11undernull(
    X|>Matrix{Float64},y|>Vector{Float64},args...; kwargs...)

comparearima11undernull(x::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-4
    ) where {Tx<:AbstractVector, Tstatfunc, Taggfunc} = compareundernull(
      x,y,statfunc,aggfunc, comparearima11undernull;
      Nsimulations, checkstationarity, stationaritytol)

nonparametricnullp(v::AbstractVector{T}, x) where T<:Real = (
  x ≥ 0 ? sum(v .> x)/length(v) : sum(v .< x)/length(v))
nonparametricnullp(v::AbstractVector, x) = nonparametricnullp(v |> skipmissing |> collect, x)
nonparametricnullp(M::AbstractMatrix, b) = (
  (v,bᵢ)->nonparametricnullp(v |> skipmissing |> collect, bᵢ)).(eachcol(M),b)


#this creates a cholesky decomposition for the covariance of a data matrix
#under the special case where one of the columns is all constant
#the purpose is to use this to generate sensible random numbers
function semicovcholl(X::TM) where TM
  K = size(X,2)

  #determine the columns that are the same
  defcols = mapslices(c->!(maximum(c)-minimum(c) + 1.0 ≈ 1.0), X, dims=1) |> vec

  #compute the full covariance matrix
  defΣ = cov(@view X[:, defcols])
  defL = cholesky(defΣ).L

  ΣcholL = zeros(eltype(TM), K, K)

  ΣcholL[defcols, defcols] .= defL

  return ΣcholL
end

function testsemicovcholl()
  Xdat = rand(Normal(),10^4,5) .+ (rand(Normal(),10^4) .* [1,2,3,4,5]') .+ collect(1:5)'
  X = deepcopy(Xdat)
  X[:,4] .= 1
  Σcheck = cov(X)
  ΣcholLcheck = deepcopy(Σcheck)
  defcols = [1,2,3,5]
  ΣcholLcheck[defcols,defcols] .= cholesky(cov(X[:,defcols])).L

  #check the computational integrity
  ΣcholL = semicovcholl(X)
  @assert (Σcheck .≈ ΣcholL * ΣcholL') |> all
  @assert (ΣcholLcheck .≈ ΣcholL) |> all

  #now check what we plan to do with it- random number generation
  X[:,4] .= Xdat[:,4]
  Σcheck = cov(X)
  ΣcholL = semicovcholl(X)
  @assert (Σcheck .≈ ΣcholL * ΣcholL') |> all

  μs = mean(X,dims=1) |> vec
  simcheck = rand(MvNormal(μs, Σcheck), 10^6)'
  sim = (ΣcholL * rand(Normal(0,1),5, 10^6).+ μs)'

  atol = 0.02
  @assert ((μ, μcheck)->isapprox(μ, μcheck; atol)
    ).(mean(sim, dims=1), mean(simcheck, dims=1)) |> all "
      mean(sim, dims=1) = $(mean(sim, dims=1))
      mean(simcheck, dims=1) = $(mean(simcheck, dims=1))"

  if !all(((σij, σcheckij)->isapprox(σij, σcheckij; atol=atol^0.5)).(cov(sim), cov(simcheck)))
    @info "cov(sim)"
    display(cov(sim))
    @info "cov(simcheck)"
    display(cov(simcheck))
    throw("failed simulation check on Σ")
  end
  @info "passed semicovcholl check"


end

function testar1(;N,K)
  testsemicovcholl()

  ϕ = [[0.8/k for k ∈ 1:K]; 0.8; ]

  #simulate a covariance matrix
  ND = max(N, min(N^2, 10^6))
  D = rand(Normal(0,1), ND,K+1) .+ [1 ./ collect(1:K);1.0;]' .* rand(Normal(0,1), ND)
  #Σtrue = reduce(vcat, 4 .* [[1/(k+r) for k ∈ 1:(K+1)]' for r ∈ 1:(K+1)])
  Σtrue  = cov(D)
    #println(Σtrue)

  μtrue = [1.0 for k ∈ 1:(K+1)]
  ε = (rand(MvNormal(Σtrue), N) .+ μtrue)'
  Xy = missings(Float64, N,K+1)
  Xy[1,:] .= ε[1]



  for t ∈ 2:N
    Xy[t,:] .= Xy[t-1, :] .* ϕ .+ ε[t, :]
  end
  global Xy = Xy


  #assign the components
  X = Xy[:,1:(end-1)]
  y = Xy[:,end]
  @info "ϕ=$ϕ, b=$(mean(ε, dims=1)), μ=$(mean(Xy, dims=1))"


  atol = 1/N^0.5
  actar1=AR1(;ϕ=[ϕ[1]], μ=[2.5], Σ=reshape([2.0], 1,1))
  xsim = simulatear1(actar1, 10^6)
  @assert isapprox(mean(xsim), 2.5; atol)  "mean(xsim): $(mean(xsim))"
  @assert isapprox(var(xsim), 2/(1-ϕ[1]^2), atol=atol^0.5) "var(xsim): $(var(xsim)), expected=$(2/(1-ϕ[1]^2))"
  println("size of xsim: $(size(xsim))")
  ar1sim = AR1(xsim)
  @info "ar1sim: $(ar1sim)"
  @assert isapprox(ar1sim.ϕ[1], ϕ[1]; atol)
  @assert isapprox(ar1sim.μ[1], 2.5; atol)
  @assert isapprox(ar1sim.Σ[1], 2; atol)


  Nsimulations=10^4
  dmeantest = 0.2
  #recall the compare func will pass the name of the relevant stat func when it is called
  aggfunc(v, s) = nonparametricnullp(v, Dict(:cor=>0.5, :dmean=>dmeantest, :beta=>1.0)[s])
  function statfunc(Xy)
    Xs = view(Xy,:, 1:(size(Xy,2)-1))
    ys = view(Xy, :, size(Xy,2))

    res = Dict(:beta=>cholesky!(Symmetric(Xs'*Xs))\(Xs'*ys),
      :cor=>cor(Xs, ys),
      :dmean=>mean(Xs, dims=1) .- mean(ys))


    return res
  end

  @info "std(Xy, dims=1): $(std(Xy, dims=1))"
  @info "expected std(Xy): $((2 ./ (1 .- ϕ.^2)).^0.5)"

  #derive the distribution for the difference in means of two AR(1) dists
  R = [0.8^abs(i-j) for i ∈ 1:N, j ∈ 1:N]
  sidmean = (1/N^2*(ones(N)' * ( R .* 2/(1-ϕ[1]^2)) * ones(N)))^0.5
  ρdmean = Σtrue[1,end] ./ (Σtrue[1,1]*Σtrue[end,end])^0.5
  #because its a difference, use a negative sign for the t-stat
  Sdmean = (2*sidmean^2-2*ρdmean*sidmean^2)^0.5
  Epdmean = cdf(Normal(), - dmeantest/Sdmean)

  #println("X[1:5,1]: $(X[1:5,1])")
  x=X[:,1]
  res = comparear1undernull(x, y, statfunc, aggfunc; Nsimulations)
  @info "compareundernull result: $res (Epdmean: $Epdmean)"

  Xtest = X#[:,1]
  @btime comparear1undernull($Xtest, $y, $statfunc, $aggfunc; Nsimulations=$Nsimulations)
  #@btime comparear1undernull($x, $y, $statfunc, $aggfunc; Nsimulations=$Nsimulations)
  #throw("simulation complete")

  X = hcat(X[:,1:2], ones(N)) |> Matrix{Float64}
  function mstatfunc(Xy)
    Xs = view(Xy,:, 1:(size(Xy,2)-1))
    ys = view(Xy, :, size(Xy,2))


    try
      return Dict(
        :beta=>cholesky!(Symmetric(Xs'*Xs))\(Xs'*ys),
        :cor=>cor(Xs,ys))
    catch err
      global X = X
      global y = y
      throw("$err:
      iter: $iter
      Xy[1:5,:]: $(Xy[1:5, :])")
    end
  end
  maggfunc(M,s) = (v->nonparametricnullp(v, Dict(:cor=>0.5, :beta=>0.3)[s])).(eachcol(M))

  res = comparear1undernull(X, y, mstatfunc, maggfunc; Nsimulations)
  @info "compareundernull matrix result: $res"


  ix = cumsum(X[:,1])
  iy = cumsum(y)
  res = comparear1undernull(ix, iy, statfunc, aggfunc; Nsimulations)
  @info "integrated compareundernull (mispecified) dmean: $res"
  res = comparearima11undernull(ix, iy, statfunc, aggfunc; Nsimulations)
  @info "integrated compareundernull dmean: $res"
  res = comparearima11undernull(X, y, mstatfunc, maggfunc; Nsimulations)
  @info "compareundernull matrix result: $res"

end


#@time testar1(N=1000, K=3)
