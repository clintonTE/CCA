using Revise

using DataFrames, Statistics, Distributions, StatsBase,
  LinearAlgebra, BenchmarkTools, Random, Optim, Zygote, UnPack

MFloat64 = Union{Missing, Float64}

Random.seed!(11)

# a simple AR(1) model
struct AR1{T<:Real, Tσ²<:Union{Nothing,Real}}
  ϕ::T
  b::T
  μ::T
  σ²::Tσ² #note- this is σ²e, not the overall σ² of the series
end

AR1(;ϕ, μ=missing, b=missing, σ²=nothing) = AR1(
  ϕ, coalesce(b, μ*(1-ϕ)), coalesce(μ, b/(1-ϕ)),σ²)

centerar1(ar1::AR1) = AR1(ar1.ϕ, 0.0, 0.0, ar1.σ²)

function estimatear1ols(x::Vector{T}) where T
  x̃ = @view x[2:end]
  Lx̃ = @view x[1:(end-1)]

  ϕ = cov(x̃, Lx̃)/(var(Lx̃))
  b = mean(x̃) - ϕ * mean(Lx̃)
  μ = b/(1-ϕ)

  ε = similar(x)
  ε[1] = x[1] - μ
  ε[2:end] = x̃ .- ϕ .* Lx̃
  σ² = var(ε)

  @assert b ≈ mean(x̃) - ϕ * mean(Lx̃)
  @assert ϕ ≈ cov(x̃ .- mean(x̃), Lx̃ .- mean(Lx̃))/var(
    Lx̃ .- mean(Lx̃))

  return AR1(;ϕ,b,σ²)
end

function estimatear1mle(x; tol=10^-10)

  #start by centering
  N=length(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  x₁ = x[1]
  #computes the next iteration of ϕ given a previous ϕ
  cσ²(ϕ) = (x₁^2*(1-ϕ^2) + sum((x̃ .- ϕ*Lx̃).^2))/N
  function objϕ(ϕ)
    σ² = cσ²(ϕ)
    return (sum((x̃ .- ϕ * Lx̃) .* Lx̃) - (σ²/(1-ϕ^2) - x₁^2)*ϕ)^2
  end

  ϕ = optimize(objϕ, -1.01, 1.01, Brent()) |> Optim.minimizer

  return AR1(;ϕ, μ, σ² = cσ²(ϕ))
end


function estimatear1mlecheck(x; tol=10^-10)

  #start by centering
  N=length(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  x₁ = x[1]

  l(ϕ, σ²) = (-N/2*log(2*π*σ²)
    + 1/2*log(1-ϕ^2)
    - x₁^2*(1-ϕ^2)/(2*σ²)
    - sum((x̃ .- ϕ .* Lx̃).^2)/(2*σ²)) * -1

  l(Φ) = l(Φ[1], Φ[2])

  ϕ₀ = cov(x̃, Lx̃)/var(Lx̃) #start by guessing the least squares solution
  ε₀= x̃ - ϕ₀*Lx̃
  σ²₀ = var(ε₀)

  ϕ, σ² = optimize(l,  [-1.0,0.0], [1.0,Inf], [ϕ₀, σ²₀], Fminbox(LBFGS()),
    autodiff = :forward) |> Optim.minimizer

  return AR1(;ϕ, μ, σ²)
end

AR1(x::AbstractVector; estimationmethod=estimatear1mle) = estimationmethod(x)

function checkmath2(x)
  if length(x) > 250
    x = x[1:250] |> deepcopy
    @warn "length(x) > 250 so truncating length to 250. This is only for the purposes
    of computing the full (not log) likelihood"
  end

  N=length(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  ϕ₀ = 0.75
  σ²₀ = 2.2
  x₁ = x̃[1]

  L(ϕ, σ²) =  (1-ϕ^2)^(1/2)/σ²^(1/2)*pdf(Normal(), x₁/(σ²/(1-ϕ^2))^(1/2))*prod(
    ((x̃ᵢ,Lx̃ᵢ)->1/σ²^(1/2)*pdf(Normal(), (x̃ᵢ-ϕ*Lx̃ᵢ)/σ²^(1/2))).(x̃,Lx̃))

  l(ϕ, σ²) = (-N/2*log(2*π*σ²)
    + 1/2*log(1-ϕ^2)
    - x₁^2*(1-ϕ^2)/(2*σ²)
    - sum((x̃ .- ϕ .* Lx̃).^2)/(2*σ²))

  @assert log(L(ϕ₀, σ²₀)) ≈ l(ϕ₀, σ²₀,) "log(L(ϕ, σ²)): $(log(L(ϕ₀, σ²₀)))
    l(ϕ, σ²): $(l(ϕ₀, σ²₀))"

  #now verify the derivitive and conditional wrt sigma^2
  lσ²(ϕ=ϕ₀, σ²=σ²₀,) = (-N/(2σ²)
    + x₁^2*(1-ϕ^2)/(2σ²^2)
    + sum((x̃ .- ϕ .* Lx̃).^2)/(2σ²^2))
  lσ²z(ϕ=ϕ₀, σ²=σ²₀,) = gradient(σ²->l(ϕ, σ²,), σ²)[1]
  lσ²n() =  (l(ϕ₀, σ²₀+1e-6,)-l(ϕ₀, σ²₀,))/1e-6
  @assert lσ²() ≈ lσ²z() "lσ²(...): $(lσ²())
    lσ²z(...): $(lσ²z())"

  cσ²(ϕ) = (x₁^2*(1-ϕ^2) + sum((x̃ .- ϕ*Lx̃).^2))/N
  @assert 1 + lσ²z(ϕ₀, cσ²(ϕ₀)) ≈ 1.0 "cσ²(ϕ₀): $(cσ²(ϕ₀))
    lσ²z(ϕ₀, cσ²(ϕ₀)): $(lσ²z(ϕ₀, cσ²(ϕ₀)))"

  #finally, verify the derivitive wrt ϕ
  lϕ(ϕ=ϕ₀, σ²=σ²₀) = (-ϕ/(1-ϕ^2) +
    x₁^2*ϕ/σ² +
    sum((x̃ .- ϕ * Lx̃) .* Lx̃)/σ²)
  lϕz(ϕ=ϕ₀, σ²=σ²₀,) = gradient(ϕ->l(ϕ, σ²), ϕ)[1]
  @assert lϕ() ≈ lϕz() "lϕ(...): $(lϕ())
    lϕz(...): $(lϕz())"

  ϕcerr(ϕ, σ²=σ²₀)= sum((x̃ .- ϕ * Lx̃) .* Lx̃) - (σ²/(1-ϕ^2) - x₁^2)*ϕ
  cϕ(σ²=σ²₀,) = optimize(ϕ->ϕcerr(ϕ)^2, -1.01, 1.01, Brent(), abs_tol=10^-14) |> Optim.minimizer
  @assert isapprox(lϕz(cϕ()), 0.0,atol=1e-5) "cϕ(ϕ₀): 2$(cϕ())
    lϕz(ϕcerr()): $(lϕz(ϕcerr()))"
  #lμ()
end

#simulates an AR1 series
function simulatear1!(ar1::AR1, x::AbstractVector)
  @unpack ϕ, b, σ², μ = ar1

  σ² === nothing && throw(" σ² must not be nothing")

  #first simulate the centered distribution
  #we will do this all in place
  σ=σ²^0.5
  rand!(Normal(0,σ), x)
  x[1] *= 1/(1-ϕ^2)^0.5 # the first entry is simulated from the unconditional distribution
  for t ∈ 2:length(x)
    x[t] += ϕ * x[t-1]
  end
  x .+= μ

  return x
end

#which approach is faster?
function simulatear1(ar1::AR1, T::Int)
  @unpack ϕ, σ², μ = ar1

  σ² === nothing && throw(" σ² must not be nothing")
  #T::Int = length(x)

  #first simulate the centered distribution
  #we will do this all in place
  σ=σ²^0.5
  x = rand(Normal(), T) .* σ
  x[1] *= 1/(1-ϕ^2)^0.5
  for t ∈ 2:T
    x[t] += ϕ * x[t-1]
  end
  x .+= μ

  return x
end

simulatear1!(ar1::AR1{T,Tσ²}, N::Int) where {T,Tσ²} = simulatear1!(ar1, Vector{T}(undef, N))

compareoncear1xy!(xar1::AR1, x::AbstractVector, yar1::AR1,  y::AbstractVector,
  comparefunc::Tcomparefunc
  ) where Tcomparefunc = comparefunc(simulatear1!(xar1,x), simulatear1!(yar1,y))


#this runs the simulations- basically a parametric bootstrap under the null
#statfunc is computed for each pass, while aggfunc is computed on the statistics
#at the end of the pass
function comparear1undernull(x::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-4
    ) where {Tx, Tstatfunc, Taggfunc}
  xar1 = AR1(x) |> centerar1
  yar1 = AR1(y) |> centerar1

  results = Dict(s => Vector{MFloat64}(undef, Nsimulations) for s ∈ keys(statfunc(x,y)))

  #this assumes stationary results are meaningless
  if checkstationarity && max(abs(xar1.ϕ), abs(yar1.ϕ)) ≥ 1.0-stationaritytol
    xar1 = AR1(NaN, NaN, xar1.μ, xar1.σ²)
    yar1 = AR1(NaN, NaN, xar1.μ, xar1.σ²)
  end

  xs = [similar(x) for t ∈ 1:Threads.nthreads()]
  ys = [similar(y) for t ∈ 1:Threads.nthreads()]

  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for n ∈ 1:Nsimulations
    t::Int = Threads.threadid()
    res = compareoncear1xy!(xar1, xs[t], yar1, ys[t], statfunc)

    #this structure allows us to compute multiple statistics
    for s ∈ keys(res)
      @assert results[s][n] ≡ missing
      results[s][n] = res[s]
    end
  end

  return (; Dict(s=>aggfunc(skipmissing(results[s]), s) for s ∈ keys(results))...)
end

function comparearima11undernull(x::Tx,y::AbstractVector,
    statfunc::Tstatfunc, aggfunc::Taggfunc; Nsimulations::Int,
    checkstationarity::Bool = true, stationaritytol::Float64=10^-4) where {Tx, Tstatfunc, Taggfunc}

    #difference everything
    dx = x[2:end] .- x[1:(end-1)]
    dy = x[2:end] .- y[1:(end-1)]

    dxar1 = AR1(dx) |> centerar1
    dyar1 = AR1(dy) |> centerar1

    #make integration funcitons
    #ix(δx) = cumsum([x[1]; δx])
    #iy(δy) = cumsum([y[1]; δy])

    #dstatfunc(δx, δy) = statfunc(δx |> ix, δy |> iy)

    results = Dict(s => Vector{MFloat64}(undef, Nsimulations) for s ∈ keys(statfunc(dx,dy)))

    #this assumes stationary results are meaningless
    if checkstationarity && max(abs(dxar1.ϕ), abs(dyar1.ϕ)) ≥ 1.0-stationaritytol
      dxar1 = AR1(NaN, NaN, dxar1.μ, dxar1.σ²)
      dyar1 = AR1(NaN, NaN, dxar1.μ, dxar1.σ²)
    end

    dxs = [similar(dx) for t ∈ 1:Threads.nthreads()]
    dys = [similar(dy) for t ∈ 1:Threads.nthreads()]

    #draw a sample and compute the aggregate function, sotring the reuslts in an array
    Threads.@threads for n ∈ 1:Nsimulations
      t::Int = Threads.threadid()
      res = compareoncear1xy!(dxar1, dxs[t], dyar1, dys[t], statfunc)

      #this structure allows us to compute multiple statistics
      for s ∈ keys(res)
        @assert results[s][n] ≡ missing
        results[s][n] = res[s]
      end
    end

    return (; Dict(s=>aggfunc(skipmissing(results[s]), s) for s ∈ keys(results))...)
end

nonparametricnullp(v::AbstractVector{T}, x) where T<:Real = (
  x ≥ 0 ? sum(v .> x)/length(v) : sum(v .< x)/length(v))
nonparametricnullp(v, x) = nonparametricnullp(v |> skipmissing |> collect, x)


function testar1(;N)
  ϕ = 0.8
  ε = rand(Normal(), N) + rand(Normal(1,1), N)
  x = missings(Float64, N)
  x[1] = ε[1]


  global x = x
  for t ∈ 2:N
    x[t] = ϕ*x[t-1] + ε[t]
  end
  @info "ϕ=$ϕ, b=$(mean(ε)), μ=$(mean(x))"

  checkmath2(x)

  #=estolsar1 = @btime AR1($x, estimationmethod=estimatear1ols)
  @info "ols conditional: $estolsar1"

  estmlear1 = @btime AR1($x, estimationmethod=estimatear1mle)
  @info "mle: $estmlear1"

  estmle2ar1 = @btime AR1($x, estimationmethod=estimatear1mlecheck)
  @info "mle2: $estmle2ar1"=#

  actar1=AR1(;ϕ, b=0.5, σ²=2)
  xsim = simulatear1(actar1, 10^6)
  @assert isapprox(mean(xsim), 2.5, atol=0.01)  "mean(xsim): $(mean(xsim))"
  @assert isapprox(var(xsim), 2/(1-ϕ^2), atol=0.01) "var(xsim): $(var(xsim))"
  ar1sim = AR1(xsim)
  @info "ar1sim: $(ar1sim)"
  @assert isapprox(ar1sim.ϕ, ϕ, atol=0.01)
  @assert isapprox(ar1sim.b, 0.5, atol=0.01)
  @assert isapprox(ar1sim.μ, 2.5, atol=0.01)
  @assert isapprox(ar1sim.σ², 2, atol=0.01)

  xsim! = simulatear1!(actar1, Vector{Float64}(undef,10^6))
  @assert isapprox(mean(xsim), 2.5, atol=0.01)  "mean(xsim): $(mean(xsim))"
  @assert isapprox(var(xsim), 2/(1-ϕ^2), atol=0.01) "var(xsim): $(var(xsim))"
  ar1sim = AR1(xsim!)
  @info "ar1sim: $(ar1sim)"
  @assert isapprox(ar1sim.ϕ, ϕ, atol=0.01)
  @assert isapprox(ar1sim.b, 0.5, atol=0.01)
  @assert isapprox(ar1sim.μ, 2.5, atol=0.01)
  @assert isapprox(ar1sim.σ², 2, atol=0.01)

  Nsimulations=10^4
  dmeantest = 0.2
  #recall the compare func will pass the name of the relevant stat func when it is called
  aggfunc(v, s) = nonparametricnullp(v, Dict(:cor=>0.1, :dmean=>dmeantest)[s])
  statfunc(x,y) = Dict(:cor=>cor(x,y), :dmean=>mean(x)-mean(y))

  @info "std(x): $(std(x))"
  @info "expected std(x): $((2/(1-ϕ^2))^0.5)"
  Sdmean = (2*2/(1-ϕ^2)/length(x))^0.5
  Epdmean = cdf(Normal(), - dmeantest/Sdmean)
  res = comparear1undernull(x, deepcopy(x), statfunc, aggfunc; Nsimulations)
  @info "compareundernull result: $res (Epdmean: $Epdmean)"
  y = deepcopy(x)

  @btime comparear1undernull($x, $y, $statfunc, $aggfunc; Nsimulations=$Nsimulations)
  throw("simulation complete")


  ix = cumsum(x .+ 0.2)
  iy = cumsum(x .+ 0.2)
  res = comparear1undernull(ix, deepcopy(ix), statfunc, aggfunc; Nsimulations)
  @info "compareundernull dmean: $res"
  res = comparearima11undernull(ix, deepcopy(ix), statfunc, aggfunc; Nsimulations)
  @info "compareundernull dmean: $res"

end

@time testar1(N=1000)
