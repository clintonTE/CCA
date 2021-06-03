using Revise

using DataFrames, Statistics, Distributions, StatsBase,
  LinearAlgebra, BenchmarkTools, Random, Optim, Zygote

Random.seed!(11)

# a simple AR(1) model
struct AR1{Tϕ<:Real, Tb<:Real, Tμ<:Real, Tσ²<:Union{Nothing,Real}}
  ϕ::Tϕ
  b::Tb
  μ::Tμ
  σ²::Tσ²
end

AR1(;ϕ, μ=missing, b=missing, σ²=nothing) = AR1(
  ϕ, coalesce(b, μ*(1-ϕ)), coalesce(μ, b/(1-ϕ)), σ²)


#runs a simple linear regression
function olsar1(x::Vector{T}) where T
  xc = @view x[2:end]
  Lx = @view x[1:(end-1)]

  ϕ = cov(xc, Lx)/(var(Lx))
  b = mean(xc) - ϕ * mean(Lx)

  @assert b ≈ mean(xc) - ϕ * mean(Lx)
  @assert ϕ ≈ cov(xc .- mean(xc), Lx .- mean(Lx))/var(
    Lx .- mean(Lx))

  return AR1(;ϕ,b)
end

function mlear1(x; tol=10^-10)

  #start by centering
  N=length(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  x₁ = x[1]
  #computes the next iteration of ϕ given a previous ϕ
  SS(ϕ) = (x₁^2*(1-ϕ^2) + sum((x̃ .- ϕ*Lx̃).^2))
  function objϕ(ϕ)
    σ² = SS(ϕ)/N
    return (sum((x̃ .- ϕ * Lx̃) .* Lx̃) - (σ²/(1-ϕ^2) - x₁^2)*ϕ)^2
  end

  ϕ = optimize(objϕ, -1.0, 1.0, Brent()) |> Optim.minimizer

  return AR1(;ϕ, μ, σ² = SS(ϕ)/N)
end


function mle2ar1(x; tol=10^-10)

  #start by centering
  N=length(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  x₁ = x[1]

  l(ϕ, σ²) = (-1/2*log(σ²/(1-ϕ^2))
    - (1-ϕ^2)/(2*σ²)*x₁^2
    - (N-1)/2*log(σ²)
    - 1/(2*σ²)*sum((x̃-ϕ*Lx̃).^2)) * -1

  l(Φ) = l(Φ[1], Φ[2])

  ϕ₀ = cov(x̃, Lx̃)/var(Lx̃) #start by guessing the least squares solution
  ε₀= x̃ - ϕ₀*Lx̃
  σ²₀ = var(ε₀)

  ϕ, σ² = optimize(l,  [-1.0,0.0], [1.0,Inf], [ϕ₀, σ²₀], Fminbox(LBFGS()),
    autodiff = :forward) |> Optim.minimizer

  return AR1(;ϕ, μ, σ²)
end


function mle3ar1(x; tol=10^-10)

  #start by centering
  N=length(x)
  x₁ = x[1]
  xc = @view(x[2:end])
  Lx = @view(x[1:(end-1)])

  #=l(ϕ, σ², μ) = (-N/2*log(σ²)
    + 1/2*log(1-ϕ^2)
    - (x₁-μ)^2*(1-ϕ^2)/(2*σ²)
    - sum((xc .- ϕ .* Lx .- μ*(1-ϕ)).^2)/(2*σ²))*-1=#

  x₁ = x[1]
  #computes the next iteration of ϕ given a previous ϕ
  cμ(ϕ) = (x₁*(1+ϕ)+sum(xc-Lx .* ϕ))/(N*(1-ϕ)+2ϕ)
  cσ²(ϕ,μ) = ((x₁-μ)^2*(1-ϕ^2) + sum((xc .- ϕ*Lx .- μ*(1-ϕ)).^2))/N
  function objϕ(ϕ)
    μ = cμ(ϕ)
    σ² = cσ²(ϕ,μ)
    #return (sum((x .- ϕ * Lx) .* Lx) - (σ²/(1-ϕ^2) - x₁^2))^2
    #the below is just proportional to the derivitive of the likelihood wrt ϕ
    err = -ϕ*σ² + (x₁-μ)^2*ϕ*(1-ϕ^2) + sum((xc .- ϕ * Lx .- μ*(1-ϕ)) .* (Lx .- μ))*(1-ϕ^2)
    return err*err
  end

  ϕ = optimize(objϕ, -1.0, 1.0, Brent(), abs_tol=tol) |> Optim.minimizer
  μ = cμ(ϕ)
  σ² = cσ²(ϕ,μ)
  return AR1(;ϕ, μ, σ²)
end

function mle4ar1(x; tol=10^-10)

  #start by centering
  N=length(x)

  x₁ = x[1]
  xc = @view(x[2:end])
  Lx = @view(x[1:(end-1)])

  l(ϕ, σ², μ) = (-N/2*log(2*π*σ²)
    + 1/2*log(1-ϕ^2)
    - (x₁-μ)^2*(1-ϕ^2)/(2*σ²)
    - sum((xc .- ϕ .* Lx .- μ*(1-ϕ)).^2)/(2*σ²))*-1

  l(Φ) = l(Φ[1], Φ[2], Φ[3])

  ϕ₀ = cov(xc, Lx)/var(Lx) #start by guessing the least squares solution
  b₀ = mean(xc)-ϕ₀*mean(Lx)
  μ₀ = b₀/(1-ϕ₀)
  ε₀= xc .- ϕ₀ .* Lx .- b₀
  σ²₀ = var(ε₀)

  ϕ, σ², μ = optimize(l,  [-1.0,0.0,-Inf], [1.0,Inf,Inf], [ϕ₀, σ²₀, μ₀], Fminbox(LBFGS()),
    autodiff = :forward) |> Optim.minimizer

  return AR1(;ϕ, μ, σ²)
end



#below doesn't work
#=
function spm(x)
  μ = mean(x)
  x̃ = @view(x[2:end]) .- μ
  Lx̃ = @view(x[1:(end-1)]) .- μ

  r = SARIMA(x,) |> StateSpaceModels.fit! |> results

  ϕ=r[1]

=#
function checkmath(x)
  N=length(x)
  x₁ = x[1]
  xc = @view(x[2:end])
  Lx = @view(x[1:(end-1)])

  ϕ₀ = 0.7
  μ₀ = 2.7
  σ²₀ = 0.6

  L(ϕ, σ², μ) =  (1-ϕ^2)^(1/2)/σ²^(1/2)*pdf(Normal(), (x₁-μ)/(σ²/(1-ϕ^2))^(1/2))*prod(
    [1/σ²^(1/2)*pdf(Normal(), (x[n]-ϕ*x[n-1]-μ*(1-ϕ))/σ²^(1/2)) for n ∈ 2:N])

  l(ϕ, σ², μ) = (-N/2*log(2*π*σ²)
    + 1/2*log(1-ϕ^2)
    - (x₁-μ)^2*(1-ϕ^2)/(2*σ²)
    - sum((xc .- ϕ .* Lx .- μ*(1-ϕ)).^2)/(2*σ²))

  @assert log(L(ϕ₀, σ²₀, μ₀)) ≈ l(ϕ₀, σ²₀, μ₀) "log(L(ϕ, σ², μ)): $(log(L(ϕ₀, σ²₀, μ₀)))
    l(ϕ, σ², μ): $(l(ϕ₀, σ²₀, μ₀))"

  #check the derivitive wrt to mu as well as the conditional
  lμ(ϕ,σ²=σ²₀, μ=μ₀) = (x₁-μ)*(1-ϕ^2)/σ²+sum((xc .- ϕ .* Lx .- μ*(1-ϕ))*(1-ϕ))/σ²
  lμz(ϕ, μ=μ₀) = gradient(μ->l(ϕ, σ²₀, μ), μ)[1]
  @assert lμ(ϕ₀) ≈ lμz(ϕ₀) "lμ(ϕ): $(lμ(ϕ₀))
    lμz(ϕ): $(lμz(ϕ₀)))"

  cμ(ϕ) = (x₁*(1+ϕ)+sum(xc-Lx .* ϕ))/(N*(1-ϕ)+2ϕ)
  @assert 1 + lμz(ϕ₀, cμ(ϕ₀)) ≈ 1.0 "cμ(ϕ₀): 2$(cμ(ϕ₀))
    lμz(ϕ₀, cμ(ϕ₀)): $(lμz(ϕ₀, cμ(ϕ₀)))"

  #now verify the derivitive and conditional wrt sigma^2
  lσ²(ϕ=ϕ₀, σ²=σ²₀, μ=μ₀) = (-N/(2σ²)
    + (x₁-μ)^2*(1-ϕ^2)/(2σ²^2)
    + sum((xc .- ϕ .* Lx .- μ*(1-ϕ)).^2)/(2σ²^2))
  lσ²z(ϕ=ϕ₀, σ²=σ²₀, μ=μ₀) = gradient(σ²->l(ϕ, σ², μ), σ²)[1]
  lσ²n() =  (l(ϕ₀, σ²₀+1e-6, μ₀)-l(ϕ₀, σ²₀, μ₀))/1e-6
  @assert lσ²() ≈ lσ²z() "lσ²(...): $(lσ²())
    lσ²z(...): $(lσ²z())"

  cσ²(ϕ,μ=μ₀) = ((x₁-μ)^2*(1-ϕ^2) + sum((xc .- ϕ*Lx .- μ*(1-ϕ)).^2))/N
  @assert 1 + lσ²z(ϕ₀, cσ²(ϕ₀)) ≈ 1.0 "cσ²(ϕ₀): 2$(cσ²(ϕ₀))
    lσ²z(ϕ₀, cσ²(ϕ₀)): $(lσ²z(ϕ₀, cσ²(ϕ₀)))"

  #finally, verify the derivitive wrt ϕ
  lϕ(ϕ=ϕ₀, σ²=σ²₀, μ=μ₀) = (-ϕ/(1-ϕ^2) +
    (x₁-μ)^2*ϕ/σ² +
    sum((xc .- ϕ * Lx .- μ*(1-ϕ)) .* (Lx .- μ))/σ²)
  lϕz(ϕ=ϕ₀, σ²=σ²₀, μ=μ₀) = gradient(ϕ->l(ϕ, σ², μ), ϕ)[1]
  @assert lϕ() ≈ lϕz() "lϕ(...): $(lϕ())
    lϕz(...): $(lϕz())"

  cϕ(σ²=σ²₀,μ=μ₀) = optimize(ϕ->lϕ(ϕ)^2, -1.0, 1.0, Brent(), abs_tol=10^-14) |> Optim.minimizer
  @assert isapprox(lϕz(cϕ()), 0.0,atol=1e-5) "cϕ(ϕ₀): 2$(cϕ())
    lϕz(cϕ()): $(lϕz(cϕ()))"
  #lμ()
end

function testar1(;N)
  ϕ = 0.8
  ε = rand(Normal(), N) + rand(N)
  x = missings(Float64, N)
  x[1] = ε[1]


  global x = x
  for t ∈ 2:N
    x[t] = ϕ*x[t-1] + ε[t]
  end
  @info "ϕ=$ϕ, b=$(mean(ε)), μ=$(mean(x))"

  checkmath(x)

  estolsar1 = @btime olsar1($x)
  @info "ols conditional: $estolsar1"

  estmlear1 = @btime mlear1($x)
  @info "mle: $estmlear1"

  estmle2ar1 = @btime mle2ar1($x)
  @info "mle2: $estmle2ar1"

  #=estmle3ar1 = @btime mle3ar1($x)
  @info "mle3: $estmle3ar1"

  estmle4ar1 = @btime mle4ar1($x)
  @info "mle4: $estmle4ar1"=#
end

@time testar1(N=250)
