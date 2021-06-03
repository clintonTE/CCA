using Revise, Random, LinearAlgebra, Optim, StatsBase,
  Statistics, Distributions, COSMO, BenchmarkTools, Plots

Random.seed!(11)

#compute the minvol portfolio given a return target, expected returns z, and Σ
function minvolw(μ, z, Σ, Σinv = Σ\I)
  K= length(z)
  some1s = ones(K)

  A = some1s' * Σinv * some1s
  B = some1s' * Σinv * z
  C = z' * Σinv * z

  λ =(μ * B - C) /(B^2 - A*C)
  γ = (μ * A - B) /(A*C -B^2)


  w = λ * Σinv * some1s .+ γ * Σinv * z

  return w
end

#same as before but acquires the solution numerically
function nminvolw(μ, z, Σ)
  K= length(z)
  some1s = ones(K)

  m = COSMO.Model()
  c1 = COSMO.Constraint(ones(1,K), -1.0, COSMO.ZeroSet)
  c2 = COSMO.Constraint(z', -μ, COSMO.ZeroSet)
  COSMO.assemble!(m, Σ, zeros(K), [c1; c2])
  res = COSMO.optimize!(m)


  return res.x

  #return w
end


function gengraphs(;N=10^5, K=15, Δ=0.001)
  Erm = 0.08 #expected market reutn
  rf = 0.02 #risk free rate
  σm = 0.1 #market standard deviation

  ivol = rand(K)*0.4 .+ 0.1 #stock ivols
  beta = rand(Normal(1,0.2), K) #stock betas

  rm = rand(Normal(Erm, σm), N) #the market returns
  noise = reduce(hcat, (σe->rand(Normal(0, σe), N)).(ivol))
  r = rm .* beta' .+ noise

  z = mean(r, dims=1) |> vec
  Σ = cov(r)

  #sanity check
  Σinv = Σ\I
  w = minvolw(0.1, z, Σ, Σinv)
  nw = nminvolw(0.1, z, Σ)
  @assert w ≈ nw

  minvolσ(μ) = minvolw(μ, z, Σ, Σinv) |> (w)->(w'*Σ*w)^0.5
  nminvolσ(μ) = nminvolw(μ, z, Σ) |> (w)->(w'*Σ*w)^0.5

  #compute the minimum variance frontier
  μs = 0.01:Δ:0.2 |> collect
  σs = similar(μs)
  Threads.@threads for i ∈ 1:length(μs)
    σs[i] = minvolσ(μs[i])
  end

  #compute global min var portfolio
  some1s = ones(K)
  A = some1s' * Σinv * some1s
  B = some1s' * Σinv * z

  wg = Σinv * some1s / A
  σg = (wg'*Σ*wg)^(0.5)
  rg = z' * wg

  #compute the tangency portfolio
  wt = Σinv * (z .- some1s .* rf) ./ (B-rf*A)
  rt = z'*wt
  σt = (wt'*Σ*wt)^(0.5)

  pth = "C:\\Users\\Clinton\\Dropbox\\AAtawork\\ChernovInvestments\\Investments\\2020\\Exams\\final"
  p = plot(σs, μs,
    yaxis = ("E(r)", (0.0,0.16), :none),
    yticks=(0.0:0.02:0.16),
    xaxis=("vol", (0.0,0.4),),
    xticks=(0.0:0.05:0.4),
    title="Minimum Variance Frontier",
    legend=:none,
    linestyle=:dot,
    linecolor=:black)

  pans = p |> deepcopy

  plot!(pans, [0.0; σt; σs[μs .> rt]], [rf; rt; μs[μs .> rt]],
    linestyle=:dash, seriestype=:path, linecolor=:black, linewidth=4,
    annotations=(σt*1.5, rt*1.2,
    Plots.text("min var frontier\n (rf asset, no borrowing)", :left, 7,)))
  plot!(pans, [0.0; σt; σt*10], [rf; rt; (rt-rf)*10],
    linestyle=:solid, seriestype=:path, linecolor=:black, linewidth=2,
    annotations=(σt*1.1, rt*1.6,
    Plots.text("min varfrontier\n (rf asset, borrowing)", :left, 7,)))

  #scatter!(pans, [0.0], [rf], markersize=6)
  scatter!(pans, [σt], [rt], markersize=6,
    annotations=(σt+0.005, rt-0.0025, Plots.text("tangency portfolio", :left, 7,)))
  scatter!(pans, [σg], [rg],markersize=6,
    annotations=(σg+0.005, rg-0.0025, Plots.text("global min var portfolio", :left, 7,)))
  savefig(p, "$pth\\problem3.pdf")
  savefig(pans, "$pth\\problemans3.pdf")

end


@time gengraphs()
