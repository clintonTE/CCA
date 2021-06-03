using Revise,LinearAlgebra, Random, StatsBase, Statistics, Distributions

Random.seed!(1111)

printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end

function testex2(;
  Ξ = ([4,8], [4,-8], [8,4], [-4, -8], [-4, 8], [-8, -4]),
  dims=(N=6, T=10^5, K=2))

  #set up the factor structure and dividend realizations
  Θ0 = reduce(hcat, [rand(Normal(), dims.T) for k in 1:dims.K])
  @assert size(Θ0) ≡ (dims.T, dims.K)
  D̄ = 1.0
  ϵ = [rand(Normal(), dims.T) for n ∈ 1:dims.N]
  Dvec = [D̄ .+ Θ0*Ξ[n] .+ ϵ[n] for n ∈ 1:dims.N]
  D = reduce(hcat, Dvec)
  Σ = cov(D)

  #select a baseline and a bunch of comp portfolios
  w1 = [1/6 for n ∈ 1:6]
  walt = [rand(6) for i ∈ 1:10^5]
  walt ./= sum.(walt)

  varbase = w1'*Σ*w1
  @info "baseline: w1 $w1: $(w1'*Σ*w1)"

  vars = (w-> w'*Σ*w).(walt)
  for (i,w) ∈ enumerate(walt)
    if vars[i] > varbase
      #@info "w$i $w: $(vars[i])"
      nothing
    else
      @warn "w$i $w: $(vars[i])"
    end
  end

  printmln(Σ)
end

testex2()
