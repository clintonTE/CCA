using Revise
using DataFrames, BenchmarkTools, LinearAlgebra, Distributions, StatsBase, Statistics, Random
Random.seed!(11)

#set up a siimilar regression situation
function crosssection(N,K; loud::Bool=true)

  #generate our data
  ϵ = rand(Laplace(0,1), N)
  Ξ = hcat((k->rand(Laplace(0,2),N)).(1:K)...)
  þ = 1 ./ (þᵢ->abs(þᵢ) ≤ 10^-5 ? 10.0^-5 * sign(þᵢ) : þᵢ).(rand(Chisq(1), N)./4)
  #þ .*= rand(Bernoulli(),N)
  X = Ξ .* þ
  #abstotX = sum(abs.(X), dims=1)

  function normalizemat!(X)
    X .-= mean(X, dims=1)
    abstot = sum(abs.(X), dims=1)
    #abstot₊ = sum((x->x≥0 ? x : 0.).(X), dims=1)
    #abstot₋ = sum((x->x<0 ? -x : 0.).(X), dims=1)
    #X .= hcat(((t)-> (x->x≥0 ? x/(abstot₊[t[1]]) : x/(abstot₋[t[1]])).(t[2])).(enumerate(eachcol(X)))...)
    X ./= abstot .* 2
  end
  #normalize the unobservable actual
  normalizemat!(X)

  #normalize the observable
  normalizemat!(Ξ)


  #σx = std(X, dims=1)
  #X ./= σx
  loud && display(abstotX)
  #isplay(X)
  X = hcat(ones(N),X)
  β = [50.0; (Float64).(1:K) .* 300  ...;]

  y = X*β  .+ ϵ

  #display(β)
  function err(Xqr,y,b)
    u = y.-X*Xb
    Λ = u .* u
    RinvQp = Xqr.R\I * Xqr.Q'
    Σ = RinvQp * Diagonal(Λ) * (RinvQp')
    return diag(Σ).^0.5
  end

  #errchk(X,y,b) = diag(((X'*X)\I) * X'* Diagonal((y.-X*Xb) .* (y.-X*Xb)) * X * ((X'*X)\I)).^0.5

  loud && @info "Test non-confounded version (β=$β)"
  Xqr = qr(X)
  Xrinv = Xqr.R\I
  Xb =  Xrinv * Xqr.Q' * y
  eX = (y.-X*Xb)
  errX = err(Xqr,y,Xb)#diag(eX' * eX ./ (N-K-1) .* (Xqr.R'*Xqr.R)\I) .^ 0.5
  #errXchk = errchk(X,y,Xb)

  loud && @info "Test observed version"
  Ξ1 =hcat(ones(N),Ξ)
  Ξqr = qr(Ξ1)
  Ξrinv = Ξqr.R\I
  Ξb =  Ξrinv * Ξqr.Q' * y
  errΞ = err(Ξqr,y,Ξb)
  loud && @info "Non-confounded version: $(Xb). ErrX = $errX"

  ΣΞbabs = sum(abs.(Ξ1 .* Ξb'), dims=1)
  loud && @info "Results: Raw b: $Ξb, raw error: $errΞ \n"
  loud && @info "Abs total: $(ΣΞbabs)"

  #=@info "Test confounded version"
  Xqr = qr(X)
  Xrinv = Xqr.R\I
  Xb =  Xrinv * Xqr.Q' * y
  @info "Xb: $(Xb) (vs non-confounded of $(Ξb))"
  ΣXb = vec(sum(, dims=1))=#
  return (Ξb=Ξb, ΣΞbabs=ΣΞbabs)
end


tups = vcat((crosssection(2_000,2, loud=false) for i in 1:1000)...)
#display(tups)
Ξb, As = vcat((t->t.Ξb').(tups)...), vcat((t->t.ΣΞbabs).(tups)...)
display(Ξb)
@info "mean raw: $(mean(Ξb, dims=1))"
@info "mean: $(mean(As, dims=1))"
@info "std: $(std(As, dims=1))"

#crosssection(2_000,2)
