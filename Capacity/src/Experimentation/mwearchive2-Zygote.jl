using Revise






#Question- do residuals of the within model equal residuals of a FI model?
#=using FixedEffectModels, DataFrames
function testwithinresid( N::Int=100, K=10; NG = N ÷ 2K)
  #form the test data
  X = rand(N,K)
  X .*= collect(1:K)'
  Gval = (i->i% NG).(1:N)

  for (i,r) ∈ enumerate(eachrow(X))
    r .+= Gval[i]
  end
  β = sqrt.(collect(1:K))
  y = X*β .+ K*rand(N)

  df = DataFrame(X)
  xnames = propertynames(df)
  df.y = y
  numericnames = propertynames(df)
  df.G = (g->Symbol(:G,g)).(Gval)
  df.i = (i->Symbol(:i,i)).(1:N)

  ffixed = @eval @formula y ~ $(Meta.parse(join(xnames, "+"))) + fe(G)

  dfwi = deepcopy(df)
  dfwis = groupby(df, :G)
  demean(v::AbstractVector{T}) where T<:Real = v .- mean(v)
  demean(v::Any) = v
  dfwi = combine(dfwis, [numericnames...; :i] .=> demean)

  #fix the names
  rename!(dfwi,(n->replace(n, "_demean"=>"")).(names(dfwi)))

  sort!(df, :i)
  sort!(dfwi, :i)


  display(dfwi)

  fwithin = @eval @formula y ~ $(Meta.parse(join(xnames, "+"))) + 0

  rfixed = reg(df, ffixed, save=:residuals)
  rwithin = reg(dfwi, fwithin, save=:residuals)

  df.resid = residuals(rfixed, df)
  dfwi.resid = residuals(rwithin, dfwi)
  println(rfixed)
  println(rwithin)

  @info "Residual difference: $(sum(abs.(df.resid .- dfwi.resid)))"
  return nothing
end

testwithinresid(100,2)=#

#=using CUDA, LinearAlgebra,SparseArrays
function hybridmatrixmultmwe()
  #generate test data
  dX = rand(Float32, 1000,10) |> cu
  S = (x->ifelse(x>0.99f0,1f0,0f0)).(rand(Float32, 1000,1000)) |> SparseMatrixCSC
  dS = S |> CUSPARSE.CuSparseMatrixCSR
  dShyb = dS |> CUSPARSE.switch2hyb

  println("sum of product (CSR):", sum(dS*dX))
  println("sum of product (hybrid):", sum(dShyb*dX))
end
hybridmatrixmultmwe()=#

#=using Revise, CUDA, LinearAlgebra, BenchmarkTools
simplefunc(Aᵢ,Bᵢ) = CUDA.abs(Aᵢ - Bᵢ)
function mweptr(resultsalloc, As, Bs)

  #use the pre-allocation
  results = resultsalloc.results
  Threads.@threads for i ∈ 1:length(As)
    results[i] .= simplefunc.(As[i], Bs[i])
  end

  return reduce(vcat,results)
end

function mweptr(N)
  As::Vector{CuMatrix} = [rand(rand(1:N), 3) |> cu for i ∈ 1:N]
  Bs::Vector{CuMatrix} = [rand(size(As[i],1), 3) |> cu for i ∈ 1:N]
  results::Vector{CuMatrix} = [CUDA.zeros(size(As[i],1), 3) |> cu for i ∈ 1:N]

  resultssum = 0f0
  resultsalloc = (results=results,)
  for i ∈ 1:10
    resultssum += mweptr(resultsalloc, As, Bs) |> sum
  end

  @info "Completed test with sum $(resultssum)"
end

@btime mweptr(500)=#

#=using Revise, LinearAlgebra, CUDA, BenchmarkTools


CUDA.allowscalar(false)
function mwesum(N)
  cpuv= rand(Float32, N)
  print("\nStandard cpu sum: ")
  @btime sum($cpuv)

  cuv=CUDA.cu(cpuv)
  print("Standard cuda sum: ")
  @btime CUDA.@sync sum($cuv)

  print("Summing using allocating a vector of 1s and dot: ")
  onesum(v) = dot(cuv, CUDA.ones(N))
  @btime CUDA.@sync $onesum($cuv)

  cuvones = CUDA.ones(N)
  print("Summing just using dot and a pre-allocated vector of 1s: ")
  onesum2(v) = dot(cuv, cuvones)
  @btime CUDA.@sync $onesum2($cuv)
end

sleep(0.5)
mwesum(256^2)=#


###MWE zygote
using Zygote, CUDA
CUDA.allowscalar(false)

function zygotedevice2host(M::TM, ::Type{T}=eltype(TM)
  ) where {T<:Real, TM<:CuArray}
  buf = Zygote.Buffer(Array{T}(undef, size(M)))
  copyto!(buf, M)
  return copy(buf)
end

function mwe(::Type{T}=Float32; N=5) where T
  x = CUDA.rand(N)

  #verify functions work
  M=CUDA.rand(T, N,N^2)
  v=CUDA.rand(T, N)

  function cuones(::Type{T}, N)
    #v1 = CUDA.rand(T, N) .* T(0.0) .+ T(1.0)
    v1 = ones(N) |> CuVector{T}
  end

  brokenfunctions = (
    shouldwork=()->CUDA.sum(M),
    ones=()->CUDA.ones(N),
    cuones=()->cuones(T, N),
    zeros=()->CUDA.zeros(N),
    fill=()->CUDA.fill(2f0, N),
    cumprod=()->CUDA.cumprod(M, dims=2),
    device2host=()-> rand(Float64, N,N) |> CuMatrix{T} |> zygotedevice2host |>sum,
    hcat=()->CUDA.hcat(M,M),
    hcat2=()->CUDA.hcat(M,cuones(T, N)),
    reduce=()->CUDA.reduce(+,M),
    reduce2=()->CUDA.reduce(hcat,[M,M])
    )

  testcuda(f,x) = CUDA.sum(f() .+ x)
  for n ∈ propertynames(brokenfunctions)
    print("testing $n: ")
    f = brokenfunctions[n]
    #check function works
    @assert testcuda(f,x) > 0
    ∇testcuda(x) = gradient((x)->testcuda(f,x), x)
    try
      ∇testcuda(x) = gradient((x)->testcuda(f,x), x)
      grad = ∇testcuda(x)
      @info "passed with val $grad"
    catch err
      @warn " function $n is broken ($err)"
    end
  end
end
mwe(Float64)

#below copy/paste from https://github.com/FluxML/Zygote.jl/issues/730
#=function mweoriginal(N=5)
  x = CUDA.rand(N)

  M=CUDA.rand(N,N^2)

  brokenfunctions = (
    shouldwork=()->CUDA.sum(M),
    ones=()->CUDA.ones(N),
    zeros=()->CUDA.zeros(N),
    fill=()->CUDA.fill(2f0, N),
    cumprod=()->CUDA.cumprod(M, dims=2),
    device2host=()-> CUDA.rand(N,N) |> Matrix |>sum,
    reduce=()->CUDA.reduce(+,M),
    reduce2=()->CUDA.reduce(hcat,[M,M])
    )

  testcuda(f,x) = CUDA.sum(f() .+ x)
  for n ∈ propertynames(brokenfunctions)
    print("testing $n: ")
    f = brokenfunctions[n]
    #check function works
    @assert testcuda(f,x) > 0
    try
      ∇testcuda(x) = gradient((x)->testcuda(f,x), x)
      grad = ∇testcuda(x)
      @info "passed with val $grad"
    catch err
      @warn " function $n is broken ($err)"
    end
  end
end=#

#mweoriginal()



#=using Revise, Optim
mutable struct AdditionalState
  fevalctr
end

function f(x, Ψ)
  Ψ.fevalctr += 1

  return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

function g!(G, x)
  G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
  G[2] = 200.0 * (x[2] - x[1]^2)
end

function optimcallback(opt,Ψ)

  @info "fevalctr: $(Ψ.fevalctr)"
  return Ψ.fevalctr > 3
end

function mwe()
  Ψ=AdditionalState(0)

  #prep for the optimization and set parameters
  obj(x)=f(x,Ψ)
  callback(opt) = optimcallback(opt,Ψ)
  optimoptions = Optim.Options(iterations = 10, callback=callback)
  x = zeros(2)
  lower = [-100.,-100.]
  upper = [Inf, Inf]

  res=optimize(obj, g!, lower, upper, x, Fminbox(LBFGS()), optimoptions)
  println(res)
end

mwe()=#

###CUDA tests of prodmat
#=using Revise, CUDA, Zygote, LinearAlgebra, BenchmarkTools, Flux, Transducers, BangBang

prodmatbaseline(M::CuMatrix) = CUDA.hcat((i->CUDA.prod(M[:,1:i], dims=2)).(1:size(M,2))...)
prodmatbaseline(M::Matrix) = hcat((i->prod(M[:,1:i], dims=2)).(1:size(M,2))...)

function prodmattrans(M::TM,
  ::Type{TV} = CuVector{eltype(TM)},
  ::Type{T}=eltype(TM)) where {TM<:CuMatrix,TV<:CuVector, T}

  Nrow,Ncol = size(M)
  cumprodcol(c1)::TV = c1
  cumprodcol(c1,c2)::TV = c1 .* c2

  return foldl(hcat, Scan(cumprodcol), eachcol(M))
  #return reshape(eachcol(M) |> Scan(cumprodcol) |> Cat() |> collect, Nrow, Ncol)
end

function prodmat(M::TM,
  ::Type{T}=eltype(TM)) where {TM<:CuMatrix,T}

  Nrow,Ncol = size(M)


  cumprodcol(c1) = c1
  cumprodcol(c1,c2) = c1 .* c2

  cols = eachcol(M) |> collect
  accumcols = Iterators.accumulate(cumprodcol,cols) |> TCat(1)
  return foldl(hcat, Iterators.accumulate(cumprodcol,cols))
end


function mwe(N=10)

  #verify functions work
  M=CUDA.rand(Float64, N,N^2)
  Mcpu = M |> Matrix

  brokenfunctions = (
    #cumprod=(M)->CUDA.cumprod(M, dims=2),
    #prodmattrans=(M)->prodmattrans(M),
    #prodmat=(M)->prodmat(M),
    prodmatbaseline=(M)->prodmatbaseline(M),
    )

  testcuda(f,M) = CUDA.sum(f(M))
  for n ∈ propertynames(brokenfunctions)
    print("testing $n performance: ")
    f = brokenfunctions[n]
    #check function works
    testval = testcuda(f,M)
    verval = sum(cumprod(M, dims=2))
    @assert (testval ≈ verval) || (testval - N ≈ verval)
    try
      ∇testcuda() = gradient(()->testcuda(f,M), Flux.params(M))
      grad = @btime $∇testcuda()
      @btime $testcuda($f,$M)
      @info "passed with val $(sum(grad[M]))"
    catch err
      @warn " function $n is broken ($err)- performance of call-"
      @btime $testcuda($f,$M)
    end
  end
end

mwe()=#

###regression tests
#=using Revise, CUDA, Zygote, LinearAlgebra, BenchmarkTools, Flux, Transducers, BangBang


reg(X::AbstractMatrix{T},y::Vector{T}, ::Val{:glm}) where T = coef(lm(X, y))
reg(X::AbstractMatrix{T},y::AbstractArray{T}, ::Val{:cholesky}) where T = cholesky!(X'*X)\ (X' * y)
reg(X::AbstractMatrix{T},y::Vector{T}, ::Val{:qr}) where T = qr(X)\y
reg(X::AbstractMatrix{T},y::Vector{T}, ::Val{:svd}) where T = svd(X)\y
reg(X::CuMatrix{T}, y::CuArray{T}, ::Val{:cholesky}) where T = hyperreg(X,y) |> CUDA.vec
reg(X::AbstractMatrix{T}, y::AbstractArray{T}, ::Val{:choleskyzygote}) where T = cholesky(X'*X)\ (X' * y)
reg(X::CuMatrix{T}, y::CuVector{T}, ::Val{:choleskyzygote}) where T = CUDA.cholesky(X'*X)\ (X' * y)
reg(X::AbstractMatrix{T}, y::AbstractArray{T}, ::Val{:pinv}) where T = pinv(X'*X)*(X'*y)
function reg(X::CuMatrix{T}, y::CuArray{T}, ::Val{:pinv}) where T
  #@warn "No reasonable pinv function for CuArrays- falling back to cholesky"
  return reg(X, y, Val{:choleskyzygote}())
end

function mwe(N=10)

  #verify functions work
  X=CUDA.rand(Float64, N*10^3, N)
  y=CUDA.rand(Float64, size(X,1))
  Xcpu = X |> Matrix
  ycpu = y |> Vector

  brokenfunctions = (
    #cumprod=(M)->CUDA.cumprod(M, dims=2),
    reg=(X)->reg(X,y, Val{:choleskyzygote}()),
    #prodmat=(M)->prodmat(M),
    )

  testcuda(f,X) = CUDA.sum(f(X))
  for n ∈ propertynames(brokenfunctions)
    print("testing $n performance: ")
    f = brokenfunctions[n]
    #check function works
    testval = testcuda(f,X)
    verval = (cholesky(Xcpu'*Xcpu)\(Xcpu'*ycpu)) |> sum
    @assert testval ≈ verval
    try
      ∇testcuda() = gradient(()->testcuda(f,X), Flux.params(X))
      grad = @btime $∇testcuda()
      @btime $testcuda($f,$X)
      @info "passed with val $(sum(grad[X]))"
    catch err
      @warn " function $n is broken ($err)- performance of call-"
      @btime $testcuda($f,$X)
    end
  end
end

mwe()=#

#=
###zygote cumprod cpu tests
using Revise, Zygote, LinearAlgebra, BenchmarkTools, Transducers, BangBang, LoopVectorization, CUDA

#prodmatbaseline(M::CuMatrix) = CUDA.hcat((i->CUDA.prod(M[:,1:i], dims=2)).(1:size(M,2))...)
prodmatbaseline(M::Matrix) = hcat((i->prod(M[:,1:i], dims=2)).(1:size(M,2))...)

function prodmattrans(M::TM,
  ::Type{TV} = Vector{eltype(TM)},
  ::Type{T}=eltype(TM)) where {TM<:Matrix,TV<:Vector, T}

  Nrow,Ncol = size(M)
  cumprodcol(c1)::TV = c1
  cumprodcol(c1,c2)::TV = c1 .* c2

  return foldl(hcat, Scan(cumprodcol), eachcol(M))
  #return reshape(eachcol(M) |> Scan(cumprodcol) |> Cat() |> collect, Nrow, Ncol)
end

function prodmat(M::TM,
  ::Type{T}=eltype(TM)) where {TM<:Matrix,T}

  Nrow,Ncol = size(M)
  cumprodcol(c1) = c1
  cumprodcol(c1,c2) = c1 .* c2

  cols = eachcol(M) |> collect
  return foldl(hcat, Iterators.accumulate(cumprodcol,cols))
end

function modprodmat(M::TM, threshold=0.1, ::Type{T}=eltype(TM)) where {TM <: Matrix, T}
  #Mt = M' |> Matrix

  modprod(x) = ifelse(abs(x) > threshold, x, threshold*sign(x))
  function modprod(x,y)
    xy = x*y
    ifelse(abs(xy) > threshold, xy, threshold*sign(xy))
  end

  modprodrow(r) = foldl(modprod, r, init=T(1.0))
  #prodpart(Mpart) = foldl(vcat, modprodrow.(eachrow(Mpart)))
  prodpart(Mpart) = map(r->modprodrow(Mpart[r,:]), 1:size(M,1))
  #cols = (i->prodpart(view(M, :, 1:i))).(1:size(M,2))
  cols = map(i->prodpart(M[:,1:i]), 1:size(M,2))#(i->prodpart(view(M, :, 1:i))).(1:size(M,2))


  return hcat(cols...)
end


#NOTE NOTE NOTE - integrate the below?
#conditiongrowth(g, ::Val{:nozero}, nozerotol=PARAM[:iternozerotol]) = (
#  g .+ nozerotol .* sign.(g) .* exp.(-g .* g))
function boundedcumprod(M::TM, scale::TV = rand(size(M,1)), ::Type{T}=eltype(TM);
   threshold::T=0.1, init=ones(T,size(M,1))) where {TM<:Matrix,TV<:Vector, T}

  #NOTE: more idiomatically, this is x*y > threshold ? x*y : threshold*sign(x*y)
  thresholds = threshold .* scale
  function checkthreshold(x, y, threshold)
    xy = x * y
    xy + threshold * sign(xy) * exp(-xy^2/(2*threshold))
    #ifelse(abs(xy) > threshold, xy, threshold*sign(xy))
  end

  cumprodcol(c1,c2) = vmap(checkthreshold, c1, c2, thresholds)

  cols = eachcol(M) |> collect
  return foldl(hcat, Iterators.accumulate(cumprodcol,cols, init=init))
end
boundedcumprod(M::TM, ::Type{T}=eltype(TM)) where {TM<:CuMatrix,T} = (
  boundedcumprod(M |> Matrix{T}) |> TM)


function mwe(N=25)

  #verify functions work
  M=rand(Float64, N,N^2)# |> CuMatrix{Float64}

  brokenfunctions = (
    #cumprod=(M)->cumprod(M, dims=2),
    #prodmattrans=(M)->prodmattrans(M),
    #modprodmat=(M)->modprodmat(M),
    boundedcumprod=(M)->boundedcumprod(M),
    #prodmatbaseline=(M)->prodmatbaseline(M),
    )

  #(boundedcumprod(M) ≈ modprodmat(M)) || error("boundedcumprod(M) ≈/ modprodmat(M)
  #  Δ: $(mean(abs.(boundedcumprod(M) .- modprodmat(M)), dims=1))
  #  boundedcumprod: $(boundedcumprod(M)[:,1])")
  test(f,M) = sum(f(M))
  for n ∈ propertynames(brokenfunctions)
    sleep(0.1)
    print("testing $n performance: ")
    f = brokenfunctions[n]
    #check function works
    testval = test(f,M)
    #verval = sum(modprodmat(M))
    #(testval ≈ verval) || N > 10 || error("testval ≈/ verval: testval=$testval, verval=$verval")
    try
      ∇test(M) = gradient((M)->test(f,M), M)
      grad = @btime $∇test($M)
      #@btime $test($f,$M)
      println("passed with val $(sum(grad[1]))")
    catch err
      @warn " function $n is broken ($err)- performance of call-"
      @btime $test($f,$M)
    end
  end
end

mwe()
=#

###mcmc tests
#=
using Revise, Distributions, Random, Plots, UnicodePlots, LinearAlgebra

function testacceptreject(N=10^7, ; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()


  σstar = 1/(1/σ1^2+1/σ2^2)
  μstar = σstar*(μ1/σ1^2+μ2/σ2^2)
  @info "μstar: $μstar, σstar^0.5: $(σstar^0.5)"

  #pdf distribution (unkown)
  p(x) = pdf(Normal(μstar, σstar^0.5),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1-x)^2/(2σ1^2)-(μ2-x)^2/(2σ2^2))

  q(x) = pdf(Normal(μ2, σ2),x)
  Mq(x) = exp(-(μ2-x)^2/(2σ2^2))

  horiz = rand(Normal(μ2, σ2), N)
  vertu = rand(Uniform(), N)

  propps = (propp).(horiz)
  Mqs = Mq.(horiz)
  @assert all(Mqs .≥ propps)

  ratios = propps ./ Mqs
  keeps = vertu .≤ ratios

  @info "We drew $(sum(keeps)) samples out of a maximum $N"

  ardraws = horiz[keeps]
  @info "Drawn mean: $(mean(ardraws)), stdev = $(std(ardraws))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end

function testacceptreject2(N=10^7, K=7; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σ2)^0.5: $(diag(Σ2).^0.5)
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-x)'*Λ2*(μ2.-x)/2)

  q(x) = pdf(MvNormal(μ2, Σ2),x)
  Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)

  horiz = rand(MvNormal(fill(μ2,K), Σ2), N)' |> Matrix
  vertu = rand(Uniform(), N)

  t0 = time()
    propps = (propp).(eachrow(horiz))
    Mqs = Mq.(eachrow(horiz))
    @assert all(Mqs .≥ propps)

    ratios = propps ./ Mqs
    keeps = (vertu .≤ ratios) |> vec
  Δt = time()-t0

  Nkeeps = sum(keeps)

  @info "We drew $(Nkeeps) samples out of a maximum $N ($(N/Nkeeps) draws per sample,"*
    "$(Δt/Nkeeps) seconds)"

  ardraws = horiz[keeps, :]
  @info "Drawn mean: $(mean(ardraws, dims=1)), stdev = $(std(ardraws, dims=1))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end

#TODO- finish the below. Use a counter to compute rejected draws
#fill a matrix with accepted draws- compare performance with previous
function testacceptreject3(N=10^6, K=7; μ1=1.0, μ2=-2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"
    #diag(Σ2)^0.5: $(diag(Σ2).^0.5)

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-x)'*Λ2*(μ2.-x)/2)

  q(x) = pdf(MvNormal(μ2, Σ2),x)
  Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)

  σk=1 ./ diag(Λ2).^0.5 #Schur complements
  Σk=[Σ2[k:k, 1:K .≠ k] for k in 1:K] #given Σ, selects row vector Σ_k,-k
  function drawqcomponent(x,k)
    #@info "2nd part: $((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))"
    μ = μ2 - ((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))[]
    rand(Normal(μ,σk[k]))
  end

  horiz = rand(MvNormal(fill(μ2,K), Σ2), N)' |> Matrix
  vertu = rand(Uniform(), N)

  ctr::Int = 0
  draws::Matrix{Float64} = Matrix{Float64}(undef, N, K)
  seed::Vector{Float64} = ones(K)
  x::Vector{Float64} = seed |> deepcopy
  t0 = time()
  for n ∈ 1:N
    for k ∈ 1:K #for each x component
      while true
        ctr += 1
        x[k] = drawqcomponent(x,k)
        cut = propp(x)/Mq(x)
        (rand(Uniform())<cut) && break #accept if true
      end
    end
    draws[n,:] .= x
  end

  Δt = time()-t0

  @info "We drew $(N) samples after $ctr attempts ($(round(ctr/N, digits=2)) draws per sample).
    Time per sample = $(round(Δt/N, digits=6))"

  @info "Drawn mean: $(mean(draws, dims=1)), stdev = $(std(draws, dims=1))"
  drawupperhalf = draws[(end÷2):end, :]
  @info "Drawn 2nd half mean: $(mean(drawupperhalf, dims=1)), 2nd half stdev = $(std(drawupperhalf, dims=1))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end

function testMHC(N=10^6, K=7; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"
    #diag(Σ2)^0.5: $(diag(Σ2).^0.5)

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-x)'*Λ2*(μ2.-x)/2)

  q(x) = pdf(MvNormal(μ2, Σ2),x)
  Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)

  σk=1 ./ diag(Λ2).^0.5 #Schur complements
  Σk=[Σ2[k:k, 1:K .≠ k] for k in 1:K] #given Σ, selects row vector Σ_k,-k
  function drawqcomponent(x,k)
    #@info "2nd part: $((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))"
    μ = μ2 - ((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))[]
    d=Normal(μ,σk[k])
    (x)->pdf(d,x), rand(d)
  end

  horiz = rand(MvNormal(fill(μ2,K), Σ2), N)' |> Matrix
  vertu = rand(Uniform(), N)

  ctr::Int = 0
  draws::Matrix{Float64} = Matrix{Float64}(undef, N, K)
  #seed::Vector{Float64} = ones(K)
  x::Vector{Float64} = ones(K)
  xcandidate = x |> deepcopy
  acceptancedraws = rand(Uniform(),N,K)
  t0 = time()
  for n ∈ 1:N
    for k ∈ 1:K #for each x component
      Mf,xcandidate[k] = drawqcomponent(x,k)
      #cut = (propp(xcandidate)/Mq(xcandidate))/(propp(x)/Mq(x)) #NOTE- this can be used as a test
      cut = (propp(xcandidate)/Mf(xcandidate[k]))/(propp(x)/Mf(x[k]))
      x[k] = ifelse(acceptancedraws[n,k]<cut, xcandidate[k], x[k])
      xcandidate[k] = x[k]
    end
    draws[n,:] .= x
  end

  Δt = time()-t0

  @info "We drew $(N) samples after $ctr attempts ($(round(ctr/N, digits=2)) draws per sample).
    Time per sample = $(round(Δt/N, digits=6))"

  @info "Drawn mean: $(mean(draws, dims=1)), stdev = $(std(draws, dims=1))"
  drawupperhalf = draws[(end÷2):end, :]
  @info "Drawn 2nd half mean: $(mean(drawupperhalf, dims=1)), 2nd half stdev = $(std(drawupperhalf, dims=1))"
  @info "uniques: $((c->length(unique(c))).(eachcol(draws)))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end

function testMH(N=10^6, K=7; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"
    #diag(Σ2)^0.5: $(diag(Σ2).^0.5)

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-x)'*Λ2*(μ2.-x)/2)

  q(x) = pdf(MvNormal(μ2, Σ2),x)
  Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)

  ctr::Int = 0
  draws::Matrix{Float64} = Matrix{Float64}(undef, N, K)
  drawqs::Matrix{Float64} = rand(MvNormal(fill(μ2,K), Σ2), N)' |> Matrix
  #seed::Vector{Float64} = ones(K)
  x::Vector{Float64} = ones(K)
  xcandidate = x |> deepcopy
  acceptancedraws = rand(Uniform(),N)
  t0 = time()
  for n ∈ 1:N
    xcandidate .= drawqs[n,:]
    cut = (propp(xcandidate)/Mq(xcandidate))/(propp(x)/Mq(x))
    x .= ifelse(acceptancedraws[n]<cut, xcandidate, x)
    draws[n,:] .= x
  end

  Δt = time()-t0

  @info "We drew $(N) samples after $ctr attempts ($(round(ctr/N, digits=2)) draws per sample).
    Time per sample = $(round(Δt/N, digits=6))"

  @info "Drawn mean: $(mean(draws, dims=1)), stdev = $(std(draws, dims=1))"
  drawupperhalf = draws[(end÷2):end, :]
  @info "Drawn 2nd half mean: $(mean(drawupperhalf, dims=1)), 2nd half stdev = $(std(drawupperhalf, dims=1))"
  @info "uniques: $((c->length(unique(c))).(eachcol(draws)))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end


function testacceptreject2abs(N=10^7, K=7; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σ2)^0.5: $(diag(Σ2).^0.5)
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-abs.(x))'*Λ2*(μ2.-abs.(x))/2)

  #q(x) = pdf(MvNormal(μ2, Σ2),x)
  #Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)
  q(x) = rand(Uniform(μ2-5*σ2, μ2+5*σ2), K)
  Mq(x)=1.0

  horiz = rand(Uniform(μ2-5*σ2, μ2+5*σ2), K, N)' |> Matrix
  vertu = rand(Uniform(), N)

  t0 = time()
    propps = (propp).(eachrow(horiz))
    Mqs = Mq.(eachrow(horiz))
    @assert all(Mqs .≥ propps)

    ratios = propps ./ Mqs
    keeps = (vertu .≤ ratios) |> vec
  Δt = time()-t0

  Nkeeps = sum(keeps)

  @info "We drew $(Nkeeps) samples out of a maximum $N ($(N/Nkeeps) draws per sample,"*
    "$(Δt/Nkeeps) seconds)"

  ardraws = horiz[keeps, :]
  @info "Drawn mean: $(mean(ardraws, dims=1)), stdev = $(std(ardraws, dims=1))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end

function testMHCabs(N=10^6, K=7; μ1=-1.0, μ2=2.0, σ1 = 2.0, σ2=3.0)
  unicodeplots()

  X = rand(Normal(0,σ2),N,K)# .+ rand(N)
  Σ2 = cov(X)
  #@info ""
  Λ2 = inv(Σ2)
  Σstar = (diagm(fill(1/σ1^2,K)) .+ Λ2)\I
  μstar = Σstar*(μ1/σ1^2 .+ Λ2*fill(μ2,K))
  @info "μstar: $μstar
    diag(Σstar) ^0.5: $(diag(Σstar).^0.5)"
    #diag(Σ2)^0.5: $(diag(Σ2).^0.5)

  #pdf distribution (unkown)
  p(x) = pdf(MvNormal(μstar, Σstar),x)
  #propp - known to a constant of proportionality, unsamplable
  propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-abs.(x))'*Λ2*(μ2.-abs.(x))/2)
  #propp(x) = exp(-(μ1.-x)'*(μ1.-x)/(2σ1^2)-(μ2.-x)'*Λ2*(μ2.-x)/2)

  q(x) = pdf(MvNormal(μ2, Σ2),x)
  Mq(x) = exp(-(μ2 .- x)'*Λ2*(μ2 .- x)/2)

  σk=1 ./ diag(Λ2).^0.5 #Schur complements
  Σk=[Σ2[k:k, 1:K .≠ k] for k in 1:K] #given Σ, selects row vector Σ_k,-k
  ΣMkinv = (k->Σ2[1:K .≠ k, 1:K .≠ k]\I).(1:K)
  function drawqcomponent(x,k)
    #@info "2nd part: $((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))"
    #μ = μ2 - ((Σk[k]./Σ2[k,k]) * (x[1:K .≠ k] .- μ2))[]
    μ = μ2 - ((Σk[k]*ΣMkinv[k]) * (x[1:K .≠ k] .- μ2))[]
    #d=Normal(μ,σk[k])
    d = Normal(-40.,3.)
    (x)->pdf(d,x), rand(d)
  end

  horiz = rand(MvNormal(fill(μ2,K), Σ2), N)' |> Matrix
  vertu = rand(Uniform(), N)

  ctr::Int = 0
  draws::Matrix{Float64} = Matrix{Float64}(undef, N, K)
  #seed::Vector{Float64} = ones(K)
  x::Vector{Float64} = ones(K) .* -39.0
  xcandidate = x |> deepcopy
  acceptancedraws = rand(Uniform(),N,K)
  t0 = time()
  for n ∈ 1:N
    for k ∈ 1:K #for each x component
      Mf,xcandidate[k] = drawqcomponent(x,k)
      #cut = (propp(xcandidate)/Mq(xcandidate))/(propp(x)/Mq(x)) #NOTE- this can be used as a test
      cut = (propp(xcandidate)/Mf(xcandidate[k]))/(propp(x)/Mf(x[k]))
      x[k] = ifelse(acceptancedraws[n,k]<cut, xcandidate[k], x[k])
      xcandidate[k] = x[k]
    end
    draws[n,:] .= x
  end

  Δt = time()-t0

  @info "We drew $(N) samples after $ctr attempts ($(round(ctr/N, digits=2)) draws per sample).
    Time per sample = $(round(Δt/N, digits=6))"

  @info "Drawn mean: $(mean(draws, dims=1)), stdev = $(std(draws, dims=1))"
  drawupperhalf = draws[(end÷2):end, :]
  @info "Drawn 2nd half mean: $(mean(drawupperhalf, dims=1)), 2nd half stdev = $(std(drawupperhalf, dims=1))"
  @info "uniques: $((c->length(unique(c))).(eachcol(draws)))"
  #Plots.histogram(ardraws, normalize=true, bins=-10.:0.01:10.0) |> display
  #Plots.histogram(rand(Normal(μstar, σstar), sum(keeps)), normalize=true, bins=-10.:0.01:10.0) |> display
end


#NOTE: OK
#testacceptreject(10^6) #1D only
#testacceptreject2(10^6,1) #slow but looks good- might be a nice check for few vars
#testMHC(10^5,3) #trust this a little more, looks great!

testacceptreject2abs(10^7,3)
testMHCabs(10^7,3)

#rejected
#testacceptreject3(10^6,10) #slow but looks good- might be a nice check for few vars
#testMH(10^4,10) #slow anbd not very good
=#


###regression zygote
#=using Zygote, CUDA, LinearAlgebra, Distributions, StatsBase, Statistics, GLM
CUDA.allowscalar(false)

zygotedevice2host(M::TM, ::Type{T}=eltype(TM)
  ) where {T<:Real, TM<:CuArray} = copyto!(Zygote.Buffer(Array{T}(undef, size(M))), M) |> copy

function regtest(::Type{TM}, ::Type{TV}, ::Type{T}; N=100, K=5) where {TM, TV, T}

  y = rand(N^2) |> TV
  A = rand(N,K)
  idx = rand([1:N...;],N*N)



  function runreg(A, f)
    X = A[idx, :] |> TM
    b=f(X,y)
    return sum((X*b .- y).^2)
  end

  regchol(X,y) = cholesky(X'*X)\(X'*y)
  regsvd(X,y) = svd(X)\y
  regqr(X,y) = qr(X'*X)\(X'*y)
  regpinv(X,y) = pinv(X'*X)*(X'*y)
  reglm(X,y) = coef(lm(X,y))

  regfunctions = (
    shouldwork=(A)->runreg(A, regchol),
    #svd=(A)->runreg(A, regsvd),
    #qr=(A)->runreg(A, regqr),
    pinv=(A)->runreg(A, regpinv),
    lm=(A)->runreg(A, reglm),
    )

  for n ∈ propertynames(regfunctions)
    print("\ntesting $n: ")
    f = regfunctions[n]
    #check function works


    try
      @assert f(A) > 0
      ∇test(A) = gradient((A)->f(A), A)
      grad = ∇test(A)[1]

      @info("passed with $(sum(grad .≠ 0.0)) non-zero gradient values")
    catch err
      errmessage = "$err"
      if length(errmessage) > 1000
        errmessage=errmessage[1:1000]
      end
      @warn " function $n is broken ($errmessage)"
    end
  end
end

#regtest(CuMatrix{Float64}, CuVector{Float64}, Float64, N=10, K=5)
regtest(Matrix{Float64}, Vector{Float64}, Float64, N=10, K=5)=#
