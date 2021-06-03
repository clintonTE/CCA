module Tutorial
using Revise #useful when revising code frequently

#these are all fairly well established packages
using Distributions, LinearAlgebra, Random, Revise, CUDA, GLM, BenchmarkTools, UnPack, CSV
sleep(0.1)


#a simple ols wrapper
#use parametric types to show high performance code
struct OLS{TX, Ty, Tcoef, T}
  X::TX
  y::Ty
  coef::Tcoef

  #inner construction method
  function OLS(X::TX, y::Ty, #main inputs
    ::Type{T}=eltype(Ty); #types for generic code- not necessary
      addintercept::Bool = true) where {TX,Ty, Tcoef ,T}
    K = size(X,2) + addintercept #set the number of coefficients
    N = size(X,1)
    @assert N ≥ K #checks like these are easy in Julia- if this is not true, it will throw an error

    coef = Vector{Union{Missing,T}}(undef, K) #preallocate
    if addintercept
      #augment a vector of ones to X
      some1s = ones(T, N) |> Ty
      Xout = [some1s X]
      @assert typeof(Xout) == TX
      return new{TX,Ty,typeof(coef), T}([X some1s],y,coef) #create the object
    else
      return new{TX,Ty,typeof(coef), T}(X,y,coef) #create the object
    end
  end
end


function simulatedata(N,K)
  ε = rand(Normal(), N)
  β = 1:K |> collect #true coefficient values
  X = rand(Normal(), N, K-1)
  y = X * β[1:(end-1)] .+ β[end] .+ ε #note the '.' implies a vectorized broadcast operation

  return (X=X,y=y,β=β)
end

function estimatebasic!(ols::OLS)
  ols.coef .= inv(ols.X'*ols.X)*(ols.X'*ols.y) #elementwise assign
  return nothing
end

rms(test,ver) = ((test .- ver).^2 |> mean)^0.5

#Basic regression
function benchbasic(N,K)
  d = simulatedata(N,K)
  ols = OLS(d.X,d.y)

  @info "Running regresion using basic estimation methedology"
  @btime estimatebasic!($ols)
  println("Coef error: $(rms(ols.coef, d.β))")
end

benchbasic(10^6,250)

#can we do better?
function estimatecholesky!(ols::OLS)
  ols.coef .= cholesky!(ols.X'*ols.X)\(ols.X'*ols.y) #elementwise assign
  return nothing
end

#cholesky
function benchcholesky(N,K)
  d = simulatedata(N,K)
  ols = OLS(d.X,d.y)

  @info "Running regresion estimates with cholesky methedology"
  @btime estimatecholesky!($ols)
  println("Coef error: $(rms(ols.coef, d.β))")
end

benchcholesky(10^6,250)
#use the GPU!
function estimatecuda!(ols::OLS{TX,<:Any, <:Any, <:Any}) where TX
  ols.coef .= cholesky!(ols.X'*ols.X)\(ols.X'*ols.y) #elementwise assign
  return nothing
end

function benchcuda(N,K)
  d = simulatedata(N,K)

  ols = OLS(d.X |> CuMatrix{Float64},d.y |> CuVector{Float64})

  @info "Running regresion estimates with CUDA methedology"
  @btime estimatecuda!($ols)
  println("Coef error: $(rms(ols.coef, d.β))")
end

benchcuda(10^6,250)

end #module end
