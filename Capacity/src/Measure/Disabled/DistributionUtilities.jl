############related to gamma function
#see https://math.stackexchange.com/questions/1441753/approximating-the-digamma-function
#and https://en.wikipedia.org/wiki/Digamma_function
function cudaapproxdigamma(x::Real)
  adj = 0.0f0
  ψ = x
  #the polynomial is only accurate for large x
  for i ∈ 1:3
    adj -= 1f0/ψ
    ψ += 1f0
  end

  ψ = (CUDA.log(ψ) - 1f0 / 2f0 / ψ - 1f0 / 12f0 / (ψ * ψ) +
    1f0/120f0 * CUDA.pow(ψ, -4) - 1f0/252f0 * CUDA.pow(ψ, -6) +
    1f0/240f0 * CUDA.pow(ψ, -8) - 5f0/660f0 * CUDA.pow(ψ, -10) +
    691f0/32760f0 * CUDA.pow(ψ, -12) - 1f0/12f0 * CUDA.pow(ψ, -14)) + adj

  return ψ
end


#runs a customized gamma log pdf that uses a gpu array when available
#adapted from StatsFuns
#https://github.com/JuliaStats/StatsFuns.jl/blob/master/src/distrs/gamma.jl
#and CUDA code
#GPU vectorized log gamma pdf
function vecgammalogpdf(k::TV, θ::TV, x::TV)::TV where TV<:CUDA.CuVector
  Γs_part = cudagammalogpdf_part.(k,θ,x)

  return Γs_part .- cudaveclgamma(k)
end

#does everything except compute the log gamma function
#adapted from StatsFuns
#https://github.com/JuliaStats/StatsFuns.jl/blob/master/src/distrs/gamma.jl
cudagammalogpdf_part(k,θ,x) = - k * CUDA.log(θ) + (k - 1f0) * CUDA.log(x) - x / θ

#vectorized log gamma function
cudaveclgamma(k) = (CUDA.lgamma).(k)

#define the adjoint
Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
Zygote.refresh()


#lightly optimized cpu version for testing
function vecgammalogpdf(k::Any, θ::Any, x::TV)::TV where TV<:AbstractVector

  Γs::TV = -(SpecialFunctions.loggamma).(k) .- k .* log.(θ) .+ (k .- 1f0) .* log.(x) .- x ./ θ

  return Γs
end

slowvecgammalogpdf(k, θ, x) = ((kᵢ, θᵢ, xᵢ)->logpdf(Gamma(kᵢ,θᵢ), xᵢ)).(k, θ, x)

#test function for vecgammalogpdf
#test it
function testgammagrad(N::Int)
  τ = abs.(randn(Float32, N)) .+ .001f0
  absμ = abs.(randn(Float32, N)) .+ .001f0

  θ = 1.0f0 ./ τ
  k = absμ .* τ
  v = rand(Float32, N)

  cuθ = θ |> gpu
  cuk = k |> gpu
  cuv = v |> gpu

  @info "testing function values (N=$N) cpu: $(sum(vecgammalogpdf(k, θ, v))) " *
    "gpu: $(sum(vecgammalogpdf(cuk, cuθ, cuv)))"
  @info "Additional check of value: $(sum(slowvecgammalogpdf(k, θ, v)))"

  @info "cpu version gradient N=$N"
  if N ≤ 10
    gs = gradient((k,θ)->sum(vecgammalogpdf(k, θ, v)),k,θ)
    @info "gs: $gs"
  else
    gs = @btime gradient((k,θ)->sum(vecgammalogpdf(k, θ, $v)),$k,$θ)
  end

  @info "gpu version gradient N=$N"
  if N ≤ 10
    cugs = gradient((cuk,cuθ)->sum(vecgammalogpdf(cuk, cuθ, cuv)),cuk,cuθ)
    @info "cugs: $cugs"
  else
    cugs = @btime gradient((cuk,cuθ)->sum(vecgammalogpdf(cuk, cuθ, $cuv)),$cuk,$cuθ)
  end
end


###log normal dpf
#sortof a cuda kernal
using Revise, BenchmarkTools, Distributions, Zygote, Flux, CUDA, CUDA
function cudanormlogprop(μᵢ::Real, logτᵢ::Real, xᵢ::Real)
  μMxᵢ = μᵢ-xᵢ
  return -1f0 * μMxᵢ*μMxᵢ*CUDA.exp(logτᵢ) + logτᵢ
end

function vecnormlogprop(μ::TV, logτ::TV, x::TV)::TV where TV<:AbstractVector
  uMx = μ .- x
  return -1f0 .* uMx .* uMx .* exp.(logτ) .+ logτ
end


#Zygote.@adjoint cudaveclgamma(k) = cudaveclgamma(k), y->(y .* cudaapproxdigamma.(k), )
#Zygote.refresh()
function vecnormlogprop(μ::TV, logτ::TV, x::TV)::TV where TV<:CUDA.CuVector
  return cudanormlogprop.(μ,logτ,x)
end
cudanormlogpdffromprop(x::Real) = 0.5f0*x - 0.5f0*log(2f0*3.1415927f0)
slowvecnormlogpdf(μ, σ, x) = ((μᵢ, σᵢ, xᵢ)->logpdf(Normal(μᵢ,σᵢ), xᵢ)).(μ, σ, x)

function vecnormlogpdf(μ::TV, logτ::TV, x::TV)::TV where TV<:CUDA.CuVector
  return cudanormlogpdffromprop.(vecnormlogprop(μ, logτ, x))
end

function vecnormlogpdf(μ::TV, logτ::TV, x::TV)::TV where TV<:AbstractVector
  return 0.5f0 .* vecnormlogprop(μ, logτ, x) .- 0.5f0 .* log.(2f0 .* 3.1415927f0)
end

function testnormgrad(N::Int)
  τ = log.(abs.(randn(Float32, N)) .+ .001f0)
  μ = abs.(randn(Float32, N)) .+ .001f0

  σ = (1.0f0 ./ exp.(τ)).^0.5f0
  v = rand(Float32, N)

  cuμ = μ |> gpu
  cuτ = τ |> gpu
  cuσ = σ |> gpu
  cuv = v |> gpu

  @info "testing prop values (N=$N) cpu: $(sum(vecnormlogprop(μ, τ, v))) " *
    "gpu: $(sum(vecnormlogprop(cuμ, cuτ, cuv)))"
  @info "corrected cpu values: $(sum(vecnormlogpdf(μ, τ, v)))"
  @info "corrected gpu values: $(sum(vecnormlogpdf(cuμ, cuτ, cuv)))"
  @info "Additional check of value: $(sum(slowvecnormlogpdf(μ, σ, v)))"

  @info "cpu version gradient N=$N"
  if N ≤ 10
    gs = gradient((μ, τ)->sum(vecnormlogprop(μ, τ, v)),μ, τ)
    @info "gs: $gs"
  else
    gs = @btime gradient((μ,τ)->sum(vecnormlogprop(μ, τ, $v)),$μ,$τ)
  end

  @info "gpu version gradient N=$N"
  if N ≤ 10
    cugs = gradient((cuμ,cuτ)->sum(vecnormlogprop(cuμ, cuτ, cuv)),cuμ,cuτ)
    @info "cugs: $cugs"
  else
    cugs = @btime gradient((cuμ,cuτ)->sum(vecnormlogprop(cuμ, cuτ, $cuv)),$cuμ,$cuτ)
  end
end
