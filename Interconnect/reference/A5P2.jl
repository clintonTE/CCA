module A5P2


#Use this to print
#=
weave(Pkg.dir("$(pwd())\\A5P2.jl"),
  informat="script",
  out_path = "$(pwd())\\A5P2.html",
  doctype = "md2html")
=#
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt,
  HCubature, StaticArrays, CTNumerical, Optim


const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "m/d/y"
const SP500_NAME = "GSPC"
const TBILL_NAME = "t-bill"
const TERM_NAME = "termStructure"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const MAX_YEAR = 2016
const IN_SAMPLE_CUTOFF_YEAR = 1975

const DEFAULT_NUM_SIMS = 10000
const PlotContainer = Union{Plot,Gadfly.Compose.Context}
const CSRegSpec = Tuple{Symbol, Symbol, Int}
const CPRegSpec = Tuple{Symbol, Int}

const MIN_VALS_FOR_EST = 10
const Ι = Complex(0.0,1.0)

abstract type AbstractΘState end
abstract type AbstractΦFunction end

#convenience function for integration in the pure julia implementation
function I1D(f::Function, a::Float64, b::Float64)

  wrapFunc(x::SVector{1,Float64})::SVector{1,Float64} = SVector{1,Float64}(f(x[1]))
  ans::SVector{1,Float64}, err::Float64 =
    hcubature(wrapFunc, SVector{1,Float64}(a), SVector{1,Float64}(b))

  return ans[1]
end


function preProcessSP500()::Void
  dateFormat::String = DATE_FORMAT_STR
  sp500DF::DataFrame = readtable("$DATA_PATH\\$SP500_NAME.csv")

  #fix the dates
  rename!(sp500DF, :Date, :DATE)
  sp500DF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(sp500DF[:,:DATE])::DataVector{Date}
  sort!(sp500DF, cols = [:DATE])
  sp500DF[:,:returns] = similar(sp500DF[:,:Adj_Close])
  sp500DF[2:end,:returns] = log.(sp500DF[2:end,:Adj_Close] ./  sp500DF[1:(end-1),:Adj_Close])

  sp500DF[:,:dayOfQuarter] = zeros(Int,size(sp500DF,1))
  for i::Int ∈ 2:size(sp500DF,1)
    if Dates.firstdayofquarter(sp500DF[i,:DATE]) ≠ Dates.firstdayofquarter(sp500DF[i-1,:DATE])
      sp500DF[i,:dayOfQuarter] = 1
    else
      sp500DF[i,:dayOfQuarter] = sp500DF[i-1,:dayOfQuarter] + 1
    end
  end
  #showcols(sp500DF)

  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls", "w")
  serialize(stream, sp500DF)
  close(stream)

  return nothing
end

struct ΘBSM <: AbstractΘState
  S₀::Float64
  τ::Float64
  X::Float64
  r::Float64
  q::Float64
  σ::Float64
  σ²::Float64

  Φ::Function
end

#constructors for the BSM parameter set
ΘBSM(S₀::Float64, τ::Float64, X::Float64, r::Float64, q::Float64, σ::Float64)::ΘBSM =
  ΘBSM(S₀, τ, X, r, q, σ, σ*σ, ΦBSM)

ΘBSM(;S₀::Float64=0, τ::Float64=0, X::Float64=0, r::Float64=0, q::Float64=0, σ::Float64=0)::ΘBSM =
  ΘBSM(S₀, τ, X, r, q, σ)

#characteristic function of the black scholes
ΦBSM(Θ::ΘBSM, ω::T where T<:Union{Complex{Float64}, Float64})::Complex{Float64}  =
  exp(im*(ω * log(Θ.S₀) + ω * Θ.τ * (Θ.r - Θ.q - 0.5*Θ.σ²)) - 0.5 * Θ.σ²  *ω*ω * Θ.τ)

#parameter set for the BATES model
struct ΘBates <: AbstractΘState
  S₀::Float64
  τ::Float64
  X::Float64
  r::Float64
  q::Float64

  σ::Float64
  ρ::Float64
  v₀::Float64

  λ::Float64
  κ::Float64
  θ::Float64
  μⱼ::Float64
  σⱼ::Float64

  σ²::Float64
  σⱼ²::Float64
  Φ::Function
end

#constructors for the Bates parameter set
ΘBates(S₀::Float64, τ::Float64, X::Float64, r::Float64, q::Float64,
  σ::Float64, ρ::Float64, v₀::Float64, λ::Float64,
 κ::Float64, θ::Float64, μⱼ::Float64, σⱼ::Float64)::ΘBates =
  ΘBates(S₀, τ, X, r, q, σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ, σ*σ, σⱼ*σⱼ, ΦBates)

ΘBates(;S₀::Float64=0, τ::Float64=0, X::Float64=0, r::Float64=0, q::Float64=0,
  σ::Float64=0, ρ::Float64=0, v₀::Float64=0, λ::Float64=0,
  κ::Float64=0, θ::Float64=0, μⱼ::Float64=0, σⱼ::Float64=0)::ΘBates =
  ΘBates(S₀, τ, X, r, q, σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ)

#this is the characteristic function for Bates
function ΦBates(Θ::ΘBates, ω::T)::Complex{Float64} where T<:Union{Complex{Float64}, Float64}
  ι::Complex{Float64} = Ι

  #caclulate helper variables
  ρσiω::Complex{Float64} = Θ.ρ*Θ.σ*ι*ω
  d::Complex{Float64} = √((ρσiω-Θ.κ)^2+Θ.σ²*(ι*ω+ω*ω))
  g::Complex{Float64} = (Θ.κ-ρσiω-d)/(Θ.κ-ρσiω+d)

  #calculate next level of helper variables
  A::Complex{Float64} = ι * ω * log(Θ.S₀) + ι*ω*(Θ.r-Θ.q) * Θ.τ
  B::Complex{Float64} = Θ.θ*Θ.κ/Θ.σ² * ((Θ.κ-ρσiω-d)*Θ.τ-2.0*log((1.0-g*e^(-d*Θ.τ))/(1.0-g)))
  C::Complex{Float64} = Θ.v₀/Θ.σ² * (Θ.κ-ρσiω-d)*(1.0-e^(-d*Θ.τ)) / (1.0-g*e^(-d*Θ.τ))
  D::Complex{Float64} = -Θ.λ*Θ.μⱼ*ι*ω*Θ.τ + Θ.λ*Θ.τ*((1.0+Θ.μⱼ)^(ι*ω)*e^(0.5*Θ.σⱼ²*ι*ω*(ι*ω-1.0))-1.0)

  return e^(A+B+C+D)
end

function ΦtestBates(Θ::ΘBates, ω::T)::Complex{Float64} where T<:Union{Complex{Float64}, Float64}
  ι::Complex{Float64} = Ι

  #caclulate helper variables
  ρσiω::Complex{Float64} = Θ.ρ*Θ.σ*ι*ω
  d::Complex{Float64} = √((ρσiω-Θ.κ)^2+Θ.σ²*(ι*ω+ω*ω))
  g::Complex{Float64} = (Θ.κ-ρσiω-d)/(Θ.κ-ρσiω+d)

  #calculate next level of helper variables
  A::Complex{Float64} = ι * ω * log(Θ.S₀) + ι*ω*(Θ.r-Θ.q) * Θ.τ
  B::Complex{Float64} = Θ.θ*Θ.κ/Θ.σ² * ((Θ.κ-ρσiω-d)*Θ.τ-2.0*log((1.0-g*e^(-d*Θ.τ))/(1.0-g)))
  C::Complex{Float64} = Θ.v₀/Θ.σ² * (Θ.κ-ρσiω-d)*(1.0-e^(-d*Θ.τ)) / (1.0-g*e^(-d*Θ.τ))
  D::Complex{Float64} = -Θ.λ*Θ.μⱼ*ι*ω*Θ.τ + Θ.λ*Θ.τ*((1.0+Θ.μⱼ)^(ι*ω)*e^(0.5*Θ.σⱼ²*ι*ω*(ι*ω-1.0))-1.0)

  println("""Bates Intermediate Results:
     A: $A
     B: $B
     C: $C
     D: $D
     d: $d
     g: $g
     Answer: $(e^(A+B+C+D))""")

  return e^(A+B+C+D)
end

struct ΘHeston <: AbstractΘState
  S₀::Float64
  τ::Float64
  X::Float64
  r::Float64
  q::Float64

  σ::Float64
  ρ::Float64
  v₀::Float64

  κ::Float64
  θ::Float64

  σ²::Float64
  Φ::Function
end

#constructors for the Bates parameter set
ΘHeston(S₀::Float64, τ::Float64, X::Float64, r::Float64, q::Float64,
  σ::Float64, ρ::Float64, v₀::Float64, κ::Float64, θ::Float64)::ΘHeston =
  ΘHeston(S₀, τ, X, r, q, σ, ρ, v₀, κ, θ, σ*σ, ΦHeston)

ΘHeston(;S₀::Float64=0, τ::Float64=0, X::Float64=0, r::Float64=0, q::Float64=0,
  σ::Float64=0, ρ::Float64=0, v₀::Float64=0,
  κ::Float64=0, θ::Float64=0)::ΘHeston =
  ΘHeston(S₀, τ, X, r, q, σ, ρ, v₀, κ, θ)

#this is the characteristic function for Bates
function ΦHeston(Θ::ΘHeston, ω::T)::Complex{Float64} where T<:Union{Complex{Float64}, Float64}
  ι::Complex{Float64} = Ι

  #caclulate helper variables
  ρσiω::Complex{Float64} = Θ.ρ*Θ.σ*ι*ω
  d::Complex{Float64} = √((ρσiω-Θ.κ)^2+Θ.σ²*(ι*ω+ω*ω))
  g::Complex{Float64} = (Θ.κ-ρσiω-d)/(Θ.κ-ρσiω+d)

  #calculate next level of helper variables
  A::Complex{Float64} = ι * ω * log(Θ.S₀) + ι*ω*(Θ.r-Θ.q) * Θ.τ
  B::Complex{Float64} = Θ.θ*Θ.κ/Θ.σ² * ((Θ.κ-ρσiω-d)*Θ.τ-2.0*log((1.0-g*e^(-d*Θ.τ))/(1.0-g)))
  C::Complex{Float64} = Θ.v₀/Θ.σ² * (Θ.κ-ρσiω-d)*(1.0-e^(-d*Θ.τ)) / (1.0-g*e^(-d*Θ.τ))

  return e^(A+B+C)
end

function priceFromΦ(Θ::T)::Float64 where {T<:AbstractΘState}
  ι::Complex{Float64} = Ι

  #Integrand for Π1 (transformed to integrate from 0 to 1)
  function Φ1(ω::Float64)::Float64
    Ω::Float64 = ω/(1.0-ω)
    return real(e^(-ι * Ω * log(Θ.X)) * Θ.Φ(Θ,Ω - ι) / (ι * Ω * Θ.Φ(Θ,-ι))) / (1.0 - ω)^2
  end

  #Integrand for Π2 (transformed to integrate from 0 to 1)
  function Φ2(ω::Float64)::Float64
    Ω::Float64 = ω/(1.0-ω)
    return real(e^(-ι * Ω * log(Θ.X)) * Θ.Φ(Θ,Ω)/(ι * Ω)) / (1.0 - ω)^2
  end

  #do the integrations

  #implement with pure julia code
  Π1::Float64 = 0.5 + 1.0/π * I1D(Φ1, 0.0, 1.0)
  Π2::Float64 = 0.5 + 1.0/π * I1D(Φ2, 0.0, 1.0)

  #use the c implementation for a small performance boost
  #Π1::Float64 = 0.5 + 1.0/π * (hquadrature(Φ1, 0.0, 1.0))[1]
  #Π2::Float64 = 0.5 + 1.0/π * (hquadrature(Φ2, 0.0, 1.0))[1]
  #hquadrature(f, a, b)

  return Θ.S₀ * e^(-Θ.τ * Θ.q) * Π1 - Θ.X * e^(-Θ.τ * Θ.r) * Π2
end



#gets the price of a normal black scholes model
function priceBSMCall(S₀::Float64, τ::Float64, X::Float64, r::Float64, q::Float64, σ::Float64)::Float64
  d1::Float64 = log(S₀/X) + (r-2.0*q+σ^2/2.0)
  d1= d1 / (τ^0.5 * σ)

  return S₀*exp(-q*τ)*CTΦ(d1)-X*exp(-r*τ)*CTΦ(d1-σ*τ^0.5)
end

#implied volatility
function getImpliedσ(price::Float64, S₀::Float64, τ::Float64, X::Float64,
  r::Float64, q::Float64)::Float64


  #need to define special functions here for the solver


  #Algorithms: :LD_LBFGS, :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQP, :LD_TNEWTON, :LD_VAR2/:LD_VAR1,
  #  :LD_MMA, :GN_ESCH, :GN_DIRECT, GN_DIRECT_L, :GD_STOGO
  # GN_ISRES, :LD_LBFGS, :LN_SBPLX,
  #GN_CRS2_LM
  #:GN_CRS2_LM***, :GN_ISRES*** (slow), :LD_MMA**, LN_SBPLX**, LD_TNEWTON** best:62,848


  function obj(σ::Real)::Real
    d1 = (log(S₀/X) + τ*(r-q+σ^2/2.0))/ (σ*τ^0.5)

    calcPrice = S₀*exp(-q*τ)*CTΦ(d1) - X*exp(-r*τ)*CTΦ(d1-σ*τ^0.5)
    return (price-calcPrice)^2
  end

  #works well:LN_NELDERMEAD , LN_SBPLX, LD_SLSQP, LN_BOBYQA
  mod = Model(solver=NLoptSolver(algorithm=:LN_BOBYQA, maxtime=1))
  #JuMP.register(mod, :priceFromVol, 1, priceFromVol, autodiff=true)
  JuMP.register(mod, :CTΦ, 1, CTΦ, autodiff=true)
  JuMP.register(mod, :obj, 1, obj, autodiff=true)

  #constrain vol to 0 and 1000%
  @variable(mod, 1.0 >= σ >= .0001, start = 0.2)

  #set the objective function
  @NLobjective(mod, Min, obj(σ))

  σBest::Float64  = 0.0
  objBest::Float64 = -9999.0

  #try
  output = solve(mod)
  σBest  = getvalue(σ)

  return σBest

end

#helper for getting the implied volatiltiy from a parameter obhect
getImpliedσ(price::Float64, Θ::T where T<:AbstractΘState) =
  getImpliedσ(price, Θ.S₀, Θ.τ, Θ.X, Θ.r, Θ.q)

#gets the implied volatility directly
function getImpliedσ(Θ::T where T<:AbstractΘState)::Float64
  price::Float64 = priceFromΦ(Θ)
  #println("price: $price")
  return getImpliedσ(price, Θ)
end


#script function for running the GARCH problem
function A5P2Script(;refreshData::Bool=true, doGraphs::Bool = true)
  if !isfile("$DATA_PATH\\$SP500_NAME.jls") || refreshData
    preProcessSP500()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls")
  sp500DF::DataFrame = deserialize(stream)
  close(stream)

  returns::Vector{Float64} = sp500DF[2:end,:returns]

  S₀::Float64 = 100.0
  τ::Float64 = 0.33
  X::Float64 = 120.0
  r::Float64 = 0.04
  q::Float64 = 0.0
  σ::Float64 = 0.15

  testBSMOption::ΘBSM = ΘBSM(S₀=S₀, τ=τ, X=X, r=r, q=q, σ=σ)
  testBSMPrice::Float64 = priceFromΦ(testBSMOption)
  println("Price of BSM option: $testBSMPrice")

  testImpliedσ::Float64 = getImpliedσ(testBSMPrice, testBSMOption)
  println("Implied vol: $testImpliedσ")

  #test the parameters for the Bates model
  ρ::Float64 = -0.43
  v₀::Float64 = 0.15^2

  λ::Float64 = 1.764
  κ::Float64 = 0.011
  θ::Float64 = 0.09
  μⱼ::Float64 = -0.0301
  σⱼ::Float64 = 0.11

  testBatesOption::ΘBates = ΘBates(S₀=S₀, τ=τ, X=X, r=r, q=q,
    σ=σ, ρ=ρ, v₀=v₀, λ=λ, κ=κ, θ=θ, μⱼ=μⱼ, σⱼ=σⱼ)

  testBatesPrice::Float64 = priceFromΦ(testBatesOption)
  println("Price of Bates option: $testBatesPrice")

  testImpliedσ = getImpliedσ(testBatesOption)
  println("Implied vol: $testImpliedσ")

  ΦtestBates(testBatesOption, 1.0)

  testHestonOption::ΘHeston = ΘHeston(S₀=S₀, τ=τ, X=X, r=r, q=q,
    σ=σ, ρ=ρ, v₀=v₀, κ=κ, θ=θ)

  testHestonPrice::Float64 = priceFromΦ(testHestonOption)
  println("Price of Heston option: $testHestonPrice")

  testImpliedσ = getImpliedσ(testHestonOption)
  println("Implied vol: $testImpliedσ")

  #this sets the ranges
  if doGraphs
    τs::Vector{Float64} = collect(0.01:0.01:5.0)
    Xs::Vector{Float64} = collect(50:0.1:150)

    Nτ::Int = length(τs)
    NX::Int = length(Xs)

    #setup the option parameters
    BatesVsτ::Vector{ΘBates} = ((t::Float64)->ΘBates(S₀, t, X, r, q,
      σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ)).(τs)
    BatesVsX::Vector{ΘBates} = ((x::Float64)->ΘBates(S₀, τ, x, r, q,
      σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ)).(Xs)

    HestonVsτ::Vector{ΘHeston} = ((t::Float64)->ΘHeston(S₀, t, X, r, q,
      σ, ρ, v₀, κ, θ)).(τs)
    HestonVsX::Vector{ΘHeston} = ((x::Float64)->ΘHeston(S₀, τ, x, r, q,
      σ, ρ, v₀, κ, θ)).(Xs)

    #pre-allocate the implied volatilities
    BatesVsτIV::Vector{Float64} = Vector{Float64}(length(τs))
    BatesVsXIV::Vector{Float64} = Vector{Float64}(length(Xs))
    HestonVsτIV::Vector{Float64} = Vector{Float64}(length(τs))
    HestonVsXIV::Vector{Float64} = Vector{Float64}(length(Xs))

    println("Getting $Nτ τ estimates")
    Threads.@threads for i::Int ∈ 1:Nτ
      BatesVsτIV[i] = getImpliedσ(BatesVsτ[i])
      HestonVsτIV[i] = getImpliedσ(HestonVsτ[i])
    end

    println("Getting $NX X estimates")
    Threads.@threads for i::Int ∈ 1:NX
      BatesVsXIV[i] = getImpliedσ(BatesVsX[i])
      HestonVsXIV[i] = getImpliedσ(HestonVsX[i])
    end

    #form the data so that it can be graphed
    τDF::DataFrame = DataFrame(τ=τs, σ=BatesVsτIV)
    τDF[:,:model] = :Bates
    τDF = [τDF; DataFrame(DataFrame(τ=τs, σ=HestonVsτIV, model=[:Heston for i∈1:Nτ]))]

    XDF::DataFrame = DataFrame(X=Xs, σ=BatesVsXIV)
    XDF[:,:model] = :Bates
    XDF = [XDF; DataFrame(DataFrame(X=Xs, σ=HestonVsXIV, model=[:Heston for i∈1:NX]))]

    #get the graphs
    Gadfly.push_theme(:default)
    plots::Vector{PlotContainer} = Vector{PlotContainer}()
    plotNames::Vector{String} = Vector{String}()

    push!(plotNames, "TimeSmile")
    push!(plots, plot(τDF, x=:τ, y=:σ, color=:model,
      Geom.line, style(key_position=:right),
      Guide.ylabel("Implied σ"), Guide.xlabel("Time to Maturity"), Coord.Cartesian(ymin=0.0,ymax=1.0, xmin=0.0, xmax=5.0),
      Guide.title("Volatility Smile by Time (μⱼ=-0.0301, θ=.09, σⱼ=0.11)")))

    push!(plotNames, "MoneySmile")
    push!(plots, plot(XDF, x=:X, y=:σ, color=:model,
       Geom.line, style(key_position=:right),
      Guide.ylabel("Implied σ"), Guide.xlabel("X Given S₀=100"), Coord.Cartesian(ymin=0.0,ymax=1.0, xmin=50, xmax=150),
      Guide.title("Volatility Smile by Moneyness (μⱼ=-0.0301, θ=.09, σⱼ=0.11))")))

    ####Now repeat for a wider smile
    μⱼ = -0.20
    σ = 0.4
    σⱼ = 0.25
    v₀ = 0.2^2

    #get the new option models
    BatesVsτ = ((t::Float64)->ΘBates(S₀, t, X, r, q,
      σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ)).(τs)
    BatesVsX = ((x::Float64)->ΘBates(S₀, τ, x, r, q,
      σ, ρ, v₀, λ, κ, θ, μⱼ, σⱼ)).(Xs)

    HestonVsτ = ((t::Float64)->ΘHeston(S₀, t, X, r, q,
      σ, ρ, v₀, κ, θ)).(τs)
    HestonVsX = ((x::Float64)->ΘHeston(S₀, τ, x, r, q,
      σ, ρ, v₀, κ, θ)).(Xs)


    println("Getting $Nτ τ estimates")
    Threads.@threads for i::Int ∈ 1:Nτ
      BatesVsτIV[i] = getImpliedσ(BatesVsτ[i])
      HestonVsτIV[i] = getImpliedσ(HestonVsτ[i])
    end

    println("Getting $NX X estimates")
    Threads.@threads for i::Int ∈ 1:NX
      BatesVsXIV[i] = getImpliedσ(BatesVsX[i])
      HestonVsXIV[i] = getImpliedσ(HestonVsX[i])
    end

    #form the data so that it can be graphed
    τDF = DataFrame(τ=τs, σ=BatesVsτIV)
    τDF[:,:model] = :Bates
    τDF = [τDF; DataFrame(DataFrame(τ=τs, σ=HestonVsτIV, model=[:Heston for i∈1:Nτ]))]

    XDF = DataFrame(X=Xs, σ=BatesVsXIV)
    XDF[:,:model] = :Bates
    XDF = [XDF; DataFrame(DataFrame(X=Xs, σ=HestonVsXIV, model=[:Heston for i∈1:NX]))]


    push!(plotNames, "TimeSmile Enhanced")
    push!(plots, plot(τDF, x=:τ, y=:σ, color=:model,
      Geom.line, #Theme(line_width = 0.25pt),
      Guide.ylabel("Implied σ"), Guide.xlabel("Time to Maturity"), Coord.Cartesian(ymin=0.0,ymax=1.0, xmin=0.0, xmax=5.0),
      Guide.title("Volatility Smile by Time (μⱼ=-0.15, σ=.4, σⱼ=0.25, v₀ = 0.2²)")))

    push!(plotNames, "MoneySmile Enhanced")
    push!(plots, plot(XDF, x=:X, y=:σ, color=:model,
       Geom.line, #Theme(line_width = 0.25pt),
      Guide.ylabel("Implied σ"), Guide.xlabel("X Given S₀=100"), Coord.Cartesian(ymin=0.0,ymax=1.0, xmin=50, xmax=150),
      Guide.title("Volatility Smile by Moneyness (μⱼ=-0.15, σ=.4, σⱼ=0.25, v₀ = 0.2²))")))

      for i ∈ 1:length(plots)
        draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 7.5inch, 9inch),plots[i])
      end
  end




end


@time begin
  A5P2Script(refreshData=false, doGraphs=true)
end

end
