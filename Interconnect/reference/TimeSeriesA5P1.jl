module TimeSeriesA5P1

#Use this to print
#=
weave(Pkg.dir("$(pwd())\\FamaMacbeth.jl"),
  informat="script",
  out_path = "$(pwd())\\FamaMacbethOUT.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD, ForwardDiff,
  Gadfly, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  Measures, Formatting, ParallelAccelerator

importall CTIO, CTReg, CTStat

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const YAHOO_DATE_FORMAT = "m/d/yyyy"
const FRED_DATE_FORMAT = "m/d/yyyy"
const PPI_DATE_FORMAT = "m/d/yyyy"
const SP500_NAME = "GSPC"
const DBV_NAME = "DBV"
const FRED_NAME = "CPI"
const PPI_NAME = "PPI"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

n2s(x::T where T<:Real) = num2Str(x, 4, Ints=true)


#holds parameters for a measurement error simulation
mutable struct AR1Measurement
  c::Vector{T} where T<:Real
  ϕ::T where T<:Real
  dist::T where T<:MultivariateDistribution
  ARMA11::ARMA
end

#main constructor, also creates an ARMA11
AR1Measurement(ϕ::Float64 ;c::Vector{Float64}=zeros(2), Σ::Matrix{Float64}=eye(2),
  dist::T where T<:MultivariateDistribution = MvNormal(2,Σ))::AR1Measurement =
    AR1Measurement(c,ϕ,dist,
    ARMA(ϕ=[ϕ], ϑ=[-ϕ^2/2],σ²=ones(2)'*Σ*ones(2)))

mutable struct AR1MeasurementRet
  X::Vector{T} where T<:Real
  Y::Vector{T} where T<:Real
  ε::Matrix{T} where T<:Real
end

#simulates returns
function AR1MeasurementRet(Θ::AR1Measurement, N::Int;
    x₀::T=0.0,  Y::Vector{T} = Vector{typeof(x₀)}(N), X::Vector{T} = similar(Y),
    ε::Matrix{Float64} = rand(Θ.dist, N))::AR1MeasurementRet where T<:Real

  X[1] = Θ.ϕ*x₀ + ε[1,1]
  Y[1] = X[1] + ε[2,1]

  for i ∈ 2:N
      X[i] = Θ.ϕ*X[i-1] + ε[1,i]
      Y[i] = X[i] + ε[2,i]
  end

  return AR1MeasurementRet(X,Y,ε)
end

#calculate the analytical auto correlations
ARMeasurementAC(Θ::AR1Measurement, k::Int)::Float64 = Θ.ϕ^k*(2-Θ.ϕ^2.)/(1.0-Θ.ϕ^2)
ARMeasurementAC(Θ::AR1Measurement, r::Range)::Vector{Float64} =
  ((k::Int)->ARMeasurementAC.(Θ,k)).(r)

#problem specific graphs
function visualizeA5P1(Θ::AR1Measurement,
  ACRange::Range)::Void

  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)

  push!(plotNames, "ACF Graphs")
  push!(plots, plot(
    layer(x=collect(ACRange), y=ARMeasurementAC(Θ,ACRange),  Geom.line, Geom.point,
      Theme(default_color=COLORS[1])),
    #Guide.manual_color_key("Legend", ["Actual", "AR1"], COLORS[1:2]),
    Guide.xlabel("Lag"),Guide.ylabel("Autocorrelation"),
    Guide.title("AutoCorrelation")))

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing
end

#script for the problem
function A5P1Script(;refreshData::Bool = true)::Void

  ϕ::Float64 = 0.85 #Process parameter
  Σ::Matrix{Float64} = eye(2) #covariance matrix for shocks
  ACRange::Range = 1:25 #range of autocorrelations
  N::Int = 50_000

  Θ::AR1Measurement = AR1Measurement(ϕ, c=zeros(2), Σ=Σ, dist = MvNormal(Σ))
  returns::AR1MeasurementRet = AR1MeasurementRet(Θ,N)
  println("Analytical AC: $(ARMeasurementAC(Θ,ACRange))")
  #println("Y: $(returns.ε)")
  println("Mean ε: $(mean(returns.ε))")
  println("Σ: $(cov(returns.X,returns.Y))")
  println("σ^2: $(var(returns.Y))")
  println("Analytical: $((2-Θ.ϕ^2)/(1-Θ.ϕ^2))")

  #Algorithms: :LD_LBFGS* , :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQPx , :LD_TNEWTON*, :LD_VAR2/:LD_VAR1 ,
  #  :LD_MMA , :GN_ESCH, :GN_DIRECT , GN_DIRECT_L ,
  # GN_ISRES, :LN_NELDERMEAD*, :LN_SBPLX, :LD_LBFGS , :GD_STOGO
  #GN_CRS2_LM

  startingARMA = ARMA(c=0.0, ϕ=[0.85], ϑ=[-0.35])

  est::ARMAEstimate =
      ARMAEstimate(ARMA(1,1), returns.Y,
        alg = :LD_TNEWTON, verbose = true, bound=5.0, boundσ=5.0, maxTime=120., suppressC=false)
  println("σ²: $(est.Θ.σ²)")
  println("c: $(est.Θ.c)")
  println("ϕ: $(est.Θ.ϕ)")
  println("ϑ: $(est.Θ.ϑ)")
  println("VParam: $((2est.Θ.ϑ[1]^2+4est.Θ.ϕ[1]*est.Θ.ϑ[1]+2)/(1-est.Θ.ϕ[1]^2))")
  println("Unc: $(var(est.YHat))")



  #visualizeA5P1(Θ,ACRange)


  #visualizeA4P2!()


  return nothing
end


@time begin
  srand(11)
  A5P1Script(refreshData = true)
end

end
