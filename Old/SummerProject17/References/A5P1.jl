module A5P1


#Use this to print
#=
weave(Pkg.dir("$(pwd())\\A5P1.jl"),
  informat="script",
  out_path = "$(pwd())\\A5P1.html",
  doctype = "md2html")
=#
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, ForwardDiff, Roots

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
const DEFAULT_LOGμ_TOL = 10.0^-10.0
const PlotContainer = Union{Plot,Gadfly.Compose.Context}
const CSRegSpec = Tuple{Symbol, Symbol, Int}
const CPRegSpec = Tuple{Symbol, Int}
const TOLERANCE = 10.0^-8.0

const MIN_VALS_FOR_EST = 10

abstract type AbstractARCH end


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

#holds a GARCH parameter object
mutable struct ΘGARCH <: AbstractARCH
  λ::Real
  β::Real
  α::Real
  μ::Real
  valid::Bool

end

#constructor
ΘGARCH(λ::Real, β::Real, α::Real, μ::Real)::ΘGARCH = ΘGARCH(λ, β, α, μ, true)

mutable struct ΘEGARCH <: AbstractARCH
  λ::Real
  β::Real
  α::Real
  μ::Real
  γ::Real
  valid::Bool

end

#constructor
ΘEGARCH(λ::Real, β::Real, α::Real, μ::Real, γ::Real)::ΘEGARCH = ΘEGARCH(λ, β, α, μ, γ, true)



function estimateGARCH(r::Vector{Float64}, verbose::Bool = false)::ΘGARCH

  #for convenience
  T::Int = length(r)

  function obj(λ::Real, β::Real, α::Real, μ::Real)

    #σ²::Vector{Real} = Vector{Real}(T-1)
    σ²tM1::Real = λ/(1.0-β)
    #l::Float64 = -(r[1]-μ)*(r[1]-μ)/(2.0*σ²₀)-0.5*log(σ²₀) #initial likelihood
    l::Real = 0.0

    #calculate the likelihood and get the next volatiltiy in the sequence
    for i::Int ∈ 1:T
      u²::Real = (r[i]-μ)*(r[i]-μ)
      l += -u²/σ²tM1 - 0.5*log(σ²tM1)
      σ²tM1 = λ+β*σ²tM1+α*u²
    end

    #println("λ=$λ, β=$β, α=$α, μ=$μ")

    return l
  end

  meanRet = mean(r)

  #set initial meta parameters
  minVal = 10.0^-10.0
  maxVal = 1.0-minVal

  #Algorithms: :LD_LBFGS, :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQP, :LD_TNEWTON, :LD_VAR2/:LD_VAR1,
  #  :LD_MMA, :GN_ESCH, :GN_DIRECT, GN_DIRECT_L, :GD_STOGO
  # GN_ISRES, :LN_NELDERMEAD, :LN_SBPLX,
  #GN_CRS2_LM
  #:GN_CRS2_LM***, :GN_ISRES*** (slow), :LD_MMA**, LN_SBPLX**, LD_TNEWTON** best:62,848
  #MMA: λ=2.991604487710108e-6, β=0.8691809689019686, α=0.2447041114406467, μ=0.0005059861641042039
  #NEWTON: 3.7136249890797323e-6, β=0.8641416287811661, α=0.23581482405030202, μ=0.0005248887137138578
  #CRS2 λ=1.8304983359585568e-6, β=0.9091365756370056, α=0.1664779062477393, μ=0.00048124017464920696
  #SBPLX λ=1.7187359905806954e-6, β=0.9122999784465299, α=0.16106473023211287, μ=0.00029331108088162115
  #LR Vol>0 #:LN_COBYLA(61450), LD_MMA(61144), ISRES(62578),


  mod = Model(solver=NLoptSolver(algorithm=:LD_TNEWTON_PRECOND,maxtime=10))

  @variable(mod, minVal <= λ <= maxVal)
  @variable(mod, minVal <= β <= maxVal)
  @variable(mod, minVal <= α <= maxVal)
  @variable(mod, minVal <= μ <= maxVal)
  #@NLconstraint(mod, 1.0-α-β>=0.0)

  setvalue(λ, 0.5)
  setvalue(β, 0.5)
  setvalue(α, 0.5)
  setvalue(μ, mean(r)) #guess the unconditional mean

  JuMP.register(mod, :obj, 4, obj, autodiff=true)
  @NLobjective(mod, Max, obj(λ, β, α, μ))

  ΘBest::ΘGARCH  = ΘGARCH(0.0, 0.0, 0.0, 0.0, false)
  objBest::Float64 = -9999.0
  try
    output = solve(mod)
    ΘBest  = ΘGARCH(getvalue(λ), getvalue(β), getvalue(α), getvalue(μ))
    objBest= obj(getvalue(λ), getvalue(β), getvalue(α), getvalue(μ))
  end


  if verbose
    impliedAnnualLRVol::Float64 = (ΘBest.λ/(1.0 - ΘBest.β)*255)^0.5
    historicalLRVol::Float64 = std(r)*255.0^0.5

    impliedReturn::Float64 = exp(ΘBest.μ*255)-1.0
    historicalReturn::Float64 = exp(mean(r)*255)-1.0

    println("λ=$(ΘBest.λ), β=$(ΘBest.β), α=$(ΘBest.α), μ=$(ΘBest.μ)")
    println("Vol: Implied=$impliedAnnualLRVol, Hist=$historicalLRVol")
    println("Return: Implied=$impliedReturn, Hist=$historicalReturn")
    println("Objective: $objBest")
    println("For Reference: x̄=$(mean(r))")
  end

  return ΘBest

end


function estimateEGARCH(r::Vector{Float64}, verbose::Bool = false)::ΘEGARCH

  T::Int = length(r)


  function obj(λ::Real, β::Real, α::Real, μ::Real, γ::Real)

    #σ²::Vector{Real} = Vector{Real}(T-1)
    EAbsZ::Float64 = (2.0/π)^0.5
    lσ²tM1::Real = (λ-γ*EAbsZ)/(1.0-β)
    #lσ²tM1::Real = λ/(1.0-β)


    l::Real = 0.0

    #calculate the likelihood and get the next volatiltiy in the sequence
    for i::Int ∈ 1:T
      u²::Real = (r[i]-μ)*(r[i]-μ)
      σtM1::Real = (exp(lσ²tM1))^0.5
      l += -u²/(σtM1*σtM1) - 0.5*lσ²tM1
      lσ²tM1 = λ+β*lσ²tM1+α*u²+γ*(abs((r[i]-μ)/σtM1)-EAbsZ)
    end

    return l
  end

  meanRet = mean(r)

  #set initial meta parameters
  minVal = -10.0
  maxVal = 100.0

  #Algorithms: :LD_LBFGS, :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQP, :LD_TNEWTON, :LD_VAR2/:LD_VAR1,
  #  :LD_MMA, :GN_ESCH, :GN_DIRECT, GN_DIRECT_L, :GD_STOGO
  # GN_ISRES, :LN_NELDERMEAD, :LN_SBPLX,
  #GN_CRS2_LM
  #:LD_LBFGS(62774), LD_VAR2(62774), LN_NELDERMEAD(61880), :LN_SBPLX (62791), :LD_TNEWTON(62800)


  mod = Model(solver=NLoptSolver(algorithm=:LD_TNEWTON_PRECOND, maxtime=20))

  @variable(mod, minVal <= λ <= maxVal)
  @variable(mod, minVal <= β <= maxVal)
  @variable(mod, minVal <= α <= maxVal)
  @variable(mod, 0.0 <= μ <= maxVal)
  @variable(mod, minVal <= γ <= maxVal)
  #@NLconstraint(mod, 1.0-α-β>=0.0)

  setvalue(λ, 0.5)
  setvalue(β, 0.5)
  setvalue(α, 0.5)
  setvalue(μ, mean(r)) #guess the unconditional mean
  setvalue(γ, 0.5) #guess the unconditional mean

  JuMP.register(mod, :obj, 5, obj, autodiff=true)
  @NLobjective(mod, Max, obj(λ, β, α, μ, γ))

  ΘBest::ΘEGARCH  = ΘEGARCH(getvalue(λ), getvalue(β), getvalue(α), getvalue(μ), getvalue(γ), false)
  objBest::Float64 = -9999.0
  try
    output = solve(mod)
    ΘBest  = ΘEGARCH(getvalue(λ), getvalue(β), getvalue(α), getvalue(μ), getvalue(γ))
    objBest= obj(getvalue(λ), getvalue(β), getvalue(α), getvalue(μ), getvalue(γ))
  end


  if verbose
    impliedAnnualLRVol::Float64 = (exp((ΘBest.λ-ΘBest.γ*(2/π)^0.5)/(1.0-ΘBest.β))*255)^0.5
    #impliedAnnualLRVol::Float64 = (exp(ΘBest.λ)/(1-ΘBest.β)*255)^0.5
    historicalLRVol::Float64 = std(r)*255.0^0.5

    impliedReturn::Float64 = exp(ΘBest.μ*255)-1.0
    historicalReturn::Float64 = exp(mean(r)*255)-1.0

    println("λ=$(ΘBest.λ), β=$(ΘBest.β), α=$(ΘBest.α), μ=$(ΘBest.μ), γ=$(ΘBest.γ)")
    println("Vol: Implied=$impliedAnnualLRVol, Hist=$historicalLRVol")
    println("Return: Implied=$impliedReturn, Hist=$historicalReturn")
    println("Objective: $objBest")
    println("For Reference: x̄=$(mean(r))")
  end

  return ΘBest

end

function estimateGARCHByQ(sp500DF::DataFrame, estimateModel::Function)#::T where T <: AbstractARCH

  dateDict::Dict = Dict(sp500DF[i,:DATE]=>i for i::Int ∈ 1:size(sp500DF,1))

  minValsForEst::Int = MIN_VALS_FOR_EST
  r::Vector{Float64} = sp500DF[2:end,:returns]

  #get a list of target indices
  quarterlyDates::Vector{Date} = sp500DF[sp500DF[:,:dayOfQuarter].==1,:DATE]
  terminalIndicesToEst::Vector{Int} = ((d::Date)->dateDict[d]).(quarterlyDates)
  #terminalIndicesToEst = terminalIndicesToEst[terminalIndicesToEst[:].>minValsForEst]

  #set up a container of contianers
  estimated::Vector{Vector{T where T<:AbstractARCH}} =
    [Vector{T where T<:AbstractARCH}() for i::Int ∈ 1:length(terminalIndicesToEst)]

  println("beginning $(length(estimated)) simulations")
  fails::Int = 0
  @fastmath for i ∈ 1:length(terminalIndicesToEst)

  #estimate the grabbed indices
    estimated[i] = [estimateModel(r[1:(terminalIndicesToEst[i]-1)], false)]
    if estimated[i][1].valid == false
      fails += 1
    end
  end
  println("simulations complete. $fails failures out of $(length(estimated)) simulations.")

  return [estimated[i][1] for i::Int ∈ 1:length(terminalIndicesToEst)]

end

#script function for running the GARCH problem
function A5P1Script(;refreshData::Bool=true, refreshGARCH::Bool=true, refreshEGARCH::Bool=true,)
  if !isfile("$DATA_PATH\\$SP500_NAME.jls") || refreshData
    preProcessSP500()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls")
  sp500DF::DataFrame = deserialize(stream)
  close(stream)

  r::Vector{Float64} = sp500DF[2:end,:returns]
  quarterlyDates::Vector{Date} = sp500DF[sp500DF[:,:dayOfQuarter].==1,:DATE]

  #estimate the GARCH model and print the output

  estimateGARCH(r, true)
  estimateEGARCH(r, true)

  #this setup is so that we don't ahve to solve the model each quarter repeatedly
  if refreshGARCH
    estimatedGARCH::Vector{ΘGARCH} = estimateGARCHByQ(sp500DF, estimateGARCH)

    #cache the results on disk
    stream = open("$OUTPUT_PATH\\GARCH_RESULTS.jls", "w")
    serialize(stream, estimatedGARCH)
    close(stream)
  else

    #access the resutls
    stream = open("$OUTPUT_PATH\\GARCH_RESULTS.jls")
    estimatedGARCH = deserialize(stream)
    close(stream)
  end

  if refreshEGARCH
    estimatedEGARCH::Vector{ΘEGARCH} = estimateGARCHByQ(sp500DF, estimateEGARCH)

    #cache the results on disk
    stream = open("$OUTPUT_PATH\\EGARCH_RESULTS.jls", "w")
    serialize(stream, estimatedEGARCH)
    close(stream)
  else

    #access the resutls
    stream = open("$OUTPUT_PATH\\EGARCH_RESULTS.jls")
    estimatedEGARCH = deserialize(stream)
    close(stream)
  end

  ##############Now work on the outputs. First make them into a dataframe
  garchDF::DataFrame = DataFrame(DATE=quarterlyDates, λ=[g.λ for g ∈ estimatedGARCH], β=[g.β for g ∈ estimatedGARCH],
    α=[g.α for g ∈ estimatedGARCH], μ=[g.μ for g ∈ estimatedGARCH],
    μLR=[exp(255*g.μ)-1.0 for g ∈ estimatedEGARCH], σ=[(g.λ/(1.0 - g.β)*255)^0.5 for g ∈ estimatedGARCH])

  egarchDF::DataFrame  = DataFrame(DATE=quarterlyDates, λ=[g.λ for g ∈ estimatedEGARCH], β=[g.β for g ∈ estimatedEGARCH],
    α=[g.α for g ∈ estimatedEGARCH], μ=[g.μ for g ∈ estimatedEGARCH], γ=[g.γ for g ∈ estimatedEGARCH],
    μLR=[exp(255*g.μ)-1.0 for g ∈ estimatedEGARCH], σ=[(exp((g.λ-g.γ*(2/π)^0.5)/(1.0-g.β))*255)^0.5 for g ∈ estimatedEGARCH])

  nQuarters::Int = length(quarterlyDates)

  for i::Int∈1:nQuarters
    if estimatedGARCH[i].valid == false
      garchDF[i,[:λ, :β, :α, :μ, :σ]] = NA
    end

    if estimatedEGARCH[i].valid == false
      egarchDF[i,[:λ, :β, :α, :μ,:γ, :σ]] = NA
    end
  end

  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  push!(plotNames, "GarchParameters")
  push!(plots, vstack(plot(x=garchDF[:,:DATE], y=garchDF[:,:λ],
    Geom.point, #Theme(line_width = 0.25pt),
    Guide.ylabel("Value"), Guide.xlabel("Date"),
    Guide.title("GARCH λ Param")),
    plot(x=garchDF[:,:DATE], y=garchDF[:,:β],
       Geom.point, #Theme(line_width = 0.25pt),
      Guide.ylabel("Value"), Guide.xlabel("Date"),
      Guide.title("GARCH β Param")),
    plot(x=garchDF[:,:DATE], y=garchDF[:,:α],
      Geom.point, #Theme(line_width = 0.25pt),
      Guide.ylabel("Value"), Guide.xlabel("Date"),
      Guide.title("GARCH α Param"))))

  egarchDF = egarchDF[((!isna).(egarchDF[:,:σ])),:]

  push!(plotNames, "EGarchParameters")
  push!(plots, vstack(plot(x=egarchDF[:,:DATE], y=egarchDF[:,:λ],
    Geom.point, #Theme(line_width = 0.25pt),
    Guide.ylabel("Value"), Guide.xlabel("Date"),
    Guide.title("GARCH λ Param")),
    plot(x=egarchDF[:,:DATE], y=egarchDF[:,:β],
       Geom.point, #Theme(line_width = 0.25pt),
      Guide.ylabel("Value"), Guide.xlabel("Date"),
      Guide.title("GARCH β Param")),
    plot(x=egarchDF[:,:DATE], y=egarchDF[:,:α],
      Geom.point, #Theme(line_width = 0.25pt),
      Guide.ylabel("Value"), Guide.xlabel("Date"),
      Guide.title("GARCH α Param"))))

  push!(plotNames, "GarchVol")
  push!(plots, vstack(plot(x=garchDF[:,:DATE], y=garchDF[:,:μLR],
    Geom.point, #Theme(line_width = 0.25pt),
    Guide.ylabel("Value"), Guide.xlabel("Date"), Coord.Cartesian(ymin=0.0,ymax=0.5),
    Guide.title("GARCH μ Param")),
    plot(x=garchDF[:,:DATE], y=garchDF[:,:σ],
       Geom.point, #Theme(line_width = 0.25pt),
      Guide.ylabel("Value"), Guide.xlabel("Date"),
      Guide.title("GARCH σ Param"))))

  push!(plotNames, "EGarchVol")
  push!(plots, vstack(plot(x=egarchDF[:,:DATE], y=egarchDF[:,:μLR],
    Geom.point, #Theme(line_width = 0.25pt),
    Guide.ylabel("Value"), Guide.xlabel("Date"), Coord.Cartesian(ymin=0.0,ymax=0.5),
    Guide.title("GARCH μ Param")),
  plot(x=egarchDF[:,:DATE], y=egarchDF[:,:σ],
     Geom.point, #Theme(line_width = 0.25pt),
    Guide.ylabel("Value"), Guide.xlabel("Date"),
    Guide.title("GARCH σ Param"))))


  for i ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 7.5inch, 9inch),plots[i])
  end


  #graph the GARCH model

end


@time begin
  #A5P1Script(refreshData=false, refreshGARCH=false, refreshEGARCH=false)
end

end
