module  A2

#this is useful if we have mdoules
#=
weave(Pkg.dir("$(pwd())\\A2.jl"),
  informat="script",
  out_path = "$(pwd())\\A2.html",
  doctype = "md2html")
=#
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, ForwardDiff

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "m/d/y"
const SP500_NAME = "GSPC"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"

const DEFAULT_ITER = 12000
const DEFAULT_BURN = 2000
const DO_NOT_WRITE = false
const PlotContainer = Union{Plot,Gadfly.Compose.Context}


function preProcessSP500()::Void
  dateFormat::String = DATE_FORMAT_STR
  sp500DF::DataFrame = readtable("$DATA_PATH\\$SP500_NAME.csv")

  #fix the dates
  rename!(sp500DF, :Date, :DATE)
  sp500DF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(sp500DF[:,:DATE])::DataVector{Date}
  sort!(sp500DF, cols = [:DATE])
  sp500DF[:,:returns] = similar(sp500DF[:,:Adj_Close])
  sp500DF[2:end,:returns] = log.(sp500DF[2:end,:Adj_Close] ./  sp500DF[1:(end-1),:Adj_Close])
  #showcols(sp500DF)

  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls", "w")
  serialize(stream, sp500DF)
  close(stream)

  return nothing
end

#contains all prior parameters
struct GibbsPriors
  priorNames::Vector{Symbol}
  priorValues::Vector{Float64}
  index::Dict
end

#constructor which forms the index
GibbsPriors(priorNames::Vector{Symbol}, priorValues::Vector{Float64})::GibbsPriors =
            GibbsPriors(priorNames, priorValues,
              Dict(priorNames[i]=>priorValues[i] for i::Int ∈ 1:length(priorNames)))

#make the priors
function GibbsPriors(priorTuples::Vector{Tuple{Symbol,Float64}})::GibbsPriors
  priorNames = Vector{Symbol}(length(priorTuples))
  priorValues = Vector{Float64}(length(priorTuples))

  for i::Int ∈ 1:length(priorTuples)
    priorNames[i], priorValues[i] = priorTuples[i]
  end

  return GibbsPriors(priorNames, priorValues)
end

#this contains parameters and all distributions
mutable struct GibbsProblem
  names::Vector{Symbol}
  values::Vector{Vector{Float64}}
  posteriors!::Vector{Function}
  priors::GibbsPriors
  Y::Vector{Float64} #the data
  paramDistribution::Matrix{Vector{Float64}}
  recordParams::Vector{Bool}

  index::Dict
  targetIterations::Int
  burnIn::Int
  iter::Int
  T::Int
end

#constructor for a Gibbs Problem
#Requires all symbols, posterior distribution functions, initial values, and simulation parameters
function GibbsProblem(names::Vector{Symbol},
            posteriors!::Vector{Function},
            priors::GibbsPriors,
            Y::Vector{Float64};
            initialValues::Vector{Vector{Float64}}=1.0/length(names),
            targetIterations::Int=DEFAULT_ITER,
            burnIn::Int = DEFAULT_BURN,
            recordParams::Vector{Bool} = trues(length(names)))::GibbsProblem

  #this is the space where we will store the parameters
  paramDistribution::Matrix{Vector{Float64}} =
    Matrix{Vector{Float64}}(length(initialValues), targetIterations - burnIn)

  #pre-allocate the space
  for c::Int ∈ 1:size(paramDistribution, 2)
    for r::Int ∈ 1:length(initialValues)
      if recordParams[r]
        paramDistribution[r,c] = similar(initialValues[r])
      end
    end
  end

  index::Dict =   Dict(names[i] => i for i::Int ∈ 1:length(initialValues))

  return GibbsProblem(names, initialValues, posteriors!, priors, Y,
  paramDistribution, #pre-allocated space to store realizations
  recordParams,
  index, #an index of the variable names
  targetIterations,
  burnIn,
  1, #start the iteration counter
  length(Y))
end

#realizes a complete set of parameters
function IterateProblem!(p::GibbsProblem)::Void

  for i::Int ∈ 1:length(p.names)
    p.posteriors![i](p)

    if p.iter > p.burnIn && p.recordParams[i]
      p.paramDistribution[i, p.iter-p.burnIn] .= p.values[i]
    end
  end

  p.iter += 1
  return nothing
end

#performans the Gibbs sampling exercise
function getParameterDistribution!(p::GibbsProblem)::Void
  @fastmath for i::Int ∈ 1:p.targetIterations
    IterateProblem!(p)
  end

  return nothing
end

getSeriesLabels(N::Int, prefix::Symbol = :Y)::Vector{Symbol} = [Symbol(prefix,i) for i::Int ∈ 1:N]

#easy access of values
get(p::GibbsProblem, s::Symbol)::Vector{Float64} = p.values[p.index[s]]
get(priors::GibbsPriors, s::Symbol)::Float64 = priors.index[s]

function getResult(p::GibbsProblem, s::Symbol)::Vector{Float64}

  valIndex::Int = p.index[s]
  result = Vector{Float64}(p.targetIterations - p.burnIn)

  for i::Int ∈ 1:length(result)
    result[i] = p.paramDistribution[valIndex, i][1]
  end

  return result
end

getResults(p::GibbsProblem, names::Vector{Symbol} = p.names)::Dict =
  Dict(s => getResult(p,s) for s::Symbol ∈ names)

#helper function to draw from the posterior mean of a normal distribution
function μPosterior!(p::GibbsProblem)::Void
  ##Priors
  θ::Float64 = get(p.priors, :θ)
  δ²::Float64 = get(p.priors, :δ²)

  #parameters
  σ²::Float64 = get(p,:σ²)[1]
  Z::Vector{Float64} = get(p, :Z)
  ξ::Vector{Float64} = get(p, :ξ)
  pμ::Vector{Float64} = get(p,:μ)

  Δ²::Float64 = 1.0/(p.T/σ²+1.0/δ²)

  #intermediate variables
  YₜMinusξₜZₜOverσ²::Float64 = sum(p.Y .-  (ξ .* Z))/σ²

  #form the distribution
  pμ[1] = rand(Normal(Δ²*(YₜMinusξₜZₜOverσ²+θ/δ²), Δ²^0.5))


  return nothing
end

function μₛPosterior!(p::GibbsProblem)::Void
  ##Priors
  θₛ::Float64 = get(p.priors, :θₛ)
  δ²ₛ::Float64 = get(p.priors, :δ²ₛ)

  #parameters
  σ²ₛ::Float64 = get(p,:σ²ₛ)[1]
  ξ::Vector{Float64} = get(p, :ξ)
  pμₛ::Vector{Float64} = get(p,:μₛ)

  Δ²::Float64 = 1.0/(p.T/σ²ₛ+1.0/δ²ₛ)

  #intermediate variables
  ξOverσ²::Float64 = sum(ξ)/σ²ₛ

  #form the distribution
  pμₛ[1] = rand(Normal(Δ²*(ξOverσ²+θₛ/δ²ₛ), Δ²^0.5))


  return nothing
end

#helper function to draw from the psoterior variance of a normal distribution
function σ²Posterior!(p::GibbsProblem)::Void
  ##Priors
  α::Float64 = get(p.priors, :α)
  β::Float64 = get(p.priors, :β)

  #parameters
  μ::Float64 = get(p, :μ)[1]
  Z::Vector{Float64} = get(p, :Z)
  ξ::Vector{Float64} = get(p, :ξ)
  pσ²::Vector{Float64} = get(p, :σ²)

  #intermediate variables
  YtMinusμ²::Float64 = sum((p.Y .- μ .- (ξ .* Z)).^2.0)

  #form the distribution
  pσ²[1] = rand(InverseGamma(p.T/2.+α,0.5*(YtMinusμ²+2.0*β)))
  #draw
  return nothing
end

#helper function to draw from the psoterior variance of a normal distribution
function σ²ₛPosterior!(p::GibbsProblem)::Void
  ##Priors
  αₛ::Float64 = get(p.priors, :αₛ)
  βₛ::Float64 = get(p.priors, :βₛ)

  #parameters
  μₛ::Float64 = get(p, :μₛ)[1]
  ξ::Vector{Float64} = get(p, :ξ)
  pσ²ₛ::Vector{Float64} = get(p, :σ²ₛ)

  #intermediate variables
  YtMinusμ²ₛ::Float64 = sum((ξ .- μₛ).^2.0)

  #form the distribution
  pσ²ₛ[1] = rand(InverseGamma(p.T/2.+αₛ,0.5*(YtMinusμ²ₛ+2.0*βₛ)))
  #draw
  return nothing
end

#helper function to draw λ
function λPosterior!(p::GibbsProblem)::Void
  ##Priors
  γ::Float64 = get(p.priors, :γ)
  η::Float64 = get(p.priors, :η)

  #parameters
  Z::Vector{Float64} = get(p, :Z)
  pλ::Vector{Float64} = get(p, :λ)

  #intermediate variables
  ΣZ::Float64 = sum(Z)

  #form the distribution and draw
  pλ[1] = rand(Beta(ΣZ+γ,p.T-ΣZ+η))

  return nothing
end

#helper function to draw ξₜ
function ξPosterior!(p::GibbsProblem)::Void
  #parameters
  μ::Float64 = get(p,:μ)[1]
  μₛ::Float64 = get(p,:μₛ)[1]
  σ²::Float64 = get(p,:σ²)[1]
  σ²ₛ::Float64 = get(p,:σ²ₛ)[1]

  Z::Vector{Float64} = get(p,:Z)

  pξ::Vector{Float64} = get(p,:ξ)

  for t::Int ∈ 1:p.T
    Δ²::Float64 = (1.0 / (Z[t] / σ² + 1.0 / σ²ₛ))

    pξ[t] =  rand(Normal(Δ² * ((p.Y[t] - μ) * Z[t] / σ² + μₛ/σ²ₛ), Δ²^0.5))
  end


  return nothing
end

#helper function to draw Zₜ
function ZPosterior!(p::GibbsProblem)::Void

  μ::Float64 = get(p,:μ)[1]
  σ²::Float64 = get(p,:σ²)[1]
  λ::Float64 = get(p,:λ)[1]
  ξ::Vector{Float64} = get(p,:ξ)

  pZ::Vector{Float64} = get(p,:Z)

  for t::Int ∈ 1:p.T
    propZ1::Float64 = λ * exp(-0.5 * (p.Y[t] - μ - ξ[t])^2.0 / σ²)
    propZ0::Float64 = (1.0 - λ) * exp(-0.5 * (p.Y[t] - μ)^2.0 / σ²)
    pZ[t] = rand(Bernoulli(propZ1/(propZ0+propZ1)))
  end

  return nothing
end



function problemA2Script(;refreshData::Bool = true, doNotWrite::Bool = DO_NOT_WRITE)::Void
  if !isfile("$DATA_PATH\\$SP500_NAME.jls") || refreshData
    preProcessSP500()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls")
  sp500DF::DataFrame = deserialize(stream)
  close(stream)
  Y::Vector{Float64} = sp500DF[2:end, :returns]
  T::Int = length(Y)

  #################Priors
  priors = GibbsPriors([
  (:θ, mean(Y)), (:δ², var(Y)), # μ~N(θ,δ²) #use std deviation as upper bound for δ
  (:α, 0.001), (:β, 0.001), #σ²~IG(α,β)
  (:θₛ, 0.0 ),  (:δ²ₛ, 1.0),  #μₛ~N(θₛ,δ²ₛ)
  (:αₛ, 0.001 ), (:βₛ, 0.001), #σ²ₛ~IG(αₛ,βₛ)
  (:γ, 1.0 ), (:η, 1.0)]) #λ~β(γ,η) (use uniform distribution)

  #Form posterior distributions
  names::Vector{Symbol} = Vector{Symbol}()
  posteriors!::Vector{Function} = Vector{Function}()
  initialValues::Vector{Vector{Float64}} =  Vector{Float64}()

  #assign variable labels, posterior functions, and initial values
  push!(names, :μ)
  push!(posteriors!, μPosterior!)
  push!(initialValues, [0.0])

  push!(names, :μₛ)
  push!(posteriors!, μₛPosterior!)
  #push!(posteriors!, μₛPosterior)
  push!(initialValues, [0.0])

  push!(names, :σ²)
  push!(posteriors!, σ²Posterior!)
  push!(initialValues, [0.01])

  push!(names, :σ²ₛ)
  push!(posteriors!, σ²ₛPosterior!)
  push!(initialValues, [0.01])

  push!(names, :λ)
  push!(posteriors!, λPosterior!)
  push!(initialValues, [0.01])

  push!(names, :ξ)
  push!(posteriors!, ξPosterior!)
  push!(initialValues, zeros(T))

  push!(names, :Z)
  push!(posteriors!, ZPosterior!)
  push!(initialValues, zeros(T))


  println("""Beginning Gibbs Problem
    targetIterations=$DEFAULT_ITER,
    burnIn = $DEFAULT_BURN
    length Of Values = $(length(initialValues))
    length Of posteriors! = $(length(posteriors!))
    length Of Names = $(length(names))

  """)

  p::GibbsProblem = GibbsProblem(names, posteriors!, priors, Y,
              initialValues=initialValues,
              targetIterations=DEFAULT_ITER,
              burnIn = DEFAULT_BURN,
              recordParams = [true, true, true, true, true, false, false])

  @time getParameterDistribution!(p)

  results::Dict = getResults(p, [:μ, :σ², :μₛ, :σ²ₛ, :λ])


  μMean::Float64 = mean(results[:μ])
  σ²Mean::Float64 = mean(results[:σ²])
  μₛMean::Float64 = mean(results[:μₛ])
  σ²ₛMean::Float64 = mean(results[:σ²ₛ])
  λMean::Float64 = mean(results[:λ])

  println("""
  Results (Iterations=$(p.targetIterations), Burnin=$(p.burnIn)): (Mean, annualized)
    μ: $(round(exp(μMean)-1.0,6))
    σ: $(round(σ²Mean^0.5,6))
    μₛ: $(round(exp(μₛMean)-1.0,6))
    σₛ: $(round(σ²ₛMean^0.5,6))
    λ: $(round(λMean,6))

    μ (Ann.): $(round(exp(255.0*μMean)-1.0,6))
    σ (Ann.): $(round((255.0 * σ²Mean)^0.5,6))
  """)

  if !doNotWrite
    gstream::GZipStream = gzopen("$DATA_PATH\\A2Out_$(p.targetIterations)_$(p.burnIn).jls.gz", "w")
    serialize(gstream, results)
    close(gstream)
  end

  return nothing
end

function getA2Graphs(;targetIterations::Int = DEFAULT_ITER, burnIn::Int = DEFAULT_BURN)::Void

  gstream::GZipStream = gzopen("$DATA_PATH\\A2Out_$(targetIterations)_$(burnIn).jls.gz")
  results::Dict =  deserialize(gstream)
  close(gstream)

  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls")
  sp500DF::DataFrame = deserialize(stream)
  close(stream)

  numPoints::Int = targetIterations-burnIn

  push!(results, :σ=>results[:σ²].^0.5)
  push!(results, :σₛ=>results[:σ²ₛ].^0.5)

  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()


  push!(plotNames, "IterationPlots1")
  push!(plots, vstack(
    plot(x=1:numPoints, y=results[:μ], Geom.line, Theme(line_width = 0.25pt),
      Guide.ylabel("μ"), Guide.xlabel("Iteration"),
      Guide.xticks(ticks=collect(0:(numPoints÷5):(numPoints))),
      Guide.title("Draws of μ")),
    plot(x=1:numPoints, y=(results[:σ]), Geom.line, Theme(line_width = 0.2pt),
      Guide.ylabel("σ"), Guide.xlabel("Iteration"),
      Guide.xticks(ticks=collect(0:(numPoints÷5):(numPoints))),
      Guide.title("Draws of σ"))))

  push!(plotNames, "IterationPlots2")
  push!(plots, vstack(
    plot(x=1:numPoints, y=results[:μₛ], Geom.line, Theme(line_width = 0.2pt),
      Guide.ylabel("μₛ"), Guide.xlabel("Iteration"),
      Guide.xticks(ticks=collect(0:(numPoints÷5):(numPoints))),
      Guide.title("Draws of μₛ")),
    plot(x=1:numPoints, y=(results[:σₛ]), Geom.line, Theme(line_width = 0.2pt),
      Guide.ylabel("σₛ"), Guide.xlabel("Iteration"),
      Guide.xticks(ticks=collect(0:(numPoints÷5):(numPoints))),
      Guide.title("Draws of σₛ")),
    plot(x=1:numPoints, y=(results[:λ]), Geom.line, Theme(line_width = 0.2pt),
      Guide.ylabel("λ"), Guide.xlabel("Iteration"),
      Guide.xticks(ticks=collect(0:(numPoints÷5):(numPoints))),
      Guide.title("Draws of λ"))))

  push!(plotNames, "HistogramPlots1")
  push!(plots, vstack(
    plot(x=results[:μ], Geom.histogram, Theme(line_width = 0.2pt),
      Guide.ylabel("μ"), Guide.xlabel("Iteration"),
      Guide.title("Distribution of μ")),
    plot(x=(results[:σ]), Geom.histogram, Theme(line_width = 0.2pt),
      Guide.ylabel("σ"), Guide.xlabel("Iteration"),
      Guide.title("Distribution of :σ"))))

  push!(plotNames, "HistogramPlots2")
  push!(plots, vstack(
    plot(x=results[:μₛ], Geom.histogram, Theme(line_width = 0.2pt),
      Guide.ylabel("μₛ"), Guide.xlabel("Iteration"),
      Guide.title("Distribution of μₛ")),
    plot(x=results[:σₛ], Geom.histogram, Theme(line_width = 0.2pt),
      Guide.ylabel("σₛ"), Guide.xlabel("Iteration"),
      Guide.title("Distribution of σₛ")),
    plot(x=results[:λ], Geom.histogram, Theme(line_width = 0.2pt),
      Guide.ylabel("λ"), Guide.xlabel("Iteration"),
      Guide.title("Distribution of λ"))))


  for i ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 7.5inch, 9inch),plots[i])
  end

  ########################make the summary data table
  summaryDecimals::Int = 6
  summaryCols::Vector{Symbol} = [:μ, :σ, :μₛ,:σₛ, :λ]
  summaryRowNames::Vector{String} =
    ["mean", "median", "Ann. Median", "Std Dev", "5th Percentile",
      "95th Percentile", "N"]
  summaryNumRows::Int = length(summaryRowNames)

  annualizedMedians::Vector{Float64} = [exp(median(results[:μ])*255)-1.0,
    median(results[:σ])*255^0.5,
    exp(median(results[:μₛ])*255)-1.0,
    median(results[:σₛ])*255^0.5,
    round(median(results[:λ])*255,1)]

  #make the content matrix
  summaryContent::Vector{Vector{String}} =
    [Vector{String}() for i∈ 1:summaryNumRows]

  summaryContent[1] = string.([
    round(mean(results[s]),summaryDecimals) for s::Symbol ∈ summaryCols])

  summaryContent[2] = string.([
    round(median(results[s]),summaryDecimals) for s::Symbol ∈ summaryCols])

  summaryContent[3] = string.([
    round(annualizedMedians[i],summaryDecimals) for i::Int ∈ 1:length(summaryCols)])

  summaryContent[4] = string.([
    round(std(results[s]),summaryDecimals) for s::Symbol ∈ summaryCols])

  summaryContent[5] = string.([
    round(quantile(results[s], 0.05),summaryDecimals) for s::Symbol ∈ summaryCols])

  summaryContent[6] = string.([
    round(quantile(results[s], 0.95),summaryDecimals) for s::Symbol ∈ summaryCols])

  summaryContent[7] = ["$numPoints"]

  tableSummary::String = texTable( "Summary of Data",
    """See Tex File""", #caption
    [["\$\\mu\$", "\$\\sigma\$", "\$\\mu_{s}\$","\$\\sigma_{s}\$", "\$\\lambda\$"]], #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    summaryRowNames, #descRowNames
    summaryContent, #descContent
    Vector{String}(),
    columnSepPt = -5,
    widthDescContent = [[ones(Int, length(summaryCols)) for i∈ 1:(length(summaryContent)-1)];
      [[length(summaryCols)]]]::Vector{Vector{Int}})

  writeTables2File([tableSummary],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "SummaryTable.tex")


  return nothing
end

#enable these to run the program
problemA2Script(refreshData=false)
#getA2Graphs()
end
