module TimeSeriesA1P2

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

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTIO, CTReg, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  ParallelAccelerator, Measures, Formatting

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const YAHOO_DATE_FORMAT = "m/d/yyyy"
const SP500_NAME = "GSPC"
const DBV_NAME = "DBV"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
const PlotContainer = Union{Plot,Gadfly.Compose.Context} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]


function getNormalData(N::Int; μ::Float64=0., σ²::Float64=1.)::Vector{Float64}
  dist::Distribution = Normal(μ,σ²^0.5)

  return rand(dist, N)
end

function getJumpData(N::Int; μ::Float64=0., σ²::Float64=1.,
  p::Float64=0.5, μⱼ::Float64 = 0.0, σ²ⱼ::Float64 = 1.)::Vector{Float64}

  normDist::Distribution = Normal(μ,σ²^0.5)
  bernDist::Distribution = Bernoulli(p)
  jumpDist::Distribution = Normal(μⱼ,σ²ⱼ^0.5)

  println("kurtosis: ", kurtosis(rand(normDist, N) .+ rand(bernDist, N) .* rand(jumpDist, N)))

  return rand(normDist, N) .+ rand(bernDist, N) .* rand(jumpDist, N)
end

jumpVar(; μ::Float64=0., σ²::Float64=1.,
  p::Float64=0.5, μⱼ::Float64 = 0.0, σ²ⱼ::Float64 = 1.) = σ²+p*(-1*(p-1)*μⱼ^2+σ²ⱼ)

jumpSkew(; μ::Float64=0., σ²::Float64=1.,
  p::Float64=0.5, μⱼ::Float64 = 0.0, σ²ⱼ::Float64 = 1.) =
    ((-1+p)*p*μⱼ*((-1+2*p)*μⱼ^2-3*σ²ⱼ))/(σ²+p*(-1*(-1+p)*μⱼ^2+σ²ⱼ))^(1.5)

jumpKurt(; μ::Float64=0., σ²::Float64=1.,
  p::Float64=0.5, μⱼ::Float64 = 0.0, σ²ⱼ::Float64 = 1.) =
    ((-1+p)*p*((-1-6*(-1+p)*p)*μⱼ^4+6*(-1+2*p)*μⱼ^2*σ²ⱼ-3*σ²ⱼ^2))/(σ²+p*(-1*(p-1)*μⱼ^2+σ²ⱼ))^2


formNum(x::T where T<:Real)::String = format(NUM_FORMAT, x)

#Calculates the summary statistics
function summaryA1P1(normalReturns::Vector{Float64}, jumpReturns::Vector{Float64})
  α::Float64 = 0.05

  #label the table
  colNames::Vector{Vector{String}} =
    [["mean", "\$\\sigma\$", "skewness", "excess kurtosis"]]
  numCols = length(colNames[end])

  rowNames::Vector{String} = ["Normal", "Normal Jumps"]
  numRows::Int = length(rowNames)

  NNormal::Int = length(normalReturns)
  NJumps::Int = length(jumpReturns)

  data::Vector{Vector{Float64}} = [normalReturns, jumpReturns]

  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #run the tests and record the results
  for r::Int ∈ 1:numRows
    content[r][1] = formNum(mean(data[r]))
    content[r][2] = formNum(std(data[r]))
    content[r][3] = formNum(skewness(data[r]))
    content[r][4] = formNum(kurtosis(data[r]))
  end

  #write the table
  summaryTable::String = texTable( "Normal and Normal Jump Distributions",
    """See Tex File""", #caption
    colNames, #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    rowNames, #descRowNames
    content, #descContent
    Vector{String}(), #notes
  )

  return summaryTable
end

#plots the descriptive graphs
function visualizeA1P1(normalReturns::Vector{Float64}, jumpReturns::Vector{Float64})

  #convenience containers
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  DF::DataFrame = DataFrame(N = 1:length(normalReturns),
    normalReturns = normalReturns, jumpReturns = jumpReturns)
  outDF::DataFrame = melt(DF, :N) #put data in long form

  push!(plotNames, "Histogram P1")
  push!(plots, plot(outDF, x=:value, color=:variable,
      Guide.title("Histograms"),
      Guide.xlabel("Log Return"),Guide.ylabel("Density"),
      Coord.Cartesian(xmin=-0.2, xmax=0.2),
      Geom.histogram(bincount=20, density=true, position=:dodge)))

  push!(plotNames, "Time Series P1")
  push!(plots, vstack(
        plot(DF, y=:normalReturns, x=:N, Geom.line,
          Guide.xlabel("Year"),Guide.ylabel("Return"), Guide.title("Normal")),
        plot(DF, y=:jumpReturns, x=:N, Geom.line,
          Guide.xlabel("Year"),Guide.ylabel("Return"), Guide.title("Normal with Normal Jumps"))))


  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
  end

end



function A1P1Script(; refreshData::Bool = false)
  N::Int = 600
  normalReturns = getNormalData(N, μ=0.008, σ²=0.063^2)
  jumpReturns = getJumpData(N, μ=0.012, σ²=0.05^2, p=0.15, μⱼ=-0.03, σ²ⱼ=0.1^2)
  visualizeA1P1(normalReturns, jumpReturns)
  summaryTable::String = summaryA1P1(normalReturns, jumpReturns)

  println("Jump Var: $(jumpVar(μ=0.012, σ²=0.05^2, p=0.15, μⱼ=-0.03, σ²ⱼ=0.1^2))")
  println("Jump Skew: $(jumpSkew(μ=0.012, σ²=0.05^2, p=0.15, μⱼ=-0.03, σ²ⱼ=0.1^2))")
  println("Jump Kurt: $(jumpKurt(μ=0.012, σ²=0.05^2, p=0.15, μⱼ=-0.03, σ²ⱼ=0.1^2))")

  writeTables2File([summaryTable],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "A1P1 Summary.tex")

end

@time begin
#  A1P1Script()
end


end
