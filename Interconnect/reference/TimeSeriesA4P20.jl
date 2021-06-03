module TimeSeriesA4P2

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
  ParallelAccelerator, Measures, Formatting

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

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

#preprocess PPI data
function preProcessPPI()::Void

  ppiDF::DataFrame = readtable("$DATA_PATH\\$PPI_NAME.csv")

  #read and pre-process the data
  dateFormat::String = PPI_DATE_FORMAT
  if :Date ∈ names(ppiDF)
    rename!(ppiDF, :Date, :DATE)
  elseif :_Date ∈ names(ppiDF)
    rename!(ppiDF, :_Date, :DATE)
  elseif :_DATE ∈ names(ppiDF)
    rename!(ppiDF, :_DATE, :DATE)
  end

  ppiDF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(ppiDF[:,:DATE])::DataVector{Date}
  sort!(ppiDF, cols = [:DATE])

  #difference

  ppiDF[:,:Δppi] = similar(ppiDF[:,:ppi])
  ppiDF[:,:Δppi] .= NA
  ppiDF[2:end,:Δppi] = ppiDF[2:end,:ppi] .- ppiDF[1:(end-1),:ppi]

  #log and log difference
  ppiDF[:,:lppi] = log.(ppiDF[:,:ppi])
  ppiDF[:,:Δlppi] = similar(ppiDF[:,:lppi])
  ppiDF[:,:Δlppi] .= NA
  ppiDF[2:end,:Δlppi] = ppiDF[2:end,:lppi] .- ppiDF[1:(end-1),:lppi]

  stream::IOStream = open("$DATA_PATH\\$PPI_NAME.jls", "w")
  serialize(stream, ppiDF)
  close(stream)

  return nothing
end



#to get the df
function getDF(pathStr::String, refreshData::Bool, preProcessFunc::Function)::DataFrame

  #re-process the data if needed or desired
  if !isfile(pathStr) || refreshData
    preProcessFunc()
  end

  #read the binary file
  stream::IOStream = open(pathStr)
  DF::DataFrame = deserialize(stream)
  close(stream)

  return DF
end


#problem specific graphs
function visualizeA4P2(Θ::Vector{ARMA}, ppiDF::DataFrame)::Void

  ACFRange::Range = 1:24
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  outDF = deepcopy(ppiDF)
  completecases!(outDF)

  push!(plotNames, "Summary Graphs")
  push!(plots, Vector{Plot}())
  for series::Symbol ∈ [:ppi, :Δppi, :lppi, :Δlppi]
    push!(plots[end],
      plot(outDF, x=:DATE, y=series,  Geom.line, Geom.point,
        Guide.xlabel("Date"),Guide.ylabel("$series"),
        Guide.title("Summary $series")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "ACF Graphs")
  push!(plots, Vector{Plot}())
  for series::Symbol ∈ [:ppi, :Δppi, :lppi, :Δlppi]
    push!(plots[end],
      plot(x=collect(ACFRange), y=autoρ(Vector{Float64}(outDF[:,series]),ACFRange),
        Geom.line, Geom.point,
        Guide.xlabel("Date"),Guide.ylabel("$series ρ"),
        Guide.title("ACF $series")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "PACF Graphs")
  push!(plots, Vector{Plot}())
  for series::Symbol ∈ [:ppi, :Δppi, :lppi, :Δlppi]
    push!(plots[end],
      plot(x=collect(ACFRange), y=pAutoρ(Vector{Float64}(outDF[:,series]),ACFRange),
        Geom.line, Geom.point,
        Guide.xlabel("Date"),Guide.ylabel("$series ∂ρ"),
        Guide.title("PACF $series")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))


  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing
end

#script for the problem
function A4P2Script(;refreshData::Bool = true)::Void

  targetSeries::Symbol = :Δlppi

  #select the models
  testModel::ARMA =  ARMA(4,4)#ARMA(c=0.1,ϕ=[0.1,-0.2], ϑ=[0.1,0.2])
  models::Vector{ARMA} = [ARMA(1,1), ARMA(4,1), ARMA(4,4), ARMA(8,4), ARMA(8,8)]
  ppiDF::DataFrame = getDF("$DATA_PATH\\$PPI_NAME.jls", refreshData, preProcessPPI)
  #visualizeA4P2(models, ppiDF)

  #get the data in an array
  Y::Vector{Float64} = Vector{Float64}(dropna(ppiDF[:,:Δlppi]))
  generalizeARMA!(testModel)
  narrowARMA!(testModel)
  println("test log like: $(ARMALogLike(testModel, Y))")diag
  ARMAEstimate(testModel, Y)

  return nothing
end


@time begin
  A4P2Script(refreshData = true)
end

end
