module TimeSeriesA5P2

#Use this to print
#=
using Weave
codeName = "TimeSeriesA5P2"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD, ForwardDiff,
  Gadfly, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  Measures, Formatting

importall CTIO, CTReg, CTStat, CTVAR

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const YAHOO_DATE_FORMAT = "m/d/yyyy"
const FRED_DATE_FORMAT = "m/d/yyyy"
const PPI_DATE_FORMAT = "m/d/yyyy"
const NIPA_DATE_FORMAT = "m/d/yyyy"
const SP500_NAME = "GSPC"
const DBV_NAME = "DBV"
const FRED_NAME = "CPI"
const PPI_NAME = "PPI"
const NIPA_NAME = "NIPA_GDP_UNRATE"

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
function preProcessNIPA(dropTop::Bool = true)::Void

  nipaDF::DataFrame = readtable("$DATA_PATH\\$NIPA_NAME.csv")

  #read and pre-process the data
  dateFormat::String = NIPA_DATE_FORMAT
  if :Date ∈ names(nipaDF)
    rename!(nipaDF, :Date, :DATE)
  elseif :_Date ∈ names(nipaDF)
    rename!(nipaDF, :_Date, :DATE)
  elseif :_DATE ∈ names(nipaDF)
    rename!(nipaDF, :_DATE, :DATE)
  end

  nipaDF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(nipaDF[:,:DATE])::DataVector{Date}
  sort!(nipaDF, cols = [:DATE])

  #difference
  #println("got here")
  #get GDP growth
  nipaDF[:,:gdp] = Vector{Float64}(nipaDF[:,:gdp]) #convert to float from int
  nipaDF[:,:Δlgdp] = similar(nipaDF[:,:gdp])
  nipaDF[2:end,:Δlgdp] = (log).(nipaDF[2:end,:gdp]) .- (log).(nipaDF[1:(end-1),:gdp])
  nipaDF[:,:unrate] /= 100.

  #get change in employement
  nipaDF[:,:Δunrate] = similar(nipaDF[:,:unrate])
  nipaDF[2:end,:Δunrate] = (nipaDF[2:end,:unrate]) .- (nipaDF[1:(end-1),:unrate])


  stream::IOStream = open("$DATA_PATH\\$NIPA_NAME.jls", "w")
  serialize(stream, nipaDF[2:end,:])
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

#Calculates the summary statistics
function summaryA5P2(nipaDF::DataFrame, LHS::Vector{Symbol}, RHS::Vector{Symbol}, est::VAREstimate;
    title="Estimate")
  #label the table
  K::Int = length(LHS)
  colNames::Vector{Vector{String}} =
    [["\$\\phi_0\$"; (String).(RHS); "N"; "\$R^2\$"; "F"]]
  numCols = length(colNames[end])

  rowNames::Vector{String} = ((i::Int)->isodd(i)?(String(LHS[(i+1)÷2])):"").(1:(2K))
  numRows::Int = length(rowNames)
  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #run the tests and record the results
  for k::Int ∈ 1:K
    r::Int = 2(k-1) + 1 #pointer to the beginning of each double row
    content[r][1] = n2s(est.var.ϕ₀[k])
    content[r+1][1] = est.coeffErrors[k] ≠ nothing ? "($(n2s(est.coeffErrors[k][1])))" : ""
    for c::Int ∈ 2:(numCols-3)
      content[r][c] = n2s(est.var.ϕ₁[k, c-1])
      content[r+1][c] = est.coeffErrors[k] ≠ nothing ? "($(n2s(est.coeffErrors[k][c])))" : ""
    end

    content[r][end-2] = est.N[k] ≠ nothing ? "$(n2s(est.N[k]))" : ""
    content[r+1][end-2] = ""
    content[r][end-1] = est.MSE[k] ≠ nothing ? "$(n2s(1.-est.MSE[k]/est.MST[k]))" : ""
    content[r+1][end-1] = ""
    content[r][end] = est.MSE[k] ≠ nothing ? "$(n2s(est.MSR[k]/est.MSE[k]))" : ""
    content[r+1][end] = ""

  end



  #write the table
  summaryTable::String = texTable( title,
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
#problem specific graphs
function visualizeA5P2(nipaDF::DataFrame, predictions::Vector{VARPrediction})::Void

  ACFRange::Range = 1:24
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  outDF = deepcopy(nipaDF)
  completecases!(outDF)

  push!(plotNames, "NIPA Summary Graphs")
  push!(plots, Vector{Plot}())
  for series::Symbol ∈ [:Δlgdp, :unrate]
    push!(plots[end],
      plot(outDF, x=:DATE, y=series,  Geom.line, Geom.point,
        Guide.xlabel("Date"),Guide.ylabel("$series"),
        Guide.title("Summary $series")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "UnrestΔlGDP")
  push!(plots, plot(
    layer(x=predictions[1].dates, y=predictions[1].Y[1,:],  Geom.line,
      Theme(default_color=COLORS[1])),
    layer(x=predictions[1].dates, y=predictions[1].YHat[1,:],  Geom.line,
      Theme(default_color=COLORS[2])),
    layer(x=predictions[2].dates, y=predictions[2].YHat[1,:],  Geom.line,
      Theme(default_color=COLORS[3])),
    Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction", "t+4 Prediction"], COLORS[1:3]),
    Guide.xlabel("Date"),Guide.ylabel("ΔlGDP"),
    Guide.title("Unrestricted VAR Prediction of GDP")))

    push!(plotNames, "UnrestUNRATE")
    push!(plots, plot(
      layer(x=predictions[1].dates, y=predictions[1].Y[2,:],  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(x=predictions[1].dates, y=predictions[1].YHat[2,:],  Geom.line,
        Theme(default_color=COLORS[2])),
      layer(x=predictions[2].dates, y=predictions[2].YHat[2,:],  Geom.line,
        Theme(default_color=COLORS[3])),
      Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction", "t+4 Prediction"], COLORS[1:3]),
      Guide.xlabel("Date"),Guide.ylabel("Unemployment"),
      Guide.title("Unrestricted VAR Prediction of Unemployment")))

    push!(plotNames, "RestrictedΔlGDP")
    push!(plots, plot(
      layer(x=predictions[3].dates, y=predictions[3].Y[1,:],  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(x=predictions[3].dates, y=predictions[3].YHat[1,:],  Geom.line,
        Theme(default_color=COLORS[2])),
      layer(x=predictions[4].dates, y=predictions[4].YHat[1,:],  Geom.line,
        Theme(default_color=COLORS[3])),
      Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction", "t+4 Prediction"], COLORS[1:3]),
      Guide.xlabel("Date"),Guide.ylabel("ΔlGDP"),
      Guide.title("Restricted VAR Prediction of GDP")))


  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing
end

function lagDF!(DF::DataFrame, target::Symbol, lags::Int = 1; name= Symbol("$(target)L$lags"))::DataFrame
  DF[:,name] = similar(DF[:,target])
  DF[:,name] .= NA
  DF[(lags+1):end,name] = DF[1:(end-lags),target]

  return DF
end

#script for the problem
function A5P2Script(;refreshData::Bool = true)::Void


  LHS::Vector{Symbol} = [:Δlgdp, :unrate, :ΔlgdpL1, :unrateL1]
  RHS::Vector{Symbol} = [:ΔlgdpL1, :unrateL1, :ΔlgdpL2, :unrateL2]
  nipaDF::DataFrame = getDF("$DATA_PATH\\$NIPA_NAME.jls", refreshData, preProcessNIPA)

  for s ∈ LHS[1:2]
    lagDF!(nipaDF,s, 1)
    lagDF!(nipaDF,s, 2)
  end

  ######first, estimate the unrestricted model
  est::VAREstimate = VAREstimate(nipaDF, LHS, RHS)
  tables::Vector{String} = Vector{String}()
  push!(tables, summaryA5P2(nipaDF, LHS, RHS, est, title="Unrestricted VAR 2"))

  #get the predicted values
  validXY::Vector{Bool} = completecases(nipaDF[RHS]) .& completecases(nipaDF[LHS])
  X::Matrix{Float64} = Matrix{Float64}(nipaDF[validXY, RHS])'
  Y::Matrix{Float64} = Matrix{Float64}(nipaDF[validXY, LHS])'
  dates::Vector{Date} = nipaDF[validXY, :DATE]

  #container for holding all predictions
  predictions::Vector{VARPrediction} = Vector{VARPrediction}()

  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 1))
  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 4))

  ##################Now do the restricted model
  LHS = [:Δlgdp, :Δunrate, :ΔlgdpL1]
  RHS = [:ΔlgdpL1, :ΔunrateL1, :ΔlgdpL2]
  lagDF!(nipaDF,:Δunrate, 1)

  #estimate the VAR
  est = VAREstimate(nipaDF, LHS, RHS)
  push!(tables, summaryA5P2(nipaDF, LHS, RHS, est, title="Restricted VAR 2"))

  #get the predicted values
  validXY = completecases(nipaDF[RHS]) .& completecases(nipaDF[LHS])
  X = Matrix{Float64}(nipaDF[validXY, RHS])'
  Y = Matrix{Float64}(nipaDF[validXY, LHS])'
  dates= nipaDF[validXY, :DATE]

  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 1))
  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 4))

  #make the grpahs
  visualizeA5P2(nipaDF, predictions)

  #write out the tables
  writeTables2File(tables,
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "TimeSeries A5P2.tex")

  return nothing
end


@time begin
# Uncomment to run
#  A5P2Script(refreshData = true)
end

end
