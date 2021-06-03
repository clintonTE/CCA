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


#does common tasks on a Yahoo dataframe
function processYahooDF!(df::DataFrame)::DataFrame

  #read and pre-process the data
  dateFormat::String = YAHOO_DATE_FORMAT
  if :Date ∈ names(df)
    rename!(df, :Date, :DATE)
  elseif :_Date ∈ names(df)
    rename!(df, :_Date, :DATE)
  end

  df[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(df[:,:DATE])::DataVector{Date}
  sort!(df, cols = [:DATE])
  df[:,:returns] = similar(df[:,:Adj_Close])
  df[2:end,:returns] = log.(df[2:end,:Adj_Close] ./  df[1:(end-1),:Adj_Close])

  return df
end

#pre-process the s&p500 data
function preProcessSP500()::Void

  #read and pre-process the data
  sp500DF::DataFrame = readtable("$DATA_PATH\\$SP500_NAME.csv")
  sp500DF = processYahooDF!(sp500DF)

  #write to JLS binary file
  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls", "w")
  serialize(stream, sp500DF)
  close(stream)

  return nothing
end

function preProcessDBV()::Void

  #read and pre-process the data
  dbvDF::DataFrame = readtable("$DATA_PATH\\$DBV_NAME.csv")
  dbvDF = processYahooDF!(dbvDF)

  #write to JLS binary file
  stream::IOStream = open("$DATA_PATH\\$DBV_NAME.jls", "w")
  serialize(stream, dbvDF)
  close(stream)

  return nothing
end

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

#plots the descriptive graphs
function visualizeA1P2(df::DataFrame)

  #convenience containers
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  outDF::DataFrame = meltdf(df, :DATE) #put data in long form
  completecases!(outDF)

  push!(plotNames, "Histogram")
  push!(plots, plot(outDF, x=:value, color=:variable,
      Guide.title("Histograms"),
      Guide.xlabel("Log Return"),Guide.ylabel("Density"),
      Coord.Cartesian(xmin=-0.12, xmax=0.12),
      Geom.histogram(bincount=20, density=true, position=:dodge)))

  push!(plotNames, "Time Series")
  push!(plots, vstack(
        plot(outDF[outDF[:,:variable].==:sp500,:], y=:value, x=:DATE, Geom.line,
          Guide.xlabel("Year"),Guide.ylabel("Log Return"), Guide.title("SP500")),
        plot(outDF[outDF[:,:variable].==:dbv,:], y=:value, x=:DATE, Geom.line,
          Guide.xlabel("Year"),Guide.ylabel("Log Return"), Guide.title("DBV"))))


  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
  end

end

#calculates the test statistic and T value of a data vector
skewnessStat(data::Vector{Float64})::Float64 = skewness(data)/√(6./length(data))
kurtosisStat(data::Vector{Float64})::Float64 = kurtosis(data)/√(24./length(data))
JBStat(data::Vector{Float64})::Float64 =
  skewness(data)^2/(6./length(data)) + kurtosis(data)^2/(24./length(data))


mutable struct UnivariateTest
  statistic::Function
  distribution::Distribution
  α::Float64
  H₀::Float64
  testType::Symbol
end

#convenience constructor for univariate test
UnivariateTest(statistic::Function, distribution::Distribution, α::Float64;
  H₀::Float64=0.0, testType::Symbol=:upper)::UnivariateTest =
  UnivariateTest(statistic, distribution, α, H₀, testType)

formNum(x::T where T<:Real)::String = format(NUM_FORMAT, x)


function testUni(U::UnivariateTest, data::Vector{Float64}; verbose::Bool = true)::NTuple{2,Float64}
  θ::Float64 = U.statistic(data .- U.H₀)

  p::Float64 = cdf(U.distribution, θ)

  #adjust the p value for the specific test type
  if U.testType == :upper
    p = 1.0 - p
  elseif  U.testType == :lower
  elseif U.testType == :twoTailed
    p = p > 0.5 ? (1-p)*2 : 2*p
  else
    warn("Test type not found in univariate test")
  end


  return (θ, p)
end

#script for doing the hypothesis test
function hypothesisA1P2(df::DataFrame)
  α::Float64 = 0.05

  #label the table
  colNames::Vector{Vector{String}} =
    [["SP500", "DBV"],
    ["statistic", "p-value", "reject \$H_0\$?", "statistic", "p-value", "reject \$H_0\$?"]]
  numCols = length(colNames[end])

  rowNames::Vector{String} = ["Skewness", "Kurtosis", "Normality"]
  numRows::Int = length(rowNames)

  #form the data
  sp500::Vector{Float64} = dropna(df[:,:dbv])
  dbv::Vector{Float64} = dropna(df[:,:dbv])

  NSP500::Int = length(sp500)
  NDBV::Int = length(dbv)

  #set the test parameters
  sp500Tests::Vector{UnivariateTest} =
    [UnivariateTest(skewnessStat, TDist(NSP500-1), α, testType=:twoTailed),
    UnivariateTest(kurtosisStat, TDist(NSP500-1), α, testType=:twoTailed),
    UnivariateTest(JBStat, TDist(NSP500), α, testType=:upper)]

  dbvTests::Vector{UnivariateTest} =
    [UnivariateTest(skewnessStat, TDist(NDBV-1), α, testType=:twoTailed),
    UnivariateTest(kurtosisStat, TDist(NDBV-1), α, testType=:twoTailed),
    UnivariateTest(JBStat, Chisq(2), α, testType=:upper)]

  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #run the tests and record the results
  for r::Int ∈ 1:numRows
    sp500θ::Float64, sp500P::Float64 = testUni(sp500Tests[r], sp500)
    content[r][1] = formNum(sp500θ)
    content[r][2] = formNum(sp500P)
    content[r][3] = sp500P < α ? "X" : ""

    dbvθ::Float64, dbvP::Float64 = testUni(dbvTests[r], dbv)
    content[r][4] = formNum(dbvθ)
    content[r][5] = formNum(dbvP)
    content[r][6] = dbvP < α ? "X" : ""
  end

  #write the table
  hypothesisTable::String = texTable( "Summary of ST and LT Interest Yields",
    """See Tex File""", #caption
    colNames, #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    rowNames, #descRowNames
    content, #descContent
    Vector{String}(), #notes
    widthColNames = [[3, 3], ones(Int,6)],
    alignmentColNames = [["c", "c"], ["r" for i::Int ∈ 1:length(6)]]
  )

  return hypothesisTable
end

#script for doing the regressions
function regressionA1P2(df::DataFrame, SE::Function; caption::String="TBD", title::String="TBD")


  #model specification
  XSpec::Symbol = parse("sp500")
  XName::Vector{Symbol} = [:intercept; :sp500]
  YSpec::Symbol = :dbv
  colNames::Vector{String} = ["DBV (1)"]

  #run the regression
  reg::Vector{CTLM} = [CTLM(df, XSpec, YSpec, XNames=XName, YName = :dbv)]

  #get some descriptive stats
  descRowNames::Vector{String} = ["N", "\$R^{2}\$"]
  descContent::Vector{Vector{String}} = [["$(reg[1].N)"], ["$(round(getR(reg[1])^2.0,2))"]]

  TableText::String = texTable(reg, SE, [:sp500; :intercept],
    titleCaption = title,
    #colNames = [["($i)" for i∈1:length(models)]],
    colNames = [colNames],
    contentRowNames = ["SP500", "Intercept"],
    descRowNames = descRowNames,
    descContent = descContent,
    stars=true,
    decimalDigits = 4,
    caption=caption,
    starStrings = OVERRIDE_STAR_STRINGS)

end


function A1P2Script(; refreshData::Bool = false)
  sp500DF::DataFrame = getDF("$DATA_PATH\\$SP500_NAME.jls", refreshData, preProcessSP500)
  dbvDF::DataFrame = getDF("$DATA_PATH\\$DBV_NAME.jls", refreshData, preProcessDBV)

  #merge the dataframes
  rename!(sp500DF, :returns, :sp500)
  rename!(dbvDF, :returns, :dbv)
  df::DataFrame = join(sp500DF[:,[:DATE, :sp500]],dbvDF[:,[:DATE, :dbv]], on=[:DATE], kind=:inner)

  #make the graphs
  visualizeA1P2(df)

  #do the hypothesis tast
  hypothesisTable::String = hypothesisA1P2(df)
  homoskedasticTable::String = regressionA1P2(df, getHomoskedΣ!, title="DBV vs S\\\&P500",
    caption="Homoskedastic Standard Errors")
  modWhiteTable::String = regressionA1P2(df, getModWhiteΣ!, title="DBV vs S\\\&P500",
    caption="Modified White Standard Errors")
  whiteTable::String = regressionA1P2(df, getWhiteΣ!, title="DBV vs S\\\&P500",
    caption="White Standard Errors")
  #=whiteTableSlow::String = regressionA1P2(df, getWhiteΣSlow, title="DBV vs S\\\&P500",
    caption="White Standard Errors (Test)")=#


  #record the tables
  writeTables2File([hypothesisTable, homoskedasticTable, modWhiteTable, whiteTable],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "OutputTables.tex")

end

@time begin
  A1P2Script()
end


end
