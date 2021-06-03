module TimeSeriesA6

#Use this to print
#=
using Weave
codeName = "TimeSeriesA6"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD, ForwardDiff,
  Gadfly, JuMP, NLopt, HCubature, StaticArrays, Measures, Formatting, CSV


importall CTIO, CTReg, CTStat, CTVAR, CTNumerical

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
const A6_DATA_NAME = "MktRet_DP_TermSpread"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

#preprocess rate data
function preProcessA6()::Void
  rateDF::DataFrame = CSV.read("$DATA_PATH\\$A6_DATA_NAME.csv")

  #showcols(rateDF)
  rename!(rateDF, [:MktExRet, :Mkt_DP,	:y10minFedFunds], [:re, :dp, :y10e])


  stream::IOStream = open("$DATA_PATH\\$A6_DATA_NAME.jls", "w")
  serialize(stream, rateDF)
  close(stream)

  return nothing
end

#Calculates the summary statistics
function summaryA6(rateDF, cols::Vector{Symbol} = [:re, :dp, :y10e],
      rowNames::Vector{String} = ["mean", "\$\\sigma\$", "\$\\rho_{t-1}\$", "Half Life", "N"]
    )
  #label the table
  colNames::Vector{Vector{String}} =
    [((s::Symbol)->replace(String(s),"_","")).(cols)]
  numCols = length(colNames[end])

  numRows::Int = length(rowNames)
  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #get the statistics and record the results
  for c::Int ∈ 1:numCols
    ρ₁::Float64 = autoρ(rateDF[:, cols[c]], 1)

    content[1][c] = n2s(mean(rateDF[:, cols[c]]))
    content[2][c] = n2s(std(rateDF[:, cols[c]]))
    content[3][c] = n2s(ρ₁)
    content[4][c] = n2s(log(0.5)/log(ρ₁))
    content[5][c] = n2s(length(rateDF[:, cols[c]]))
  end

  #write the table
  summaryTable::String = texTable( "Summary of Data Series",
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

#VAR table
function VARTableA6(rateDF::DataFrame, LHS::Vector{Symbol}, RHS::Vector{Symbol}, est::VAREstimate;
    title="Estimate")
  #label the table

  K::Int = length(LHS)
  colNames::Vector{Vector{String}} =
    [["\$\\phi_0\$"; (String).(RHS); "N"; "\$R^2\$";"\$R^2 (Adj)\$"; "F"]]
  numCols = length(colNames[end])

  rowNames::Vector{String} = ((i::Int)->isodd(i)?(String(LHS[(i+1)÷2])):"").(1:(2K))
  numRows::Int = length(rowNames)
  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #run the tests and record the results
  for k::Int ∈ 1:K
    r::Int = 2(k-1) + 1 #pointer to the beginning of each double row

    #need these for the adjusted R2 and F statistics
    MSE::Float64 = est.SSE[k] / (est.N[k]-K-1.)
    MSR::Float64 = est.SSR[k] / K
    MST::Float64 = est.SST[k] /(est.N[k] - 1.)


    content[r][1] = n2s(est.var.ϕ₀[k])
    content[r+1][1] = est.coeffErrors[k] ≠ nothing ? "($(n2s(est.coeffErrors[k][1])))" : ""
    for c::Int ∈ 2:(numCols-4)
      content[r][c] = n2s(est.var.ϕ₁[k, c-1])
      content[r+1][c] = est.coeffErrors[k] ≠ nothing ? "($(n2s(est.coeffErrors[k][c])))" : ""
    end

    content[r][end-3] = est.N[k] ≠ nothing ? "$(n2s(est.N[k]))" : ""
    content[r+1][end-3] = ""
    content[r][end-2] = est.SSE[k] ≠ nothing ? "$(n2s(est.SSR[k]/est.SST[k]))" : ""
    content[r+1][end-2] = ""
    content[r][end-1] = est.SSE[k] ≠ nothing ? "$(n2s(1.-MSE/MST))" : ""
    content[r+1][end-1] = ""
    content[r][end] = est.SSE[k] ≠ nothing ? "$(n2s(MSR/MSE))" : ""
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

#make the problem-specific graphs
function visualizeA6(rateDF::DataFrame, predictions::Vector{VARPrediction},
      responses::Vector{Matrix{Float64}})::Void

  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  outDF = deepcopy(rateDF)
  dropmissing!(outDF)


  push!(plotNames, "summaryRatesAIO")
  push!(plots, plot(
    layer(outDF, x=:DATE, y=:y10e,  Geom.line, Theme(default_color=COLORS[1])),
    layer(outDF, x=:DATE, y=:dp,  Geom.line, Theme(default_color=COLORS[2])),
    layer(outDF, x=:DATE, y=:re,  Geom.line, Theme(default_color=COLORS[3])),
    Guide.manual_color_key("Legend", ["Excess Yield", "dp",
      "Excess Return"], COLORS[1:3]),
    Guide.xlabel("Date"),Guide.ylabel("Rate"),
    Coord.Cartesian(xmin=1950., xmax=2018.),
    Guide.title("Excess 10-year yield, dividend yield, and excess market return")))


  push!(plotNames, "10Y Excess Yield Predictions")
  push!(plots, plot(
    layer(x=predictions[1].dates, y=predictions[1].Y[1,:],  Geom.line,
      Theme(default_color=COLORS[1])),
    layer(x=predictions[1].dates, y=predictions[1].YHat[1,:],  Geom.line,
      Theme(default_color=COLORS[2])),
    layer(x=predictions[2].dates, y=predictions[2].YHat[1,:],  Geom.line,
      Theme(default_color=COLORS[3])),
    layer(x=predictions[3].dates, y=predictions[3].YHat[1,:],  Geom.line,
      Theme(default_color=COLORS[4])),
    Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction",
      "t+4 Prediction", "t+20 Prediction"], COLORS[1:4]),
    Coord.Cartesian(xmin=1950., xmax=2018.),
    Guide.xlabel("Date"),Guide.ylabel("Excess Yield"),
    Guide.title("Excess yield predictions")))

  push!(plotNames, "Dividend yield predictions")
  push!(plots, plot(
    layer(x=predictions[1].dates, y=predictions[1].Y[2,:],  Geom.line,
      Theme(default_color=COLORS[1])),
    layer(x=predictions[1].dates, y=predictions[1].YHat[2,:],  Geom.line,
      Theme(default_color=COLORS[2])),
    layer(x=predictions[2].dates, y=predictions[2].YHat[2,:],  Geom.line,
      Theme(default_color=COLORS[3])),
    layer(x=predictions[3].dates, y=predictions[3].YHat[2,:],  Geom.line,
      Theme(default_color=COLORS[4])),
    Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction",
      "t+4 Prediction", "t+20 Prediction"], COLORS[1:4]),
    Guide.xlabel("Date"),Guide.ylabel("Dividend Yield"),
    Coord.Cartesian(xmin=1950., xmax=2018.),
    Guide.title("Dividend yield predictions")))

  push!(plotNames, "Excess Market Predictions")
  push!(plots, plot(
    layer(x=predictions[1].dates, y=predictions[1].Y[3,:],  Geom.line,
      Theme(default_color=COLORS[1])),
    layer(x=predictions[1].dates, y=predictions[1].YHat[3,:],  Geom.line,
      Theme(default_color=COLORS[2])),
    layer(x=predictions[2].dates, y=predictions[2].YHat[3,:],  Geom.line,
      Theme(default_color=COLORS[3])),
    layer(x=predictions[3].dates, y=predictions[3].YHat[3,:],  Geom.line,
      Theme(default_color=COLORS[4])),
    Guide.manual_color_key("Legend", ["Actual", "t+1 Prediction",
      "t+4 Prediction", "t+20 Prediction"], COLORS[1:4]),
    Guide.xlabel("Date"),Guide.ylabel("Rate of Return"),
    Coord.Cartesian(xmin=1950., xmax=2018.),
    Guide.title("Excess market return predictions")))

    #Impulse responses
    xLags::Vector{Int} = collect(0:(size(responses[1],2)-1))

    push!(plotNames, "Yield Impulse Response")
    push!(plots, plot(
      layer(x=xLags, y=responses[1][1,:],  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(x=xLags, y=responses[1][2,:],  Geom.line,
        Theme(default_color=COLORS[2])),
      layer(x=xLags, y=responses[1][3,:],  Geom.line,
        Theme(default_color=COLORS[3])),
      Guide.manual_color_key("Legend", ["Excess Yield", "dp",
        "Excess Return"], COLORS[1:3]),
      Guide.xlabel("Lag"),Guide.ylabel("Rate"), Guide.title("Yield Impulse Response")))

    push!(plotNames, "Dividend Yield Impulse Response")
    push!(plots, plot(
      layer(x=xLags, y=responses[2][1,:],  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(x=xLags, y=responses[2][2,:],  Geom.line,
        Theme(default_color=COLORS[2])),
      layer(x=xLags, y=responses[2][3,:],  Geom.line,
        Theme(default_color=COLORS[3])),
      Guide.manual_color_key("Legend", ["Excess Yield", "dp",
        "Excess Return"], COLORS[1:3]),
      Guide.xlabel("Lag"),Guide.ylabel("Rate"), Guide.title("Dividend Yield Impulse Response")))

    push!(plotNames, "Market Return Impulse Response")
    push!(plots, plot(
      layer(x=xLags, y=responses[3][1,:],  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(x=xLags, y=responses[3][2,:],  Geom.line,
        Theme(default_color=COLORS[2])),
      layer(x=xLags, y=responses[3][3,:],  Geom.line,
        Theme(default_color=COLORS[3])),
      Guide.manual_color_key("Legend", ["Excess Yield", "dp",
        "Excess Return"], COLORS[1:3]),
      Guide.xlabel("Lag"),Guide.ylabel("Rate"), Guide.title("Excess Return Impulse")))

    push!(plotNames, "DRₜ-EDRₜ and dpₜ-Edpₜ")
    push!(plots, plot(
      layer(outDF, x=:DATE, y=:rDREDR,  Geom.line,
        Theme(default_color=COLORS[1])),
      layer(outDF, x=:DATE, y=:ldpEdp,  Geom.line,
        Theme(default_color=COLORS[2])),
      Guide.manual_color_key("Legend", ["DRₜ-EDRₜ", "dpₜ"], COLORS[1:2]),
      Guide.xlabel("Date"),Guide.ylabel("Rate"), Guide.title("DRₜ-EDRₜ and dpₜ-Edpₜ")))

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

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

function A6StationarityTest(var::VAR; verbose::Bool = true)::Void

  #get the eigenvectors
  λFact::Base.LinAlg.Eigen = eigfact(var.ϕ₁)
  λ::Vector{T where T<:Number} = λFact[:values]

  stationaryTests::BitArray = (abs2).(λ) .< 1.

  if verbose
    for i::Int ∈ 1:length(stationaryTests)
      println("$(stationaryTests[i]?"Passed":"Failed"): λ=$(n2s(λ[i])) |λ|=$(n2s(abs(λ[i])))")
    end
  end

  println("The VAR is $(prod(stationaryTests)?"Stationary":"Not Stationary")")

  return nothing
end


function CSA6(var::VAR, Y::Matrix{Float64}; ρ::Float64=0.96)
  K::Int = size(Y,1)
  T::Int = size(Y,2)
  I::Matrix{Float64} = eye(K)

  pre::Matrix{Float64} = ρ*(I .- ρ*var.ϕ₁)\I
  μ::Vector{Float64} = (I .- var.ϕ₁)\I*var.ϕ₀

  M::Matrix{Float64} = pre * (Y .- μ)

  return M
end


#script for the problem
function A6Script(;refreshData::Bool = true)::Void

  lags::Int = 20
  LHS::Vector{Symbol} = [:y10e, :dp, :re]
  RHS::Vector{Symbol} = [:y10eL1, :dpL1, :reL1]

  #get the dataframe
  rateDF::DataFrame = getDF("$DATA_PATH\\$A6_DATA_NAME.jls", refreshData, preProcessA6)
  T::Int = size(rateDF, 1)

  #make the summary tables
  tables::Vector{String} = Vector{String}()
  push!(tables, summaryA6(rateDF))

  lagDF!(rateDF, LHS) #make the lags

  #esitmate the VAR
  est = VAREstimate(rateDF, LHS, RHS, errors=getNeweyWestFunc(40))
  push!(tables, VARTableA6(rateDF, LHS, RHS, est, title="VAR 1 for Rates"))

  #make the predictions
  validXY::Vector{Bool} = completecases(rateDF[RHS]) .& completecases(rateDF[LHS])
  X::Matrix{Float64} = Matrix{Float64}(rateDF[validXY, RHS])'
  Y::Matrix{Float64} = Matrix{Float64}(rateDF[validXY, LHS])'
  dates::Vector{Float64} = rateDF[validXY, :DATE]

  #container for holding all predictions
  predictions::Vector{VARPrediction} = Vector{VARPrediction}()
  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 1))
  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 4))
  push!(predictions, VARPrediction(est.var, dates, X, Y, τ = 20))

  responses::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}()
  push!(responses, propagateImpulse(est.var, [1.,0.,0.]))
  push!(responses, propagateImpulse(est.var, [0.,1.,0.]))
  push!(responses, propagateImpulse(est.var, [0.,0.,1.]))

  #get the DR-EDR quantity and store in the dataframe
  rateDF[:, :rDREDR] = Vector{Union{Float64,Missing}}(T)
  rateDF[:, :rDREDR] = missing
  rateDF[validXY, :rDREDR] = CSA6(est.var, Y)[3,:]
  rateDF[:,:ldpEdp] = (log).(rateDF[:,:dp])

  rateDF[:,:ldpEdp] -= mean(skipmissing(rateDF[:,:ldpEdp]))

  #make the graphs

  visualizeA6(rateDF, predictions, responses)

  A6StationarityTest(est.var)

  println("Predicted return volatility: ", ((diag(varΣ(est.var)).*4).^0.5)[3])
  println("Historical Vol of Return: $(std(Y[3,:])*2)")
  println("Conditional volatility of return: $((est.var.L*est.var.L')[3,3]^0.5*2)")

  #write out the tables
  writeTables2File(tables,
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "TimeSeries A6.tex")



  return nothing
end

@time begin
  A6Script(refreshData = true)
end

end
