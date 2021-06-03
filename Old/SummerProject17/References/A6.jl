module A6


#Use this to print
#=
weave(Pkg.dir("$(pwd())\\A6.jl"),
  informat="script",
  out_path = "$(pwd())\\A6.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  ParallelAccelerator, Measures



const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const ALT_DATE_FORMAT_STR = "yyyymmdd"
const DATE_FORMAT_STR = "m/d/y"
const SP500_NAME = "GSPC"
const TBILL_NAME = "t-bill"
const TERM_NAME = "termStructure"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const OPTION_STD_NAME = "VIXoptionsStd"
const OPTION_NAME = "VIXoptions"

const PlotContainer = Union{Plot,Gadfly.Compose.Context}
const CSRegSpec = Tuple{Symbol, Symbol, Int}
const CPRegSpec = Tuple{Symbol, Int}

const MIN_VALS_FOR_EST = 10
const Ι = Complex(0.0,1.0)
const STANDARD_DAYS = [30, 60, 91, 122]


#pre-process the s&p500 data
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

#a simple routine for making an easy to read dataframe binary
function preProcessTBills()::Void
  dateFormat::String = DATE_FORMAT_STR
  tBillDF::DataFrame = readtable("$DATA_PATH\\$TBILL_NAME.csv")

  #fix the dates
  tBillDF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(tBillDF[:,:DATE])::DataVector{Date}
  sort!(tBillDF, cols = [:DATE])
  tBillDF[:,:TB3MS] .*= .01
  tBillDF[:,:rates]= tBillDF[:,:TB3MS] ./12.0 #(1.0+tBillDF[:,:TB3MS]).^(1.0/12.0) .- 1.0

  stream::IOStream = open("$DATA_PATH\\$TBILL_NAME.jls", "w")
  serialize(stream, tBillDF)
  close(stream)

  return nothing
end

#pre-process standard options data
function preProcessVIXOptionsStd(;verbose::Bool = false)::Void
  dateFormat::String = ALT_DATE_FORMAT_STR
  stdDF::DataFrame = readtable("$DATA_PATH\\$OPTION_STD_NAME.csv")

  stdDays::SVector{length(STANDARD_DAYS), Int} = SVector{length(STANDARD_DAYS), Int}(STANDARD_DAYS)

  #adjust date type
  stdDF[:, :date] = ((i::Int)->Date("$i", dateFormat)::Date).(stdDF[:,:date])::DataVector{Date}
  rename!(stdDF, :date, :DATE)


  stdDF = stdDF[completecases(stdDF[:,[:DATE, :days, :forward_price, :delta, :cp_flag]]),:]

  if verbose
    println("total record")
    showcols(stdDF)

    for d::Int ∈ unique(stdDF[:,:days])
      println("records for $d")
      showcols(stdDF[stdDF[:,:days].==d,:])
    end
  end

  maximumTau::Int = maximum(STANDARD_DAYS)
  stdDF = stdDF[stdDF[:,:days].≤maximumTau,:]

  #make this of type char
  delete!(stdDF, [:secid, :ticker, :index_flag])
  stdDF[:, :cp_flag] = ((s::String)->s[1]).(stdDF[:,:cp_flag])::DataVector{Char}
    sort!(stdDF, cols = [:DATE, :cp_flag, :days])

  #showcols(stdDF)

  stream::IOStream = open("$DATA_PATH\\$OPTION_STD_NAME.jls", "w")
  serialize(stream, stdDF)
  close(stream)

  return nothing

end

#pre-process standard options data
function preProcessVIXOptions()::Void
  #constants that will be used later
  dateFormat::String = ALT_DATE_FORMAT_STR
  optionDF::DataFrame = readtable("$DATA_PATH\\$OPTION_NAME.csv")


  #adjust date type
  optionDF[:, :date] = ((i::Int)->Date("$i", dateFormat)::Date).(optionDF[:,:date])::DataVector{Date}
  rename!(optionDF, :date, :DATE)

  optionDF[:, :exdate] = ((i::Int)->Date("$i", dateFormat)::Date).(optionDF[:,:exdate])::DataVector{Date}
  rename!(optionDF, :exdate, :exDATE)

  maximumTau::Int = maximum(STANDARD_DAYS)
  optionDF[:,:days] = (Dates.value).(optionDF[:, :exDATE])  .- (Dates.value).(optionDF[:, :DATE])
  optionDF = optionDF[(!isna).(optionDF[:,:days]), :]
  optionDF = optionDF[optionDF[:,:days] .≤ maximumTau, :]
  optionDF[:, :τ] = optionDF[:,:days] ./ 365.0

  #sort and re-scale
  optionDF[:, :strike_price] = optionDF[:, :strike_price] ./ 1000.0
  sort!(optionDF, cols = [:DATE, :exDATE, :cp_flag, :strike_price])

  optionDF[:, :price] = (optionDF[:, :best_offer] .+ optionDF[:, :best_bid]) ./ 2.0

  #filter out bad data
  optionDF[:, :cp_flag] = ((s::String)->s[1]).(optionDF[:,:cp_flag])::DataVector{Char}
  optionDF = optionDF[completecases(optionDF[:,[:price, :delta, :cp_flag]]),:]
  delete!(optionDF, [:best_offer, :best_bid, :secid, :ticker, :index_flag, :issuer, :exercise_style])



  #showcols(optionDF)

  stream::IOStream = open("$DATA_PATH\\$OPTION_NAME.jls", "w")
  serialize(stream, optionDF)
  close(stream)

  return nothing

end


#this gets a mapping of a given day to a standard day, and the weight for interpolation
function getInterpolationDayIndex()::Tuple{Vector{Int}, Vector{Float64}}
  standardDays::Vector{Int} = STANDARD_DAYS
  maxStandardDays = maximum(standardDays)

  #allocate the return containers, which include the mapping and the weights
  standardDayIndex::Vector{Int} = Vector{Int}(maximum(standardDays))
  weights::Vector{Float64} = Vector{Float64}(maxStandardDays)

  #iterate through all of the standardized days and map them to a standard day
  standardCtr::Int = 1
  for i::Int = 1:maxStandardDays
    if i>standardDays[standardCtr]
      standardCtr += 1
    end
    standardDayIndex[i] = standardCtr

    #calculate the weight
    weights[i] = standardCtr==1? 0.0 :
      (i-standardDays[standardCtr-1]) / (standardDays[standardCtr] - standardDays[standardCtr-1])

  end

  return (standardDayIndex, weights)
end

#merge the optionDF and stdDF
function mergeOptionDFs!(optionDF::DataFrame, stdDF::DataFrame)::DataFrame

  #calculate information about the standard time periods
  standardDays::Vector{Int} = STANDARD_DAYS
  minStandardDays = minimum(standardDays)
  maxStandardDays = maximum(standardDays)
  standardDayIndex::Vector{Int}, weights::Vector{Float64} = getInterpolationDayIndex();

  #index up the standard maturities (will help with merge)
  Dateτ::DataVector{Tuple{Date, Int, Char}} =
    ((d::Date, f::Int, c::Char)->(d,f,c)).(stdDF[:,:DATE],stdDF[:,:days], stdDF[:,:cp_flag])

  stdDict::Dict = Dict(Dateτ[i]=>i for i::Int ∈ 1:size(stdDF,1))
  FStd::DataVector{Float64} =  stdDF[:, :forward_price]

  #for each date and ex date pair
  N::Int =size(optionDF,1)

  optionDF[:,:F] = DataVector{Float64}(Vector{Float64}(N))
  optionDF[:,:F] .= NA
  optionDF[:,:invalid] = false


  #interpolate the forward price
  by(optionDF, :DATE) do dateDF::SubDataFrame
    by(dateDF, :exDATE) do exDateDF::SubDataFrame
      by(exDateDF, :cp_flag) do subDF::SubDataFrame

        period::Int = min(subDF[1,:days], maxStandardDays)
        current::Date = subDF[1,:DATE]
        cp::Char = subDF[1,:cp_flag]

        #now calculate the interpolated value
        if !haskey(stdDict, (current, standardDays[standardDayIndex[period]], cp))
          subDF[:,:invalid] = true
        elseif period ∈ standardDays #if no need for interpolation
          subDF[:,:F] = FStd[stdDict[(current, period, cp)]]
        elseif standardDayIndex[period] == 1
          subDF[:,:F] = FStd[stdDict[(current, standardDays[1], cp)]]
        else
          #take the weighted average of the nearest values
          subDF[:,:F] =
            FStd[stdDict[(current, standardDays[standardDayIndex[period]], cp)]] * weights[period] +
            FStd[stdDict[(current, standardDays[standardDayIndex[period]-1], cp)]] * (1.0-weights[period])
        end
      end
    end
  end

  optionDF = optionDF[(!isna).(optionDF[:,:F]),:]
  optionDF = optionDF[!optionDF[:,:invalid], :]
  delete!(optionDF, :invalid)
  N = size(optionDF,1)

  optionDF[:, :M] = optionDF[:, :F] ./ optionDF[:, :strike_price]

  #now get the high/low metric for the VIX

  medianST::Float64 = median(stdDF[stdDF[:,:days].==minStandardDays,:forward_price])
  medianLT::Float64 = median(stdDF[stdDF[:,:days].==maxStandardDays,:forward_price])
  medianTot::Float64 = median(stdDF[(stdDF[:,:days].==minStandardDays) .| (
    stdDF[:,:days].==maxStandardDays),:forward_price])


  optionDF[:,:LT_VIX] = Vector{Symbol}(N)
  optionDF[:,:LT_VIX] = NA

  optionDF[:,:ST_VIX] = Vector{Symbol}(N)
  optionDF[:,:ST_VIX] = NA

  optionDF[:,:LT_VIX_BC] = Vector{Symbol}(N)
  optionDF[:,:LT_VIX_BC] = NA

  optionDF[:,:ST_VIX_BC] = Vector{Symbol}(N)
  optionDF[:,:ST_VIX_BC] = NA

  optionDF[:,:VIX_BC] = Vector{Symbol}(N)
  optionDF[:,:VIX_BC] = NA

  by(optionDF, :DATE) do dateDF::SubDataFrame
    by(dateDF, :cp_flag) do subDF::SubDataFrame
      #appropriate standard indices

      #create a view for calculating if the VIX is high or low
      stdSubDF =
        view(stdDF, (stdDF[:,:DATE] .== subDF[1,:DATE]) .& (stdDF[:,:cp_flag] .== subDF[1,:cp_flag]))

      #compute if high or low
      maxDays::Int = maximum(stdSubDF[:,:days])
      minDays::Int = minimum(stdSubDF[:,:days])
      LTVIX::Float64 = (stdSubDF[stdSubDF[:,:days].==maxDays,:forward_price])[1]
      STVIX::Float64 = (stdSubDF[stdSubDF[:,:days].==minDays,:forward_price])[1]

      #update the dataframe
      subDF[:,:LT_VIX] = LTVIX>medianLT?(:H):(:L)
      subDF[:,:ST_VIX] = STVIX>medianST?(:H):(:L)

      subDF[:,:LT_VIX_BC] = LTVIX>medianTot?(:H):(:L)
      subDF[:,:ST_VIX_BC] = STVIX>medianTot?(:H):(:L)
      subDF[:,:VIX_BC] = STVIX>LTVIX*exp(-subDF[1,:ρ]*(maxDays-minDays)/365.0)?(:B):(:C)

    end
  end

  return optionDF
end

#this function merges in the risk-free rate
function mergeOptionDFρ!(optionDF::DataFrame, tBillDF::DataFrame)::DataFrame

  rates::Vector{Float64} = tBillDF[:,:TB3MS]
  rateDates::DataVector{Date} = tBillDF[:,:DATE]
  rateDict::Dict = Dict(rateDates[i]=>rates[i] for i::Int ∈ 1:length(rates))

  #allocate the interest rate column
  optionDF[:, :ρ] = DataVector{Float64}(Vector{Float64}(size(optionDF,1)))
  optionDF[:, :ρ] .= NA
  #for each date, look up the interest rate
  by(optionDF, :DATE) do subDF
    dt::Date = Dates.firstdayofmonth(subDF[1, :DATE])
    if haskey(rateDict, dt)
      subDF[:, :ρ] = rateDict[dt]
    end
  end

  #get rid of invalid interest rate data
  optionDF = optionDF[(!isna).(optionDF[:,:ρ]), :]
  optionDF[:,:ρ] = (log).(1.0 .+ optionDF[:,:ρ])


  return optionDF

end

#convert all puts to calls
function conformToCalls!(optionDF::DataFrame)::DataFrame
  putDF::SubDataFrame = view(optionDF, optionDF[:,:cp_flag].=='P')

  #use put call parity to convert all puts to calls
  putDF[:, :price] += ((exp).(-1.0 .* putDF[:,:τ] .* putDF[:,:ρ])) .*
    (putDF[:,:F] .- putDF[:,:strike_price])


  delete!(optionDF, :cp_flag)

  return optionDF
end

#black model with no dividends
function getImpliedBlackσ(price::Float64, F::Float64, τ::Float64, X::Float64,
  r::Float64, startσ::Float64 = 0.15)::Float64
  #Algorithms: :LD_LBFGS, :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQP, :LD_TNEWTON, :LD_VAR2/:LD_VAR1,
  #  :LD_MMA, :GN_ESCH, :GN_DIRECT, GN_DIRECT_L, :GD_STOGO
  # GN_ISRES, :LD_LBFGS, :LN_SBPLX,
  #GN_CRS2_LM
  #:GN_CRS2_LM***, :GN_ISRES*** (slow), :LD_MMA**, LN_SBPLX**, LD_TNEWTON** best:62,848

  B::Float64 = exp(-1.0*r*τ)

  function obj(σ::Real)::Real
    d1::Real = (log(F/X) + τ*(σ^2/2.0))/ (σ*τ^0.5)

    calcPrice::Real = (F*cdf(Normal(),d1) - X*cdf(Normal(),d1-σ*τ^0.5))*B
    return (price-calcPrice)*(price-calcPrice)
  end

  #println("test obj: $(obj(0.3))")

  #works well:LN_NELDERMEAD , LN_SBPLX, LD_SLSQP, LN_BOBYQA
  mod = Model(solver=NLoptSolver(algorithm=:LD_SLSQP, maxtime=1))
  JuMP.register(mod, :obj, 1, obj, autodiff=true)

  #constrain vol to 0 and 1000%
  @variable(mod, 3.0 >= σ >= 10.^-6.0, start = startσ)

  #set the objective function
  @NLobjective(mod, Min, obj(σ))

  σBest::Float64  = -99.0
  objBest::Float64 = -9999.0

  try
    solve(mod)
    σBest = getvalue(σ)
  catch
    σBest = -99.0
  end

  if σBest ≤   2.0*10.0 ^ -6.0
    σBest = -99.0
  end

  return σBest
end

#calcualte the IV for all options in the dataframe
function calculateIVs!(optionDF)::DataFrame
  N::Int = size(optionDF,1)

  Cs::Vector{Float64} = optionDF[:, :price]
  Fs::Vector{Float64} = optionDF[:, :F]
  startσs::Vector{Float64} = optionDF[:, :impl_volatility]
  τs::Vector{Float64} = optionDF[:, :τ]
  Xs::Vector{Float64} = optionDF[:, :strike_price]
  ρs::Vector{Float64} = optionDF[:, :ρ]

  #pre-allocate the answer space
  σs::Vector{Float64} = Vector{Float64}(N)
  #  σs .= NA

  @acc begin for i::Int  ∈  1:N
    σs[i] = getImpliedBlackσ(Cs[i], Fs[i], τs[i], Xs[i], ρs[i], startσs[i])
  end end

  optionDF[:,:σ] = DataVector{Float64}(Vector{Float64}(N))
  optionDF[:,:σ] .= ((f::Float64)->f<0.0?NA:f).(σs)
  println("Solution failures: $(length(optionDF[isna.(optionDF[:,:σ]),:σ]))
    out of $(length(optionDF[:,:σ]))")


  return optionDF
end



#makes the volatility smile graph
function makeSmileGraphs(optionDF::DataFrame; label::Symbol = :full,
    binScaleM::Float64 = 25.0, binScaleDays::Float64 = 0.1,
    width::Measures.Length = 4inch, height::Measures.Length= 4inch)::Void

  ####first smile plot
  plotDF::DataFrame = deepcopy(optionDF)
  delete!(plotDF, [:DATE, :exDATE, :strike_price, :F, :ρ, :impl_volatility, :τ, :price,
    :ST_VIX, :LT_VIX, :ST_VIX_BC, :LT_VIX_BC, :VIX_BC])
  completecases!(plotDF)

  plotDFNoBin::DataFrame = deepcopy(plotDF)

  println("$label Records: $(length(plotDF[:,:σ]))")


  plotDF[:,:M] = ((f::Float64)->round(f*binScaleM)/binScaleM).(plotDF[:,:M])
  plotDF[:,:days] = ((i::Int)->Int(round(i*binScaleDays)/binScaleDays)).(plotDF[:,:days])
  plotDF=aggregate(plotDF, [:days, :M], mean)

  #get the graphs
  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  push!(plotNames, "multiSmile")
  push!(plots, plot(plotDF, x=:M, y=:σ_mean, color=:days,
    Geom.line,
    Guide.ylabel("σ"), Guide.xlabel("F/X"),
    Coord.Cartesian(ymin=0.0,ymax=2.5, xmin=0.5, xmax=2.0),
    Guide.title("Volatility Smile  by Time and Moneyness $(label)")))

  push!(plotNames, "multiSmileDot")
  push!(plots, plot(plotDF, x=:M, y=:σ_mean, color=:days, shape=:days,
    Geom.point,
    Guide.ylabel("σ"), Guide.xlabel("F/X"),
    Coord.Cartesian(ymin=0.0,ymax=2.5, xmin=0.5, xmax=2.0),
    Guide.title("Volatility Smile  by Time and Moneyness $(label)")))


  for i ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i])_$label.svg", width, height), plots[i])
  end

  return nothing
end

#this encapsulates the processing of the option DF
function processOptionData!(optionDF::DataFrame, stdDF::DataFrame, tBillDF::DataFrame)
  optionDF = mergeOptionDFρ!(optionDF, tBillDF) #merge in interest rates
  optionDF = mergeOptionDFs!(optionDF, stdDF) #get the merged data
  optionDF = conformToCalls!(optionDF) #convert puts to calls

  println("beginning implied volatility calculation")
  optionDF =  calculateIVs!(optionDF) #calculate the implied volatilities

  #write this out to a file
  stream::IOStream = open("$DATA_PATH\\$OPTION_NAME-merged.jls", "w")
  serialize(stream, optionDF)
  close(stream)
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

function makeAllSmileGraphs(optionDF::DataFrame)
  println("\nTotal Records: $(length(optionDF[:,:ST_VIX]))")
  makeSmileGraphs(optionDF, label = :smooth, binScaleM = 10.0, binScaleDays=0.05, width=7inch, height=8.5inch)
  makeSmileGraphs(optionDF, label = :fine, binScaleM = 1000.0, binScaleDays=0.025, width=7inch, height=8.5inch)

  makeSmileGraphs(optionDF[(optionDF[:,:ST_VIX] .== :H) .& (optionDF[:,:LT_VIX] .== :H), :],
    label = :HH, binScaleM = 10.0, binScaleDays=0.05)
  makeSmileGraphs(optionDF[(optionDF[:,:ST_VIX] .== :L) .& (optionDF[:,:LT_VIX] .== :H), :],
    label = :LH, binScaleM = 10.0, binScaleDays=0.05)
  makeSmileGraphs(optionDF[(optionDF[:,:ST_VIX] .== :H) .& (optionDF[:,:LT_VIX] .== :L), :],
    label = :HL, binScaleM = 10.0, binScaleDays=0.05)
  makeSmileGraphs(optionDF[(optionDF[:,:ST_VIX] .== :L) .& (optionDF[:,:LT_VIX] .== :L), :],
    label = :LL, binScaleM = 10.0, binScaleDays=0.05)

    makeSmileGraphs(optionDF[(optionDF[:,:VIX_BC] .== :C), :],
      label = :C, binScaleM = 10.0, binScaleDays=0.05, width=3.7inch, height=8.5inch)
    makeSmileGraphs(optionDF[(optionDF[:,:VIX_BC] .== :B), :],
      label = :B, binScaleM = 10.0, binScaleDays=0.05, width=3.7inch, height=8.5inch)



end

function makeSP500VIXGraphs(stdDF::DataFrame, sp500DF::DataFrame; label::Symbol = :full,
      width::Measures.Length = 7inch, height::Measures.Length= 3.5inch)::Void

    standardDays::Vector{Int} = STANDARD_DAYS

    minStandardDays = minimum(standardDays)
    maxStandardDays = maximum(standardDays)

    #aggregate the forward price
    plotDF::DataFrame = deepcopy(stdDF)
    delete!(plotDF, [:impl_volatility, :delta, :cp_flag])
    completecases!(plotDF)
    plotDF=aggregate(plotDF, [:DATE, :days], mean)

    Gadfly.push_theme(:default)
    plots::Vector{PlotContainer} = Vector{PlotContainer}()
    plotNames::Vector{String} = Vector{String}()

    sp500DF = sp500DF[2:end, :]
    sp500DF = sp500DF[sp500DF[:,:DATE].≥minimum(plotDF[:,:DATE]), :]

    push!(plotNames, "VIX")
    push!(plots, plot(plotDF, x=:DATE, y=:forward_price_mean, color=:days,
      Geom.line,
      Guide.ylabel("Price"), Guide.xlabel("Date"), Theme(line_width = 0.25pt),
      Coord.Cartesian(ymin=0.0,ymax=75, xmin=Date(2007,1,1), xmax=Date(2017,1,1)),
      Guide.title("Time Series of VIX Forwards by Tenure")))

    push!(plotNames, "SP500")
    push!(plots, plot(sp500DF[2:end,:], x=:DATE, y=:returns,
      Geom.line,
      Guide.ylabel("Return"), Guide.xlabel("Date"),
      Coord.Cartesian(ymin=-0.15,ymax=0.15, xmin=Date(2007,1,1), xmax=Date(2018,1,1)),
      Guide.title("Time Series of SP500 Returns")))


    for i ∈ 1:length(plots)
      draw(SVG("$OUTPUT_PATH\\$(plotNames[i])_$label.svg", width, height), plots[i])
    end

    return nothing
end

#script for getting and processing the data
function A6Script(;refreshData::Bool=true, reprocessOptions=true)::Void

  #read the binary file
  tBillDF::DataFrame = getDF("$DATA_PATH\\$TBILL_NAME.jls", refreshData, preProcessTBills)
  sp500DF::DataFrame = getDF("$DATA_PATH\\$SP500_NAME.jls", refreshData, preProcessSP500)
  stdDF::DataFrame = getDF("$DATA_PATH\\$OPTION_STD_NAME.jls", refreshData, preProcessVIXOptionsStd)
  optionDF::DataFrame = getDF("$DATA_PATH\\$OPTION_NAME.jls", refreshData, preProcessVIXOptions)

  #we need to create an anonymous function to use the getDF interface
  #mergeOptionDFs!(deepcopy(optionDF),deepcopy(stdDF))
  optionDF = getDF("$DATA_PATH\\$OPTION_NAME-merged.jls", reprocessOptions,
    ()->processOptionData!(optionDF, stdDF,tBillDF))

  #makeAllSmileGraphs(optionDF)
  makeSP500VIXGraphs(stdDF,sp500DF)

  #print the summary stats
  preProcessVIXOptionsStd(verbose = true)

  #write out the processed option file for testing purposes
  writetable("$OUTPUT_PATH\\$OPTION_NAME.csv", optionDF)
  return nothing
end



@time begin
  #A6Script(refreshData=false, reprocessOptions=false)
end

end
