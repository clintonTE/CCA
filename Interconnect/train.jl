module train

#Use this to print
#=
using Weave
codeName = "train"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

#NOTE: To prepare the data, conform missing data to "". That means eliminating
#missing codes "N/A", "missing", and "--"

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, CSV, StatsBase, GZip, ForwardDiff,# Query,
  Gadfly, NLopt, StaticArrays, Measures, Formatting, Dates, Serialization


using CTCopy

#############parameters
const VOL_PERIODS = 1

#path constants
const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT = "m/d/yyyy"
const CHINA_DS_NAME = "ChinaDS1.1"
const CHINA_DS_NAME_L = "$(CHINA_DS_NAME)_L"
const CHINA_DS_NAME_W = "$(CHINA_DS_NAME)_W"
const R_TEST_PATH = "C:\\Users\\Clinton Tepper\\Dropbox\\Projects\\InterconnectRTester\\trainTester\\trainTester"
const ACCEPTABLE_RANGE = #use this to check if values are in acceptable ranges
    Dict(:value=>(-10., 10.), :marketCap=>(0.1*10. ^ 9.,Inf))
const DROP_OUTSIDE_RANGE = true
const MIN_SZSECOM_MARKETCAP = 6. * 10. ^ 9.
const KEEP_MISSING_MARKETCAP = false

const MIN_SAMPLES_IN_INDUSTRY = 5
#const MIN_SAMPLES_INTERACTED = 50 #obsolete
const MIN_COMPANY_QUARTERS = 1 #threshold for dropping companies
const SET_MISSING_IF_BELOW = true
const KEEP_LARGEST_WEIGHT_IF_DUP = true # keep the largest weight if multiple shareclasses

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
#const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs

const DATE_SYMBOLS =  Symbol("12/31/2010"), Symbol("3/31/2011"), Symbol("6/30/2011"),
  Symbol("9/30/2011"), Symbol("12/31/2011"), Symbol("3/31/2012"), Symbol("6/30/2012"),
  Symbol("9/30/2012"), Symbol("12/31/2012"), Symbol("3/31/2013"), Symbol("6/30/2013"),
  Symbol("9/30/2013"), Symbol("12/31/2013"), Symbol("3/31/2014"), Symbol("6/30/2014"),
  Symbol("9/30/2014"), Symbol("12/31/2014"), Symbol("3/30/2015"), Symbol("6/30/2015"),
  Symbol("9/30/2015"), Symbol("12/31/2015"), Symbol("3/30/2016"), Symbol("6/30/2016"),
  Symbol("9/30/2016"), Symbol("12/31/2016"), Symbol("3/30/2017"), Symbol("6/30/2017"),
  Symbol("9/30/2017"), Symbol("12/31/2017"), Symbol("3/31/2018")

const LAUNCH_DATE_SH = Date("11/17/2014",  DATE_FORMAT)
const LAUNCH_DATE_SZ = Date("12/5/2016",  DATE_FORMAT)
const FOCAL_DATES = (Date("6/30/2013",  DATE_FORMAT), Date("3/31/2018",  DATE_FORMAT))
const FOCAL_DATES_DD1 = (Date("6/30/2013",  DATE_FORMAT), Date("3/30/2016",  DATE_FORMAT))
const FOCAL_DATES_DD2 = (Date("6/30/2015",  DATE_FORMAT),Date("3/31/2018",  DATE_FORMAT))
const OUTCOMES = ["NORMALIZED_ACCRUALS_CF_METHOD"]

const MFloat64 = Union{Missing, Float64}
const MInt = Union{Missing, Int}
const MBool = Union{Missing, Bool}
const MSymbol = Union{Missing, Symbol}

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

#preprocess rate data
function preProcessChinaDSW()::Nothing
  wideDF::DataFrame = CSV.read("$DATA_PATH\\$CHINA_DS_NAME.csv", strings=:raw)

  stream::IOStream = open("$DATA_PATH\\$CHINA_DS_NAME_W.jls", "w")
  serialize(stream, wideDF)
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

function summaryVTable(wideDFV::DataFrame,
    decimals::Int = DECIMALS)::String

  numColumns = 8

  #setup the table headers
  columnNames::Vector{Vector{String}} = [["Shanghai Volume \$\\times 10^{-12}\$", "Shenzhen Volume \$\\times 10^{-12}\$"], ["Shr Vol", "N.B. Buys",
      "N.B. Sells", "Total N.B.", "Shr Vol", "N.B. Buys", "N.B. Sells", "Total N.B."]]
  columnSizes::Vector{Vector{Int}} = [ones(Int, numColumns ÷ 4) .* 4, ones(Int, numColumns)]

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Mean", "Median", "Std Deviation",
    "10\\% Quantile", "90\\% Quantile", "N>0 (Quarters)"]


  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(undef,numColumns)).(1:length(descRowNames))

  #get and scale the summary variables
  series::Vector{Vector{MFloat64}} = [
      wideDFV[:SHCOMP_Index],
      wideDFV[:C1DBTO_Index],
      wideDFV[:C1DSTO_Index],
      wideDFV[:C1DBTO_Index] .+ wideDFV[:C1DSTO_Index],
      wideDFV[:SZCOMP_Index],
      wideDFV[:C2DBTO_Index],
      wideDFV[:C2DSTO_Index],
      wideDFV[:C2DBTO_Index] .+ wideDFV[:C2DSTO_Index]]

  #remove missing
  for i::Int ∈ 1:length(series)
    series[i] = series[i][((x::MFloat64)->!ismissing(x)).(series[i])]
    series[i] .= series[i] .* 10^-9.
  end



  for c::Int ∈ 1:numColumns
    descContent[1][c] = num2Str(mean(series[c]), decimals)
    descContent[2][c] = num2Str(median(series[c]), decimals)
    descContent[3][c] = num2Str(std(series[c]), decimals)
    descContent[4][c] = num2Str(quantile(Vector{Float64}(series[c]),0.1), decimals)
    descContent[5][c] = num2Str(quantile(Vector{Float64}(series[c]),0.9), decimals)
    descContent[6][c] = num2Str(length(series[c]), Ints=true)
  end


  summaryVTable::String = texTable( "Summary of Volume Data",
      """See Tex File""", #caption
      columnNames, #colNames
      Vector{String}(),#contentRowNames
      Vector{Matrix{String}}(), #content
      descRowNames, #descRowNames
      descContent, #descContent
      Vector{String}(),
      widthColNames = columnSizes,
      alignmentColNames = [["c", "c"], ["r", "r", "r", "r|", "r" , "r", "r", "r"]],
      lineSpacer = "\\\\"
    )

    writeTables2File([summaryVTable],
        HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
        outName = "trainVTables.tex")

  return summaryVTable
end

#computes policy outcomes for the different exchanges
function getExchangeVolumeData(longDFV::DataFrame)::DataFrame

  wideDFV::DataFrame = unstack(longDFV, :date, :security, :value)
  N::Int = size(wideDFV,1)

  #sum buy and sell orders

  wideDFV[:C1Vol] = Vector{MFloat64}(wideDFV[:C1DBTO_Index] .+ wideDFV[:C1DSTO_Index])
  wideDFV[:C2Vol] = Vector{MFloat64}(wideDFV[:C2DBTO_Index] .+ wideDFV[:C2DSTO_Index])

  wideDFV[:lC1Vol] = Vector{MFloat64}(wideDFV[:C1DBTO_Index] .- wideDFV[:C1DSTO_Index])
  wideDFV[:lC2Vol] = Vector{MFloat64}(wideDFV[:C2DBTO_Index] .- wideDFV[:C2DSTO_Index])

  #container for the policy variable
  wideDFV[:C1Effect] = Vector{MFloat64}(missings(N))
  wideDFV[:C2Effect] = Vector{MFloat64}(missings(N))

  wideDFV[:lC1Effect] = Vector{MFloat64}(missings(N))
  wideDFV[:lC2Effect] = Vector{MFloat64}(missings(N))

  sort!(wideDFV, :date)

  #this loop computes running totals of volue and sums them
  #the formula is connect period buys + sells / total period volume
  for i::Int ∈ VOL_PERIODS:N
    wideDFV[i, :C1Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:C1Vol])) /
      sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:SHCOMP_Index]))

    wideDFV[i, :C2Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:C2Vol])) /
      sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:SZCOMP_Index]))

    wideDFV[i, :lC1Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:lC1Vol])) /
      sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:SHCOMP_Index]))

    wideDFV[i, :lC2Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:lC2Vol])) /
      sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:SZCOMP_Index]))
  end


  return wideDFV
end

#runs any checks they are required on the wideDF post data-quality
function wideDiagnostics(wideDF::DataFrame)::Nothing

  #checks to see if we have cross-listings or duplicates
  countDF = by(wideDF, :security, nrow)
  NSecurities::Int = size(countDF,1)
  securities = countDF[:security]
  countIndex::Dict{Symbol, Int} = Dict(countDF[i,:security]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))

  #check for overloaded securities
  for i::Int ∈ 1:NSecurities
    if countIndex[securities[i]] > 1
      warn("Overloaded security: $(securities[i]) ($(countIndex[securities[i]])")
    end
  end

  return nothing
end

function graphTrends(longDF::DataFrame)::Nothing
  dateFormat::String = "yyyy-mm-dd"

  #get the names of all columns to be NOT unstacked
  #nonDateCols::Vector{Symbol} = setdiff(names(longDF), [variable, value])
  #wideDF::DataFrame = unstack(longDF, nonDateCols, variable, value)

  #set up the plotting environment
  Gadfly.push_theme(:default)
  plotNames::Vector{String} = Vector{String}()
  plots::Vector{PlotContainer} = Vector{PlotContainer}()


  #copy this to transform the data
  outDF::DataFrame = deepcopy(longDF[[:oldDate, :exchange, :value]])
  dropmissing!(outDF)
  aggDF::DataFrame = aggregate(outDF, [:exchange, :oldDate], mean)


  #showcols(aggDF)

  #D10::DateTime = Date(2010,12,31)
  #D11::DateTime = Date(2011,12,31)
  #D12::DateTime = Date(2012,12,31)
  D13::DateTime = Date(2013,12,31)
  D14::DateTime = Date(2014,12,31)
  D15::DateTime = Date(2015,12,31)
  D16::DateTime = Date(2016,12,31)
  D17::DateTime = Date(2017,12,31)
  majorDates::Vector{DateTime}= [#=D10, D11, D12, =#D13, D14, D15, D16, D17]
  #majorDates:::Vector{Int} = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017]

  push!(plotNames, "Parallel Trends")
  push!(plots, plot(aggDF, x=:oldDate, y=:value_mean, color=:exchange,
      Stat.xticks(ticks=majorDates), Geom.line, Guide.xlabel("Year"),
      Guide.ylabel("Average Firm Cashflow Accrual"),
      Scale.x_continuous(labels=x::DateTime->"$(Dates.year(x))"),
      Coord.Cartesian(xmin=Date(2013,6,30), xmax=D17), Guide.colorkey(title="Exchange:"),
      style(key_position=:bottom, guide_title_position=:left),
      xintercept=[LAUNCH_DATE_SH, LAUNCH_DATE_SZ], Geom.vline(style=[:solid, [1mm,1mm]])))

  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing


end

#conforms the industries
function conformIndustries!(longDF::DataFrame;
    variable::Symbol=:date, value::Symbol=:value,
    threshold::Int = MIN_SAMPLES_IN_INDUSTRY, newCol::Symbol = :gics,
    setMissing::Bool = SET_MISSING_IF_BELOW,
    runWideDiagnostics=true)::DataFrame

  dateFormat::String = "yyyy-mm-dd"

  #get the names of all columns to be NOT unstacked
  nonDateCols::Vector{Symbol} = setdiff(names(longDF), [variable, value])
  wideDF::DataFrame = unstack(longDF, nonDateCols, variable, value)

  #convert to a symbol
  wideDF[newCol] = Vector{MSymbol}((Symbol).(wideDF[:subIndustry]))



  #start at the lowest level of industry classification
  countDF::DataFrame = by(wideDF, newCol, nrow)

  countIndex::Dict{Symbol, Int} = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))

  #For each industry, examine the minimum count for subindustries. If its less than desirable
  #then conform the subindustries to industries
  by(wideDF, :industry) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[newCol] = Symbol(df[1,:industry])
    end
  end

  #repeat for industry groups and industries
  countDF = by(wideDF, newCol, nrow)
  countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
  by(wideDF, :industryGroup) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[newCol] = Symbol(df[1,:industryGroup])
    end
  end

  #and again one more time to conform some industry groups to sectors
  countDF = by(wideDF, newCol, nrow)
  countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
  by(wideDF, :sector) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[newCol] = Symbol(df[1,:sector])
    end
  end

  #finally, if the threshold is sitll not met, set to missing (if the option is selected)
  if setMissing
    countDF = by(wideDF, newCol, nrow)
    countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
    by(wideDF, :sector) do df::SubDataFrame

      subGICS::Vector{Symbol} = unique(df[newCol])
      minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

      #conform iff the min count is below the desired threshold
      if minCount < threshold
        df[newCol] = missing
      end
    end
  end

  #reform the long dataframe
  if runWideDiagnostics
    wideDiagnostics(wideDF)
  end

  longDF = melt(wideDF, [nonDateCols; newCol])
  longDF[:date] = (s::Symbol->Date(string(s),dateFormat)).(longDF[:variable])
  delete!(longDF, :variable)
  longDF = longDF[false .== (ismissing).(longDF[:value]),:] # need to refilter missing variables
  return longDF
end


#conveience function to drop bad rows
function pruneLong!(longDF::DataFrame)::DataFrame

  #make a tracking column and store some convient quantities
  dateFormat::String = DATE_FORMAT
  N::Int = size(longDF,1)
  longDF[:keep] = Vector{Bool}(trues(N))
  numQuarters::Int =
      sum((s::Symbol->Date(string(s),dateFormat) ∈ minimum(FOCAL_DATES):Day(1):maximum(
      FOCAL_DATES)).(DATE_SYMBOLS))

  println(unique(longDF[:date]))
  #now process the firm-level data


  by(longDF, :security) do df::SubDataFrame
    nSub::Int = size(df,1)

    #println("$nSub")
    #println("nSub: $nSub numQuarters:$numQuarters")
    #check for missing values, inconsistent lengths, and stocks switching exchanges
    df[(ismissing).(df[:value]) .|
      (ismissing).(df[:exchange]) .|
      (ismissing).(df[:marketCap]),:keep] = false


    #drop a stock if it has less than two years of data, if its too small,
    # or its listed on multiple exchanges
    if ((sum(df[:keep] .== true) < MIN_COMPANY_QUARTERS) ||
        (length(unique(df[:exchange])) ≠ 1) ||
        sum((!ismissing).(df[:marketCap])) == 0 || #require at least one marketcap
        (minimum(skipmissing(df[:marketCap])) < 5.) )

      df[:keep] = false

    end
  end



  longDF = longDF[longDF[:keep],:]
  ##finally need to check for  stocks listed multiple times

  by(longDF, :description) do dfName::SubDataFrame #for each company name

    by(dfName, :date) do dfUniqueDt #for each date
      if size(dfUniqueDt,1) > 1 #if there is a dup
        if KEEP_LARGEST_WEIGHT_IF_DUP #if true, we keep the largest weighted name
          maxWeight::Float64 = maximum(dfUniqueDt[:weight])
          dfUniqueDt[dfUniqueDt[:weight] .≠ maxWeight,:keep] = false
        else
          dfUniqueDt[:keep] = false
        end
      end
    end
  end

  longDF = longDF[longDF[:keep] .== true,:]


  #finally, do the winsorization
  keepSeries::Vector{Bool} = longDF[:keep]
  for s::Symbol ∈ [:value, :marketCap]
    series::Vector{MFloat64} = longDF[s]

    #μ::Float64 = mean(skipmissing(series))
    #σ::Float64 = std(skipmissing(series))
    lower::Float64 = ACCEPTABLE_RANGE[s][1]
    upper::Float64 = ACCEPTABLE_RANGE[s][2]
    for i::Int ∈ 1:length(series)
      if !ismissing(series[i])
        if series[i] < lower
          if DROP_OUTSIDE_RANGE
            keepSeries[i] = false
          else
            series[i] = lower
          end
        elseif series[i] > upper
          if DROP_OUTSIDE_RANGE
            keepSeries[i] = false
          else
            series[i] = upper
          end
        end
      end
    end

    longDF[s] = series
  end

  #use an anonymous function to execute the SZCECOM cutoff
  longDF[(longDF[:index] .== :SZSECOM) .& ((f::MFloat64->
    ismissing(f) ?  false : (f < MIN_SZSECOM_MARKETCAP)
    ).(longDF[:marketCap])), :keep] = false

  longDF[:keep] = keepSeries
  longDF = longDF[longDF[:keep] .== true,:]
  delete!(longDF, :keep)


  return longDF
end

#processing activities for the long dataframe
function processLong!(wideDF::DataFrame)::DataFrame
  dateFormat::String = DATE_FORMAT

  #use all columns except the dates as ID columns
  #showcols(wideDF)
  nonDateCols::Vector{Symbol} = setdiff(names(wideDF), DATE_SYMBOLS)
  longDF::DataFrame = melt(wideDF, nonDateCols)

  #convert the newly melted columns to dates
  longDF[:date] = (s::Symbol->Date(string(s),dateFormat)).(longDF[:variable])
  delete!(longDF, :variable)
  longDF = longDF[(d::Date->d ∈ minimum(FOCAL_DATES):Day(1):maximum(FOCAL_DATES)).(longDF[:date]),:] #filter to the correct dates
  #println("names: $(names(longDF))")
  println(unique(longDF[:date]))
  #clean up names
  longDF[:security] =
      (s::String->Symbol(replace(strip(s), " "=>"_"))).(longDF[:security])

  #drop and collect non-outcome rows. Start with a flag for outcome variables
  println("typeof: $(typeof(longDF[:command]))")
  longDF[:Outcome] =(
      (s::String->s ∈ OUTCOMES ? true : false).(longDF[:command]))
  longDFV::DataFrame =
      longDF[.!(longDF[:Outcome]),[:date, :security, :value]]
  longDF = longDF[longDF[:Outcome],:]
  delete!(longDF, :Outcome) #no longer need this


  #flag for deleting rows with bad data
  longDF = longDF[(!).((ismissing).(longDF[:exchange])), :]
  longDF[:exchange] = Vector{MSymbol}((Symbol).(longDF[:exchange]))

  #delete marked rows
  println("Pruning Step 1. Starting rows: $(size(longDF,1))")
  longDF = pruneLong!(longDF)
  println("Pruning complete. Rows: $(size(longDF,1))")
  N::Int = size(longDF,1) #udpate size

  #conform industry codes so each industry gets some variation
  #a broader classificaiton is used for the interacted version
  longDF = conformIndustries!(longDF, runWideDiagnostics=false)
  #longDF = conformIndustries!(longDF, threshold=MIN_SAMPLES_INTERACTED, newCol = :gicsInteracted)


  N = size(longDF,1) #udpate size
  #get wide data for the input, and join to the long DF
  wideDFV::DataFrame = getExchangeVolumeData(longDFV)
  serialize(open("$DATA_PATH\\$(CHINA_DS_NAME)_VW.jls", "w"), wideDFV) #save a copy for later
  longDF = join(longDF,
      wideDFV[[:date, :C1Effect, :C2Effect, :lC1Effect, :lC2Effect]],
      on = :date, kind=:inner)

  summaryVTable(wideDFV)
  #now assign the treatment variables
  #begin by handling missing values and allocating space
  exchanges::Vector{MSymbol} = longDF[:exchange]

  #work with vectors for performance reasons
  C1Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x) ? 0.0 : Float64(x)).(longDF[:C1Effect])
  C2Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x) ? 0.0 : Float64(x)).(longDF[:C2Effect])
  CEffect::Vector{MFloat64} = Vector{MFloat64}(missings(N))
  lC1Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x) ? 0.0 : Float64(x)).(longDF[:lC1Effect])
  lC2Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x) ? 0.0 : Float64(x)).(longDF[:lC2Effect])
  lCEffect::Vector{MFloat64} = Vector{MFloat64}(missings(N))
  postTreatment = Vector{MInt}(missings(N))

  for i::Int ∈ 1:N
    if exchanges[i]==:Shanghai
      CEffect[i] = C1Effect[i]
      lCEffect[i] = lC1Effect[i]
    elseif exchanges[i]==:Shenzhen
      CEffect[i] = C2Effect[i]
      lCEffect[i] = lC2Effect[i]
    end

    postTreatment[i] = Int(abs(CEffect[i]) > 0.0)
  end

  longDF[:CEffect] = CEffect #focal policy variable
  longDF[:lCEffect] = lCEffect #focal policy variable
  longDF[:postTreatment] = postTreatment
  longDF[:gicsXDate] =
    Vector{MSymbol}(((s1::MSymbol,s2::Date)->ismissing(s1) ? missing : Symbol(
    s1,s2)).(longDF[:gics],longDF[:date]))

  longDF[:lMarketCap] = Vector{MFloat64}((log).(longDF[:marketCap]))
  longDF[:oldDate] = Vector{Date}(deepcopy(longDF[:date]))

  longDF[:DD1] = (d::Date->Int(
      d ∈ minimum(FOCAL_DATES_DD1):Day(1):maximum(FOCAL_DATES_DD1))).(longDF[:date])
  longDF[:DD2] = (d::Date->Int(
      d ∈ minimum(FOCAL_DATES_DD2):Day(1):maximum(FOCAL_DATES_DD2))).(longDF[:date])

  categorical!(longDF, :security)
  categorical!(longDF, :exchange)
  categorical!(longDF, :gics)
  #  categorical!(longDF, :gicsInteracted)
  categorical!(longDF, :date)
  categorical!(longDF, :gicsXDate)

  #println("unique date= $(length(unique(longDF[:date])))")
  #println("unique industry= $(length(unique(longDF[:gics])))")
  #println("unique interacted = $(length(unique(longDF[:gicsXDate])))")

  println("\nLength: $(size(longDF,1)), Length Unique: $(size(unique(longDF[:security]),1))")
  return longDF
end

function processDF(refreshData::Bool = true)
  wideDF::DataFrame = getDF("$DATA_PATH\\$CHINA_DS_NAME_W.jls", refreshData, preProcessChinaDSW)
  longDF::DataFrame = processLong!(wideDF)

  stream::IOStream = open("$DATA_PATH\\$CHINA_DS_NAME_L.jls", "w")
  serialize(stream, longDF)
  close(stream)
end

function discrete(longDF::DataFrame ;
    YSym::Symbol = :value, focalSpec::String = "postTreatment",
    title::String = "Discrete Specifications",
    focalName=:integrated)::String

    decimals::Int = DECIMALS

  #showcols(longDF[[:lMarketCap, :value, :gics, :ey, :exchange, :postTreatment, :date]])


  #println("got here")
  #setup the containers
  YSpecs::Vector{Symbol} = Vector{Symbol}() # holds the LHS
  XSpecs::Vector{CTExpr} = Vector{CTExpr}() #holds the RHS
  XNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() # holds the names of the x variables
  withinSpec::Vector{CTExpr} = Vector{CTExpr}() #holds the within transformation variable
  clusteredSpec::Vector{CTExpr} = Vector{CTExpr}() #holds the clustering variable
  models::Vector{CTLM} = Vector{CTLM}() #holds the regression models

  #showall(describe(longDF))
  #println(longDF[1:5, [:security, :postTreatment, :exchange, :date]])

  #############specifications here

  push!(XSpecs, Meta.parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :exchange)
  push!(clusteredSpec, :gics) #tag1

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap])
  push!(YSpecs, YSym)
  push!(withinSpec, :exchange)
  push!(clusteredSpec, :gics) #tag1

  #=push!(XSpecs, Meta.parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, Meta.parse("$focalSpec + exchange + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + exchange + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap ])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + exchange"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + exchange"))
  push!(XNames, [:intercept; focalName; :lMarketCap])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gics)=#


  ###################end new regression specs

  NSpecs::Int = length(XSpecs)
  columnNames::Vector{Vector{String}} = [((i::Int)->"($i)").(1:NSpecs)]

  println("Running $title with focal var $focalSpec")
  for i::Int ∈ 1:NSpecs
      println("spec: $(XSpecs[i])")
      push!(models, CTLM(longDF, XSpecs[i],  YSpecs[i], withinSym = withinSpec[i],
          clusterSym = clusteredSpec[i], XNames=XNames[i], YName = YSpecs[i]))
  end

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Date F.E.",
      "Firm F.E.", "Exchange F.E.", "Industry F.E.", "Industry-Date F.E",
      "\$R^{2}\$ (within)", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} =
      ((i::Int)->Vector{String}(undef, length(models))).(1:length(descRowNames))

  #this builds the descriptive rows. There is Probably a better way to do this,
  #but its fairly specific to the project.
  for i ∈ 1:length(XSpecs)
    descContent[1][i] = "$(withinSpec[i]==:date ||
        occursin("date", string(XSpecs[i])) ? "X" : "")"
    descContent[2][i] = "$(withinSpec[i]==:security ||
        occursin("security", string(XSpecs[i])) ? "X" : "")"
    descContent[3][i] = "$(withinSpec[i]==:exchange ||
        occursin("exchange", string(XSpecs[i])) ? "X" : "")"
    descContent[4][i] = "$(withinSpec[i]==:gics ? "X" : "")"
    descContent[5][i] = "$(withinSpec[i]==:gicsXDate ? "X" : "")"
    descContent[6][i] =
    "$(num2Str(getR²(models[i], adjusted=false),  decimals))"
    descContent[7][i] = "$(models[i].N)"
  end

  tableText::String = texTable(models,
      getClustered!#=getNeweyWestFunc(5) (lm::CTLM)->getNeweyWestSlow(lm, 5)=#,
      [focalName, :lMarketCap],
      titleCaption = title,
      colNames = columnNames,
      contentRowNames = ["Integrated", "lMarketCap"],
      descRowNames = descRowNames,
      descContent = descContent,
      decimalDigits = 4,
      stars=true,
      starStrings = OVERRIDE_STAR_STRINGS,
      #clearMem = USE_AGGRESSIVE_GC,
      caption = "to be written")

  return tableText
end

function policy(longDF::DataFrame ;
    YSym::Symbol = :value, focalSpec::String = "CEffect",
    title::String = "Policy Specifications",
    focalName=:IMeasure)::String

    decimals::Int = DECIMALS

  #showcols(longDF[[:lMarketCap, :value, :gics, :ey, :exchange, :postTreatment, :date]])


  #println("got here")
  #setup the containers
  YSpecs::Vector{Symbol} = Vector{Symbol}() # holds the LHS
  XSpecs::Vector{CTExpr} = Vector{CTExpr}() #holds the RHS
  XNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() # holds the names of the x variables
  withinSpec::Vector{CTExpr} = Vector{CTExpr}() #holds the within transformation variable
  clusteredSpec::Vector{CTExpr} = Vector{CTExpr}() #holds the clustering variable
  models::Vector{CTLM} = Vector{CTLM}() #holds the regression models

  #showall(describe(longDF))
  #println(longDF[1:5, [:security, :postTreatment, :exchange, :date]])

  #############specifications here

  push!(XSpecs, Meta.parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :exchange)
  push!(clusteredSpec, :gics) #tag1

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap])
  push!(YSpecs, YSym)
  push!(withinSpec, :exchange)
  push!(clusteredSpec, :gics) #tag1

  #=push!(XSpecs, Meta.parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, Meta.parse("$focalSpec + exchange + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + exchange + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + exchange"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gics)

  push!(XSpecs, Meta.parse("$focalSpec + lMarketCap + exchange"))
  push!(XNames, [:intercept; focalName; :lMarketCap])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gics)=#


  ###################end new regression specs

  NSpecs::Int = length(XSpecs)
  columnNames::Vector{Vector{String}} = [((i::Int)->"($i)").(1:NSpecs)]

  println("Running $title with focal var $focalSpec")
  for i::Int ∈ 1:NSpecs
      println("spec: $(XSpecs[i])")
      push!(models, CTLM(longDF, XSpecs[i],  YSpecs[i], withinSym = withinSpec[i],
          clusterSym = clusteredSpec[i], XNames=XNames[i], YName = YSpecs[i]))
  end

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Date F.E.",
      "Firm F.E.", "Exchange F.E.", "Industry F.E.", "Industry-Date F.E",
      "\$R^{2}\$ (within)", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} =
      ((i::Int)->Vector{String}(undef, length(models))).(1:length(descRowNames))

  #this builds the descriptive rows. There is Probably a better way to do this,
  #but its fairly specific to the project.
  for i ∈ 1:length(XSpecs)
    descContent[1][i] = "$(withinSpec[i]==:date ||
        occursin("date", string(XSpecs[i])) ? "X" : "")"
    descContent[2][i] = "$(withinSpec[i]==:security ||
        occursin("security", string(XSpecs[i])) ? "X" : "")"
    descContent[3][i] = "$(withinSpec[i]==:exchange ||
        occursin("exchange", string(XSpecs[i])) ? "X" : "")"
    descContent[4][i] = "$(withinSpec[i]==:gics ? "X" : "")"
    descContent[5][i] = "$(withinSpec[i]==:gicsXDate ? "X" : "")"
    descContent[6][i] =
    "$(num2Str(getR²(models[i], adjusted=false),  decimals))"
    descContent[7][i] = "$(models[i].N)"
  end

  tableText::String = texTable(models,
      getClustered!#=getNeweyWestFunc(5) (lm::CTLM)->getNeweyWestSlow(lm, 5)=#,
      [focalName, :lMarketCap],
      titleCaption = title,
      colNames = columnNames,
      contentRowNames = ["$focalName", "lMarketCap"],
      descRowNames = descRowNames,
      descContent = descContent,
      decimalDigits = 4,
      stars=true,
      starStrings = OVERRIDE_STAR_STRINGS,
      #clearMem = USE_AGGRESSIVE_GC,
      caption = "to be written")

  return tableText
end

function summary(longDF::DataFrame, decimals::Int = DECIMALS)::String
  numColumns = 6

  #setup the table headers
  columnNames::Vector{Vector{String}} = [["Accrual Ratio", "Market Cap \$\\times 10^{-9}\$"],
      ["Shanghai", "Shenzhen", "Total", "Shanghai", "Shenzhen", "Total"]]
  columnSizes::Vector{Vector{Int}} = [ones(Int, numColumns ÷ 3) .* 3, ones(Int, numColumns)]

  #make views for conveneince
  longDFSH::SubDataFrame = view(longDF, longDF[:exchange] .== :Shanghai)
  longDFSZ::SubDataFrame = view(longDF, longDF[:exchange] .== :Shenzhen)

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Mean", "Median", "Std Deviation",
    "10\\% Quantile", "90\\% Quantile", "Unique Firms", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(undef, numColumns)).(1:length(descRowNames))

  #get and scale the summary variables
  series::Vector{Vector{MFloat64}} = [
      longDFSH[:value], longDFSZ[:value], longDF[:value],
        longDFSH[:marketCap].*10. ^ -9.,
        longDFSZ[:marketCap].*10. ^ -9.,
        longDF[:marketCap].*10. ^ -9.]

  #remove missing
  for i::Int ∈ 1:length(series)
    series[i] = series[i][((x::MFloat64)->!ismissing(x)).(series[i])]
  end


  #get the unique number of firms
  descContent[6] = [num2Str(length(unique(longDFSH[:security])), Ints=true),
      num2Str(length(unique(longDFSZ[:security])), Ints=true),
      num2Str(length(unique(longDF[:security])), Ints=true), "", "", ""]

  for c::Int ∈ 1:numColumns
    descContent[1][c] = num2Str(mean(series[c]), decimals)
    descContent[2][c] = num2Str(median(series[c]), decimals)
    descContent[3][c] = num2Str(std(series[c]), decimals)
    descContent[4][c] = num2Str(quantile(Vector{Float64}(series[c]),0.1), decimals)
    descContent[5][c] = num2Str(quantile(Vector{Float64}(series[c]),0.9), decimals)
    descContent[7][c] = num2Str(length(series[c]), Ints=true)
  end


  summaryTable::String = texTable( "Summary of Data",
      """See Tex File""", #caption
      columnNames, #colNames
      Vector{String}(),#contentRowNames
      Vector{Matrix{String}}(), #content
      descRowNames, #descRowNames
      descContent, #descContent
      Vector{String}(),
      widthColNames = columnSizes,
      alignmentColNames = [["c", "c"], ["r" for i::Int ∈ 1:numColumns]],
      lineSpacer = "\\\\"
    )



  return summaryTable
end



function runRegressions(longDF::DataFrame; runDiscrete::Bool = true, runPolicy::Bool = true,
    runSummary::Bool = true)

  longDF1::DataFrame = longDF[longDF[:DD1] .== 1, :]
  longDF2::DataFrame = longDF[longDF[:DD2] .== 1, :]
  #run the regression with the discrete specification
  discreteTable1::String = runDiscrete ? discrete(longDF1, title="Shanghai Interconnect Discrete") : ""
  discreteTable2::String = runDiscrete ? discrete(longDF2, title="Shenzhen Interconnect Discrete") : ""
  discreteTable::String = runDiscrete ? discrete(longDF, title="Combined Interconnect Discrete") : ""

  policyTableVol1::String = runPolicy ? policy(longDF1, focalSpec="CEffect",
    title="Shanghai Interconnect Volume Policy", focalName=:IMeasureVol) : ""
  policyTableVol2::String = runPolicy ? policy(longDF2, focalSpec="CEffect",
    title="Shenzhen Interconnect Volume Policy", focalName=:IMeasureVol) : ""
  policyTableVol::String = runPolicy ? policy(longDF, focalSpec="CEffect",
    title="Combined Interconnect Volume Policy", focalName=:IMeasureVol) : ""

  policyTableNet1::String = runPolicy ? policy(longDF1, focalSpec="lCEffect",
    title="Shanghai Interconnect Net Flow Policy", focalName=:IMeasureNet) : ""
  policyTableNet2::String = runPolicy ? policy(longDF2, focalSpec="lCEffect",
    title="Shenzhen Interconnect Net Flow Policy", focalName=:IMeasureNet) : ""
  policyTableNet::String = runPolicy ? policy(longDF, focalSpec="lCEffect",
    title="Combined Interconnect Net Flow Policy", focalName=:IMeasureNet) : ""

  #calculate the summary statistics
  summaryTable::String = runSummary ? summary(longDF) : ""
  graphTrends(longDF)

  #write the processed data to a CSV for any tests in R
  CSV.write("$R_TEST_PATH\\longDF.csv", longDF)

  #print tables
  writeTables2File([summaryTable,
      discreteTable1, discreteTable2, discreteTable,
      policyTableVol1, policyTableVol2, policyTableVol,
      policyTableNet1, policyTableNet2, policyTableNet],
      HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
      outName = "trainTables.tex")
end
#main script
function mainScript(;refreshData::Bool = true, processData::Bool =true,
    runDiscrete::Bool = true, runPolicy::Bool = true, runSummary::Bool = false)::Nothing

  longDF::DataFrame = getDF(
      "$DATA_PATH\\$CHINA_DS_NAME_L.jls", processData, ()->processDF(refreshData))

  # run the regressions
  runRegressions(longDF, runDiscrete = runDiscrete, runPolicy = runPolicy, runSummary = runSummary)

  return nothing
end


# Uncomment to run, commented to print the code appendix

@time begin #run main
  mainScript(refreshData = true, processData=true, runDiscrete=true,
      runPolicy=true, runSummary=true)
end

end #end module
