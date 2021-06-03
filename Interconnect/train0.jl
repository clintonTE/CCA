module train

#Use this to print
#=
using Weave
codeName = "TimeSeriesA6"
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

using DataFrames, Distributions, CSV, StatsBase, GZip, ForwardDiff, Query,
  Gadfly, NLopt, HCubature, StaticArrays, Measures, Formatting


importall CTIO, CTReg, CTStat, CTVAR, CTNumerical

#############parameters
const VOL_PERIODS = 1

#path constants
const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT = "m/d/yyyy"
const CHINA_DS_NAME = "ChinaDS9"
const CHINA_DS_NAME_L = "$(CHINA_DS_NAME)_L"
const CHINA_DS_NAME_W = "$(CHINA_DS_NAME)_W"
const R_TEST_PATH = "C:\\Users\\Clinton Tepper\\Dropbox\\Projects\\InterconnectRTester\\trainTester\\trainTester"
const ACCEPTABLE_RANGE =
    Dict(:value=>(-0.5, 0.5), :marketCap=>(100*10.^6.,Inf), :ey =>(-200., 200.))


const MIN_SAMPLES_IN_INDUSTRY = 5
const MIN_SAMPLES_INTERACTED = 50 #obsolete
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

const DATE_SYMBOLS =  Symbol("12/31/2010"), Symbol("3/31/2011"), Symbol("6/30/2011"),
  Symbol("9/30/2011"), Symbol("12/31/2011"), Symbol("3/31/2012"), Symbol("6/30/2012"),
  Symbol("9/30/2012"), Symbol("12/31/2012"), Symbol("3/31/2013"), Symbol("6/30/2013"),
  Symbol("9/30/2013"), Symbol("12/31/2013"), Symbol("3/31/2014"), Symbol("6/30/2014"),
  Symbol("9/30/2014"), Symbol("12/31/2014"), Symbol("3/31/2015"), Symbol("6/30/2015"),
  Symbol("9/30/2015"), Symbol("12/31/2015"), Symbol("3/31/2016"), Symbol("6/30/2016"),
  Symbol("9/30/2016"), Symbol("12/31/2016"), Symbol("3/31/2017"), Symbol("6/30/2017"),
  Symbol("9/30/2017"), Symbol("12/31/2017"), Symbol("3/31/2018")

const FOCAL_DATES = Date("12/31/2010",  DATE_FORMAT):Dates.Month(3):Date("3/31/2018",  DATE_FORMAT)
const OUTCOMES = ["NORMALIZED_ACCRUALS_CF_METHOD"]

const MFloat64 = Union{Missing, Float64}
const MInt = Union{Missing, Int}
const MBool = Union{Missing, Bool}
const MSymbol = Union{Missing, Symbol}

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

#preprocess rate data
function preProcessChinaDSW()::Void
  wideDF::DataFrame = CSV.read("$DATA_PATH\\$CHINA_DS_NAME.csv", weakrefstrings=false)

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

#computes policy outcomes for the different exchanges
function getExchangeVolumeData(longDFV::DataFrame)::DataFrame

  wideDFV::DataFrame = unstack(longDFV, :date, :security, :value)
  N::Int = size(wideDFV,1)

  #sum buy and sell orders

  wideDFV[:C1Vol] = Vector{MFloat64}(wideDFV[:C1DBTO_Index] .+ wideDFV[:C1DSTO_Index])
  wideDFV[:C2Vol] = Vector{MFloat64}(wideDFV[:C2DBTO_Index] .+ wideDFV[:C2DSTO_Index])

  wideDFV[:lC1Vol] = Vector{MFloat64}((log).(wideDFV[:C1DBTO_Index]) .- (log).(wideDFV[:C1DSTO_Index]))
  wideDFV[:lC2Vol] = Vector{MFloat64}((log).(wideDFV[:C2DBTO_Index]) .- (log).(wideDFV[:C2DSTO_Index]))

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

    wideDFV[i, :lC1Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:lC1Vol]))

    wideDFV[i, :lC2Effect] = sum(skipmissing(wideDFV[(i-VOL_PERIODS+1):i,:lC2Vol]))
  end


  return wideDFV
end

#runs any checks they are required on the wideDF post data-quality
function wideDiagnostics(wideDF::DataFrame)::Void

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
  wideDF[:,newCol] = Vector{MSymbol}((Symbol).(wideDF[:,:subIndustry]))



  #start at the lowest level of industry classification
  countDF::DataFrame = by(wideDF, newCol, nrow)

  countIndex::Dict{Symbol, Int} = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))

  #For each industry, examine the minimum count for subindustries. If its less than desirable
  #then conform the subindustries to industries
  by(wideDF, :industry) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[:,newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[:,newCol] = Symbol(df[1,:industry])
    end
  end

  #repeat for industry groups and industries
  countDF = by(wideDF, newCol, nrow)
  countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
  by(wideDF, :industryGroup) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[:,newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[:,newCol] = Symbol(df[1,:industryGroup])
    end
  end

  #and again one more time to conform some industry groups to sectors
  countDF = by(wideDF, newCol, nrow)
  countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
  by(wideDF, :sector) do df::SubDataFrame

    subGICS::Vector{Symbol} = unique(df[:,newCol])
    minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

    #conform iff the min count is below the desired threshold
    if minCount < threshold
      df[:,newCol] = Symbol(df[1,:sector])
    end
  end

  #finally, if the threshold is sitll not met, set to missing (if the option is selected)
  if setMissing
    countDF = by(wideDF, newCol, nrow)
    countIndex = Dict(countDF[i,newCol]=>countDF[i,:x1] for i::Int ∈ 1:size(countDF,1))
    by(wideDF, :sector) do df::SubDataFrame

      subGICS::Vector{Symbol} = unique(df[:,newCol])
      minCount = minimum(((s::Symbol)->countIndex[s]).(subGICS))

      #conform iff the min count is below the desired threshold
      if minCount < threshold
        df[:,newCol] = missing
      end
    end
  end

  #reform the long dataframe
  if runWideDiagnostics
    wideDiagnostics(wideDF)
  end

  longDF = melt(wideDF, [nonDateCols; newCol])
  longDF[:,:date] = (s::Symbol->Date(string(s),dateFormat)).(longDF[:,:variable])
  delete!(longDF, :variable)
  longDF = longDF[!(ismissing).(longDF[:,:value]),:] # need to refilter missing variables
  return longDF
end


#conveience function to drop bad rows
function pruneLong!(longDF::DataFrame)::DataFrame

  #make a tracking column and store some convient quantities
  dateFormat::String = DATE_FORMAT
  N::Int = size(longDF,1)
  longDF[:,:keep] = Vector{MBool}(trues(N))
  numQuarters::Int =
      sum((s::Symbol->Date(string(s),dateFormat) ∈ minimum(FOCAL_DATES):maximum(
      FOCAL_DATES)).(DATE_SYMBOLS))

  #now process the firm-level data
  by(longDF, :security) do df::SubDataFrame
    nSub::Int = size(df,1)

    #println("$nSub")
    #println("nSub: $nSub numQuarters:$numQuarters")
    #check for missing values, inconsistent lengths, and stocks switching exchanges
    df[(ismissing).(df[:,:value]) |
        (ismissing).(df[:exchange]) |
        (ismissing).(df[:marketCap]), :keep] = false

    #drop a stock if it has less than two years of data, if its too small,
    # or its listed on multiple exchanges
    if ((sum(df[:, :keep] .== true) < MIN_COMPANY_QUARTERS) ||
        (length(unique(df[:,:exchange])) ≠ 1) ||
        (minimum(skipmissing(df[:marketCap])) < 5.) )

      df[:,:keep] = false

    end
  end

  longDF = longDF[longDF[:keep] .== true,:]
  ##finally need to check for  stocks listed multiple times
  by(longDF, :description) do dfName::SubDataFrame #for each company name

    by(dfName, :date) do dfUniqueDt #for each date
      if size(dfUniqueDt,1) > 1 #if there is a dup
        if KEEP_LARGEST_WEIGHT_IF_DUP #if true, we keep the largest weighted name
          maxWeight::Float64 = maximum(dfUniqueDt[:,:weight])
          dfUniqueDt[dfUniqueDt[:,:weight] .≠ maxWeight,:keep] = false
        else
          dfUniqueDt[:,:keep] = false
        end
      end
    end
  end

  longDF = longDF[longDF[:keep] .== true,:]
  delete!(longDF, :keep)

  #finally, do the winsorization
  for s::Symbol ∈ [:value, :marketCap, :ey]
    series::Vector{MFloat64} = longDF[:,s]

    #μ::Float64 = mean(skipmissing(series))
    #σ::Float64 = std(skipmissing(series))
    lower::Float64 = ACCEPTABLE_RANGE[s][1]
    upper::Float64 = ACCEPTABLE_RANGE[s][2]
    for i::Int ∈ 1:length(series)
      if !ismissing(series[i])
        if series[i] < lower
          series[i] = lower
        elseif series[i] > upper
          series[i] = upper
        end
      end
    end

    longDF[:, s] = series
  end




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
  longDF[:,:date] = (s::Symbol->Date(string(s),dateFormat)).(longDF[:,:variable])
  delete!(longDF, :variable)
  longDF = longDF[(d::Date->d ∈ minimum(FOCAL_DATES):maximum(FOCAL_DATES)).(longDF[:,:date]),:] #filter to the correct dates
  #println("names: $(names(longDF))")

  #clean up names
  longDF[:,:security] =
      (s::String->Symbol(replace(strip(s), " ","_"))).(longDF[:,:security])

  #drop and collect non-outcome rows. Start with a flag for outcome variables
  longDF[:,:Outcome] =
      (s::CategoricalString->s ∈ OUTCOMES ? true : false).(longDF[:command])
  longDFV::DataFrame =
      longDF[!longDF[:Outcome],[:date, :security, :value]]
  longDF = longDF[longDF[:Outcome],:]
  delete!(longDF, :Outcome) #no longer need this


  #flag for deleting rows with bad data
  longDF = longDF[!((ismissing).(longDF[:exchange])), :]
  longDF[:,:exchange] = Vector{MSymbol}((Symbol).(longDF[:,:exchange]))

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

  #now assign the treatment variables
  #begin by handling missing values and allocating space
  exchanges::Vector{MSymbol} = longDF[:exchange]

  #work with vectors for performance reasons
  C1Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x)?0.0:Float64(x)).(longDF[:, :C1Effect])
  C2Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x)?0.0:Float64(x)).(longDF[:, :C2Effect])
  CEffect::Vector{MFloat64} = Vector{MFloat64}(missings(N))
  lC1Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x)?0.0:Float64(x)).(longDF[:, :lC1Effect])
  lC2Effect::Vector{MFloat64} = (x::MFloat64->ismissing(x)?0.0:Float64(x)).(longDF[:, :lC2Effect])
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

    postTreatment[i] = MInt(abs(CEffect[i]) > 0.0)
  end

  longDF[:, :CEffect] = CEffect #focal policy variable
  longDF[:, :lCEffect] = lCEffect #focal policy variable
  longDF[:, :postTreatment] = postTreatment
  longDF[:, :gicsXDate] =
    Vector{MSymbol}(((s1::MSymbol,s2::Date)->ismissing(s1)||ismissing(s2)?missing:Symbol(
    s1,s2)).(longDF[:gics],longDF[:date]))

  longDF[:,:lMarketCap] = Vector{MFloat64}((log).(longDF[:,:marketCap]))

  categorical!(longDF, :security)
  categorical!(longDF, :exchange)
  categorical!(longDF, :gics)
  #  categorical!(longDF, :gicsInteracted)
  categorical!(longDF, :date)
  categorical!(longDF, :gicsXDate)

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
    focalName=:integrated)::String

  #showcols(longDF[:,[:lMarketCap, :value, :gics, :ey, :exchange, :postTreatment, :date]])


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
  push!(XSpecs, parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, parse("$focalSpec + exchange + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + exchange + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, parse("$focalSpec + exchange"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gicsXDate)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + exchange"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gicsXDate)


  ###################end new regression specs

  NSpecs::Int = length(XSpecs)
  columnNames::Vector{Vector{String}} = [((i::Int)->"($i)").(1:NSpecs)]

  println("Running Discrete Spec")
  for i::Int ∈ 1:NSpecs
      println("spec: $(XSpecs[i])")
      push!(models, CTLM(longDF, XSpecs[i],  YSpecs[i], withinSym = withinSpec[i],
          clusterSym = clusteredSpec[i], XNames=XNames[i], YName = YSpecs[i]))
  end

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Date F.E.",
      "Firm F.E.", "Exchange F.E.", "Industry F.E.", "Industry-Date F.E", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} =
      ((i::Int)->Vector{String}(length(models))).(1:length(descRowNames))

  #this builds the descriptive rows. There is Probably a better way to do this,
  #but its fairly specific to the project.
  for i ∈ 1:length(XSpecs)
    descContent[1][i] = "$(withinSpec[i]==:date||
        contains(string(XSpecs[i]), "date")?"X":"")"
    descContent[2][i] = "$(withinSpec[i]==:security||
        contains(string(XSpecs[i]), "security")?"X":"")"
    descContent[3][i] = "$(withinSpec[i]==:exchange||
        contains(string(XSpecs[i]), "exchange")?"X":"")"
    descContent[4][i] = "$(withinSpec[i]==:gics?"X":"")"
    descContent[5][i] = "$(withinSpec[i]==:gicsXDate?"X":"")"
    descContent[6][i] = "$(models[i].N)"
  end

  tableText::String = texTable(models,
      getClustered!#=getNeweyWestFunc(5) (lm::CTLM)->getNeweyWestSlow(lm, 5)=#,
      [focalName, :lMarketCap, :ey],
      titleCaption = "Discrete Specifications",
      colNames = columnNames,
      contentRowNames = ["Integrated", "lMarketCap", "E/Y"],
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
    focalName=:IMeasure)::String

  #showcols(longDF[:,[:lMarketCap, :value, :gics, :ey, :exchange, :postTreatment, :date]])


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
  push!(XSpecs, parse("$focalSpec + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :security)
  push!(clusteredSpec, :security)

  push!(XSpecs, parse("$focalSpec + exchange + date"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + exchange + date"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :gics)
  push!(clusteredSpec, :gics)

  push!(XSpecs, parse("$focalSpec + exchange"))
  push!(XNames, [:intercept; focalName])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gicsXDate)

  push!(XSpecs, parse("$focalSpec + lMarketCap + ey + exchange"))
  push!(XNames, [:intercept; focalName; :lMarketCap; :ey])
  push!(YSpecs, YSym)
  push!(withinSpec, :gicsXDate)
  push!(clusteredSpec, :gicsXDate)


  ###################end new regression specs

  NSpecs::Int = length(XSpecs)
  columnNames::Vector{Vector{String}} = [((i::Int)->"($i)").(1:NSpecs)]

  println("Running Policy Spec")
  for i::Int ∈ 1:NSpecs
      println("spec: $(XSpecs[i])")
      push!(models, CTLM(longDF, XSpecs[i],  YSpecs[i], withinSym = withinSpec[i],
          clusterSym = clusteredSpec[i], XNames=XNames[i], YName = YSpecs[i]))
  end

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Date F.E.",
      "Firm F.E.", "Exchange F.E.", "Industry F.E.", "Industry-Date F.E", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} =
      ((i::Int)->Vector{String}(length(models))).(1:length(descRowNames))

  #this builds the descriptive rows. There is Probably a better way to do this,
  #but its fairly specific to the project.
  for i ∈ 1:length(XSpecs)
    descContent[1][i] = "$(withinSpec[i]==:date||
        contains(string(XSpecs[i]), "date")?"X":"")"
    descContent[2][i] = "$(withinSpec[i]==:security||
        contains(string(XSpecs[i]), "security")?"X":"")"
    descContent[3][i] = "$(withinSpec[i]==:exchange||
        contains(string(XSpecs[i]), "exchange")?"X":"")"
    descContent[4][i] = "$(withinSpec[i]==:gics?"X":"")"
    descContent[5][i] = "$(withinSpec[i]==:gicsXDate?"X":"")"
    descContent[6][i] = "$(models[i].N)"
  end

  tableText::String = texTable(models,
      getClustered!#=getNeweyWestFunc(5) (lm::CTLM)->getNeweyWestSlow(lm, 5)=#,
      [focalName, :lMarketCap, :ey],
      titleCaption = "Policy Specifications",
      colNames = columnNames,
      contentRowNames = ["$focalName", "lMarketCap", "E/Y"],
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
  numColumns = 9

  #setup the table headers
  columnNames::Vector{Vector{String}} = [["Accrual Ratio", "Market Cap \$\\times 10^{-9}\$", "Earnings Yield"],
      ["Shanghai", "Shenzhen", "Total", "Shanghai", "Shenzhen",
      "Total", "Shanghai", "Shenzhen", "Total"]]
  columnSizes::Vector{Vector{Int}} = [ones(Int, numColumns ÷ 3) .* 3, ones(Int, numColumns)]

  #make views for conveneince
  longDFSH::SubDataFrame = view(longDF, longDF[:,:exchange] .== :Shanghai)
  longDFSZ::SubDataFrame = view(longDF, longDF[:,:exchange] .== :Shenzhen)

  ####################write descriptive rows
  descRowNames::Vector{String} = ["Mean", "Median", "Std Deviation",
    "10\\% Quantile", "90\\% Quantile", "N (firm-quarters)"]

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(numColumns)).(1:length(descRowNames))

  #get and scale the summary variables
  series::Vector{Vector{MFloat64}} = [
      longDFSH[:value], longDFSZ[:value], longDF[:value],
        longDFSH[:marketCap].*10.^-9.,
        longDFSZ[:marketCap].*10.^-9.,
        longDF[:marketCap].*10.^-9.,
      longDFSH[:ey], longDFSZ[:ey], longDF[:ey]]

  #remove missing
  for i::Int ∈ 1:length(series)
    #println("$(maximum(longDFSH[:value]))")
    series[i] = series[i][((x::MFloat64)->!ismissing(x)).(series[i])]
    #println("max: $(maximum(series[i]))")
  end

  #percentiles::Vector{Function} =
#    (s::Vector{MFloat64}->ecdf(Vector{Float64}(s))).(series)

  for c::Int ∈ 1:numColumns
    descContent[1][c] = num2Str(mean(series[c]), decimals)
    descContent[2][c] = num2Str(median(series[c]), decimals)
    descContent[3][c] = num2Str(std(series[c]), decimals)
    descContent[4][c] = num2Str(quantile(Vector{Float64}(series[c]),0.1), decimals)
    descContent[5][c] = num2Str(quantile(Vector{Float64}(series[c]),0.9), decimals)
    descContent[6][c] = num2Str(length(series[c]), Ints=true)
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
      alignmentColNames = [["c", "c", "c"], ["r" for i::Int ∈ 1:numColumns]],
      lineSpacer = "\\\\"
    )

  return summaryTable
end

function runRegressions(longDF::DataFrame; runDiscrete::Bool = true, runPolicy::Bool = true,
    runSummary::Bool = true)

  #run the regression with the discrete specification
  discreteTable::String = runDiscrete?discrete(longDF):""
  lpolicyTable::String = runPolicy?policy(longDF, focalSpec="lCEffect", focalName=:IMeasureNet):""
  policyTable::String = runPolicy?policy(longDF, focalSpec="CEffect", focalName=:IMeasureVol):""

  #calculate the summary statistics
  summaryTable::String = runSummary?summary(longDF):""

  #write the processed data to a CSV for any tests in R
  CSV.write("$R_TEST_PATH\\longDF.csv", longDF)

  #print tables
  writeTables2File([summaryTable, discreteTable, lpolicyTable, policyTable],
      HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
      outName = "trainTables.tex")
end
#main script
function mainScript(;refreshData::Bool = true, processData::Bool =true,
    runDiscrete::Bool = true, runPolicy::Bool = true, runSummary::Bool = false)::Void

  longDF::DataFrame = getDF(
      "$DATA_PATH\\$CHINA_DS_NAME_L.jls", processData, ()->processDF(refreshData))

  # run the regressions
  runRegressions(longDF, runDiscrete = runDiscrete, runPolicy = runPolicy, runSummary = runSummary)

  return nothing
end

@time begin #run main
  mainScript(refreshData = true, processData=true, runDiscrete=true,
      runPolicy=true, runSummary=true)
end

end #end module
