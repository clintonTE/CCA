module HFMod

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


using DataFrames, Distributions, StatsBase, GZip, Gadfly, JuMP, NLopt, Missings
using CTIO, CTReg, CTStat, CTVAR
#importall CT

#Set this to an optimal value- I find logical cores *(3/4) works well
#BLAS.set_num_threads(Int(round((Sys.CPU_CORES*3÷4)))-1)


#if !isdefined(:constDef)
const constDef = true

const DATA_PATH = pwd() * "\\data"
const DATA_FILE_NAME = "HFPerformance_2015-2017" #pad this by a year since we are interested in flows
const YEAR_RANGE = 2015:2016
const DATA_FILE_NAME_DEAD = DATA_FILE_NAME * "_dead"
const DATE_FORMAT_STR = "yyyymmdd"
const NOTICE_PERIOD_ADJ = 0 #days to add to the notice (since data are monthly)
const MIN_PERIODS_FOR_INCLUSION = 12 #number of periods required for inclusion
const T_TEST_THRESHOLD = 2.0 #this is a parameter for identificaiton if the flows are lagged
const PERIODS_SAME_AUM_TO_DELETE = 4 #assume 4 consecutive periods of identical aum implies stale data
const LAGS_NOTICE_TO_CAPTURE = 2
const LAGS_PERFORMANCE_TO_CAPTURE = 12
const MIN_ASSETS = 10.0 #10,000,000
const START_LAGS_FROM_NOTICE = false #start the performance lags from the notice period (not the redemption date)
const RUN_LM_INTERACTIONS = false #run the interaction regressions (takes a while)
const RUN_2SLS_INTERACTIONS =  false
const RESULTS_PATH = pwd() * "\\results"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const USE_AGGRESSIVE_GC = true #this enables very agressive GC which helps memory management at cost of performance
const R_TEST_PATH = "C:\\Users\\Clinton\\Dropbox\\Projects\\Summer17RTester\\Summer17RTester"
const MIN_MONTHS_NOTICE_IV = 1 #min months notice to be included in IV analysis
const MIN_MONTHS_NOTICE_LONG = 1 #min months notice to be included in the long file
const SIG_REDEMPTION_THRESHOLDS = 0.05:0.05:0.1
const PERFORMANCE_BOUND = 0.5
const MIN_FLOW = -0.99

const MFloat64 = Union{Missing, Float64}
const MInt = Union{Missing, Int}

#const CTExpr = Union{Symbol,Expr,Void}
#const CTSym = Union{Symbol,Void}


#end

###############HELPER FUNCTIONS

#helper function to extract all positive numbers from a string
#IN: A string
#OUT: A vector of floats, possibly empty
function getPosNumsFromString(s::String)
  len::Int = length(s)
  preParseS::Vector{Char} = collect(s)

  #replace non-numerical values with delimitor Ξ
  for i ∈ 1:len
     preParseS[i] = (isnumber(s[i]) || (
      s[i]=='.' && (i>1 && isnumber(s[i-1]) || (
        i<len && isnumber(s[i+1])))))?s[i]:'Ξ'
  end

  #split the tokens, keeping only numberical strings
  sSplit::Vector{String} = split(preParseS,'Ξ',keep=false)

  #parse and return the values as an array of floats
  return broadcast((x::String)->Float64(parse(x)),sSplit)::Vector{Float64}
end

#=helper function which sorts and splits DF on symbol, extracts all positive
numeric values, averages them, and replaces the column
IN: DataFrame, column symbol, a default value, and whetehr to leave NAs
OUT: None explicit. DataFrame modified with numberic column named col,
old column is named col_str=#
function DFStringtoMeanNum!(DF::DataFrame, col::Symbol;
  default::Float64=0.0, leaveMissings::Bool=false, silent::Bool=true)

  oldCol = Symbol(col, "_str")
  rename!(DF, col, oldCol)
  DF[:,col] = default
  sort!(DF, cols=(oldCol))

  by(DF, oldCol) do DFSub::SubDataFrame
    if !ismissing(DFSub[1,oldCol])
      parsed::Vector{Float64} = #parse the string
        getPosNumsFromString(DFSub[1,oldCol])
      if length(parsed)>0 #if string wasn't all garbage
        DFSub[:,col] = mean(parsed)
      elseif !silent#Give a heads up given a garbage incentive fee
        println("WARNING: Could not parse ", col, ". Value: ", DFSub[1, oldCol],
          ", Assuming numeric value of ", default)
      end
    elseif leaveMissings
      DFSub[:,col] = missing
    end
  end
end

#takes a vector and writes a string
#=function vec2String(v::Vector{T}, delim::W = "")::String where
  {T<:Any, W<:Union{AbstractString,Char}}

  return join(String.(v), delim)
end=#

#calculates the t statistic of the Pearson Correlation (Lowry 2017)
function tStatOfCor(x::Vector{Float64},y::Vector{Float64})
  r::Float64 = cor(x,y)
  return r/sqrt((1-r^2.0)/(length(x)-2))::Float64
end

#this is a type which stores all of the notice lags
struct NoticeLagSymbols
  noticeLags::Vector{Symbol}
  noticeLLags::Vector{Symbol}
  cumNoticeLags::Vector{Symbol}
  cumNoticeLLags::Vector{Symbol}
end

#Constructor for the above type
#IN: Number of lags
#OUT: A container of symbols with the appropraite lags
function NoticeLagSymbols(lags::Int)::NoticeLagSymbols
  noticeLags::Vector{Symbol} = Vector{Symbol}(lags)
  noticeLLags::Vector{Symbol} = similar(noticeLags)
  cumNoticeLags::Vector{Symbol} = similar(noticeLags)
  cumNoticeLLags::Vector{Symbol} = similar(noticeLags)

  for i ∈ 1:lags
    noticeLags[i] = Symbol("noticeLag$i")
    noticeLLags[i] = Symbol("noticeLLag$i")
    cumNoticeLags[i] = Symbol("cumNoticeLag$i")
    cumNoticeLLags[i] = Symbol("cumNoticeLLag$i")
  end

  return NoticeLagSymbols(noticeLags,noticeLLags,cumNoticeLags,
    cumNoticeLLags)
end

#this is a type which stores all of the performance lags
struct PerformanceLagSymbols
  performanceLags::Vector{Symbol} #regular performance lags and value-weighted if aggregated
end

#Constructor for the above type
#IN: Number of lags
#OUT: A container of symbols with the appropraite lags
function PerformanceLagSymbols(lags::Int)::PerformanceLagSymbols

  performanceLags::Vector{Symbol} = Vector{Symbol}(lags)   #preallocate
  for i ∈ 1:lags
    performanceLags[i] = Symbol("performanceLag$i")
  end

  return PerformanceLagSymbols(performanceLags)
end

#this structure holds symbols for cases where flows msut meet a certain threshold
struct RedemptionThresholds
  negLFlows::Vector{Symbol}
  noticePeriod::Vector{Symbol}
  monthsPastNotice::Vector{Symbol}
end

RedemptionThresholds(rng::Range=SIG_REDEMPTION_THRESHOLDS)::RedemptionThresholds =
  RedemptionThresholds([Symbol("negLFlowsT$(Int(t*10000))") for t ∈ collect(rng)],
  [Symbol("noticePeriodT$(Int(t*10000))") for t ∈ collect(rng)],
  [Symbol("monthsPastNoticeT$(Int(t*10000))") for t ∈ collect(rng)])

contains(e::T where T<:CTExpr, s::String)::Bool  = Base.contains("$e",s)
###########Main Methods
#U for unlagged, L for lagged, C for corrected, R for corrected but with
#a lagged default

function preProcessHFData(;flowsLaggingProcedure::Char = 'U', newBinary::Bool = false,
  dropFirstMonth::Bool = false)

  nPer::Int = 0 #numer of periods
  dtFormat::DateFormat = DateFormat(DATE_FORMAT_STR) ## a date format object

  #this allows us to capture a variable amount of lags
  noticeLS::NoticeLagSymbols = NoticeLagSymbols(LAGS_NOTICE_TO_CAPTURE)
  performanceLS::PerformanceLagSymbols = PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)

  #this makes sure we have a binary of the data (Improves load times 3-4x)
  if !isfile("$DATA_PATH\\$DATA_FILE_NAME.jls.gz") || newBinary

    #extract from the zip file, package into dataframe and serialize into zip
    #Saves a lot of time not to have to do this each time
    gHandleAlive::GZipStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME.gz")
    gHandleDead::GZipStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME_DEAD.gz")
    ogStream::GZipStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME.jls.gz","w")
    serialize(ogStream, [readtable(gHandleAlive); readtable(gHandleDead)])

    close(ogStream)
    close(gHandleAlive)
    close(gHandleDead)


  end

  #load the binary
  igStream::GZipStream= gzopen("$DATA_PATH\\$DATA_FILE_NAME.jls.gz")
  HFData::DataFrame = deserialize(igStream)
  close(igStream)

  println("Initial rows: ", size(HFData,1), ", Initial columns: ", size(HFData,2))
  println("Begining pre-processing...")

  #drop some columns
  HFData = HFData[:, [:fund_id, :main_strategy,
    :sub_strategy, :fund_status, :fund_assets, :firm_assets,
    :advance_notice, :fund_assets_as_of,
    :date, :performance, :nav, :assets]]

  #Drop some null valued fields:
  HFData = dropNullsFromDF(HFData, [:date, :performance, :assets,
  :advance_notice])

  HFData = HFData[HFData[:,:assets] .≥ MIN_ASSETS, :]
  HFData[:,:main_strategy] = Symbol.(HFData[:,:main_strategy])


  #drop unlikely outliers (data quality)
  HFData[:,:performance] ./= 100.0
  HFData = HFData[(HFData[:performance].< PERFORMANCE_BOUND) & (
    HFData[:performance].> (-1.0*PERFORMANCE_BOUND)),:]

  HFData[:,:date] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:date])

  #sort the data by fund id and date
  sort!(HFData,cols = (:fund_id,:date))

  #clean out some funds with stale aum data
  by(HFData, :fund_id) do HFSub::SubDataFrame
    nPer = size(HFSub,1)

    if nPer > PERIODS_SAME_AUM_TO_DELETE #check if this procedure makes sense
      identicalAUMs::Int = 1
      i::Int =  2
      while i ≤ nPer && identicalAUMs < PERIODS_SAME_AUM_TO_DELETE
        identicalAUMs = (HFSub[i,:assets]==HFSub[i-1,:assets])?identicalAUMs+1:1
        i+=1
      end
      if identicalAUMs == PERIODS_SAME_AUM_TO_DELETE #set assets to missing to flag for deletion
        HFSub[:,:assets] = missing
      end
    end
  end

  HFData =  HFData[!((ismissing).(HFData[:,:assets])),:] #delete the identified funds

  #delete funds with gaps in their performance history, or funds with a
  #track record less than six months
  HFData[:,:mcnt]=0
  ctr::Int = 0
  by(HFData, :fund_id) do HFSub::SubDataFrame
    mcnt::Int = size(HFSub,1)
    #now check for track record breaks

    if maximum(HFSub[:,:date]) ≠
      Dates.lastdayofmonth(minimum(HFSub[:,:date])+Dates.Month(mcnt)-Dates.Month(1))
      mcnt = 0

      #=if ctr<10
        showall(HFSub[:,[:fund_id,:date,:performance,:main_strategy, :fund_assets, :fund_status]])
      end
      ctr+=1=#
    end
    HFSub[:,:mcnt] = mcnt
  end
#  println("funds dropped: $ctr")
  HFData=HFData[~(HFData[:,:mcnt].<MIN_PERIODS_FOR_INCLUSION),:]

  #get the months of advance notice
  firstMonthCorrection::Int = dropFirstMonth?:1:0
  HFData[:,:monthsNotice] = Vector{MInt}(broadcast((x::Int)->
    max(0,((x+NOTICE_PERIOD_ADJ)÷30)-(firstMonthCorrection)
    )::Int,HFData[:,:advance_notice]))####get rid of MAX

  ###The following sections create nessecary columns with placeholders
  #calculate monthly flows and construct the data set
  fillSize::Int = size(HFData,1)

  HFData[:,:flows] = Vector{MFloat64}(fill(0.0,fillSize))
  HFData[:,:redempNotice] = Vector{MFloat64}(fill(0.0,fillSize)) #notifications about flows
  HFData[:,:rollingNotice] = Vector{MFloat64}(fill(0.0,fillSize)) #notifications about flows
  HFData[:,:rollingNoticeL] = Vector{MFloat64}(fill(0.0,fillSize)) #notifications about flows
  HFData[:,:lFlows]  = Vector{MFloat64}(fill(0.0,fillSize)) #will hold log flows
  HFData[:,:redempNoticeDt] = Date(1,1,1)

  #create the appropriate columns for lagged flow notice
  for i::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE
    HFData[:,noticeLS.noticeLags[i]] = Vector{MFloat64}(fill(0.0,fillSize))
    HFData[:,noticeLS.noticeLLags[i]] = Vector{MFloat64}(fill(0.0,fillSize))
    HFData[:,noticeLS.cumNoticeLags[i]] = Vector{MFloat64}(fill(0.0,fillSize))
    HFData[:,noticeLS.cumNoticeLLags[i]] = Vector{MFloat64}(fill(0.0,fillSize))
  end

  #create columns for performance lag
  for i::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
    #HFData[:,performanceLS.performanceLags[i]] = 0.0
    HFData[:,performanceLS.performanceLags[i]] = Vector{MFloat64}(fill(missing,fillSize))
  end

  HFData[:,:isNoticePeriod] = Vector{MFloat64}(fill(0.0,fillSize))
  HFData[:, :isRedemptionPeriod] = Vector{MFloat64}(fill(0.0,fillSize))
  HFData[:,:isNoticePeriodSingle] = Vector{MFloat64}(fill(0.0,fillSize))

  ctr = 0
  ###This is the main section for manipulating data by fund. We will
  #adjust the lagging of assets and assign to appropriate lagged variables
  by(HFData, :fund_id) do HFSub::SubDataFrame
    nPer = HFSub[1,:mcnt]
    useLaggingProc::Bool = false

    #now calculate the flows in the lagged and unlagged flows
    flowsUnlagged::Vector{Float64} =Vector{Float64}(nPer-1)
    flowsLagged::Vector{Float64} =Vector{Float64}(nPer-1)
    for i ∈ 2:nPer

      #back out the flows:
      #F_unlagged=A(t)-A(t-1)*(1+r(t))
      flowsUnlagged[i-1] = HFSub[i,:assets] -
      (HFSub[i-1,:assets])*(1.0+HFSub[i,:performance])

      #F_lagged=A(t)-A(t-1)*(1+r(t-1))
      flowsLagged[i-1]  = HFSub[i,:assets] -
      (HFSub[i-1,:assets])*(1.0+HFSub[i-1,:performance])
    end

    if flowsLaggingProcedure == 'U'  #assume all unlagged
      useLaggingProc = false

    elseif flowsLaggingProcedure == 'L' #assume all lagged
      useLaggingProc = true

      #below assumes a comparison process
    elseif flowsLaggingProcedure == 'C' || flowsLaggingProcedure == 'R'

      #get the differencein performance between lagged and unlagged
      perfDelta::Vector{Float64} = HFSub[2:end,:performance] .-
        HFSub[1:(end-1),:performance]

      #check the correlations
      flowsUnlaggedT::Float64 = tStatOfCor(flowsUnlagged,perfDelta)
      flowsLaggedT::Float64 = tStatOfCor(flowsLagged,perfDelta)

      #='C' Logic: Assume unlagged. If unlagged T test for correlation against
      return deltas is <-2 (or some other threshold), switch to lagged unless
      the lagged T stat is greater than 2
      'R' Logic: Same as C but reverse.=#
      if flowsUnlaggedT > -1.0 * T_TEST_THRESHOLD && flowsLaggedT > 1.0 * T_TEST_THRESHOLD
        useLaggingProc = false
      elseif flowsUnlaggedT < -1.0 * T_TEST_THRESHOLD && flowsLaggedT < 1.0 * T_TEST_THRESHOLD
        useLaggingProc = true
      elseif flowsUnlaggedT * flowsLaggedT ≤ 0.0
        if flowsLaggingProcedure == 'C'
          useLaggingProc = false
        elseif flowsLaggingProcedure == 'R'
          useLaggingProc = true
        end
      end
    end

    #now execute the appropriate lagging procedure
    if useLaggingProc
      HFSub[end,:flows] = missing
      HFSub[1:(end-1),:flows] = flowsLagged
      HFSub[1:(end-1),:lFlows] = flowsLagged ./
        (HFSub[2:end,:assets] .-  flowsLagged)
    else

      HFSub[1,:flows] = missing
      HFSub[1,:lFlows] = missing

      HFSub[2:end,:flows] = flowsUnlagged
      HFSub[2:end,:lFlows] = flowsUnlagged ./
        (HFSub[2:end,:assets] .- flowsUnlagged)
    end

  #now we look up the notice period to determine the appropriate month
  #Only will be for negative flows
  #we will make the table two way, at the expense of redundancy
  for i::Int ∈ 1:nPer
    #note notifications of flows does not include the current month's flows

    if !ismissing(HFSub[i,:flows]) && i > HFSub[i,:monthsNotice]
      noticeStart::Int = i-HFSub[i,:monthsNotice]+1-firstMonthCorrection
      if HFSub[i,:monthsNotice] ≥ MIN_MONTHS_NOTICE_IV & (
        noticeStart ≤ i) & (noticeStart ≥ 1)
        HFSub[noticeStart : i ,:rollingNotice] += HFSub[i,:flows]
      end


        if HFSub[i,:flows] < 0.0
          HFSub[i,:isNoticePeriod] = 1.0
          HFSub[i,:isNoticePeriodSingle] += 1.0
          HFSub[i,:isRedemptionPeriod] = 1.0

          #record the notification on the notification date
          HFSub[i-HFSub[i,:monthsNotice],:redempNotice] = HFSub[i,:flows]
          HFSub[i-HFSub[i,:monthsNotice],:redempNoticeDt] = HFSub[i,:date]

          #record the redemption notificaion lags and performance lags
          if HFSub[i,:monthsNotice] > 0
            for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, HFSub[i,:monthsNotice])
              HFSub[i-l, noticeLS.noticeLags[l]] = HFSub[i,:flows]
              HFSub[i-l, noticeLS.noticeLLags[l]] = HFSub[i,:lFlows]
            end
            for l::Int ∈ 1:(min(HFSub[i,:monthsNotice],i-1))
              HFSub[i-l,:isNoticePeriod] = 1.0
              HFSub[i-l,:isNoticePeriodSingle] += 1.0
            end
          end
        end
      end
    end


    #this si where we capture the cumulative notifications of redeumptions from prior months
    for i::Int ∈ 2:nPer

      HFSub[i,:rollingNoticeL] = max(MIN_FLOW,
        ismissing(HFSub[i,:rollingNotice])?missing:(HFSub[i,:rollingNotice] / HFSub[i,:assets]))

    end


    #now assign to the performance lags
    if START_LAGS_FROM_NOTICE #if true, we start from the end of the notice period
      for i::Int ∈ 2:nPer
        if i - 1 - HFSub[i,:monthsNotice] ≥ 1 #(make sure its positive number)
          for l::Int ∈ 1:min(i-1- HFSub[i,:monthsNotice],LAGS_PERFORMANCE_TO_CAPTURE)
            HFSub[i, performanceLS.performanceLags[l]] = HFSub[i - l, :performance]
          end
        end
      end
    else
      for i::Int ∈ 2:nPer
        for l::Int ∈ 1:min(i-1,LAGS_PERFORMANCE_TO_CAPTURE)
          HFSub[i, performanceLS.performanceLags[l]] = HFSub[i - l, :performance]
        end
      end
    end

  end###


  for l::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE
    HFData[HFData[:,:monthsNotice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = missing
    HFData[HFData[:,:monthsNotice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = missing
  end

  HFData[!((ismissing).(HFData[:,:lFlows])),:lFlows] =
    ((x::Float64)->(max(x,MIN_FLOW))).(
    HFData[!((ismissing).(HFData[:,:lFlows])),:lFlows])

  #get a seperate column with net negative flows
  HFData[:,:negFlows] = Vector{MFloat64}(-1.0*HFData[:,:flows])
  HFData[(!((ismissing).(HFData[:,:flows]))) & (HFData[:,:flows].≥0.0),:negFlows] = missing

  HFData[:,:negFlowsOrZero] = -1.0*HFData[:,:flows]
  HFData[(!((ismissing).(HFData[:,:flows]))) & (HFData[:,:flows].≥0.0),:negFlowsOrZero] = 0.0

  HFData[:,:negLFlows] = Vector{MFloat64}(-1.0*HFData[:,:lFlows])
  HFData[(!((ismissing).(HFData[:,:lFlows]))) & (HFData[:,:lFlows].≥0.0),:negLFlows] = missing

  HFData[:,:negLFlowsOrZero] = Vector{MFloat64}(-1.0*HFData[:,:lFlows])
  HFData[(!((ismissing).(HFData[:,:flows]))) & (HFData[:,:flows].≥0.0),:negLFlowsOrZero] = 0.0

  HFData[:,:negRollingNoticeL] = Vector{MFloat64}(-1.0*HFData[:,:rollingNoticeL])
  HFData[(!((ismissing).(HFData[:,:negRollingNoticeL]))) & (HFData[:,:negRollingNoticeL].<0.0),
    :negRollingNoticeL] = missing

  HFData[:,:negRollingNoticeLOrZero] = Vector{MFloat64}(-1.0*HFData[:,:rollingNoticeL])
  HFData[(!((ismissing).(HFData[:,:negRollingNoticeLOrZero]))) & (
    HFData[:,:negRollingNoticeLOrZero].<0.0), :negRollingNoticeLOrZero] = 0.0

  HFData[:,:negRollingNoticeLOrZeroSingle] = copy(HFData[:,:negRollingNoticeLOrZero])
  HFData[(HFData[:,:isNoticePeriodSingle] .> 1.0)#= | (HFData[:,:monthsNotice] .!= 1)=#,
    :negRollingNoticeLOrZeroSingle] = missing

  #showall(HFData[HFData[:,:rollingNoticeL] .< -100.0,[:fund_id,:date,:lFlows,:assets, :rollingNoticeL, :negRollingNoticeLOrZero,:negRollingNoticeLOrZeroSingle]] )
  # println(maximum(HFData[:,:rollingNoticeL]))
  # println(minimum(HFData[:,:rollingNoticeL]))

  #process the notice periods
  #Premise: Data not relevant if no notice period is present
  HFData[HFData[:,:monthsNotice] .< 1,:isNoticePeriod] = missing
  HFData[HFData[:,:monthsNotice] .< 1, :isNoticePeriodSingle] = missing

  #drop overlaps from the non-overalp specification
  HFData[!((ismissing).(HFData[:,:isNoticePeriodSingle])) & (HFData[:,:isNoticePeriodSingle] .> 1.0)
    , :isNoticePeriodSingle] = missing

  fillSize = size(HFData, 1)
  HFData[:,:isNoticePeriod1Mo] = Vector{MFloat64}(fill(0.0, fillSize))

  HFData[:,:isNoticePeriod1Mo] .= HFData[:,:isNoticePeriod]

  HFData[HFData[:,:monthsNotice] .!= 1,:isNoticePeriod1Mo]  = missing
  #make the lagged flows positive numbers
  for l::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE, sVec::Vector{Symbol} ∈
      [noticeLS.cumNoticeLags, noticeLS.cumNoticeLLags, noticeLS.noticeLags, noticeLS.noticeLLags]

      HFData[!ismissing(HFData[:,sVec[l]]),sVec[l]] =
          -1.0*HFData[!ismissing(HFData[:,sVec[l]]),sVec[l]]
  end

  #showall(HFData[3000:3200, [:fund_id, :performance, :assets, :flows, :lFlows, :negLFlows, :negLLFlows, :negRollingNoticeL]])

  #convert string columns to factors
  HFTypes::Vector{Type} = eltypes(HFData) #all types
  HFNames::Vector{String} = names(HFData) #all names
  for i ∈ 1:size(HFData,2) #for all columns
    if HFTypes[i] <: Union{String,Symbol} && !(typeof(HFData[:,i]) <: CategoricalArray)
      categorical!(HFData,Symbol(HFNames[i])) #convert to factor
    end
  end

  #create a factor column for dates
  HFData[:,:fDate] = HFData[:,:date]
  categorical!(HFData,:fDate)
  #filter out years that we don't want
  #Note we are filtering here so that the endpoints are better covered
  #with respect to the flows
  HFData = HFData[(((x::Date)->Dates.year(x) ∈ collect(YEAR_RANGE)).(HFData[:,:date])),:]

  println("Processing complete.")
  println("Ending rows: ", size(HFData,1), ", Ending columns: ", size(HFData,2))

  #showall(HFData[200:300,[:fund_id,:fDate,:monthsNotice,:lFlows,:isNoticePeriod, :isRedemptionPeriod, :isNoticePeriodSingle]])


  ogStream = gzopen(
    "$DATA_PATH\\$(DATA_FILE_NAME)_$(flowsLaggingProcedure)_$(dropFirstMonth)_clean.jls.gz","w")
  serialize(ogStream, HFData)
  close(ogStream)

  return nothing
end


#This function aggregates HF data on strategyID
#We will do this in the format of the standard wide file, so that
#we can plug it into the long file creation process
function aggregateOnStrategy(HFData::DataFrame; monthsNotice::T = missing
    )::DataFrame where T<:Union{Missing,Int}


  performanceLS::PerformanceLagSymbols =
      PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)

  ###now calculate the aggregates

  #pre-allocate the aggregate dataframe
  strategies::Vector{Symbol} = unique(HFData[:,:main_strategy])
  dates::CategoricalArray =  unique(HFData[:,:fDate])

  #for each date, have a row for each strategy and a row for the total
  HFAgg::DataFrame = HFData[1:((length(strategies)+1) * length(dates)),:]

  numFill::Int = size(HFAgg, 1)
  for i ∈ 1:size(HFAgg,2)
    HFAgg[:,i] = Vector{Union{eltype(HFAgg[:,i]),Missing}}(fill(missing,numFill))
  end

  ctr::Int = 0

  for s::Symbol ∈ [strategies; :total], d::CategoricalValue ∈ dates
      ctr += 1
      HFAgg[ctr,:main_strategy] = s
      HFAgg[ctr,:fDate] = d
      HFAgg[ctr,:monthsNotice] = monthsNotice
  end


  numFill = size(HFAgg, 1)
  HFAgg[:,:var] = Vector{MFloat64}(fill(0.0,numFill))
  #HFAgg[:,:totNegLFlows] = 0.0
  HFAgg[:,:totFlowsOnNotice] = Vector{MFloat64}(fill(0.0,numFill))
  HFAgg[:,:var] = missing
  #HFAgg[:,:totNegLFlows] = missing
  HFAgg[:,:totFlowsOnNotice] = missing


  #get views for each date
  by(HFData,:fDate) do HFSubPart::SubDataFrame
    HFAggSub::SubDataFrame = view(HFAgg, HFAgg[:,:fDate] .== HFSubPart[1,:fDate])

    stratTable::Dict = Dict(HFAggSub[i,:main_strategy]=>i for i::Int ∈ 1:length(strategies))

    for s::Symbol ∈ strategies
      if haskey(stratTable,s)
        sInd = stratTable[s]
        HFSub = view(HFSubPart, HFSubPart[:,:main_strategy].==s)

        #now fill in the aggregates
        HFAggSub[sInd,:assets] = sum(skipmissing(HFSub[:,:assets]))
        HFAggSub[sInd,:flows] = sum(skipmissing(HFSub[:,:flows]))
        HFAggSub[sInd,:lFlows] = HFAggSub[sInd,:flows] / HFAggSub[sInd,:assets]

        HFAggSub[sInd,:negLFlows] = HFAggSub[sInd,:lFlows]≤0.0?-1.0*HFAggSub[sInd,:lFlows]:missing
      #  HFAggSub[sInd,:totNegLFlows] = sum(skipmissing(HFSub[:,:negFlows]))
        HFAggSub[sInd,:totFlowsOnNotice] =
          sum(skipmissing(HFSub[HFSub[:,:monthsNotice].≥1,:negFlows]))
        #use value weighted or equal weighted method to calculate performance
        HFAggSub[sInd,:performance] =
          sum(skipmissing(HFSub[:,:performance].*HFSub[:,:assets])) / HFAggSub[sInd,:assets]
        HFAggSub[sInd,:var] =
          sum(skipmissing((HFSub[:,:performance].-HFAggSub[sInd,:performance]).^2 .*
          HFSub[:,:assets])) / HFAggSub[sInd,:assets]

        end
      end
      sInd = findfirst(HFAggSub[:,:main_strategy], :total)

      #now get the totals
      #HFAggSub[sInd,:totNegLFlows] = sum(skipmissing(HFSubPart[:,:negFlows]))
      HFAggSub[sInd,:totFlowsOnNotice] =
          sum(skipmissing(HFSubPart[HFSubPart[:,:monthsNotice].≥1,:negFlows]))
      HFAggSub[sInd,:assets] = sum(skipmissing(HFSubPart[:,:assets]))
      HFAggSub[sInd,:flows] = sum(skipmissing(HFSubPart[:,:flows]))
      HFAggSub[sInd,:performance] =
          sum(skipmissing(HFSubPart[:,:assets] .* HFSubPart[:,:performance])) / HFAggSub[sInd,:assets]
      HFAggSub[sInd,:var] =
          sum(skipmissing((HFSubPart[:,:performance].-HFAggSub[sInd,:performance]).^2 .*
          HFSubPart[:,:assets])) / HFAggSub[sInd,:assets]

  end


  #now get the lags and assign a "fund" id
  strategyID::Int  = 1
  sort!(HFAgg,cols = [:main_strategy,:fDate])
  for s::Symbol ∈ strategies
      HFAggSub::SubDataFrame = view(HFAgg, HFAgg[:,:main_strategy] .== s)
      HFAggSub[:,:fund_id] = strategyID
      for l::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
        @simd for i  ∈  2:size(HFAggSub,1)
          if i-l>0
            HFAggSub[i,performanceLS.performanceLags[l]] = HFAggSub[i-l,:performance]
          end
        end
      end
      strategyID += 1
  end

  return HFAgg
end


#returns a table with the total flows
#High performance requires an initial sort
function aggregateOnDate(HFData::DataFrame, lags::Int=5,
  noFirstNoticeMonth::Bool=false)

  nDates::Int = size(unique(HFData[:,:date]),1)
  noticeLS::NoticeLagSymbols = NoticeLagSymbols(lags)

  #preallocate data frame
  HFByDate::DataFrame =
    DataFrame([Union{Date,Missing},
      Union{Int,Missing},
      Union{Float64,Missing},
      Union{Float64,Missing},
      Union{Float64,Missing},
      Union{Float64,Missing},
      Union{Float64,Missing}],
    [:date, :numFunds, :totAUM, :totFlows, :totRedemp,
    :totRollingRedemp, :totRedempNotice], nDates)
  HFByDate[:,:date] .= unique(HFData[:,:date])

  #allocate the lag columns
  numFill::Int = size(HFByDate,1)
  for s ∈ [noticeLS.noticeLags; noticeLS.cumNoticeLags]

    HFByDate[:,s] = Vector{MFloat64}(fill(0.0,numFill))
  end


  #make a lookup table
  dateTbl = Dict(HFByDate[i,:date] => i for i::Int ∈ 1:nDates)
  #might be a better way to do this with aggregates and joins,
  #although this is reasonably fast with a low number of comparisons
  #Basically we are getting totals for each symbol by a bunch of symbols by date
  by(HFData, :date) do HFSub::SubDataFrame
    dateIdx::Int = dateTbl[HFSub[1,:date]]
    applyToDate::Vector{Bool} = HFByDate[:,:date].==HFSub[1,:date]
    HFByDate[dateIdx,:numFunds]=size(unique(HFSub[:,:fund_id]),1)
    HFByDate[dateIdx,:totAUM]=sum(HFSub[!(ismissing).(HFSub[:,:assets]),:assets])
    HFByDate[dateIdx,:totFlows]=sum(HFSub[!(ismissing).(HFSub[:,:flows]),:flows])
    HFByDate[dateIdx,:totRedemp]=
      sum(HFSub[broadcast((x)->ifelse(!ismissing(x)&&x<0,true, false),HFSub[:,:flows]),:flows])
    HFByDate[dateIdx,:totRedempNotice]=
      sum(HFSub[!ismissing(HFSub[:,:redempNotice]),:redempNotice])
    HFByDate[dateIdx,:totRollingRedemp]=
      sum(HFSub[!ismissing(HFSub[:,:rollingNotice]),:rollingNotice])
    for s ∈ [noticeLS.noticeLags; noticeLS.cumNoticeLags]
      HFByDate[dateIdx,s]= sum(HFSub[!ismissing(HFSub[:,s]),s])
    end

  end

  return HFByDate

end

#this executes the linear model tests for a NAV valuation effect
function LM(HFData::DataFrame, titleCaption::String, caption::String;
    flowStr::String = "negLFlows", dependentVar=:performance,
    flowStrDesc::String = "Outflows / AUM")

    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)
        #negLLFlows
    XLMSpecs::Vector{CTExpr} = [parse("$flowStr")]
    XLMNames::Vector{Vector{Symbol}} = [[:intercept; :outflowVar]] #note only non-factor names
    FLMCluster::Vector{CTSym} = [nothing]

    #showall(HFData[1:50, [:fund_id, :main_strategy, :negLFlows, :performance]])

    #####put additional specifications here
    push!(XLMSpecs, parse("$flowStr"))
    push!(XLMNames, [:intercept; :outflowVar])
    push!(FLMCluster,:fDate)

    push!(XLMSpecs, parse("$flowStr+fDate"))
    push!(XLMNames, [:intercept; :outflowVar])
    push!(FLMCluster,:main_strategy)

    push!(XLMSpecs, parse("$flowStr+fDate"))
    push!(XLMNames, [:intercept; :outflowVar])
    push!(FLMCluster,:fund_id)

    push!(XLMSpecs, parse("$flowStr + " *
        "$(vec2String(performanceLS.performanceLags, "+"))"))
    push!(XLMNames, [:intercept; :outflowVar; performanceLS.performanceLags])
    push!(FLMCluster,nothing)

    push!(XLMSpecs, parse("$flowStr + " *
        "$(vec2String(performanceLS.performanceLags, "+"))"))
    push!(XLMNames, [:intercept; :outflowVar; performanceLS.performanceLags])
    push!(FLMCluster,:fDate)

    push!(XLMSpecs, parse("$flowStr + " *
    "$(vec2String(performanceLS.performanceLags, "+")) + fDate"))
    push!(XLMNames, [:intercept; :outflowVar; performanceLS.performanceLags])
    push!(FLMCluster,:main_strategy)

    push!(XLMSpecs, parse("$flowStr + " *
    "$(vec2String(performanceLS.performanceLags, "+")) + fDate"))
    push!(XLMNames, [:intercept; :outflowVar; performanceLS.performanceLags])
    push!(FLMCluster,:fund_id)

    modelsLM::Vector{CTLM} = Vector{CTLM}()

    println("Running LM Spec X: $flowStr Y: performance")
    #run the LM Specs
    for i ∈ 1:length(XLMSpecs)
        #println("Model Spec: $(XLMSpecs[i])")
        push!(modelsLM, CTLM(HFData, XLMSpecs[i], fixedEffectsSym = FLMCluster[i],
            dependentVar, XNames=XLMNames[i], YName = :performance))
    end
    ##########LM IO
    descRowNamesLM::Vector{String} = ["Performance Lags", "Date Fixed Effects",
        "Strategy Fixed Effects","Fund Fixed Effects", "N ('000s)"]

    #need to print the descriptive rows
    descContentLM::Vector{Vector{String}} =
        [Vector{String}(length(modelsLM))for i∈ 1:length(descRowNamesLM)]

    #this builds the descriptive rows. There is Probably a better way to do this,
    #but its fairly specific to the project.
    for i ∈ 1:length(XLMSpecs)
        descContentLM[1][i] = "$(contains(XLMSpecs[i], "performanceLag")?"X":"")"
        descContentLM[2][i] = "$(FLMCluster[i]==:fDate||contains(XLMSpecs[i], "fDate")?"X":"")"
        descContentLM[3][i] = "$(FLMCluster[i]==:main_strategy||contains(XLMSpecs[i], "main_strategy")?"X":"")"
        descContentLM[4][i] = "$(FLMCluster[i]==:fund_id?"X":"")"
        descContentLM[5][i] = "$(round(modelsLM[i].N/1000.0,1))"
    end


    TableTextLM::String = texTable(modelsLM, getModWhiteΣ!, [:outflowVar],
        titleCaption = titleCaption,
        colNames = [["($i)" for i::Int ∈ 1:length(modelsLM)]],
        contentRowNames = ["$flowStrDesc"],
        descRowNames = descRowNamesLM,
        descContent = descContentLM,
        decimalDigits = 2,
        #columnSepPt = -15,
        scaling = [100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = caption
        )


    return TableTextLM
end

#this specification runs a conventional panel using dummies for redemption and notice periods
#this executes the linear model tests for a NAV valuation effect
function PM(HFData::DataFrame, titleCaption::String, caption::String;
    focalStr= "isNoticePeriod", dependentVar=:performance, cluster=false,
    coefs::T=nothing)::String where T<:Union{Void,Vector{Float64}}

    XPMSpecs::Vector{CTExpr} = [parse("$focalStr")]
    XPMNames::Vector{Vector{Symbol}} = [[:intercept; :isNoticePeriod]] #note only non-factor names
    FPMCluster::Vector{CTSym} = [nothing]
    #####put additional specifications here

    push!(XPMSpecs, parse("$focalStr"))
    push!(XPMNames, [:intercept; :isNoticePeriod])
    push!(FPMCluster,:fDate)

    push!(XPMSpecs, parse("$focalStr + fDate"))
    push!(XPMNames, [:intercept; :isNoticePeriod])
    push!(FPMCluster,:main_strategy)

    push!(XPMSpecs, parse("$focalStr + fDate"))
    push!(XPMNames, [:intercept; :isNoticePeriod])
    push!(FPMCluster,:fund_id)

    push!(XPMSpecs, parse("$focalStr + isRedemptionPeriod"))
    push!(XPMNames, [:intercept; :isNoticePeriod; :isRedemptionPeriod])
    push!(FPMCluster,nothing)

    push!(XPMSpecs, parse("$focalStr + isRedemptionPeriod"))
    push!(XPMNames, [:intercept; :isNoticePeriod; :isRedemptionPeriod])
    push!(FPMCluster,:fDate)

    push!(XPMSpecs, parse("$focalStr + isRedemptionPeriod + fDate"))
    push!(XPMNames, [:intercept; :isNoticePeriod; :isRedemptionPeriod])
    push!(FPMCluster,:main_strategy)

    push!(XPMSpecs, parse("$focalStr + isRedemptionPeriod + fDate"))
    push!(XPMNames, [:intercept; :isNoticePeriod; :isRedemptionPeriod])
    push!(FPMCluster,:fund_id)


    #####

    if !cluster
        FPMCluster[:] .= nothing
    end

    modelsPM::Vector{CTLM} = Vector{CTLM}()
    #run the PM Specs

    println("Running PM Spec X:$focalStr Y:$dependentVar")
    for i ∈ 1:length(XPMSpecs)
        push!(modelsPM, CTLM(HFData, XPMSpecs[i], fixedEffectsSym = FPMCluster[i],
            dependentVar, XNames=XPMNames[i], YName = :performance))
    end

    ##########PM IO
    descRowNamesPM::Vector{String} = ["Date F.E.",
        "Strategy F.E.", "Fund-level F.E","N ('000s)"]

    #need to print the descriptive rows
    descContentPM::Vector{Vector{String}} =
        [Vector{String}(length(modelsPM))for i∈ 1:length(descRowNamesPM)]

    #this builds the descriptive rows. There is Probably a better way to do this,
    #but its fairly specific to the project.
    for i ∈ 1:length(XPMSpecs)
        descContentPM[1][i] = "$(FPMCluster[i]==:fDate||contains(XPMSpecs[i], "fDate")?"X":"")"
        descContentPM[2][i] = "$(FPMCluster[i]==:main_strategy||contains(XPMSpecs[i], "main_strategy")?"X":"")"
        descContentPM[3][i] = "$(FPMCluster[i]==:fund_id?"X":"")"
        descContentPM[4][i] = "$(round(modelsPM[i].N/1000.0,1))"

    end


    TableTextPM::String = texTable(modelsPM, getModWhiteΣ!, [:isNoticePeriod; :isRedemptionPeriod],
      titleCaption = titleCaption,
      colNames = [["($i)" for i∈1:length(modelsPM)]],
      contentRowNames = ["Notice Period", "Redemption Dt"],
      descRowNames = descRowNamesPM,
      descContent = descContentPM,
      decimalDigits = 3,
      #columnSepPt = -27,
      scaling = [100.0,100.0],
      #clearMem = USE_AGGRESSIVE_GC,
      caption = caption)

      #get the coefficients we want
      if coefs!=nothing
        for m::CTLM ∈ modelsPM
          push!(coefs,getTerm(m,:isNoticePeriod))
        end
      end


    return TableTextPM
end

function Summary(HFData::DataFrame)::String
    summaryDecimals::Int = 2


    #list all columns for which we will provide summary stats
    summaryCols::Vector{Symbol} =
        [:assets, :performance, :negLFlows, :monthsNotice]
    summaryNumCols::Int = length(summaryCols)

    #make an index for ease of reference
    summaryDict::Dict = Dict(summaryCols[i] => i for i::Int ∈ 1:summaryNumCols)
    summaryDat::Vector{Vector{Float64}} =
        [HFData[!(ismissing).(HFData[:,s]),s] for s::Symbol ∈ summaryCols]

    #scale for readability
    summaryDat[summaryDict[:performance]] .*= 100.0
    summaryDat[summaryDict[:negLFlows]] .*= 100.0


    #elliminate zero outflow values from the vector for summary stat purposes
    filter!((x::Float64)->x>0, summaryDat[summaryDict[:negLFlows]])
    summaryRowNames::Vector{String} =
        ["mean", "median", "\$\\sigma\$","Number of Funds",
        "Number of Points","Months (2000-2016)"]

    summaryNumRows::Int = length(summaryRowNames)
    summaryContent::Vector{Vector{String}} =
        [Vector{String}() for i∈ 1:summaryNumRows]

    #holds the width of each data cell in columns
    summaryWidth::Vector{Vector{Int}} =
        [Vector{Int}() for i∈ 1:summaryNumRows]

    #get the mean
    summaryContent[1] = string.([
        round(mean(summaryDat[i]),summaryDecimals) for i::Int ∈ 1:summaryNumCols])
    summaryWidth[1] = ones(summaryNumCols)

    summaryContent[2] = string.([
        round(median(summaryDat[i]),summaryDecimals) for i::Int ∈ 1:summaryNumCols])
    summaryWidth[2] = ones(summaryNumCols)

    #get the standard deviation
    summaryContent[3] = string.([
        round(std(summaryDat[i]),summaryDecimals) for i::Int ∈ 1:summaryNumCols])
    summaryWidth[3] = ones(summaryNumCols)

    summaryContent[4] = ["$(length(unique(HFData[:,:fund_id]))) " *
        "(Currently Active: $(length(unique(HFData[HFData[:,:fund_status].==
        "Active",:fund_id]))))"]
    summaryWidth[4] = [summaryNumCols]

    summaryContent[5] = ["$(length(HFData[:,:fund_id]))"]
    summaryWidth[5] = [summaryNumCols]

    summaryContent[6] = ["$(length(unique(HFData[:,:fDate])))"]
    summaryWidth[6] = [summaryNumCols]


    TableSummary::String = texTable("Summary of Data",
        """Includes alive and dead funds from 2000 to 2016.
        Initial Data set consisted of 1.6mn rows, reduced by preprocessing to
        approximately 380k.  Pre-processing
        filtered out funds with discontinuous data, monthly performance magnitude >50\\%, AUM < 10mn,
        track record < 12 months, 4 or more identical AUMs in a row (stale data), and missing
        relevant fields. Flows calculated via \$AUM_{t}-(1+r_t)AUM_{t-1}\$ and further
        experessed as a percentage of the relevant AUM. A t-test on
        the correlation of return and flows was used to adjust the timing of reported AUM, with
        a critical value of \$2\\sigma\$ implying an adjustment to timing (main results are robust
        to not making this adjustment.) Seven days were added to the notice period to account for
        informal notice, which was then divided by 30 to calculate months' notice.""", #caption
        [["AUM(\\\$mn)", "Performance \\%", "Redemptions \\%AUM", "Months' Notice"]], #colNames
        Vector{String}(),#contentRowNames
        Vector{Matrix{String}}(), #content
        summaryRowNames, #descRowNames
        summaryContent, #descContent
        Vector{String}(), #notes
        widthDescContent = summaryWidth,
        #columnSepPt=5
        )

    plOutRD = plot(HFData, x=:assets,
        Geom.density,
        Guide.title("Assets Historgram"),
        Guide.xlabel("AUM"),Guide.ylabel("Count"))
      draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_assets.png", 7.5inch, 6inch), plOutRD)

    return TableSummary
end

function joinDataAgg(HFData::DataFrame, HFDataAgg::DataFrame)::DataFrame
  sort!(HFDataAgg, cols=[:fDate, :main_strategy])

  dateSymTable =  Dict((HFDataAgg[i,:fDate], HFDataAgg[i,:main_strategy]) =>
      HFDataAgg[i,:performance] for i::Int∈1:size(HFDataAgg,1))

  numFill::Int = size(HFData,1)
  HFData[:,:performanceStrat] = Vector{MFloat64}(fill(0.0,numFill))
  HFData[:,:performanceRel] = Vector{MFloat64}(fill(0.0,numFill))

  HFData[:,:performanceStrat] = missing
  HFData[:,:performanceRel] = missing

  warnCtr::Int = 0
  @simd for i ∈ 1:size(HFData,1)
      if haskey(dateSymTable, (HFData[i,:fDate], HFData[i,:main_strategy]))
          #gets the strategy performance and first difference of the performance
          HFData[i,:performanceStrat] = dateSymTable[(HFData[i,:fDate], HFData[i,:main_strategy])]

          #get the benchmark relative performance and first difference
          HFData[i,:performanceRel] = HFData[i,:performance] -
              HFData[i,:performanceStrat]
      elseif warnCtr≤10
          error("Join Failure on row $i. Row Info:\n $(HFLong[i,[:fDate,
              :performance,:monthsNotice,:main_strategy]])")
          #warnCtr += 1
      end
  end

  return HFData
end

#this function gets the estimate of canceled flows
#NOTE: depending on the implementation of dataframe, it may modify HFData
function getCanceledEstimate!(HFData::DataFrame, titleCaption::String, caption::String;
  focalStr::String="isNoticePeriod", dependentVar=:performance)::Vector{String}


  HFDataModel = HFData[!(ismissing).(HFData[:,:isNoticePeriod]),:]
  dates::Vector{Date} = unique(HFData[:,:date])
  N::Int = length(dates)

  #first, construct the models for the clustering
  XSpecs::Vector{CTExpr} = [parse("$focalStr + fDate")]
  XNames::Vector{Vector{Symbol}} = [[:intercept; :isNoticePeriod]] #note only non-factor names
  FCluster::Vector{CTSym} = [nothing]
  #####put additional specifications here

  push!(XSpecs, parse("$focalStr + fDate"))
  push!(XNames, [:intercept; :isNoticePeriod])
  push!(FCluster,:main_strategy)

  push!(XSpecs, parse("$focalStr + fDate"))
  push!(XNames, [:intercept; :isNoticePeriod])
  push!(FCluster,:fund_id)

  models::Vector{CTLM} = Vector{CTLM}()

  #keep track of the variable names
  M::Int = length(XSpecs)
  μSyms::Vector{Symbol} = [Symbol("μ$i") for i::Int ∈ 1:M]
  σ2Syms::Vector{Symbol} = [Symbol("σ2$i") for i::Int ∈ 1:M]
  YSyms::Vector{Symbol} = [Symbol("Y$i") for i::Int ∈ 1:M]
  thresholdSyms::Vector{Symbol} = [Symbol("threshold$i") for i::Int ∈ 1:M]
  thresholdMinusμSyms::Vector{Symbol} = [Symbol("thresholdMinusμ$i") for i::Int ∈ 1:M]
  totSubmittedSyms::Vector{Symbol} = [Symbol("totSubmitted$i") for i::Int ∈ 1:M]
  canceledSyms::Vector{Symbol} = [Symbol("canceled$i") for i::Int ∈ 1:M]


  #run the model specs

  for i ∈ 1:M
      println("Running Spec $i for redemption calc")
      push!(models, CTLM(HFDataModel, XSpecs[i], fixedEffectsSym = FCluster[i],
          dependentVar, XNames=XNames[i], YName = :performance))

      #check to make sure we didn't drop a row
      if (models[end].N) != size(HFDataModel[:,:fDate],1)
        error("Size Mismatch spec $i- model rows: $(models[end].N), data rows: $(size(HFData[:,:fDate],1))")
      end
      HFDataModel[:,YSyms[i]] = models[end].Y
  end


  HFFlowDat::DataFrame = DataFrame([Date, Float64, Float64, Float64, Float64],
    [:fDate, :totFlowsOnNotice, :totNegFlows, :assets, :assetsModel], N)

  HFFlowDat[:,:fDate] = dates

  #write the variable column names into the dataframe
  flowRows::Int = size(HFFlowDat,1)
  for i::Int ∈ 1:M
    HFFlowDat[:,μSyms[i]]  = Vector{MFloat64}(flowRows)
    HFFlowDat[:,σ2Syms[i]] = Vector{MFloat64}(flowRows)
    HFFlowDat[:,thresholdSyms[i]]  = Vector{MFloat64}(flowRows)
    HFFlowDat[:,thresholdMinusμSyms[i]]  = Vector{MFloat64}(flowRows)
    HFFlowDat[:,totSubmittedSyms[i] ] = Vector{MFloat64}(flowRows)
    HFFlowDat[:,canceledSyms[i]] = Vector{MFloat64}(flowRows)

    HFFlowDat[:,μSyms[i]] = missing
    HFFlowDat[:,σ2Syms[i]] = missing
    HFFlowDat[:,thresholdSyms[i]]  = missing
    HFFlowDat[:,thresholdMinusμSyms[i]]  = missing
    HFFlowDat[:,totSubmittedSyms[i]]  = missing
    HFFlowDat[:,canceledSyms[i]] = missing
  end

  #now aggregate by date from the larger set
  by(HFData, :fDate) do HFSub::SubDataFrame
    r = findfirst(HFFlowDat[:,:fDate], HFSub[1,:fDate])

    #aggregate the flows
    HFFlowDat[r,:totNegFlows] = sum(skipmissing(HFSub[:,:negFlows]))
    HFFlowDat[r,:totFlowsOnNotice] =
        sum(skipmissing(HFSub[HFSub[:,:monthsNotice].≥1,:negFlows]))
    HFFlowDat[r,:assets] = sum(skipmissing(HFSub[:,:assets]))

    #this is the case where we are not demeaning the performance and variance
    #beyond the specific date
    for i::Int ∈ 1:M
      if FCluster[i] == nothing
        HFFlowDat[r,μSyms[i]] =
            sum(skipmissing(HFSub[:,dependentVar].*HFSub[:,:assets])) / HFFlowDat[r,:assets]
        HFFlowDat[r,σ2Syms[i]] =
          sum(skipmissing((HFSub[:,dependentVar].-HFFlowDat[r,μSyms[i]]).^2 .*
            HFSub[:,:assets])) / HFFlowDat[r,:assets]
      end
    end
  end

  #where we have the demeaned info
  by(HFDataModel, :fDate) do HFSub::SubDataFrame
    r = findfirst(HFFlowDat[:,:fDate], HFSub[1,:fDate])
    HFFlowDat[r,:assetsModel] = sum(skipmissing(HFSub[:,:assets]))

    #now aggregate the performance

    for i::Int ∈ 1:M
      if FCluster[i] != nothing
        HFFlowDat[r,μSyms[i]] =
            sum(skipmissing(HFSub[:,Symbol("Y$i")].*HFSub[:,:assets])) / HFFlowDat[r,:assetsModel]
        HFFlowDat[r,σ2Syms[i]] =
          sum(skipmissing((HFSub[:,Symbol("Y$i")].-HFFlowDat[r,Symbol("μ$i")]).^2 .*
          HFSub[:,:assets])) / HFFlowDat[r,:assetsModel]
      end
    end
  end

  #this is a check to make sure we aggregated correctly
  #=d::Date = HFFlowDat[5,:fDate]
  assetsChk = sum(HFData[HFData[:,:fDate].==d,:assets])
  μChk = sum(HFData[HFData[:,:fDate].==d,:performance].*
    HFData[HFData[:,:fDate].==d,:assets])/assetsChk
  σChk = sum((HFData[HFData[:,:fDate].==d,:performance]-μChk).^2.0 .*
    HFData[HFData[:,:fDate].==d,:assets])/assetsChk

  assetsChkModel = sum(HFDataModel[HFDataModel[:,:fDate].==d,:assets])
  μChkModel = sum(skipmissing(HFDataModel[HFDataModel[:,:fDate].==d,Symbol("Y3")].*
    HFDataModel[HFDataModel[:,:fDate].==d,:assets]))/assetsChkModel
  σChkModel = sum((HFDataModel[HFDataModel[:,:fDate].==d,Symbol("Y3")]-μChkModel).^2.0 .*
    HFDataModel[HFDataModel[:,:fDate].==d,:assets])/assetsChkModel

  showall(HFFlowDat[5,:])
  println("\nassets Check: ",assetsChk , " μCheck:", μChk, " σCheck: ",σChk)
  println("assetsModel Check: ",assetsChkModel ," μCheckModel:", μChkModel, " σCheckModel: ",σChkModel)=#


  #define some functions for the solver
  const sqrt2Inv::Float64 = 1/(2.0^.5)
  const sqrt2PiInv::Float64 = 1.0/(2.0*π)^0.5

  #need to define special functions here for the solver
  CTCumNorm(z) = 0.5+0.5*erf(z*sqrt2Inv)
  CTNorm(z) = sqrt2PiInv*exp(-z^2.0/2.0)

  for m ∈ 1:M

    println("Solving Values in Model $m: ")

    for r::Int ∈ 1:N
      #print(HFFlowDat[r,:fDate],", ")

        #form the objective function
      μ::Float64 = HFFlowDat[r,μSyms[m]]
      σ::Float64 = HFFlowDat[r,σ2Syms[m]]^0.5
      β::Float64 = getTerm(models[m],Symbol(focalStr))
      function obj(t)
          z = (t-μ)/σ
          ϵ = β + σ*CTNorm(z)/CTCumNorm(z)
          return ϵ*ϵ
      end

      #build the optimizer
      #Algorithms: :LD_LBFGS-fails, :LN_COBYLA-(too slow),
      #  :LD_SLSQP (works!), :LD_TNEWTON (works!), :LD_VAR2/:LD_VAR1 (fails),
      #  PRAXIS (fails/slow), :LD_MMA (fails/slow),
      mod = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
      JuMP.register(mod, :CTNorm, 1, CTNorm, autodiff=true)
      JuMP.register(mod, :CTCumNorm, 1, CTCumNorm, autodiff=true)
      @variable(mod, t)

      JuMP.register(mod, :obj, 1, obj, autodiff=true)
      @NLobjective(mod, Min, obj(t))
      output = solve(mod)

      #calculate the formulas
      HFFlowDat[r,thresholdSyms[m]] = getvalue(t)
      HFFlowDat[r,thresholdMinusμSyms[m]] = HFFlowDat[r,thresholdSyms[m]] - μ
      #below is total = notice/Φ(z)
      HFFlowDat[r,totSubmittedSyms[m]]  = HFFlowDat[r,:totFlowsOnNotice]/
          CTCumNorm((HFFlowDat[r,thresholdSyms[m]] - μ)/σ) +
          (HFFlowDat[r,:totNegFlows]- HFFlowDat[r,:totFlowsOnNotice])
      HFFlowDat[r,canceledSyms[m]] = HFFlowDat[r,totSubmittedSyms[m]] -
        HFFlowDat[r,:totNegFlows]
    end

  end


  #aggregate to the annual level
  HFFlowDat[:,:year] = (Dates.year).(HFFlowDat[:,:fDate])
  years::Vector{Float64} = unique(HFFlowDat[:,:year])
  HFFlowDatAnn = aggregate(HFFlowDat[:,[:year; :assets; :totFlowsOnNotice;
    :totNegFlows; thresholdSyms; thresholdMinusμSyms;totSubmittedSyms; canceledSyms]],
    :year, [sum,mean])

  #rename the aggregated columns
  for i∈1:M
    HFFlowDatAnn[:,thresholdSyms[i]] = HFFlowDatAnn[:,Symbol("threshold$(i)_mean")]
    HFFlowDatAnn[:,thresholdMinusμSyms[i]] = HFFlowDatAnn[:,Symbol("thresholdMinusμ$(i)_mean")]
    HFFlowDatAnn[:,totSubmittedSyms[i]] = HFFlowDatAnn[:,Symbol("totSubmitted$(i)_sum")]
    HFFlowDatAnn[:,canceledSyms[i]] = HFFlowDatAnn[:,Symbol("canceled$(i)_sum")]
  end

  HFFlowDatAnn[:,:assets] = HFFlowDatAnn[:,:assets_sum]
  HFFlowDatAnn[:,:totFlowsOnNotice] = HFFlowDatAnn[:,:totFlowsOnNotice_sum]
  HFFlowDatAnn[:,:totNegFlows] = HFFlowDatAnn[:,:totNegFlows_sum]

  HFFlowDatAnn = unique(HFFlowDatAnn[:,[:year; :assets; :totFlowsOnNotice;:totNegFlows;
    thresholdSyms; thresholdMinusμSyms; totSubmittedSyms; canceledSyms]])

  #calculate a couple basic ratios
  cancToAssetsSyms::Vector{Symbol} = [Symbol("cancToAssets$i") for i::Int ∈ 1:M]
  cancToTotalSyms::Vector{Symbol} = [Symbol("cancToTotal$i") for i::Int ∈ 1:M]

  for i::Int ∈ 1:M
    HFFlowDatAnn[:,cancToAssetsSyms[i]] =
      HFFlowDatAnn[:,canceledSyms[i]] ./ HFFlowDatAnn[:,:assets]
    HFFlowDatAnn[:,cancToTotalSyms[i]] =
      HFFlowDatAnn[:,canceledSyms[i]] ./ HFFlowDatAnn[:,totSubmittedSyms[i]]
  end

  #Make the tables

  flowTablesText::Vector{String} =  Vector{String}()
  numYears::Int = size(HFFlowDatAnn,1)

  for i::Int ∈ 1:M
    #need tp build a table
    decimals = 3

    #write the appropriate title
    modelTitle::String = " Spec ($i): \$\\beta=$(round(getTerm(models[i],Symbol(focalStr))*100,2))\$"
    modelTitle = titleCaption*modelTitle
    if FCluster[i]==:main_strategy
      modelTitle = modelTitle*" (Includes Strategy F.E.)"
    elseif FCluster[i]==:fund_id
      modelTitle = modelTitle*" (Includes Fund F.E.)"
    end

    flowContent::Vector{Vector{String}} =
        [Vector{String}() for r∈ 1:numYears]

    #the net threshold is only meaningful with absolute, not relative performance
    if dependentVar == :performance
      colNames::Vector{Vector{String}} =
        [["Avg Threshold", "Avg Threshold - \$\\mu\$", "Canceled",
        "Canceled/Planned (\\%)", "Canceled/Assets (\\%)"]]
        for r ∈ 1:numYears
          flowContent[r] = ([string(round(HFFlowDatAnn[r,thresholdSyms[i]]*100.0,decimals)),
            string(round(HFFlowDatAnn[r,thresholdMinusμSyms[i]]*100.0,decimals)),
            string(Int(round(HFFlowDatAnn[r,canceledSyms[i]]))),
            string(round(HFFlowDatAnn[r,cancToTotalSyms[i]]*100.0,decimals)),
            string(round(HFFlowDatAnn[r,cancToAssetsSyms[i]]*100.0,decimals))])
        end
      else
        colNames =
          [["Avg Threshold", "Canceled",
          "Canceled/Planned (\\%)", "Canceled/Assets (\\%)"]]

        for r ∈ 1:numYears
          flowContent[r] = ([string(round(HFFlowDatAnn[r,thresholdSyms[i]]*100.0,decimals)),
            string(Int(round(HFFlowDatAnn[r,canceledSyms[i]]))),
            string(round(HFFlowDatAnn[r,cancToTotalSyms[i]]*100.0,decimals)),
            string(round(HFFlowDatAnn[r,cancToAssetsSyms[i]]*100.0,decimals))])
        end
      end

    rowNames::Vector{String} = string.(HFFlowDatAnn[:,:year])
    numCols = length(colNames[1])

    #println(i)

    push!(flowTablesText, texTable(modelTitle,
        caption, #caption
        colNames, #colNames
        Vector{String}(),#contentRowNames
        Vector{Matrix{String}}(), #content
        rowNames, #descRowNames
        flowContent, #descContent
        Vector{String}(), #notes
        #widthDescContent = summaryWidth
        ))

  end


  return flowTablesText

end

#this is the main method for analysis
function analyzeAggregates(;flowsLaggingProcedure::Char = 'U', dropFirstMonth::Bool = false,
  verbose::Bool = true, runSummary::Bool = true, runLM::Bool = true,
  runPM::Bool = true, process=true, save=false )



  #wrapping this up to save time if the file is long
  if process
    #open the pre-processed data file
    igStream::GZipStream =
        gzopen(
          "$DATA_PATH\\$(DATA_FILE_NAME)_$(flowsLaggingProcedure)_$(dropFirstMonth)_clean.jls.gz")
    HFData::DataFrame = deserialize(igStream)
    close(igStream)

    #get the long files
    HFDataAgg::DataFrame = aggregateOnStrategy(HFData, monthsNotice=MIN_MONTHS_NOTICE_LONG) #get wide strategy aggregate

    HFData = joinDataAgg(HFData,HFDataAgg)

    titleCaption::String = string()
    caption::String = string()

    sort!(HFData,cols = (:date, :fund_id))
    HFByDate::DataFrame = aggregateOnDate(HFData, LAGS_NOTICE_TO_CAPTURE) #get date aggregate

    if save
      ogStream::GZipStream = gzopen(
        "$DATA_PATH\\$(DATA_FILE_NAME)_$(flowsLaggingProcedure)_$(dropFirstMonth)_Proc.jls.gz","w")
      serialize(ogStream, [HFData, HFDataAgg,  HFByDate])
      close(ogStream)
    end
  else # if we arn't going to process, load in the processed files
      igStream =
          gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(flowsLaggingProcedure)_$(dropFirstMonth)_Proc.jls.gz")
      HFData, HFDataAgg, HFByDate = deserialize(igStream)
      close(igStream)
  end

  ###################Plots
  Gadfly.push_theme(:dark)
  plOut = plot(HFData, x=:negLFlowsOrZero,  #=shape=:variable, size=[.1mm],=#
      Geom.histogram,
      Guide.title("Fund Flows / AUM"),
      Guide.ylabel("Density"),Guide.xlabel("OutFlows"))
    draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_Flows.png", 7.5inch, 6inch), plOut)

  ##################OLS Spec

  titleCaption = "Performance vs Redemption Size"
  caption = """ This table illustrates the relationship between redemptions atand
      performance. The key takeway is performance on the redemption date is,
      counterintuitivelly, positively correlated with the size of the redemption.
      Date fixed effects are by month. Performance lags includes the prior
      24 months of performance as a covariate. Data is from alive and dead hedge
      funds 2000-2016. Where fixed effects are present, the standard errors are
      clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6).
       Outflows were expressed as a percentage of the
      AUM. A positive value for the coefficient implies a positive relationship
      between performance and redemptions."""

  TableTextLM::String = runLM?LM(HFData, titleCaption, caption,
      flowStr = "negLFlowsOrZero", flowStrDesc = "Outflows / AUM"):""

  titleCaption = "Performance vs Cumulative Rolling Redemption Notice Size"
  caption = """This table reflects the relationship between pending notifications
          of redemptions and performance. To avoide bias in the data,
          the first day Date fixed effects are by month. Performance lags includes the prior
          24 months of performance as a covariate. Data is from alive and dead hedge
          funds 2000-2016. Where fixed effects are present, the standard errors are
          clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6).
          Outflows were expressed as a percentage of the
          AUM. A positive value for the coefficient implies a positive relationship
          between performance and redemptions."""

  TableTextLMNotice::String = runLM?LM(HFData, titleCaption, caption,
    flowStr = "negRollingNoticeLOrZero", flowStrDesc = "Rolling Notice / AUM"):""
  #########IV Note: See the postIVO backup for the IV code
  ##############Panel Approach

  coefs::Vector{Float64}=Vector{Float64}()

  titleCaption = "Effect of Notice Period on Performance"
  caption = "This specification uses a dummy for the notice period and a dummy
  for the redemption period. In each specification, funds underperform in the
  notice period. Except with fund-level fixed effects,
  hedge funds realize  further losses on the redemption day. Data is from alive and dead hedge
  funds 2000-2016. Interpret a value of 1.0 for the dummy as corresponding
  to a difference in performance of 1.0\\%. Where fixed effects are present, the standard errors are
  clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
  TableTextPM = runPM?PM(HFData,titleCaption,caption, cluster=true, coefs=coefs):""

  titleCaption = "Effect of Notice Period on Performance (BM-Relative)"
  caption = "This is the same specification as the prior table, only
  performance is net of a strategy-specific value-weighted benchmark. Results
  are generally robust. In each specification, funds underperform in the
  notice period. Except with fund-level fixed effects,
  hedge funds realize  further losses on the redemption day. Data is from alive and dead hedge
  funds 2000-2016. Interpret a value of 1.0 for the dummy as corresponding
  to a difference in performance of 1.0\\%. Where fixed effects are present, the standard errors are
  clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
  TableTextPMRel = runPM?PM(HFData,titleCaption,caption, dependentVar=:performanceRel, cluster=true):""

  titleCaption = "Effect of Notice Period on Performance (No-overlap)"
  caption = "A robustness check on underperformance in the notice period.
  Here, overlapping redemptions are dropped from the data. Again, in each
  specification, funds underperform in the notice period. Except with fund-level fixed effects,
  hedge funds realize  further losses on the redemption day. Data is
  from alive and dead hedge funds 2000-2016. Interpret a value of 1.0 for
  the dummy as corresponding to a difference in performance of 1.0\\%. Where
  fixed effects are present, the standard errors are clustered either by fund
  (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""

  TableTextPMSingle = runPM?PM(HFData,titleCaption,caption,
    focalStr="isNoticePeriodSingle", cluster=true):""

  titleCaption = "Effect of Notice Period on Performance (No-overlap, BM Relative)"
  caption = "This specification drops overlapping redemptions and uses performance
  net of the value-weighted benchmark. Except with fund-level fixed effects,
  hedge funds realize
  further losses on the redemption day. Data is
  from alive and dead hedge funds 2000-2016. Interpret a value of 1.0 for
  the dummy as corresponding to a difference in performance of 1.0\\%. Where
  fixed effects are present, the standard errors are clustered either by fund
  (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
  TableTextPMSingleRel = runPM?PM(HFData,titleCaption,caption,
      dependentVar=:performanceRel,focalStr="isNoticePeriodSingle", cluster=true):""

  titleCaption = "Effect of Notice Period on Performance (1mo Notice Only)"
  caption = "This specification drops overlapping redemptions and limits the
  data set to one months' notice. The purpose is to estimate the amount of
  canceled redemptions. Except with fund-level fixed effects,
  hedge funds realize further losses on the redemption day. Data is
  from alive and dead hedge funds 2000-2016. Interpret a value of 1.0 for
  the dummy as corresponding to a difference in performance of 1.0\\%. Where
  fixed effects are present, the standard errors are clustered either by fund
  (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
  TableTextPMSingle1Mo = runPM?PM(HFData,titleCaption,caption,
    dependentVar=:performance,focalStr="isNoticePeriod1Mo", cluster=true):""

  ###################Calculate canceled assets

  titleCaption = "Canceled Flows"
  caption = """This table shows the implied total amount of redemptions which were
    submitted to the hedge fund and canceled due to an improvement in
    performance. Data is from alive and dead hedge funds 2000-2016.  \$\\beta\$ is the impact of a
    notice period on performance as detailed in prior tables. ``Avg threshold''
    corresponds to the calculated level of performance at which an investor
    would have canceled their redemption to realize the \$\\beta\$ value, assuming
    a true \$\\beta\$ value of zero. The average difference between the threshold
    and the monthly return is also presented. ``Canceled''
    is the derived estimate of canceled redemptions based on the threshold value.
    ``Canceled/Total'' is the ratio of canceled redemptions to total submitted
    redemptions, while ``Canceled/Assets'' corresponds to canceled redemptions
    to to total fund assets. """

  flowTablesText = getCanceledEstimate!(HFData,titleCaption,caption)

  titleCaption = "Canceled Flows (Benchmark Relative)"
  caption = """This table shows the implied total amount of redemptions which were
    submitted to the hedge fund and canceled due to an improvement in
    performance. Data is from alive and dead hedge funds 2000-2016, while performance
    used is relative to a value-weighted benchmark.  \$\\beta\$ is the impact of a
    notice period on performance as detailed in prior tables. ``Avg threshold''
    corresponds to the calculated level of performance at which an investor
    would have canceled their redemption to realize the \$\\beta\$ value, assuming
    a true \$\\beta\$ value of zero. The average difference between the threshold
    and the monthly return is also presented. ``Canceled''
    is the derived estimate of canceled redemptions based on the threshold value.
    ``Canceled/Total'' is the ratio of canceled redemptions to total submitted
    redemptions, while ``Canceled/Assets'' corresponds to canceled redemptions
    to to total fund assets. """
  flowTablesText =
    [flowTablesText; getCanceledEstimate!(HFData,titleCaption,caption,dependentVar=:performanceRel)]


  ####################Do the Summary Tables and write all table to a file

  TableSummary::String = runSummary?Summary(HFData::DataFrame):""

  writeTables2File([TableSummary; TableTextLM; TableTextLMNotice;
    TableTextPM; TableTextPMRel; TableTextPMSingle; TableTextPMSingleRel;
    TableTextPMSingle1Mo; flowTablesText],
    HEADER_NAME, FOOTER_NAME, path=RESULTS_PATH,
    outName = "RegTables_$(flowsLaggingProcedure)_$(dropFirstMonth).tex") #cosmetic marker block for summary tables


end

function verify2R(;flowsLaggingProcedure::Char = 'U', dropFirstMonth::Bool=false)

  #read the seria
  igStream::GZipStream =
      gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(flowsLaggingProcedure)_$(dropFirstMonth)_clean.jls.gz")
  HFData::DataFrame = deserialize(igStream)
  close(igStream)

  #write the csv
  writetable("$R_TEST_PATH\\RHFData.csv",HFData)

end


@time begin
FLPChar = 'C'
dropFirstMonth = false

#preProcessHFData(flowsLaggingProcedure=FLPChar, newBinary=false, dropFirstMonth = dropFirstMonth)
analyzeAggregates(flowsLaggingProcedure=FLPChar, dropFirstMonth = dropFirstMonth,
  runSummary=true, runLM = false, runPM = false, process=false, save = true)
verify2R(flowsLaggingProcedure=FLPChar, dropFirstMonth =  dropFirstMonth)
end

#=@time for c ∈ ['C', 'R', 'L', 'U']
  PreProcessHFData(flowsLaggingProcedure=c, newBinary=true)
  AnalyzeAggregates(flowsLaggingProcedure=c)
end=#




end
