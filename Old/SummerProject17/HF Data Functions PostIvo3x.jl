
module HFMod

if length(findin("C:\\Users\\Clinton\\Dropbox\\Projects\\SummerProject17",LOAD_PATH)) == 0
  push!(LOAD_PATH,"C:\\Users\\Clinton\\Dropbox\\Projects\\SummerProject17")
end

#=TODO:
1) IV=#
#2) Run specifications

using DataFrames, Distributions, StatsBase, GZip, JLD, Gadfly, JuMP,CTMod,NLopt
#importall CT

#Set this to an optimal value- I find logical cores *(3/4) works well
BLAS.set_num_threads(Int(round((Sys.CPU_CORES*3÷4)))-1)


#if !isdefined(:constDef)
  const constDef = true
  #other constants here

  const DATA_PATH = pwd() * "\\data"
  const DATA_FILE_NAME = "HFPerformance_2000-2016"
  const DATA_FILE_NAME_DEAD = DATA_FILE_NAME * "_dead"
  const DATE_FORMAT_STR = "yyyymmdd"
  const NOTICE_PERIOD_ADJ = 7 #days to add to the notice (since data are monthly)
  const MIN_PERIODS_FOR_INCLUSION = 12 #number of periods required for inclusion
  const T_TEST_THRESHOLD = 2.0 #this is a parameter for identificaiton if the flows are lagged
  const PERIODS_SAME_AUM_TO_DELETE = 4 #assume 4 consecutive periods of identical aum implies stale data
  const LAGS_NOTICE_TO_CAPTURE = 2
  const LAGS_PERFORMANCE_TO_CAPTURE = 24
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
  const PERFORMANCE_BOUND = .5


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
  default::Float64=0.0, leaveNAs::Bool=false, silent::Bool=true)

  oldCol = Symbol(col, "_str")
  rename!(DF, col, oldCol)
  DF[:,col] = default
  sort!(DF, cols=(oldCol))

  by(DF, oldCol) do DFSub::SubDataFrame
    if ~isna(DFSub[1,oldCol])
      parsed::Vector{Float64} = #parse the string
        getPosNumsFromString(DFSub[1,oldCol])
      if length(parsed)>0 #if string wasn't all garbage
        DFSub[:,col] = mean(parsed)
      elseif !silent#Give a heads up given a garbage incentive fee
        println("WARNING: Could not parse ", col, ". Value: ", DFSub[1, oldCol],
          ", Assuming numeric value of ", default)
      end
    elseif leaveNAs
      DFSub[:,col] = NA
    end
  end
end

#takes a vector and writes a string
function vec2String(v::Vector{T}, delim::W = "")::String where
  {T<:Any, W<:Union{AbstractString,Char}}

  return join(String.(v), delim)
end

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

function preProcessHFData(;FlowsLaggingProcedure::Char = 'U', newBinary::Bool = false)

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
      if identicalAUMs == PERIODS_SAME_AUM_TO_DELETE #set assets to NA to flag for deletion
        HFSub[:,:assets] = NA
      end
    end
  end

  HFData =  HFData[~isna(HFData[:,:assets]),:] #delete the identified funds

  #delete funds with gaps in their performance history, or funds with a
  #track record less than six months
  HFData[:,:mcnt]=0
  by(HFData, :fund_id) do HFSub::SubDataFrame
    mcnt::Int = size(HFSub,1)
    #now check for track record breaks
    if maximum(HFSub[:,:date]) ≠
      Dates.lastdayofmonth(minimum(HFSub[:,:date])+Dates.Month(mcnt)-Dates.Month(1))
      mcnt = 0
    end
    HFSub[:,:mcnt] = mcnt
  end

  HFData=HFData[~(HFData[:,:mcnt].<MIN_PERIODS_FOR_INCLUSION),:]

  #get the months of advance notice
  HFData[:,:monthsNotice] = broadcast((x::Int)->
    ((x+NOTICE_PERIOD_ADJ)÷30)::Int,HFData[:,:advance_notice])

  ###The following sections create nessecary columns with placeholders
  #calculate monthly flows and construct the data set
  HFData[:,:flows] = 0.0 #flows
  HFData[:,:redempNotice] = 0.0 #notifications about flows
  HFData[:,:rollingNotice] = 0.0 #notifications about flows
  HFData[:,:rollingNoticeL] = 0.0 #notifications about flows
  HFData[:,:lFlows]  = 0.0 #will hold log flows
  HFData[:,:redempNoticeDt] = Date(1,1,1)

  #create the appropriate columns for lagged flow notice
  for i::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE
    HFData[:,noticeLS.noticeLags[i]] = 0.0
    HFData[:,noticeLS.noticeLLags[i]] = 0.0
    HFData[:,noticeLS.cumNoticeLags[i]] = 0.0
    HFData[:,noticeLS.cumNoticeLLags[i]] = 0.0
  end

  #create columns for performance lag
  for i::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
    HFData[:,performanceLS.performanceLags[i]] = 0.0
    HFData[:,performanceLS.performanceLags[i]] .= NA
  end

  HFData[:,:isNoticePeriod] = 0.0
  HFData[:, :isRedemptionPeriod] = 0.0
  HFData[:,:isNoticePeriodSingle] = 0.0

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

    if FlowsLaggingProcedure == 'U'  #assume all unlagged
      useLaggingProc = false

    elseif FlowsLaggingProcedure == 'L' #assume all lagged
      useLaggingProc = true

      #below assumes a comparison process
    elseif FlowsLaggingProcedure == 'C' || FlowsLaggingProcedure == 'R'

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
        if FlowsLaggingProcedure == 'C'
          useLaggingProc = false
        elseif FlowsLaggingProcedure == 'R'
          useLaggingProc = true
        end
      end
    end

    #now execute the appropriate lagging procedure
    if useLaggingProc
      HFSub[end,:flows] = NA
      HFSub[1:(end-1),:flows] = flowsLagged
      HFSub[1:(end-1),:lFlows] = flowsLagged ./
        (HFSub[2:end,:assets] .-  flowsLagged)
    else
      HFSub[1,:flows] = NA
      HFSub[1,:lFlows] = NA
      HFSub[2:end,:flows] = flowsUnlagged
      HFSub[2:end,:lFlows] = flowsUnlagged ./
        (HFSub[2:end,:assets] .- flowsUnlagged)
    end



    #now we look up the notice period to determine the appropriate month
    #Only will be for negative flows
    #we will make the table two way, at the expense of redundancy
    for i::Int ∈ 1:nPer
        if !isna(HFSub[i,:flows]) && i > HFSub[i,:monthsNotice]

            #note notifications of flows does not include the current month's flows
            if HFSub[i,:monthsNotice] ≥ MIN_MONTHS_NOTICE_IV
                HFSub[(i-HFSub[i,:monthsNotice]):(i-1),:rollingNotice] += HFSub[i,:flows]
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

        HFSub[i,:rollingNoticeL] = isna(HFSub[i,:rollingNotice])?NA:HFSub[i,:rollingNotice] / HFSub[i,:assets]

        for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, i-1)
            #we only record the lag if the total notices are net negative
            if (!isna(HFSub[i-l,:rollingNotice])) && (HFSub[i-l,:rollingNotice] < 0.0)
                HFSub[i,noticeLS.cumNoticeLags[l]] = HFSub[i-l,:rollingNotice]
                HFSub[i,noticeLS.cumNoticeLLags[l]] = HFSub[i-l,:rollingNotice] / HFSub[i-l,:assets]
            end
        end
    end

    #lagged flows values at the beginning of the series are genuinly missing
    for i::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE
        for l::Int ∈ LAGS_NOTICE_TO_CAPTURE:i
            HFSub[i,noticeLS.cumNoticeLags[l]] = NA
            HFSub[i,noticeLS.cumNoticeLLags[l]] = NA
        end
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
  end

    for l::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE
      HFData[HFData[:,:monthsNotice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = NA
      HFData[HFData[:,:monthsNotice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = NA
    end

    #get a seperate column with net negative flows
    #Also, take a log of the flows
    #HFData[:,:llFlows] = log.(1 .+ HFData[:,:lFlows] )
    HFData[:,:negFlows] = -1.0*HFData[:,:flows]
    HFData[(~isna(HFData[:,:flows])) & (HFData[:,:flows].≥0.0),:negFlows] = NA

    HFData[:,:negLFlows] = -1.0*HFData[:,:lFlows]
    HFData[(~isna(HFData[:,:lFlows])) & (HFData[:,:lFlows].≥0.0),:negLFlows] = NA
    #=HFData[:,:negLLFlows] = log.(1 .+ HFData[:,:negLFlows] )
    HFData[(~isna(HFData[:,:lFlows])) & (HFData[:,:lFlows].≥0.0),:negLLFlows] = NA=#

    HFData[:,:negRollingNoticeL] = -1.0*HFData[:,:rollingNoticeL]
    HFData[(~isna(HFData[:,:negRollingNoticeL])) & (HFData[:,:negRollingNoticeL].<0.0),
        :negRollingNoticeL] = NA


    #process the notice periods
    #Premise: Data not relevant if no notice period is present
    HFData[HFData[:,:monthsNotice] .< 1,:isNoticePeriod] = NA
    HFData[HFData[:,:monthsNotice] .< 1, :isNoticePeriodSingle] = NA

    #drop overlaps from the non-overalp specification
    HFData[HFData[:,:isNoticePeriodSingle] .> 1.0, :isNoticePeriodSingle] = NA

    HFData[:,:isNoticePeriod1Mo] = 0.0
    HFData[:,:isNoticePeriod1Mo] .= HFData[:,:isNoticePeriod]
    HFData[HFData[:,:monthsNotice] .!= 1,:isNoticePeriod1Mo]  = NA

    #make the lagged flows positive numbers
    for l::Int ∈ 1:LAGS_NOTICE_TO_CAPTURE, sVec::Vector{Symbol} ∈
        [noticeLS.cumNoticeLags, noticeLS.cumNoticeLLags, noticeLS.noticeLags, noticeLS.noticeLLags]

        HFData[~isna(HFData[:,sVec[l]]),sVec[l]] =
            -1.0*HFData[~isna(HFData[:,sVec[l]]),sVec[l]]
    end

    #showall(HFData[3000:3200, [:fund_id, :performance, :assets, :flows, :lFlows, :negLFlows, :negLLFlows, :negRollingNoticeL]])

  #convert string columns to factors
  HFTypes::Vector{Type} = eltypes(HFData) #all types
  HFNames::Vector{String} = names(HFData) #all names
  for i ∈ 1:size(HFData,2) #for all columns
    if HFTypes[i] <: Union{String,Symbol} && !(typeof(HFData[:,i]) <: PooledDataVector)
      pool!(HFData,Symbol(HFNames[i])) #convert to factor
    end
  end

  #create a factor column for dates
  HFData[:,:fDate] = HFData[:,:date]
  pool!(HFData,:fDate)


  println("Processing complete.")
  println("Ending rows: ", size(HFData,1), ", Ending columns: ", size(HFData,2))

  #showall(HFData[200:300,[:fund_id,:fDate,:monthsNotice,:lFlows,:isNoticePeriod, :isRedemptionPeriod, :isNoticePeriodSingle]])

  #save the binary
  ogStream = gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls.gz","w")
  serialize(ogStream, HFData)
  close(ogStream)

end

function processLongHFData(HFData::DataFrame; negFlowsOnly=true)::DataFrame

    println("Beginning lengthening of data. \nStarting Wide rows: $(size(HFData,1))")


    #get any variable arrays of symbols
    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)
    noticeLS::NoticeLagSymbols = NoticeLagSymbols(LAGS_NOTICE_TO_CAPTURE)
    thresholdLS::RedemptionThresholds = RedemptionThresholds(SIG_REDEMPTION_THRESHOLDS)

    #store this for convenience
    arraySyms::Vector{Symbol} = [:fund_id; :performance; performanceLS.performanceLags;
        :negLFlows; :monthsNotice; :fDate; :main_strategy]


    if negFlowsOnly #Generally only interested in redemptions. Make a view so no accidental modification
        HFSub::SubDataFrame = view(HFData[:, arraySyms], ~isna(HFData[:,:negLFlows]) & (HFData[:,:negLFlows] .> 0.0))
    else #however, we may need this if we are merging the data set later
        HFData[:,:negLFlows] = -1.0*HFData[:,:lFlows]
        HFSub = view(HFData[:, arraySyms], ~isna(HFData[:,:negLFlows]))
    end

    println("Ending wide rows: $(size(HFData,1))")
    #Stack the performance over each redemption
    HFLong::DataFrame = melt(HFSub, [:negLFlows,:monthsNotice, :fDate, :main_strategy, :fund_id])
    println("Starting long rows: $(size(HFLong,1))")


    #define placeholders for the following variables, placeholders represent the redemption date
    HFLong[:,:monthsPastNotice] = 0
    HFLong[:,:monthsPastNotice] .= HFLong[:,:monthsNotice] #months past notice, negative if before notice
    HFLong[:,:noticePeriod] = true #within the notice period
    #HFLong[:,:noticePeriod] = NA #within the notice period
    HFLong[:,:lagVal] = 1 #we will manually dummy the redemption date
    rename!(HFLong, :fDate, :fRedemptionDate)

    HFLong[:,:fDate] = Date("1111-11-11") #placeholder
    HFLong[:,:fDate] .= HFLong[:,:fRedemptionDate] #placeholder and default

    #elliminate all data outside of the notice period
    for l::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE

        #make a view so we are not contantly re-filtering
        HFSub = view(HFLong, HFLong[:,:variable] .== performanceLS.performanceLags[l])
        HFSub[:,:lagVal] = l
        HFSub[:,:monthsPastNotice] =  HFSub[:,:monthsNotice] .- l
        HFSub[~isna(HFSub[:,:monthsPastNotice]),:noticePeriod] =
            ((i::Int)->i>0?true:false).(HFSub[~isna(HFSub[:,:monthsPastNotice]),:monthsPastNotice])

        #lag the dates as needed
        HFSub[:,:fDate] = ((d::Date)-> Dates.lastdayofmonth(d-Dates.Month(l))).(
                HFSub[:,:fRedemptionDate])

    end

    #drops values with too little months notice
    HFLong = HFLong[(~isna(HFLong[:,:value])) & (
        HFLong[:,:monthsNotice].≥MIN_MONTHS_NOTICE_LONG),:]

    pool!(HFLong, :lagVal) #convert this to a factor
    pool!(HFLong, :fund_id)

    for i ∈ 1:length(thresholdLS.negLFlows)
        HFLong[:,thresholdLS.negLFlows[i]] = 0.0 #sets the column as a Float64
        HFLong[:,thresholdLS.negLFlows[i]] = NA
        HFLong[HFLong[:,:negLFlows] .> (collect(SIG_REDEMPTION_THRESHOLDS)[i]),
            thresholdLS.negLFlows[i]] = HFLong[HFLong[:,:negLFlows] .>
            collect(SIG_REDEMPTION_THRESHOLDS)[i],:negLFlows]
    end

    #make a dummy for redemption date (manually so the results are displayed)
    HFLong[:,:isRedemptionDate] = 0.0
    HFLong[HFLong[:,:variable].==:performance,:isRedemptionDate] = 1.0

    pool!(HFLong, :fRedemptionDate)
    pool!(HFLong, :fDate)

    #we now do the first differencing
    sort!(HFLong, cols = [:fund_id, :fRedemptionDate, :fDate])

    #copy these out for cacheing
    fundID::DataVector{Int} = HFLong[:,:fund_id]
    redemptionDate::DataVector{Date} = HFLong[:,:fRedemptionDate]
    performance::DataVector{Float64} =  HFLong[:,:value]
    performance1D::DataVector{Float64} = DataArray(Float64,size(HFLong,1))

    @simd for i∈2:size(HFLong,1)
        if redemptionDate[i]==redemptionDate[i-1] && fundID[i] == fundID[i-1]
            performance1D[i] = performance[i] - performance[i-1]
        end
    end
    HFLong[:,:performance1D] = performance1D

    println("Ending long rows: $(size(HFLong,1))")

    return HFLong

end

#This function aggregates HF data on strategyID
#We will do this in the format of the standard wide file, so that
#we can plug it into the long file creation process
function aggregateOnStrategy(HFData::DataFrame; monthsNotice::T = NA, )::DataFrame where T<:Union{NAtype,Int}

    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)

    ###now calculate the aggregates

    #pre-allocate the aggregate dataframe
    strategies::Vector{Symbol} = unique(HFData[:,:main_strategy])
    dates::Vector{Date} =  unique(HFData[:,:fDate])

    #for each date, have a row for each strategy and a row for the total
    HFAgg::DataFrame = HFData[1:((length(strategies)+1) * length(dates)),:]
    HFAgg[:,:] = NA #blank the rows

    ctr::Int = 0
    for s::Symbol ∈ [strategies; :total], d::Date ∈ dates
        ctr += 1
        HFAgg[ctr,:main_strategy] = s
        HFAgg[ctr,:fDate] = d
        HFAgg[ctr,:monthsNotice] = monthsNotice
    end

    HFAgg[:,:var] = 0.0
    #HFAgg[:,:totNegLFlows] = 0.0
    HFAgg[:,:totFlowsOnNotice] = 0.0
    HFAgg[:,:var] = NA
    #HFAgg[:,:totNegLFlows] = NA
    HFAgg[:,:totFlowsOnNotice] = NA


    #get views for each date
    by(HFData,:fDate) do HFSubPart::SubDataFrame
        HFAggSub::SubDataFrame = view(HFAgg, HFAgg[:,:fDate] .== HFSubPart[1,:fDate]::Date)

        stratTable::Dict = Dict(HFAggSub[i,:main_strategy]=>i for i::Int ∈ 1:length(strategies))

        for s::Symbol ∈ strategies
            if haskey(stratTable,s)
                sInd = stratTable[s]
                HFSub = view(HFSubPart, HFSubPart[:,:main_strategy].==s)

                #now fill in the aggregates
                HFAggSub[sInd,:assets] = sum(dropna(HFSub[:,:assets]))
                HFAggSub[sInd,:flows] = sum(dropna(HFSub[:,:flows]))
                HFAggSub[sInd,:lFlows] = HFAggSub[sInd,:flows] / HFAggSub[sInd,:assets]

                HFAggSub[sInd,:negLFlows] = HFAggSub[sInd,:lFlows]≤0.0?-1.0*HFAggSub[sInd,:lFlows]:NA
              #  HFAggSub[sInd,:totNegLFlows] = sum(dropna(HFSub[:,:negFlows]))
                HFAggSub[sInd,:totFlowsOnNotice] =
                    sum(dropna(HFSub[HFSub[:,:monthsNotice].≥1,:negFlows]))
                #use value weighted or equal weighted method to calculate performance
                HFAggSub[sInd,:performance] =
                    sum(dropna(HFSub[:,:performance].*HFSub[:,:assets])) / HFAggSub[sInd,:assets]
                HFAggSub[sInd,:var] =
                    sum(dropna((HFSub[:,:performance].-HFAggSub[sInd,:performance]).^2 .*
                    HFSub[:,:assets])) / HFAggSub[sInd,:assets]

            end
        end
        sInd = findfirst(HFAggSub[:,:main_strategy], :total)

        #now get the totals
        #HFAggSub[sInd,:totNegLFlows] = sum(dropna(HFSubPart[:,:negFlows]))
        HFAggSub[sInd,:totFlowsOnNotice] =
            sum(dropna(HFSubPart[HFSubPart[:,:monthsNotice].≥1,:negFlows]))
        HFAggSub[sInd,:assets] = sum(dropna(HFSubPart[:,:assets]))
        HFAggSub[sInd,:flows] = sum(dropna(HFSubPart[:,:flows]))
        HFAggSub[sInd,:performance] =
            sum(dropna(HFSubPart[:,:assets] .* HFSubPart[:,:performance])) / HFAggSub[sInd,:assets]
        HFAggSub[sInd,:var] =
            sum(dropna((HFSubPart[:,:performance].-HFAggSub[sInd,:performance]).^2 .*
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

  #=  HFAgg[:,:sigma] = HFAgg[:,:var].^.5
    showall(HFAgg[HFAgg[:,:fDate].==Date("2015-09-30"),
        [:fund_id, :main_strategy, :fDate,:flows, :assets, :performance, :sigma, :var]])
=#
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
    DataFrame([Date, Int, Float64, Float64, Float64, Float64,Float64],
    [:date, :numFunds, :totAUM, :totFlows, :totRedemp,
    :totRollingRedemp, :totRedempNotice], nDates)
  HFByDate[:,:date] .= unique(HFData[:,:date])

  #allocate the lag columns
  for s ∈ [noticeLS.noticeLags; noticeLS.cumNoticeLags]

    HFByDate[:,s] = 0.0
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
    HFByDate[dateIdx,:totAUM]=sum(HFSub[~isna(HFSub[:,:assets]),:assets])
    HFByDate[dateIdx,:totFlows]=sum(HFSub[~isna(HFSub[:,:flows]),:flows])
    HFByDate[dateIdx,:totRedemp]=
      sum(HFSub[broadcast((x)->ifelse(~isna(x)&&x<0,true, false),HFSub[:,:flows]),:flows])
    HFByDate[dateIdx,:totRedempNotice]=
      sum(HFSub[~isna(HFSub[:,:redempNotice]),:redempNotice])
    HFByDate[dateIdx,:totRollingRedemp]=
      sum(HFSub[~isna(HFSub[:,:rollingNotice]),:rollingNotice])
    for s ∈ [noticeLS.noticeLags; noticeLS.cumNoticeLags]
      HFByDate[dateIdx,s]= sum(HFSub[~isna(HFSub[:,s]),s])
    end

  end

  return HFByDate

end

#this executes the linear model tests for a NAV valuation effect
function LM(HFData::DataFrame, titleCaption::String, caption::String;
    flowStr::String = "negLFlows", flowStrDesc::String = "Outflows / AUM x 100")

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
    #run the LM Specs
    for i ∈ 1:length(XLMSpecs)
        println("Running LM Spec $i")

        push!(modelsLM, CTLM(HFData, XLMSpecs[i], fixedEffectsSym = FLMCluster[i],
            :performance, XNames=XLMNames[i], YName = :performance))
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
        colNames = [["($i)" for i∈1:length(modelsLM)]],
        contentRowNames = ["$flowStrDesc"],
        descRowNames = descRowNamesLM,
        descContent = descContentLM,
        decimalDigits = 2,
        columnSepPt = -18,
        scaling = [100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = caption)


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

    for i ∈ 1:length(XPMSpecs)
        println("Running PM Spec $i")
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
      columnSepPt = -27,
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

#this tests for a notice period negative bias to performance
function NNP(HFLong::DataFrame, titleCaption::String, caption::String;
    dependentVar::Symbol = :value)::String
    thresholdLS::RedemptionThresholds = RedemptionThresholds(SIG_REDEMPTION_THRESHOLDS)

    modelsNNP::Vector{CTLM} = Vector{CTLM}()

    XNNPSpecs::Vector{CTExpr} = Vector{CTExpr}()
    XNNPNames::Vector{Vector{Symbol}}  = Vector{Vector{Symbol}}()

    push!(XNNPSpecs, parse("noticePeriod"))
    push!(XNNPNames, [:intercept; :noticePeriod])

    push!(XNNPSpecs, parse("noticePeriod + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod])

    push!(XNNPSpecs, parse("noticePeriod + fRedemptionDate + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod])

    push!(XNNPSpecs, parse("noticePeriod + fDate + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod])

    push!(XNNPSpecs, parse("noticePeriod&negLFlows"))
    push!(XNNPNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNPSpecs, parse("noticePeriod&negLFlows + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNPSpecs, parse("noticePeriod&negLFlows + fRedemptionDate + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod_negLFlows])

    push!(XNNPSpecs, parse("noticePeriod&negLFlows + fDate + main_strategy"))
    push!(XNNPNames, [:intercept; :noticePeriod_negLFlows])



    #push the specifications (use a common name to make a clean table)
    #=for s::Symbol ∈ thresholdLS.negLFlows

        push!(XNNPSpecs, parse("noticePeriod&$s"))
        push!(XNNPNames, [:intercept; :noticePeriod])

        push!(XNNPSpecs, parse("noticePeriod&$s + fDate + main_strategy"))
        push!(XNNPNames, [:intercept; :noticePeriod])

    end=#

    #println(HFLong[1:200, [:noticePeriod, :main_strategy, :monthsNotice, :fDate, :fRedemptionDate]])
    #run the model
    for iX ∈ 1:length(XNNPSpecs)
        push!(modelsNNP, CTLM(HFLong, XNNPSpecs[iX],
            dependentVar, XNames=XNNPNames[iX], YName = :performance))
    end

    ######IO Code for performance bias analysis

    descRowNamesNNP::Vector{String} = ["Date F.E.",  "Redemp. Date F.E.","Strategy F.E.",#="Outflow \\% >", =# "N ('000s)"]

    #need to print the descriptive rows
    descContentNNP::Vector{Vector{String}} =
        [Vector{String}(length(modelsNNP)) for i∈ 1:length(descRowNamesNNP)]

    #Determine the threshold levels (searches for the threshold level in each spec)
    thresholdVec::Vector{Int} = ((i::Int)->maximum(
        ((j::Int)->contains(XNNPSpecs[i],"$(thresholdLS.negLFlows[j])")?
            Int(100.0*collect(SIG_REDEMPTION_THRESHOLDS)[j]):0).(1:length(thresholdLS.negLFlows))
        )).(1:length(XNNPSpecs))

    for i ∈ 1:length(XNNPSpecs)
        descContentNNP[1][i] = "$(contains(XNNPSpecs[i], "fDate")?"X":"")"
        descContentNNP[2][i] = "$(contains(XNNPSpecs[i], "fRedemptionDate")?"X":"")"
        descContentNNP[3][i] = "$(contains(XNNPSpecs[i], "main_strategy")?"X":"")"
        #descContentNNP[3][i] = "$(contains(XNNPSpecs[i], "negLFlows")?"X":"")"
        #descContentNNP[4][i] = "$(thresholdVec[i])\\%"
        descContentNNP[end][i] = "$(round(modelsNNP[i].N/1000.0,1))"
    end

    #rows and variable names for the regression
    rowsNNP::Vector{Symbol} = [:intercept; :noticePeriod; :noticePeriod_negLFlows]

    TableTextNNP::String = texTable(modelsNNP, getHomoskedΣ!, rowsNNP,
        titleCaption = titleCaption,
        colNames = [["($i)" for i∈1:length(modelsNNP)]],
        contentRowNames = ["intercept", "Notice Period", "Notice X Outflow"],
        descRowNames = descRowNamesNNP,
        descContent = descContentNNP,
        decimalDigits = 2,
        columnSepPt = -20,
        scaling = fill(100.0,length(rowsNNP))::Vector{Float64},
        clearMem = USE_AGGRESSIVE_GC,
        caption = caption)

    return TableTextNNP
end

#This tests for a redemption day effect
#note data other than the notice period is dropped
function RD!(HFLong::DataFrame, titleCaption::String, caption::String)::String
    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)
    thresholdLS::RedemptionThresholds =
        RedemptionThresholds(SIG_REDEMPTION_THRESHOLDS)

    for l::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
        HFLong[(HFLong[:,:monthsNotice].≤l) & (HFLong[:,:variable] .==
            performanceLS.performanceLags[l]),:value] = NA
    end

    HFLong = HFLong[~isna(HFLong[:,:value]),:]

    modelsRD::Vector{CTLM} = Vector{CTLM}()

    XRDSpecs::Vector{CTExpr} = Vector{CTExpr}()
    XRDNames::Vector{Vector{Symbol}}  = Vector{Vector{Symbol}}()


    #push the specifications (use a common name to make a clean table)
    for s::Symbol ∈ thresholdLS.negLFlows
        push!(XRDSpecs, parse("isRedemptionDate"))
        push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

        push!(XRDSpecs, parse("isRedemptionDate&$s"))
        push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    end


    #run the model
    for iX ∈ 1:length(XRDSpecs)
        push!(modelsRD, CTLM(HFLong, XRDSpecs[iX],
            :value, XNames=XRDNames[iX], YName = :performance))
    end

    ######IO Code for RD
    #only need N as a descriptive row
    #descRowNamesRD::Vector{String} = ["Date Fixed Effects",
    #    "Strategy Fixed Effects", "Seperate \$\\beta\$ on lags", "N ('000s)" , "Outflow \\% >" ]

    descRowNamesRD::Vector{String} = ["Interacted",  "Outflow \\% >", "N ('000s)"]

    #need to print the descriptive rows
    descContentRD::Vector{Vector{String}} =
        [Vector{String}(length(modelsRD)) for i∈ 1:length(descRowNamesRD)]

    #Determine the threshold levels (searches for the threshold level in each spec)
    thresholdVec::Vector{Int} = ((i::Int)->maximum(
        ((j::Int)->contains(XRDSpecs[i],"$(thresholdLS.negLFlows[j])")?
            Int(100.0*collect(SIG_REDEMPTION_THRESHOLDS)[j]):0).(1:length(thresholdLS.negLFlows))
        )).(1:length(XRDSpecs))

    for i ∈ 1:length(XRDSpecs)
        #=descContentRD[1][i] = "$(contains(XRDSpecs[i], "fDate")?"X":"")"
        descContentRD[2][i] = "$(contains(XRDSpecs[i], "main_strategy")?"X":"")"
        descContentRD[3][i] = "$(contains(XRDSpecs[i], "lagVal")?"X":"")"=#
        descContentRD[1][i] = "$(contains(XRDSpecs[i], "&")?"X":"")"

        descContentRD[2][i] = "$(thresholdVec[i])\\%"
        descContentRD[3][i] = "$(round(modelsRD[i].N/1000.0,1))"
    end

    Gadfly.push_theme(:dark)
    plOutRD = plot(HFLong[((s::Symbol)->s ∈ [:performance, :performanceLag1,
            :performanceLag2]).(HFLong[:,:variable]
             ) & ~isna(HFLong[:,:negLFlowsT1000]),:],
        y=:value,color=:variable, x=:negLFlows,  #=shape=:variable, size=[.1mm],=#
        Stat.binmean(n=25), Geom.line,
        Guide.title("Performance vs Redemption Notifications & Flows"),
        Guide.ylabel("Performance"),Guide.xlabel("Flows"))
      draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_RD.png", 7.5inch, 6inch), plOutRD)


    TableTextRD::String = texTable(modelsRD, getHomoskedΣ!, [:intercept; :negLFlows; :isRedemptionDate],
        titleCaption = titleCaption,
        colNames = [["($i)" for i∈1:length(modelsRD)]],
        contentRowNames = ["intercept", "outflows/AUM x 100", "redemption date x 100"],
        descRowNames = descRowNamesRD,
        descContent = descContentRD,
        decimalDigits = 2,
        columnSepPt = -10,
        scaling = [100.0,100.0,100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = caption)
    return TableTextRD
end

#This tests the notice period against a 0 baseline
#This tests the notice period against a 0 baseline
function NP0!(HFLong::DataFrame, titleCaption::String, caption::String;
    dependentVar=:value, dependentVar1D=:performance1D)::String

    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)
    thresholdLS::RedemptionThresholds =
        RedemptionThresholds(SIG_REDEMPTION_THRESHOLDS)

    for l::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
        HFLong[(HFLong[:,:monthsNotice].≤l) & (HFLong[:,:variable] .==
            performanceLS.performanceLags[l]),:value] = NA
    end

    HFLong = HFLong[~isna(HFLong[:,:value]),:]

    modelsNP0::Vector{CTLM} = Vector{CTLM}()

    XNP0Specs::Vector{CTExpr} = Vector{CTExpr}()
    YNP0Specs::Vector{Symbol} = Vector{CTExpr}()
    XNP0Names::Vector{Vector{Symbol}}  = Vector{Vector{Symbol}}()

    mean(HFLong[:,:value])

    push!(XNP0Specs, parse("negLFlows"))
    push!(XNP0Names, [:intercept; :negLFlows])

    push!(XNP0Specs, parse("negLFlows+main_strategy"))
    push!(XNP0Names, [:intercept; :negLFlows])

    push!(XNP0Specs, parse("negLFlows+main_strategy + fRedemptionDate"))
    push!(XNP0Names, [:intercept; :negLFlows])

    push!(XNP0Specs, parse("negLFlows+fDate + main_strategy"))
    push!(XNP0Names, [:intercept; :negLFlows])


    push!(YNP0Specs, dependentVar)
    push!(YNP0Specs, dependentVar1D)

    #set the descriptive names
    descRowNamesNP0::Vector{String} = ["Date F.E.",  "Redemption Dt F.E.", "Strategy F.E.", "N ('000s)"]

    #need to print the descriptive rows
    descContentNP0::Vector{Vector{String}} =
        [Vector{String}(length(YNP0Specs)*length(XNP0Specs))
            for i∈ 1:length(descRowNamesNP0)]

    #Determine the threshold levels (searches for the threshold level in each spec)
    thresholdVec::Vector{Int} = ((i::Int)->maximum(
        ((j::Int)->contains(XNP0Specs[i],"$(thresholdLS.negLFlows[j])")?
            Int(100.0*collect(SIG_REDEMPTION_THRESHOLDS)[j]):0).(1:length(thresholdLS.negLFlows))
        )).(1:length(XNP0Specs))

    #run the model
    for iY ∈ 1:length(YNP0Specs)
        for iX ∈ 1:length(XNP0Specs)
            push!(modelsNP0, CTLM(HFLong, XNP0Specs[iX],
                YNP0Specs[iY], XNames=XNP0Names[iX], YName = :performance))
            #fill in descriptive content
            modelCnt::Int = length(modelsNP0)
            descContentNP0[1][modelCnt] = "$(contains(XNP0Specs[iX], "fDate")?"X":"")"
            descContentNP0[2][modelCnt] = "$(contains(XNP0Specs[iX], "fRedemptionDate")?"X":"")"
            descContentNP0[3][modelCnt] = "$(contains(XNP0Specs[iX], "main_strategy")?"X":"")"
            descContentNP0[4][modelCnt] = "$(round(modelsNP0[iX].N/1000.0,1))"

        end
    end

    ######IO Code for NP0




    TableTextNP0::String = texTable(modelsNP0, getHomoskedΣ!, [:intercept; :negLFlows],
        titleCaption = titleCaption,
        colNames = [["Performance","\$\\Delta\$ Performance"],["($i)" for i∈1:length(modelsNP0)]],
        contentRowNames = ["intercept", "outflows/AUM x 100"],
        descRowNames = descRowNamesNP0,
        descContent = descContentNP0,
        decimalDigits = 2,
        columnSepPt = -15,
        scaling = [100.0,100.0],
        clearMem = USE_AGGRESSIVE_GC,
        widthColNames = [[length(modelsNP0)÷2, length(modelsNP0)-length(modelsNP0)÷2],ones(Int, length(modelsNP0))],
        caption = caption)


    plOutRD = plot(HFLong, x=:negLFlows,
        Geom.density,
        Guide.title("Outflows Histogram"),
        Guide.xlabel("Outflow"),Guide.ylabel("Count"))
      draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_OutflowLong.png", 7.5inch, 6inch), plOutRD)

    return TableTextNP0
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
        [HFData[~isna(HFData[:,s]),s] for s::Symbol ∈ summaryCols]

    #scale for readability
    summaryDat[summaryDict[:performance]] .*= 100.0
    summaryDat[summaryDict[:negLFlows]] .*= 100.0

    #elliminate zero outflow values from the vector for summary stat purposes
    filter!((x::Float64)->x>0, summaryDat[summaryDict[:negLFlows]])
    summaryRowNames::Vector{String} =
        ["mean", "median", "\$\\sigma\$","Number of Funds",
        "Number of Points","Months (2000-2017)"]

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
        """Data sourced from WRDS HFS. Includes alive and dead funds 2000-2017. Initial
        Data set consisted of 1.6mn rows, reduced by preprocessing to ~385k. Pre-processing
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
        columnSepPt=5)

    plOutRD = plot(HFData, x=:assets,
        Geom.density,
        Guide.title("Assets Historgram"),
        Guide.xlabel("AUM"),Guide.ylabel("Count"))
      draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_assets.png", 7.5inch, 6inch), plOutRD)

    return TableSummary
end

#this creates a performance benchmark
function joinLongAgg(HFLong::DataFrame, HFLongAgg::DataFrame)::DataFrame
    sort!(HFLongAgg, cols=[:fDate, :main_strategy])

    dateSymTable =  Dict((HFLongAgg[i,:fDate], HFLongAgg[i,:main_strategy])=>
        (HFLongAgg[i,:value],HFLongAgg[i,:performance1D]) for i::Int∈1:size(HFLongAgg,1))

    HFLong[:,:performanceStrat] = 0.0
    HFLong[:,:performanceStrat1D] = 0.0
    HFLong[:,:performanceRel] = 0.0
    HFLong[:,:performanceRel1D] = 0.0

    HFLong[:,:performanceStrat] = NA
    HFLong[:,:performanceStrat1D] = NA
    HFLong[:,:performanceRel] = NA
    HFLong[:,:performanceRel1D] = NA


    @simd for i ∈ 1:size(HFLong,1)
        if haskey(dateSymTable, (HFLong[i,:fDate], HFLong[i,:main_strategy]))
            #gets the strategy performance and first difference of the performance
            HFLong[i,:performanceStrat], HFLong[i,:performanceStrat1D] =
                dateSymTable[(HFLong[i,:fDate], HFLong[i,:main_strategy])]

            #get the benchmark relative performance and first difference
            HFLong[i,:performanceRel] = HFLong[i,:value] -
                HFLong[i,:performanceStrat]
            HFLong[i,:performanceRel1D] = HFLong[i,:performance1D] -
                HFLong[i,:performanceStrat1D]
        else
            error("Join Failure on row $i. Row Info:\n $(HFLong[i,[:fDate,
                :fRedemptionDate,:variable,:value,:monthsNotice,:main_strategy]])")
        end
    end

    return HFLong

end

function joinDataAgg(HFData::DataFrame, HFDataAgg::DataFrame)::DataFrame
    sort!(HFDataAgg, cols=[:fDate, :main_strategy])

    dateSymTable =  Dict((HFDataAgg[i,:fDate], HFDataAgg[i,:main_strategy]) =>
        HFDataAgg[i,:performance] for i::Int∈1:size(HFDataAgg,1))

    HFData[:,:performanceStrat] = 0.0
    HFData[:,:performanceRel] = 0.0

    HFData[:,:performanceStrat] = NA
    HFData[:,:performanceRel] = NA


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
function getCanceledEstimate!(HFData::DataFrame;
  focalStr::String="isNoticePeriod", dependentVar=:performance)::DataFrame


  HFDataModel = HFData[~isna(HFData[:,:isNoticePeriod]),:]
  dates::Vector{Date} = unique(HFData[:,:fDate])
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

  #run the model specs
  M::Int = length(XSpecs)
  for i ∈ 1:M
      println("Running Spec $i for redemption calc")
      push!(models, CTLM(HFDataModel, XSpecs[i], fixedEffectsSym = FCluster[i],
          dependentVar, XNames=XNames[i], YName = :performance))

      #check to make sure we didn't drop a row
      if (models[end].N) != size(HFDataModel[:,:fDate],1)
        error("Size Mismatch spec $i- model rows: $(models[end].N), data rows: $(size(HFData[:,:fDate],1))")
      end
      HFDataModel[:,Symbol("Y", i)] = models[end].Y
  end


  HFFlowDat::DataFrame = DataFrame([Date, Float64, Float64, Float64, Float64],
    [:fDate, :totFlowsOnNotice, :totNegFlows, :assets, :assetsModel], N)

  HFFlowDat[:,:fDate] = dates

  #write the variable column names into the dataframe
  for i::Int ∈ 1:M
    HFFlowDat[:,Symbol("μ$i")] = 0.0
    HFFlowDat[:,Symbol("σ2$i")] = 0.0
    HFFlowDat[:,Symbol("threshold$i")] = 0.0
    HFFlowDat[:,Symbol("totSubmitted$i")] = 0.0
    HFFlowDat[:,Symbol("canceled$i")] = 0.0

    HFFlowDat[:,Symbol("μ$i")] = NA
    HFFlowDat[:,Symbol("σ2$i")] = NA
    HFFlowDat[:,Symbol("threshold$i")] = NA
    HFFlowDat[:,Symbol("totSubmitted$i")] = NA
    HFFlowDat[:,Symbol("canceled$i")] = NA
  end

  #now aggregate by date from the larger set
  by(HFData, :fDate) do HFSub::SubDataFrame
    r = findfirst(HFFlowDat[:,:fDate], HFSub[1,:fDate])

    #aggregate the flows
    HFFlowDat[r,:totNegFlows] = sum(dropna(HFSub[:,:negFlows]))
    HFFlowDat[r,:totFlowsOnNotice] =
        sum(dropna(HFSub[HFSub[:,:monthsNotice].≥1,:negFlows]))
    HFFlowDat[r,:assets] = sum(dropna(HFSub[:,:assets]))

    #this is the case where we are not demeaning the performance and variance
    #beyond the specific date
    for i::Int ∈ 1:M
      if FCluster[i] == nothing
        HFFlowDat[r,Symbol("μ$i")] =
            sum(dropna(HFSub[:,:performance].*HFSub[:,:assets])) / HFFlowDat[r,:assets]
        HFFlowDat[r,Symbol("σ2$i")] =
          sum(dropna((HFSub[:,:performance].-HFFlowDat[r,Symbol("μ$i")]).^2 .*
            HFSub[:,:assets])) / HFFlowDat[r,:assets]
      end
    end
  end

  #where we have the demeaned info
  by(HFDataModel, :fDate) do HFSub::SubDataFrame
    r = findfirst(HFFlowDat[:,:fDate], HFSub[1,:fDate])
    HFFlowDat[r,:assetsModel] = sum(dropna(HFSub[:,:assets]))

    #now aggregate the performance

    for i::Int ∈ 1:M
      if FCluster[i] != nothing
        HFFlowDat[r,Symbol("μ$i")] =
            sum(dropna(HFSub[:,Symbol("Y$i")].*HFSub[:,:assets])) / HFFlowDat[r,:assetsModel]
        HFFlowDat[r,Symbol("σ2$i")] =
          sum(dropna((HFSub[:,Symbol("Y$i")].-HFFlowDat[r,Symbol("μ$i")]).^2 .*
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
  μChkModel = sum(dropna(HFDataModel[HFDataModel[:,:fDate].==d,Symbol("Y3")].*
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
      print(HFFlowDat[r,:fDate],", ")

        #form the objective function
      μ::Float64 = HFFlowDat[r,Symbol("μ$m")]
      σ::Float64 = HFFlowDat[r,Symbol("σ2$m")]^0.5
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
      JuMP.register(mod, :CTNorm, 1, CTCumNorm, autodiff=true)
      JuMP.register(mod, :CTCumNorm, 1, CTCumNorm, autodiff=true)
      @variable(mod, t)

      JuMP.register(mod, :obj, 1, obj, autodiff=true)
      @NLobjective(mod, Min, obj(t))
      output = solve(mod)

      HFFlowDat[r,Symbol("threshold$m")] = getvalue(t)
      HFFlowDat[r,Symbol("totSubmitted$m")]  = HFFlowDat[r,:totFlowsOnNotice]/
          CTCumNorm((HFFlowDat[r,Symbol("threshold$m")] - μ)/σ)
      HFFlowDat[r,Symbol("canceled$m")] = HFFlowDat[r,Symbol("totSubmitted$m")] -
        HFFlowDat[r,:totFlowsOnNotice]
    end

  end

  #aggregate to the annual level
  HFFlowDat[:,:year] = (Dates.year).(HFFlowDat[:,:fDate])
  years::vector = unique(HFFlowDat[:,:year])
  HFFlowDatAnn = aggregate(HFFlowDat[:,:year, :totFlowsOnNotice,
    :threshold1, :totSubmitted1, :canceled1,
    :threshold2, :totSubmitted2, :canceled2,
    :threshold3, :totSubmitted3, :canceled3], :year, [sum,mean])

  for i∈1:M
    HFFlowDatAnn[:,Symbol("threshold$i")] = HFFlowDatAnn[:,Symbol("threshold$(i)_mean")]
    HFFlowDatAnn[:,Symbol("threshold$i")] = HFFlowDatAnn[:,Symbol("threshold$(i)_sum")]
    HFFlowDatAnn[:,Symbol("threshold$i")] = HFFlowDatAnn[:,Symbol("threshold$(i)_sum")]
  end

  HFFlowDatAnn[:,:totFlowsOnNotice] = HFFlowDatAnn[:,:totFlowsOnNotice_sum]

  HFlowDatAnn = unique(HFFlowDatAnn[:,:year, [:totFlowsOnNotice,
    :threshold1, :totSubmitted1, :canceled1,
    :threshold2, :totSubmitted2, :canceled2,
    :threshold3, :totSubmitted3, :canceled3]])

  return HFFlowDat

end

#=
function getCanceledEstimateOld(HFDataAgg::DataFrame, β::Float64)::DataFrame

    N::Int = length(unique(HFDataAgg[:,:fDate]))

    sort!(HFDataAgg, cols=[:fDate])

    #allocate space for the result
    HFDataAgg[:,:threshold] = 0.0
    HFDataAgg[:,:totSubmitted] = 0.0
    HFDataAgg[:,:canceled] = 0.0

    HFDataAgg[:,:threshold] = NA
    HFDataAgg[:,:totSubmitted] = NA
    HFDataAgg[:,:canceled] = NA



    const sqrt2Inv::Float64 = 1/(2.0^.5)
    const sqrt2PiInv::Float64 = 1.0/(2.0*π)^0.5

    #need to define special functions here for the solver
    CTCumNorm(z) = 0.5+0.5*erf(z*sqrt2Inv)
    CTNorm(z) = sqrt2PiInv*exp(-z^2.0/2.0)



    curYear::Int = 1111
    for i::Int ∈ 1:size(HFDataAgg,1)
      if HFDataAgg[i,:main_strategy] == :total #not quite enough data for strategy level
        if curYear != Dates.year(HFDataAgg[i,:fDate])
          println("solving months in: $(Dates.year(HFDataAgg[i,:fDate]))")
        end
        curYear = Dates.year(HFDataAgg[i,:fDate])

        #form the bojective function
        μ::Float64 = HFDataAgg[i,:performance]
        σ::Float64 = HFDataAgg[i,:sigma]
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
        JuMP.register(mod, :CTNorm, 1, CTCumNorm, autodiff=true)
        JuMP.register(mod, :CTCumNorm, 1, CTCumNorm, autodiff=true)
        @variable(mod, t)

        JuMP.register(mod, :obj, 1, obj, autodiff=true)
        @NLobjective(mod, Min, obj(t))
        output = solve(mod)

        HFDataAgg[i,:threshold] = getvalue(t)
        HFDataAgg[i,:totSubmitted] = HFDataAgg[i,:totFlowsOnNotice]/
            CTCumNorm((HFDataAgg[i,:threshold]-HFDataAgg[i,:performance])/HFDataAgg[i,:sigma])
        HFDataAgg[i,:canceled] = HFDataAgg[i,:totSubmitted]-HFDataAgg[i,:totFlowsOnNotice]
      end

    end

    return HFDataAgg

end=#

#this is the main method for analysis
function analyzeAggregates(;FlowsLaggingProcedure::Char = 'U',
    verbose::Bool = true, runLM::Bool = true,
    runPM::Bool = true, runRD::Bool = true,
    runSummary::Bool=true, runNNP::Bool = true,
    runNNP1D::Bool = true, runNP0::Bool = true,
    process=true, save=false )



    #wrapping this up to save time if the file is long
    if process
        #open the pre-processed data file
        igStream::GZipStream =
            gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls.gz")
        HFData::DataFrame = deserialize(igStream)
        close(igStream)

        #get the long files
        HFLong::DataFrame = processLongHFData(HFData) #get long file
        HFDataAgg::DataFrame = aggregateOnStrategy(HFData, monthsNotice=MIN_MONTHS_NOTICE_LONG) #get wide strategy aggregate
        HFLongAgg::DataFrame = processLongHFData(HFDataAgg, negFlowsOnly = false)

        HFData = joinDataAgg(HFData,HFDataAgg)
        HFLong = joinLongAgg(HFLong, HFLongAgg)
        HFLongAgg = HFLongAgg[HFLongAgg[:,:negLFlows].>0.0,:]

        titleCaption::String = string()
        caption::String = string()

        sort!(HFData,cols = (:date, :fund_id))
        HFByDate::DataFrame = aggregateOnDate(HFData, LAGS_NOTICE_TO_CAPTURE) #get date aggregate

        if save
            ogStream::GZipStream = gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_Proc.jls.gz","w")
            serialize(ogStream, [HFData, HFLong, HFDataAgg, HFLongAgg, HFByDate])
            close(ogStream)
        end
    else # if we arn't going to process, load in the processed files
        igStream =
            gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_Proc.jls.gz")
        HFData, HFLong, HFDataAgg, HFLongAgg, HFByDate = deserialize(igStream)
        close(igStream)
    end

    ##################OLS Spec

    titleCaption = "Performance vs Redemption Size"
    caption = """ This table illustrates the relationship between redemptions atand
        performance. The key takeway is performance on the redemption date is,
        counterintuitivelly, positively correlated with the size of the redemption.
        Date fixed effects are by month. Performance lags includes the prior
        24 months of performance as a covariate. Data is from alive and dead hedge
        funds 2000-2017. Where fixed effects are present, the standard errors are
        clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6).
         Outflows were expressed as a percentage of the
        AUM. A positive value for the coefficient implies a positive relationship
        between performance and redemptions."""

    TableTextLM::String = runLM?LM(HFData, titleCaption, caption,
        flowStr = "negLFlows", flowStrDesc = "Outflows / AUM x 100"):""

    if USE_AGGRESSIVE_GC
        gc()
    end


    titleCaption = "Performance vs Cumulative Rolling Redemption Notice Size"
    caption = """This table reflects the relationship between pending notifications
            of redemptions and performance. Here, pending notifications are negatively
            correlated with performance. Date fixed effects are by month. Performance lags includes the prior
            24 months of performance as a covariate. Data is from alive and dead hedge
            funds 2000-2017. Where fixed effects are present, the standard errors are
            clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6).
            Outflows were expressed as a percentage of the
            AUM. A positive value for the coefficient implies a positive relationship
            between performance and redemptions."""

    TableTextLMNotice::String = runLM?LM(HFData, titleCaption, caption,
        flowStr = "negRollingNoticeL", flowStrDesc = "Rolling Notice / AUM x 100"):""
    if USE_AGGRESSIVE_GC
        gc()
    end

    #########IV Note: See the postIVO backup for the IV code
    ##############Panel Approach

    coefs::Vector{Float64}=Vector{Float64}()

    titleCaption = "Effect of Notice Period on Performance"
    caption = "This specification uses a dummy for the notice period and a dummy
    for the redemption period. In each specification, funds underperform in the
    notice period. Except with fund-level fixed effects,
    hedge funds realize  further losses on the redemption day. Data is from alive and dead hedge
    funds 2000-2017. Interpret a value of 1.0 for the dummy as corresponding
    to a difference in performance of 1.0\\%. Where fixed effects are present, the standard errors are
    clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
    TableTextPM = PM(HFData,titleCaption,caption, cluster=true, coefs=coefs)


    titleCaption = "Effect of Notice Period on Performance (BM-Relative)"
    caption = "This is the same specifciation as the prior table, only
    performance is net of a strategy-specific value-weighted benchmark. Results
    are generally robust. In each specification, funds underperform in the
    notice period. Except with fund-level fixed effects,
    hedge funds realize  further losses on the redemption day. Data is from alive and dead hedge
    funds 2000-2017. Interpret a value of 1.0 for the dummy as corresponding
    to a difference in performance of 1.0\\%. Where fixed effects are present, the standard errors are
    clustered either by fund (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
    TableTextPMRel = runPM?PM(HFData,titleCaption,caption, dependentVar=:performanceRel, cluster=true):""

    titleCaption = "Effect of Notice Period on Performance (No-overlap)"
    caption = "Another robustness check on underperformance in the notice period.
    Here, overlapping redemptions are dropped from the data. Again, in each
    specification, funds underperform in the notice period. Except with fund-level fixed effects,
    hedge funds realize  further losses on the redemption day. Data is
    from alive and dead hedge funds 2000-2017. Interpret a value of 1.0 for
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
    from alive and dead hedge funds 2000-2017. Interpret a value of 1.0 for
    the dummy as corresponding to a difference in performance of 1.0\\%. Where
    fixed effects are present, the standard errors are clustered either by fund
    (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
    TableTextPMSingleRel = runPM?PM(HFData,titleCaption,caption,
        dependentVar=:performanceRel,focalStr="isNoticePeriodSingle", cluster=true):""

    titleCaption = "Effect of Notice Period on Performance (1mo Notice Only)"
    caption = "This specification drops overlapping redemptions and limits the
    data set to one months' notice. The purpose is to estimate the amount of
    canceled redemptions. Except with fund-level fixed effects,
    hedge funds realize
    further losses on the redemption day. Data is
    from alive and dead hedge funds 2000-2017. Interpret a value of 1.0 for
    the dummy as corresponding to a difference in performance of 1.0\\%. Where
    fixed effects are present, the standard errors are clustered either by fund
    (4 \\& 8), strategy (3 \\& 7), or date (2 \\& 6)."""
    TableTextPMSingle1Mo = runPM?PM(HFData,titleCaption,caption,
        dependentVar=:performance,focalStr="isNoticePeriod1Mo", cluster=true):""

    begin
      #######Now the NNP Approach (comparison of notice period and non-notice period)

      titleCaption = "Performance Near Withdrawals: Impact of Notice Period"
      caption = """Compares performance prior to notice,
          with performance after notice, in both a  standard and First-Differenced
          configuration. In some specifications, the notice
          period dummy is interacted with the amount of the outflow. Threshold
          refers to the minimum outflow included, thus some versions include
          only significant outflows. Data is from alive and dead hedge funds
          2000-2017. Interpret a value of 1.0 for the non-interacted
          specification as the increase in performance (or the
          change in performance in the 1st differenced version).
          For the interacted specification, a 1.0 implies that outflows of
          10\\% corresponding to a 0.1\\% drop in monthly performance."""

      TableTextNNP::String = runNNP?NNP(HFLong, titleCaption, caption):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "First-Differenced Performance vs Outflow Notice Period"
      caption = """Compares performance prior to notice,
          with performance after notice, in both a  standard and First-Differenced
          configuration. In some specifications, the notice
          period dummy is interacted with the amount of the outflow. Threshold
          refers to the minimum outflow included, thus some versions include
          only significant outflows. Data is from alive and dead hedge funds
          2000-2017. Interpret a value of 1.0 for the non-interacted
          specification as the increase in performance (or the
          change in performance in the 1st differenced version).
          For the interacted specification, a 1.0 implies that outflows of
          10\\% corresponding to a 0.1\\% drop in monthly performance."""
      #######Same as prior NNP but with first differences
      TableTextNNP1D::String = runNNP1D?NNP(HFLong,titleCaption, caption,
          dependentVar = :performance1D):""

      if USE_AGGRESSIVE_GC
          gc()
      end


      titleCaption = "HF Strategy Outflow Performance Relative to Notice Period Performance"
      caption = """This is the same specification as the fund level spec,
          only we assume a one month notice and evaluate flows at the strategy level.
          In some specifications, the notice period dummy is interacted with the
          amount of the outflow. Threshold
          refers to the minimum outflow included, thus some versions include
          only significant outflows. Data is from alive and dead hedge funds
          2000-2017. Interpret a value of 1.0 for the non-interacted
          specification as the increase in performance (or the
          change in performance in the 1st differenced version).
          For the interacted specification, a 1.0 implies that outflows of
          10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNNPAgg::String = runNNP?NNP(HFLongAgg, titleCaption, caption):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "Net-of-BM Outflow Performance Relative to Notice Period Performance"
      caption = """This time the performance is net of value-weighted strategy performance.
          In some specifications, the notice period dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included, thus some versions include
          only significant outflows. Data is from alive and dead hedge funds
          2000-2017. Interpret a value of 1.0 for the non-interacted
          specification as the increase in performance (or the
          change in performance in the 1st differenced version).
          For the interacted specification, a 1.0 implies that outflows of
          10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNNPRel::String = runNNP?NNP(HFLong, titleCaption, caption, dependentVar=:performanceRel):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "Net-of-BM Outflow Performance Relative to Notice Period Performance"
      caption = """Finally, this specification uses the first differenced
          value-weighted strategy performance.
          In some specifications, the notice period dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included, thus some versions include
          only significant outflows. Data is from alive and dead hedge funds
          2000-2017. Interpret a value of 1.0 for the non-interacted
          specification as the increase in performance (or the
          change in performance in the 1st differenced version).
          For the interacted specification, a 1.0 implies that outflows of
          10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNNPRel1D::String = runNNP?NNP(HFLong, titleCaption, caption, dependentVar=:performanceRel1D):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      #RD Analysis
      titleCaption = "Performance Flows RD Specifications"
      caption = """RD specification uses the redemption date as the discontinuity,
          relative to performance during the redemption's notificaiton period.
          In some specifications, the redemption dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included,
          thus some versions include only significant outflows.
          Data is from alive and dead hedge funds 2000-2017. Interpret a value
          of 1.0 for the non-interacted specification as the discontinuity increasing
          performance by 1\\%. For the interacted specification, a 1.0 implies that
          outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextRD::String = runRD?RD!(HFLong, titleCaption, caption):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "Fund-Level Notice Period Performance"
      caption = """This specification examines whether performance is different
          than zero. In some specifications, the redemption dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included,
          thus some versions include only significant outflows.
          Data is from alive and dead hedge funds 2000-2017. Interpret a value
          of 1.0 for the non-interacted specification as the discontinuity increasing
          performance by 1\\%. For the interacted specification, a 1.0 implies that
          outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNP0::String = runNP0?NP0!(HFLong, titleCaption, caption):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "Fund-Level BM-Relative Notice Period Performance"
      caption = """This specification examines whether performance is different
          than zero. In some specifications, the redemption dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included,
          thus some versions include only significant outflows.
          Data is from alive and dead hedge funds 2000-2017. Interpret a value
          of 1.0 for the non-interacted specification as the discontinuity increasing
          performance by 1\\%. For the interacted specification, a 1.0 implies that
          outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNP0Rel::String = runNP0?NP0!(HFLong, titleCaption, caption,
          dependentVar=:performanceRel, dependentVar1D=:performanceRel1D):""
      if USE_AGGRESSIVE_GC
          gc()
      end

      titleCaption = "Strategy-Level Notice Period Performance"
      caption = """This specification examines whether performance is different
          than zero after flows have been aggregated to the strategy level.
          In some specifications, the redemption dummy is interacted with the
          amount of the outflow. Threshold refers to the minimum outflow included,
          thus some versions include only significant outflows.
          Data is from alive and dead hedge funds 2000-2017. Interpret a value
          of 1.0 for the non-interacted specification as the discontinuity increasing
          performance by 1\\%. For the interacted specification, a 1.0 implies that
          outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance."""
      TableTextNP0Agg::String = runNP0?NP0!(HFLongAgg, titleCaption, caption):""
      if USE_AGGRESSIVE_GC
          gc()
      end
    end
    ###################Calculate canceled assets
    HFFlowDat::DataFrame, HFFlowDat::Ann = getCanceledEstimate!(HFData)

    showall(HFFlowDat[: , [:fDate, :canceled1,:threshold1,:totSubmitted1]])


    ####################Do the Summary Tables and write all table to a file

    TableSummary::String = runSummary?Summary(HFData::DataFrame):""

    writeTables2File([TableSummary, TableTextLM, TableTextLMNotice,
        TableTextPM, TableTextPMRel, TableTextPMSingle, TableTextPMSingleRel,
        TableTextPMSingle1Mo, TableTextRD, TableTextNNP, TableTextNNPAgg,
        TableTextNNPRel, TableTextNNPRel1D, TableTextNNP1D,
        TableTextNP0, TableTextNP0Rel, TableTextNP0Agg],
        HEADER_NAME, FOOTER_NAME, path=RESULTS_PATH,
        outName = "RegTables_$FlowsLaggingProcedure.tex") #cosmetic marker block for summary tables


end

function verify2R(;FlowsLaggingProcedure::Char = 'U')

    #read the seria
    igStream::GZipStream =
        gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls.gz")
    HFData::DataFrame = deserialize(igStream)
    close(igStream)

    #write the csv
    writetable("$R_TEST_PATH\\RHFData.csv",HFData)

end




@time begin
FLPChar = 'C'
#preProcessHFData(FlowsLaggingProcedure=FLPChar, newBinary=true)
analyzeAggregates(FlowsLaggingProcedure=FLPChar, runSummary=true, runLM = false,
    runPM = true, runRD = false, runNNP=false, runNNP1D=false,
    runNP0=false, process=true, save=true)
#verify2R(FlowsLaggingProcedure=FLPChar)
end



#=@time for c ∈ ['C', 'R', 'L', 'U']
  PreProcessHFData(FlowsLaggingProcedure=c, newBinary=true)
  AnalyzeAggregates(FlowsLaggingProcedure=c)
end=#




end
