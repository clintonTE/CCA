
module HFMod

if length(findin("C:\\Users\\Clinton\\Dropbox\\Projects\\SummerProject17",LOAD_PATH)) == 0
  push!(LOAD_PATH,"C:\\Users\\Clinton\\Dropbox\\Projects\\SummerProject17")
end

#=TODO:
1) IV=#
#2) Run specifications

using DataFrames, Distributions, StatsBase, GZip, JLD, Gadfly, CTMod
#importall CT

#Set this to an optimal value- I find logical cores *(3/4) works well
BLAS.set_num_threads(Int(round((Sys.CPU_CORES*3÷4))))


#if !isdefined(:constDef)
  const constDef = true
  #other constants here

  const DATA_PATH = pwd() * "\\data"
  const DATA_FILE_NAME = "HFPerformance_2000-2017"
  const DATA_FILE_NAME_DEAD = DATA_FILE_NAME * "_dead"
  const DATE_FORMAT_STR = "yyyymmdd"
  const NOTICE_PERIOD_ADJ = 7 #days to add to the notice (since data are monthly)
  const MIN_PERIODS_FOR_INCLUSION = 12 #number of periods required for inclusion
  const T_TEST_THRESHOLD = 2.0 #this is a parameter for identificaiton if the flows are lagged
  const PERIODS_SAME_AUM_TO_DELETE = 4 #assume 4 consecutive periods of identical aum implies stale data
  const LAGS_NOTICE_TO_CAPTURE = 3
  const LAGS_PERFORMANCE_TO_CAPTURE = 11
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
  const MIN_MONTHS_NOTICE_RD = 1 #min months notice to be included in RD analysis
  const SIG_REDEMPTION_THRESHOLDS = 0.0:0.05:0.1
  const PERFORMANCE_BOUND = .25


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
  performanceLags::Vector{Symbol}
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

contains(e::T where T<:CTExpr, s::String)::Bool  = Base.contains("$e",s)
###########Main Methods
#U for unlagged, L for lagged, C for corrected, R for corrected but with
#a lagged default

function PreProcessHFData(;FlowsLaggingProcedure::Char = 'U', newBinary::Bool = false)

  nPer::Int = 0 #numer of periods
  dtFormat::DateFormat = DateFormat(DATE_FORMAT_STR) ## a date format object

  #this allows us to capture a variable amount of lags
  noticeLS::NoticeLagSymbols = NoticeLagSymbols(LAGS_NOTICE_TO_CAPTURE)
  performanceLS::PerformanceLagSymbols = PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)

  #this makes sure we have a binary of the data (Improves load times 3-4x)
  if !isfile("$DATA_PATH\\$DATA_FILE_NAME.jls") || newBinary
    #extract from the zip file
    gHandle::GZipStream = GZip.open("$DATA_PATH\\$DATA_FILE_NAME.gz")
    write("$DATA_PATH\\$DATA_FILE_NAME.csv", readlines(gHandle,chomp=false))
    close(gHandle)

    gHandle = GZip.open("$DATA_PATH\\$DATA_FILE_NAME_DEAD.gz")
    write("$DATA_PATH\\$DATA_FILE_NAME_DEAD.csv", readlines(gHandle,chomp=false))
    close(gHandle)

    #write the binary
    oStream::IOStream = open("$DATA_PATH\\$DATA_FILE_NAME.jls","w+")
    serialize(oStream,
        [readtable("$DATA_PATH\\$DATA_FILE_NAME.csv"); readtable("$DATA_PATH\\$DATA_FILE_NAME_DEAD.csv")])
    close(oStream)
  end

  #load the binary
  iStream::IOStream = open("$DATA_PATH\\$DATA_FILE_NAME.jls")
  HFData::DataFrame = deserialize(iStream)
  close(iStream)

  println("Initial rows: ", size(HFData,1), ", Initial columns: ", size(HFData,2))
  println("Begining pre-processing...")

  #drop some columns
  HFData = HFData[:, [:fund_id, :main_strategy,
    :sub_strategy, :fund_assets, :firm_assets,
    :advance_notice, :fund_assets_as_of,
    :date, :performance, :nav, :assets]]

  #Drop some null valued fields:
  HFData = dropNullsFromDF(HFData, [:date, :performance, :assets,
  :advance_notice])

  HFData = HFData[HFData[:,:assets] .≥ MIN_ASSETS, :]

  ##This can be used to constrain to only funds with an incentive fees
  #However, the incentive fee DQ is low, hence it is not included
  #=parse the incentive fees and mangement fees
  DFStringtoMeanNum!(HFData, :incentive_fee)
  DFStringtoMeanNum!(HFData, :management_fee)

  #convert from percentage to decimal
  HFData[:,:incentive_fee] ./= 100.0
  HFData[:,:management_fee] ./= 100.0

  HFData = HFData[~isna(HFData[:,:incentive_fee]) &  (HFData[:,:incentive_fee] .> 0.0),:]
  =#

  #drop unlikely outliers (data quality)
  HFData[:,:performance] ./= 100.0
  HFData = HFData[(HFData[:performance].< PERFORMANCE_BOUND) & (
    HFData[:performance].> (-1.0*PERFORMANCE_BOUND)),:]

  #broadcast an anonymous functino to do the type conversions
#  HFData[:,:inception] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:inception])
  #HFData[:,:date_added_to_db] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:date_added_to_db])
  HFData[:,:date] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:date])
 # HFData[:,:fund_assets_as_of] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:fund_assets_as_of])

  #sort the data by fund id and date
  sort!(HFData,cols = (:fund_id,:date))
  #delStaleAUMDat!(HFData)

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
  HFData[:,:months_notice] = broadcast((x::Int)->
    ((x+NOTICE_PERIOD_ADJ)÷30)::Int,HFData[:,:advance_notice])

  ###The following sections create nessecary columns with placeholders
  #calculate monthly flows and construct the data set
  HFData[:,:flows] = 0.0 #flows
  HFData[:,:redemp_notice] = 0.0 #notifications about flows
  HFData[:,:rolling_notice] = 0.0 #notifications about flows
  HFData[:,:lFlows]  = 0.0 #will hold log flows
  HFData[:,:redemp_notice_dt] = Date(1,1,1)

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
      perfDelta::Vector{Float64} = HFSub[2:end,:performance]
        .-HFSub[1:(end-1),:performance]

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
        if !isna(HFSub[i,:flows]) && i::Int > HFSub[i,:months_notice]

            #note notifications of flows does not include the current month's flows
            if HFSub[i,:months_notice] ≥ MIN_MONTHS_NOTICE_IV
                HFSub[(i-HFSub[i,:months_notice]):(i-1),:rolling_notice] += HFSub[i,:flows]
            end

            if HFSub[i,:flows] < 0.0

                #record the notification on the notification date
                HFSub[i-HFSub[i,:months_notice],:redemp_notice] = HFSub[i,:flows]
                HFSub[i-HFSub[i,:months_notice],:redemp_notice_dt] = HFSub[i,:date]

                #record the redemption notificaion lags and performance lags
                if HFSub[i,:months_notice] > 0
                    for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, HFSub[i,:months_notice])
                        HFSub[i-l, noticeLS.noticeLags[l]] = HFSub[i,:flows]
                        HFSub[i-l, noticeLS.noticeLLags[l]] = HFSub[i,:lFlows]
                    end
                end
            end
        end
    end

    #this si where we capture the cumulative notifications of redeumptions from prior months
    for i::Int ∈ 2:nPer
        for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, i-1)
            #we only record the lag if the total notices are net negative
            if (!isna(HFSub[i-l,:rolling_notice])) && (HFSub[i-l,:rolling_notice] < 0.0)
                HFSub[i,noticeLS.cumNoticeLags[l]] = HFSub[i-l,:rolling_notice]
                HFSub[i,noticeLS.cumNoticeLLags[l]] = HFSub[i-l,:rolling_notice] / HFSub[i-l,:assets]
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
        if i - 1 - HFSub[i,:months_notice] ≥ 1 #(make sure its positive number)
          for l::Int ∈ 1:min(i-1- HFSub[i,:months_notice],LAGS_PERFORMANCE_TO_CAPTURE)
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
      HFData[HFData[:,:months_notice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = NA
      HFData[HFData[:,:months_notice].<MIN_MONTHS_NOTICE_IV,noticeLS.cumNoticeLags[l]] = NA
  end

  #get a seperate column with net negative flows
  HFData[:,:negLFlows] = HFData[:,:lFlows]
  HFData[(~isna(HFData[:,:lFlows])) & (HFData[:,:lFlows].≥0.0),:negLFlows] .= NA



  #convert string columns to factors
  HFTypes::Vector{Type} = eltypes(HFData) #all types
  HFNames::Vector{String} = names(HFData) #all names
  for i ∈ 1:size(HFData,2) #for all columns
    if HFTypes[i] == String && !(typeof(HFData[:,i]) <: PooledDataVector)
      pool!(HFData,Symbol(HFNames[i])) #convert to factor
    end
  end

  #create a factor column for dates
  HFData[:,:fDate] = HFData[:,:date]
  pool!(HFData,:fDate)


  println("Max Months Notice: $(unique(HFData[:,:months_notice]))")


  println("Processing complete.")
  println("Ending rows: ", size(HFData,1), ", Ending columns: ", size(HFData,2))

  #save the binary
  oStream = open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls","w+")
  serialize(oStream, HFData)
  close(oStream)

end

getSigRedemptionThresholds(rng::Range = SIG_REDEMPTION_THRESHOLDS)::Vector{Symbol} =
    [Symbol("negLFlowsT$(Int(t*10000))") for t ∈ collect(rng)]


function processLongHFData(;FlowsLaggingProcedure::Char = 'U')
    iStream::IOStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls")
    HFData::DataFrame = deserialize(iStream)
    close(iStream)

    println("Beginning lengthening of data. \nStarting Wide rows: $(size(HFData,1))")

    #get any variable arrays of symbols
    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)
    noticeLS::NoticeLagSymbols = NoticeLagSymbols(LAGS_NOTICE_TO_CAPTURE)
    thresholdSyms::Vector{Symbol} = getSigRedemptionThresholds()
    println(thresholdSyms)

    #store this for convenience
    arraySyms::Vector{Symbol} = [:performance; performanceLS.performanceLags;
        :negLFlows; :months_notice; :fDate; :main_strategy]

    #Only interested in redemptions
    HFData = HFData[~isna(HFData[:,:negLFlows]) & (HFData[:,:negLFlows] .< 0.0), arraySyms]

    println("Ending wide rows: $(size(HFData,1))")
    #Stack the performance over each redemption
    HFLong::DataFrame = melt(HFData, [:negLFlows,:months_notice, :fDate, :main_strategy])
    println("Starting long rows: $(size(HFLong,1))")

    #elliminate all data outside of the notice period
    HFLong[:,:lagVal] = 1 #we will manually dummy the redemption date
    for l::Int ∈ 1:LAGS_PERFORMANCE_TO_CAPTURE
        HFLong[(HFLong[:,:months_notice].≤l) & (HFLong[:,:variable] .==
            performanceLS.performanceLags[l]),:value] = NA
        HFLong[HFLong[:,:variable] .== performanceLS.performanceLags[l],:lagVal] = l
    end

    pool!(HFLong, :lagVal) #convert this to a factor

    HFLong = HFLong[~isna(HFLong[:,:value]),:]
    HFLong = HFLong[HFLong[:,:months_notice].≥MIN_MONTHS_NOTICE_RD,:]

    for i::Int ∈ 1:length(thresholdSyms)
        HFLong[:,thresholdSyms[i]] = 0.0 #sets the column as a Float64
        HFLong[:,thresholdSyms[i]] = NA
        HFLong[HFLong[:,:negLFlows] .< -1.0*(collect(SIG_REDEMPTION_THRESHOLDS)[i]),
            thresholdSyms[i]] = HFLong[HFLong[:,:negLFlows] .<
            -1.0*collect(SIG_REDEMPTION_THRESHOLDS)[i],:negLFlows]
    end

    #make a dummy for redemption date (manually so the results are displayed)
    HFLong[:,:isRedemptionDate] = 0
    HFLong[HFLong[:,:variable].==:performance,:isRedemptionDate] = 1
    #pool!(HFLong, :isRedemptionDate)

    println("Ending long rows: $(size(HFLong,1))")

    #showall(HFLong[1:1000:end,:])

    #save the binary
    oStream = open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_long.jls","w+")
    serialize(oStream, HFLong)
    close(oStream)

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
      sum(HFSub[~isna(HFSub[:,:redemp_notice]),:redemp_notice])
    HFByDate[dateIdx,:totRollingRedemp]=
      sum(HFSub[~isna(HFSub[:,:rolling_notice]),:rolling_notice])
    for s ∈ [noticeLS.noticeLags; noticeLS.cumNoticeLags]
      HFByDate[dateIdx,s]= sum(HFSub[~isna(HFSub[:,s]),s])
    end

  end

  return HFByDate

end

#this is the main method for analysis
function AnalyzeAggregates(;FlowsLaggingProcedure::Char = 'U',
    verbose::Bool = true)

    #open the pre-processed data file
    iStream::IOStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls")
    HFData::DataFrame = deserialize(iStream)
    close(iStream)

    sort!(HFData,cols = (:date, :fund_id))
    HFByDate::DataFrame = aggregateOnDate(HFData, LAGS_NOTICE_TO_CAPTURE)
    performanceLS::PerformanceLagSymbols =
        PerformanceLagSymbols(LAGS_PERFORMANCE_TO_CAPTURE)

    XLMSpecs::Vector{CTExpr} = [parse("negLFlows")]
    XLMNames::Vector{Vector{Symbol}} = [[:intercept; :negLFlows]] #note only non-factor names
    #####put additional specifications here
    push!(XLMSpecs, parse("negLFlows + fDate"))
    push!(XLMNames, [:intercept; :negLFlows])

    push!(XLMSpecs, parse("negLFlows + main_strategy"))
    push!(XLMNames, [:intercept; :negLFlows])

    push!(XLMSpecs, parse("negLFlows + " *
        "$(vec2String(performanceLS.performanceLags, "+"))"))
    push!(XLMNames, [:intercept; :negLFlows; performanceLS.performanceLags])

    push!(XLMSpecs, parse("negLFlows + " *
    "$(vec2String(performanceLS.performanceLags, "+")) + fDate"))
    push!(XLMNames, [:intercept; :negLFlows; performanceLS.performanceLags])

    push!(XLMSpecs, parse("negLFlows + " *
    "$(vec2String(performanceLS.performanceLags, "+")) + main_strategy"))
    push!(XLMNames, [:intercept; :negLFlows; performanceLS.performanceLags])

    push!(XLMSpecs, parse("negLFlows + " *
    "$(vec2String(performanceLS.performanceLags, "+")) + fDate + main_strategy"))
    push!(XLMNames, [:intercept; :negLFlows; performanceLS.performanceLags])


    if RUN_LM_INTERACTIONS
        push!(XLMSpecs, parse("negLFlows + fDate * main_strategy"))
        push!(XLMNames, [:intercept; :negLFlows])
        push!(XLMSpecs, parse("negLFlows + "  *
          "$(vec2String(performanceLS.performanceLags, "+")) + " *
          "fDate * main_strategy"))
        push!(XLMNames, [:intercept; :negLFlows; performanceLS.performanceLags])
    end
    #####

    modelsLM::Vector{CTLM} = Vector{CTLM}()
    #run the LM Specs
    for i ∈ 1:length(XLMSpecs)

        if verbose
            println("Begin Spec RHS: $(XLMSpecs[i])")
        end


        push!(modelsLM, CTLM(HFData, XLMSpecs[i],
            :performance, XNames=XLMNames[i], YName = :performance))

        if verbose
            println("Linear Model X Specification: $(XLMSpecs[i])")
            println("Names: $(XLMNames[i])")
            println("β: $(modelsLM[i].β[1:min(10,end)])")
            println("σ: $(sqrt.(diag(getHomoskedΣ!(modelsLM[i])))[1:min(10,end)])")
            println("Data points: $(size(modelsLM[i].X,1))\n")
        end

    end
    ##########LM IO
    descRowNamesLM::Vector{String} = ["Performance Lags", "Date Fixed Effects",
        "Strategy Fixed Effects", "F.E. Interactions", "N ('000s)"]

    #need to print the descriptive rows
    descContentLM::Vector{Vector{String}} =
        #[fill("",length(modelsLM)) for i∈ 1:length(descRowNamesLM)]
        [Vector{String}(length(modelsLM))for i∈ 1:length(descRowNamesLM)]

    #this builds the descriptive rows. There is Probably a better way to do this,
    #but its fairly specific to the project.
    for i ∈ 1:length(XLMSpecs)
        descContentLM[1][i] = "$(contains(XLMSpecs[i], "performanceLag")?"X":"")"
        descContentLM[2][i] = "$(contains(XLMSpecs[i], "fDate")?"X":"")"
        descContentLM[3][i] = "$(contains(XLMSpecs[i], "main_strategy")?"X":"")"
        descContentLM[4][i] = "$(contains(XLMSpecs[i], "*")?"X":"")"
        descContentLM[5][i] = "$(round(modelsLM[i].N/1000.0,1))"
    end


    TableTextLM::String = texTable(modelsLM, getHomoskedΣ!, [:intercept; :negLFlows],
        titleCaption = "Performance Flows OLS Specifications",
        colNames = [["($i)" for i∈1:length(modelsLM)]],
        contentRowNames = ["intercept", "outflows/AUM x 100"],
        descRowNames = descRowNamesLM,
        descContent = descContentLM,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [100.0,100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = """Ordinary linear regression of outflows on performance, using a
                variety of different covariate. Date fixed effects are by month.
                Performance lags includes the prior 11 months of performance
                as an exogeneous covariate. Data is from alive and dead hedge
                funds 2000-2017. Interpret the focal variable as outflows of 10\\%
                corresponding to a 0.1\\% drop in monthly performance.""")

    if USE_AGGRESSIVE_GC
        gc()
    end

    #####################2SLS code###################
    X2SLSSpecs::Vector{CTExpr} = [parse("negLFlows+0")]
    W2SLSSpecs::Vector{CTExpr} = [Symbol("")]
    Z2SLSSpecs::Vector{CTExpr} = [parse("cumNoticeLLag1+0")]

    W2SLSNames::Vector{Vector{Symbol}} = [[:intercept]]
    X2SLSNames::Vector{Vector{Symbol}}  = [[:negLFlows]]
    Z2SLSNames::Vector{Vector{Symbol}}  = [[:cumNoticeLLag1]]

    #####put additional 2SLS specifications here
    push!(W2SLSSpecs, parse("fDate"))
    push!(W2SLSNames, [:intercept])

    push!(W2SLSSpecs, parse("main_strategy"))
    push!(W2SLSNames, [:intercept])

    push!(W2SLSSpecs, parse("$(vec2String(performanceLS.performanceLags, "+"))"))
    push!(W2SLSNames, [:intercept; performanceLS.performanceLags])

    push!(W2SLSSpecs, parse("$(vec2String(performanceLS.performanceLags, "+"))" *
        " + fDate"))
    push!(W2SLSNames, [:intercept; performanceLS.performanceLags])

    push!(W2SLSSpecs, parse("$(vec2String(performanceLS.performanceLags, "+"))" *
        " + main_strategy"))
    push!(W2SLSNames, [:intercept; performanceLS.performanceLags])

    push!(W2SLSSpecs, parse("$(vec2String(performanceLS.performanceLags, "+"))" *
        " + fDate + main_strategy"))
    push!(W2SLSNames, [:intercept; performanceLS.performanceLags])


    if RUN_2SLS_INTERACTIONS #needs too much memory. Would need a mem-mapped array or simialr solution to do this.
        push!(W2SLSSpecs, parse("fDate * main_strategy"))
        push!(W2SLSNames, [:intercept])

        push!(W2SLSSpecs, parse("$(vec2String(performanceLS.performanceLags, "+"))" *
          " + fDate * main_strategy"))
        push!(W2SLSNames, [:intercept; performanceLS.performanceLags])
    end

    ################# run the 2SLS models


    models2SLS::Vector{CT2SLS} = Vector{CT2SLS}()

    for iX ∈ 1:length(X2SLSSpecs), iW ∈ 1:length(W2SLSSpecs), iZ ∈ 1:length(Z2SLSSpecs)
        if verbose #print the specification to run
            println("Spec X: $(X2SLSSpecs[iX])
            Spec W: $(W2SLSSpecs[iW])
            Spec Z: $(Z2SLSSpecs[iZ])")
        end

        #build the regression model
        push!(models2SLS, CT2SLS(HFData, X2SLSSpecs[iX], W2SLSSpecs[iW],
            :performance, Z2SLSSpecs[iZ], XNames=X2SLSNames[iX],
            WNames = W2SLSNames[iW], YName = :performance, ZNames = Z2SLSNames[iZ]))

        if verbose
            println("IV Model Specification:
                X: $(X2SLSSpecs[iX])
                W: $(W2SLSSpecs[iW])
                Z: $(Z2SLSSpecs[iZ])")
            println("WNames: $(W2SLSNames[iW])")
            println("β: $(models2SLS[i].δ2[1:min(10,end)])")
            println("σ: $(sqrt.(diag(getHomoskedΣ!(models2SLS[i])))[1:min(10,end)])")
            println("Data points: $(size(models2SLS[i].X,1))\n")
        end



    end
##################IO Code for 2SLS



    descRowNames2SLS::Vector{String} = ["Performance Lags", "Date Fixed Effects",
        "Strategy Fixed Effects", "N ('000s)"]

    #need to print the descriptive rows
    descContent2SLS::Vector{Vector{String}} =
        #[fill("",length(modelsLM)) for i∈ 1:length(descRowNamesLM)]
        [Vector{String}(length(models2SLS))for i∈ 1:length(descRowNames2SLS)]

    #this builds the descriptive rows. There is Probably a better way to do this,
    #but its fairly specific to the project.
    for i ∈ 1:length(W2SLSSpecs)
        descContent2SLS[1][i] = "$(contains(W2SLSSpecs[i], "performanceLag")?"X":"")"
        descContent2SLS[2][i] = "$(contains(W2SLSSpecs[i], "fDate")?"X":"")"
        descContent2SLS[3][i] = "$(contains(W2SLSSpecs[i], "main_strategy")?"X":"")"
        descContent2SLS[4][i] = "$(round(models2SLS[i].N/1000.0,1))"
    end
    #note this only works for a single focal var
    stage12SLS::Vector{CTLM} = [get1stStage(m)[1] for m ∈ models2SLS]
    TableText1stStage::String = texTable(stage12SLS, getHomoskedΣ!, [:intercept; Z2SLSNames[1][1]],
        titleCaption = "Performance Flows IV 1st Stage (Focal: Outflows)",
        colNames = [["($i)" for i∈1:length(models2SLS)]],
        contentRowNames = ["intercept", "Net Notice/AUM x100"],
        descRowNames = descRowNames2SLS,
        descContent = descContent2SLS,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [100.0,100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = """1st stage of 2SLS specification uses cumulative negative flows
            one-month prior as the instrumental varaible. Date fixed effects are by month.
                Performance lags includes the prior 11 months of performance
                as an exogenous covariate. Data is from alive and dead hedge
                funds 2000-2017. Interpret a 1.0 value for the focal variable as a
                net outflow notification of 10\\% corresponding to a 0.1\\% realized
                outflow at the redemption date.""")


    TableText2SLS::String = texTable(models2SLS, getHomoskedΣ!, [:intercept; :negLFlows],
        titleCaption = "Performance Flows 2SLS Specifications",
        colNames = [["Z=1 Month Notification of Flows"],["($i)" for i∈1:length(models2SLS)]],
        contentRowNames = ["intercept", "outflows/AUM x 100"],
        descRowNames = descRowNames2SLS,
        descContent = descContent2SLS,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [100.0,100.0],
        widthColNames = [[length(W2SLSSpecs)],ones(Int,length(W2SLSSpecs))],
        clearMem = USE_AGGRESSIVE_GC,
        caption = """2SLS specification uses cumulative negative flows one-month
            prior as the instrumental varaible. Date fixed effects are by month.
            Performance lags includes the prior 11 months of performance
            as an exogenous covariate. Data is from alive and dead hedge
            funds 2000-2017. Interpret a value of 1.0 for the focal variable as
            outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance.
            LATE does not apply as instrument is not randomly selected.""")

    if USE_AGGRESSIVE_GC
        gc()
    end
    #######Now the RD Approach

    #read in the RD long file
    iStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_long.jls")
    HFLong::DataFrame = deserialize(iStream)
    close(iStream)

    thresholdSyms::Vector{Symbol} = getSigRedemptionThresholds()
    modelsRD::Vector{CTLM} = Vector{CTLM}()

    XRDSpecs::Vector{CTExpr} = Vector{CTExpr}()
    XRDNames::Vector{Vector{Symbol}}  = Vector{Vector{Symbol}}()


    #push the specifications (use a common name to make a clean table)
    for s::Symbol ∈ thresholdSyms
        push!(XRDSpecs, parse("$s+isRedemptionDate"))
        push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

        push!(XRDSpecs, parse("$s+isRedemptionDate&$s"))
        push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    end
    #=push!(XRDSpecs, parse("negLFlows+isRedemptionDate&negLFlows"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate+ isRedemptionDate&negLFlows"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate+ isRedemptionDate&negLFlows+lagVal"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate + isRedemptionDate&negLFlows +" *
        "isRedemptionDate&main_strategy+lagVal+main_strategy"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate+ isRedemptionDate&negLFlows  +" *
        "isRedemptionDate&main_strategy + fDate+main_strategy"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate + isRedemptionDate&negLFlows + fDate+lagVal"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])

    push!(XRDSpecs, parse("negLFlows+isRedemptionDate+isRedemptionDate&negLFlows+" *
        "isRedemptionDate&main_strategy+fDate+lagVal + main_strategy"))
    push!(XRDNames, [:intercept; :negLFlows; :isRedemptionDate])=#

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
        ((j::Int)->contains(XRDSpecs[i],"$(thresholdSyms[j])")?
            Int(100.0*collect(SIG_REDEMPTION_THRESHOLDS)[j]):0).(1:length(thresholdSyms))
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
             ) & (HFLong[:,:value].<.5) & (HFLong[:,:value].>-.5),:],
        y=:value,color=:variable, x=:negLFlows,  #=shape=:variable, size=[.1mm],=#
        Stat.binmean(n=25), Geom.line,
        Guide.title("Performance vs Redemption Notifications & Flows"),
        Guide.ylabel("Performance"),Guide.xlabel("Flows"))
      draw(PNG("$RESULTS_PATH\\$(DATA_FILE_NAME)_RD.png", 7.5inch, 6inch), plOutRD)


    TableTextRD::String = texTable(modelsRD, getHomoskedΣ!, [:intercept; :negLFlows; :isRedemptionDate],
        titleCaption = "Performance Flows RD Specifications",
        colNames = [["($i)" for i∈1:length(modelsRD)]],
        contentRowNames = ["intercept", "outflows/AUM x 100", "redemption date"],
        descRowNames = descRowNamesRD,
        descContent = descContentRD,
        decimalDigits = 2,
        columnSepPt = -20,
        scaling = [100.0,100.0,100.0],
        clearMem = USE_AGGRESSIVE_GC,
        caption = """RD specification uses the redemption date as the discontinuity,
            relative to performance during the redemption's notificaiton period.
            Date fixed effects are by month. Strategy is interacted with
            the redemption date dummy. Data is from alive and dead hedge
            funds 2000-2017. Interpret a value of 1.0 for the focal variable as
            outflows of 10\\% corresponding to a 0.1\\% drop in monthly performance.""")

    writeTables2File([TableTextLM, TableText1stStage, TableText2SLS, TableTextRD],
        HEADER_NAME, FOOTER_NAME, path=RESULTS_PATH,
        outName = "RegTables_$FlowsLaggingProcedure.tex")

    #=plOut::Plot = plot(melt(HFByDate[:,[:date,
    #:totFlows,  :totRedemp, :noticeLag1,
    :totRedempNotice,:noticeLag1,:noticeLag2]],:date),
        x=:date,y=:value,color=:variable,Geom.line,
        Guide.title("Total Fund Flows over Time"),
        Guide.xlabel("Date"),Guide.ylabel("Flows"))
      draw(PNG("$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_flows.png", 7.5inch, 6inch), plOut)
      =#

end

function verify2R(;FlowsLaggingProcedure::Char = 'U')

    #read the seria
    iStream::IOStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls")
    HFData::DataFrame = deserialize(iStream)
    close(iStream)

    #read in the RD long file
    iStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_long.jls")
    HFLong::DataFrame = deserialize(iStream)
    close(iStream)


    #write the csv
    writetable("$R_TEST_PATH\\RHFData.csv",HFData)
    writetable("$R_TEST_PATH\\RHFLong.csv",HFLong)
end




@time begin
PreProcessHFData(FlowsLaggingProcedure='C', newBinary=false)
processLongHFData(FlowsLaggingProcedure='C')
AnalyzeAggregates(FlowsLaggingProcedure='C', verbose=false)
verify2R(FlowsLaggingProcedure='C')
end

##Full Pass:

#=@time for c ∈ ['C', 'R', 'L', 'U']
  PreProcessHFData(FlowsLaggingProcedure=c, newBinary=true)
  AnalyzeAggregates(FlowsLaggingProcedure=c)
end=#




end
