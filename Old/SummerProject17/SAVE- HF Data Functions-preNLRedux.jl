
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
  const DATA_FILE_NAME = "HFPerformance_2015-2017"
  const DATA_FILE_NAME_DEAD = DATA_FILE_NAME * "_dead"
  const DATE_FORMAT_STR = "yyyymmdd"
  const NOTICE_PERIOD_ADJ = 7 #days to add to the notice (since data are monthly)
  const MIN_PERIODS_FOR_INCLUSION = 6 #number of periods required for inclusion
  const T_TEST_THRESHOLD = 2.0 #this is a parameter for identificaiton if the flows are lagged
  const PERIODS_SAME_AUM_TO_DELETE = 4 #assume 4 consecutive periods of identical aum implies stale data
  const LAGS_NOTICE_TO_CAPTURE = 3
  const LAGS_PERFORMANCE_TO_CAPTURE = 12
  const MIN_ASSETS = 10.0 #10,000,000
  const START_LAGS_FROM_NOTICE = true #start the performance lags from the notice period (not the redemption date)
  const RUN_LM_INTERACTIONS = false #run the interaction regressions (takes a while)
  const RUN_2SLS_INTERACTIONS =  false
  const RESULTS_PATH = pwd() * "\\results"
  const FOOTER_NAME = "footer.tex"
  const HEADER_NAME = "header.tex"
  const USE_AGGRESSIVE_GC = false #this enables very agressive GC which helps memory management at cost of performance
  const R_TEST_PATH = "C:\\Users\\Clinton\\Dropbox\\Projects\\Summer17RTester\\Summer17RTester"


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
  noticeLagsNo1st::Vector{Symbol}
  noticeLLags::Vector{Symbol}
  noticeLLagsNo1st::Vector{Symbol}
end

#Constructor for the above type
#IN: Number of lags
#OUT: A container of symbols with the appropraite lags
function NoticeLagSymbols(lags::Int)::NoticeLagSymbols
  noticeLags::Vector{Symbol} = Vector{Symbol}(lags)
  noticeLagsNo1st::Vector{Symbol} = similar(noticeLags)
  noticeLLags::Vector{Symbol} = similar(noticeLags)
  noticeLLagsNo1st::Vector{Symbol} = similar(noticeLags)

  for i ∈ 1:lags
    noticeLags[i] = Symbol("noticeLag$i")
    noticeLagsNo1st[i] = Symbol("noticeLag$(i)No1st")
    noticeLLags[i] = Symbol("noticeLLag$i")
    noticeLLagsNo1st[i] = Symbol("noticeLLag$(i)No1st")
  end

  return NoticeLagSymbols(noticeLags,noticeLagsNo1st, noticeLLags,
    noticeLLagsNo1st)
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
  HFData = HFData[:, [:fund_id, :inception, :main_strategy,
    :sub_strategy, :leverage, :fund_assets, :firm_assets,
    :returns_denomination, :management_fee, :incentive_fee, :high_watermark,
    :hurdle_rate, :subscriptions, :redemptions,
    :advance_notice, :lockup, :fund_assets_as_of,
    :fund_assets_denomin, :sales_fee, :other_fees,
    :date_added_to_db, :fund_status, :ucitsiii, :date, :performance,
    :nav, :assets]]

  #Drop some null valued fields:
  HFData = dropNullsFromDF(HFData, [:date, :performance, :assets,
  :returns_denomination, :management_fee,:advance_notice])

  #only want active funds
  #HFData = HFData[HFData[:,:fund_status].=="Active",:]

  delete!(HFData,:fund_status) #no need for this anymore
  HFData = HFData[HFData[:,:assets] .≥ MIN_ASSETS, :]


  #parse the incentive fees and mangement fees
  DFStringtoMeanNum!(HFData, :incentive_fee)
  DFStringtoMeanNum!(HFData, :management_fee)

  #convert from percentage to decimal
  HFData[:,:incentive_fee] ./= 100.0
  HFData[:,:management_fee] ./= 100.0
  HFData[:,:performance] ./= 100.0

  #broadcast an anonymous functino to do the type conversions
  HFData[:,:inception] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:inception])
  HFData[:,:date_added_to_db] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:date_added_to_db])
  HFData[:,:date] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:date])
  HFData[:,:fund_assets_as_of] = ((x::Int)->Date("$x",dtFormat)).(HFData[:,:fund_assets_as_of])

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
    HFData[:,noticeLS.noticeLagsNo1st[i]] = 0.0
    HFData[:,noticeLS.noticeLLags[i]] = 0.0
    HFData[:,noticeLS.noticeLLagsNo1st[i]] = 0.0
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
      if !isna(HFSub[i,:flows]) && HFSub[i,:flows] < 0.0 &&
        i::Int > HFSub[i,:months_notice]
        HFSub[i-HFSub[i,:months_notice],:redemp_notice] = HFSub[i,:flows]
        HFSub[(i-HFSub[i,:months_notice]):(i),:rolling_notice] += HFSub[i,:flows]
        HFSub[i-HFSub[i,:months_notice],:redemp_notice_dt] = HFSub[i,:date]

        #record the redemption notificaion lags and performance lags
        if HFSub[i,:months_notice] > 0
          for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, HFSub[i,:months_notice])
            HFSub[i-l, noticeLS.noticeLags[l]] = HFSub[i,:flows]
            HFSub[i-l, noticeLS.noticeLLags[l]] = HFSub[i,:lFlows]
          end
        end

        #record the redemption notificaion lags, excluding the first months' notice
        if HFSub[i,:months_notice] > 1
          for l::Int ∈ 1:min(LAGS_NOTICE_TO_CAPTURE, HFSub[i,:months_notice]-1)
            HFSub[i-l, noticeLS.noticeLagsNo1st[l]] = HFSub[i,:flows]
            HFSub[i-l, noticeLS.noticeLLagsNo1st[l]] = HFSub[i,:lFlows]
          end
        end


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


  println("Processing complete.")
  println("Ending rows: ", size(HFData,1), ", Ending columns: ", size(HFData,2))

  #save the binary
  oStream = open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls","w+")
  serialize(oStream, HFData)
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
  for s ∈ [noticeLS.noticeLags; noticeLS.noticeLagsNo1st]

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
    for s ∈ [noticeLS.noticeLags; noticeLS.noticeLagsNo1st]
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
        caption = "Performance Flows OLS Specifications",
        colNames = [["($i)" for i∈1:length(modelsLM)]],
        contentRowNames = ["intercept", "outflows/AUM x 100"],
        descRowNames = descRowNamesLM,
        descContent = descContentLM,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [1000.0,1000.0],
        clearMem = USE_AGGRESSIVE_GC)

    if USE_AGGRESSIVE_GC
        gc()
    end

    #####################2SLS code###################
    X2SLSSpecs::Vector{CTExpr} = [parse("negLFlows+0")]
    W2SLSSpecs::Vector{CTExpr} = [Symbol("")]
    Z2SLSSpecs::Vector{CTExpr} = [parse("noticeLag1+0")]

    W2SLSNames::Vector{Vector{Symbol}} = [[:intercept]]
    X2SLSNames::Vector{Vector{Symbol}}  = [[:negLFlows]]
    Z2SLSNames::Vector{Vector{Symbol}}  = [[:noticeLag1]]

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
        println(models2SLS[end].ZNames)
        println(models2SLS[end].Π1)

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
    TableText1stStage::String = texTable(stage12SLS, getHomoskedΣ!, [:intercept; :noticeLag1],
        caption = "Performance Flows IV 1st Stage (Focal: Outflows)",
        colNames = [["($i)" for i∈1:length(models2SLS)]],
        contentRowNames = ["intercept", "Notifications (frac of fund x1000)"],
        descRowNames = descRowNames2SLS,
        descContent = descContent2SLS,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [1000.0,1000.0],
        clearMem = USE_AGGRESSIVE_GC)


    TableText2SLS::String = texTable(models2SLS, getHomoskedΣ!, [:intercept; :negLFlows],
        caption = "Performance Flows 2SLS Specifications",
        colNames = [["Z=1 Month Notification of Flows"],["($i)" for i∈1:length(models2SLS)]],
        contentRowNames = ["intercept", "outflows/AUM x 1000"],
        descRowNames = descRowNames2SLS,
        descContent = descContent2SLS,
        decimalDigits = 3,
        columnSepPt = -20,
        scaling = [1.0,1000.0],
        widthColNames = [[length(W2SLSSpecs)],ones(Int,length(W2SLSSpecs))],
        clearMem = USE_AGGRESSIVE_GC)


    writeTables2File([TableTextLM, TableText1stStage, TableText2SLS], HEADER_NAME, FOOTER_NAME, path=RESULTS_PATH,
        outName = "RegTables_$FlowsLaggingProcedure.tex")
    #=plOut::Plot = plot(melt(HFByDate[:,[:date,
    #:totFlows,  :totRedemp, :noticeLag1,
    :totRedempNotice,:noticeLag1,:noticeLag2]],:date),
        x=:date,y=:value,color=:variable,Geom.line,
        Guide.title("Total Fund Flows over Time"),
        Guide.xlabel("Date"),Guide.ylabel("Flows"))
      draw(PNG("$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_flows.png", 7.5inch, 6inch), plOut)
      =#
    close(iStream)
end

function verify2R(;FlowsLaggingProcedure::Char = 'U')

    #read the seria
    iStream::IOStream =
        open("$DATA_PATH\\$(DATA_FILE_NAME)_$(FlowsLaggingProcedure)_clean.jls")
    HFData::DataFrame = deserialize(iStream)
    close(iStream)

    #write the csv
    writetable("$R_TEST_PATH\\RHFData.csv",HFData)
end




@time begin
#PreProcessHFData(FlowsLaggingProcedure='C', newBinary=false)
AnalyzeAggregates(FlowsLaggingProcedure='C', verbose=false)
verify2R(FlowsLaggingProcedure='C')
end

##Full Pass:

#=@time for c ∈ ['C', 'R', 'L', 'U']
  PreProcessHFData(FlowsLaggingProcedure=c, newBinary=true)
  AnalyzeAggregates(FlowsLaggingProcedure=c)
end=#




end
