
#MAIN ENTRY POINT
#final columns for the decompressed file

#merges crsp and compustat
# note refresh merge will be automatically run if any of the other two are run
function loadccm(;
  datapath::String=DATA_PATH,
  ccmname::String = CCM_NAME,
  incsvstream::Function = IN_CSV_STREAM,
  csvextension::String = CSV_EXTENSION)

  local ccm::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  ccm = incsvstream("$datapath\\$ccmname.$csvextension") |> CSV.File |> DataFrame

  return ccm

end


function mergecompccm!(comp::DataFrame;
    testoutput::Bool = TEST_OUTPUT,
    months2stale::Int =  MONTHS_2_STALE_LONG)

  #first
  local ccm::DataFrame = loadccm()
  ccm = preprocessccm!(ccm)

  testoutput && CSV.write("output\\ccm.csv", ccm)

  select!(ccm, Not([:cusip]))
  #println("Comp rows before merge: $(size(comp, 1))")
  comp = innerjoin(comp, ccm, on=[:gvkey])
  comp.keep = trues(size(comp,1))

  comp.keep = (r::DataFrameRow->
    r.adate ∈ r.linkeffdate:Day(1):r.linkenddate).(eachrow(comp))

  comp = comp[comp.keep,:]
  select!(comp, Not(:keep))
  println("Comp rows after merge, post-filter: $(size(comp, 1))")

  #Note: we still will have overlapping lineffdates and linkenddates, which will be corrected
  #with the below method
  comp = reconcilecompccmintervals!(comp)

  #  lagwithin!(comp, :fdate, :gvkey, :fyear)
  sort!(comp, [:lpermno, :fdate])
  issorted(comp, [:lpermno, :adate]) || error("I thought comp was sorted by [:lpermno, :adate]")
  period2stale::Month = Month(months2stale)
  lagwithin2!(comp, [:fdate, :fyrmonth], :lpermno, date=:fdate, maxnotstale = period2stale)
  lagwithin2!(comp, [:Lfdate, :Lfyrmonth], :lpermno, date=:fdate, maxnotstale = period2stale)
  for r ∈ eachrow(comp)
    ismissing(r.Lfdate) && continue
    (r.Lfdate < r.fdate -  period2stale) && error("failed date check $r")
  end
  @info "past date validatation check"


  #=
  NOTE: we accept dup fyears due to changing fiscal years
  function testfordups(df::AbstractDataFrame, groupfields::Vector{Symbol})
    for sdf ∈ groupby(df, groupfields)
      if size(sdf,1) > 1
        println(sdf[!, [:gvkey, :lpermno, :fdate,
          :Lfdate, :adate, :Nadate, :linkeffdate, :linkenddate, :begindate, :enddate, :fyear]])
        error("duplicate fyear found")
      end
    end
  @info "passed dedup fyear test"
  end
  testfordups(comp, [:lpermno, :fyear])=#

  return comp
end


#formats the ccm table in a reasonable way
function preprocessccm!(ccm::DataFrame;
  retainedcolsccm::Vector{Symbol} = RETAINED_COLS_CCM,
  lastdatestring::String = LAST_DATE_STRING,
  testoutput::Bool = TEST_OUTPUT)

  #names ot lower case
  cleannames::Vector{String} = (s::String->lowercase(s)).(names(ccm))
  rename!(ccm, ((oldn, newn)->oldn=>newn).(names(ccm), cleannames))

  #fix the end date
  ccm.linkenddt .= (s::MString->
    coalesce(s,"")=="E" ? lastdatestring : s).(ccm.linkenddt)

  #below didn't work for some reason
  #filter!(r::DataFrameRow->((!ismissing(r.gvkey)) && (!ismissing(r.lpermno))), ccm)
  ccm = ccm[(r::DataFrameRow->(
    (!ismissing(r.gvkey)) && (!ismissing(r.lpermno)))).(eachrow(ccm)),:]
  ccm.linkprim = (s::String->length(s)==1 ? s[1] : error("invalid linkprim in ccm")).(ccm.linkprim)
  filter!(r::DataFrameRow->r.linkprim =='P' || r.linkprim =='C', ccm)

  #parse the dates
  ccm.linkeffdate = (s->parseccm(Date, s)).(ccm.linkdt)
  ccm.linkenddate = (s->parseccm(Date, s)).(ccm.linkenddt)

  #dedup if possible
  ccm.tokeep = trues(size(ccm,1))
  ccm.dateranges = ((linkeffdate::Date,linkenddate::Date)->
    linkeffdate:Day(1):linkenddate).(ccm.linkeffdate,ccm.linkenddate)

  gccm::GroupedDataFrame = groupby(ccm, [:gvkey, :lpermno])
  @mpar for i ∈ 1:length(gccm)
    sccm::SubDataFrame = gccm[i]
    Nsccm::Int = size(sccm,1)

    #construct the union of all of the date sets
    if Nsccm >  1
      local periodrange::StepRange{Date,Day} = minimum(
        sccm.linkeffdate):Day(1):maximum(sccm.linkenddate)
      local dateunionlength::Int = length(unique([sccm.dateranges...;]))

      if dateunionlength == length(periodrange)
        sccm.linkeffdate[1] = minimum(periodrange)
        sccm.linkenddate[1] = maximum(periodrange)
        sccm.tokeep[2:end] .= false
      end

    end
  end

  println("ccm size pre-dedup $(size(ccm))")
  ccm = ccm[ccm.tokeep,:]
  println("ccm size post-dedup $(size(ccm))")


  if testoutput
    testoutput && CSV.write("output\\ccm_preproc.csv", ccm)
  end
  #rename!(ccm, :indfmt=>:indfmtold)
  #ccm.indfmt = (s::MString->parsecomp(Symbol, s)).(ccm.indfmtold)

  #only need some of the fields
  ccm = select!(ccm, retainedcolsccm)

  return ccm
end

#this function reconciles the link intervals with the announcement intervals
#NOTE: A weird situation can arise if gvkey changes. The the linkenddate finishes too soon,
#we end up with a gap . Considering the laternatives, it seems like the lesser evil
#to extending the record beyond the linkenddate
function reconcilecompccmintervals!(comp::DataFrame)
  #subtract one day since the effective interval is inclusive of the
  #first date while the begin-end date is exclusve of the first date
  comp.begindate .= (max).(comp.begindate, comp.linkeffdate .- Day(1))
  comp.enddate .= (min).(comp.enddate, comp.linkenddate)



  #this gets rid of one day records
  #(this assumes the previous effective record is the current record)
  #not many of these ~140
  #println("records before interval reconciliation: $(size(comp,1))")
  #println(comp[comp.begindate .≥ comp.enddate,
  #  [:gvkey, :adateq, :begindate, :enddate, :linkeffdate, :linkenddate]])
  comp = comp[comp.begindate .< comp.enddate, :]
  println("records after interval reconciliation: $(size(comp,1))")

  return comp
end


function parseccm(::Type{T}, s::String;
    #=parsedmissings = CCM_PARSED_MISSINGS=#) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  #(!ismissing(v)) && (v ∈ parsedmissings) && (v=missing)

  return v
end

#the date case
parseccm(::Type{Date}, s::String;
    ccmdateformat::DateFormat = CCM_DATE_FORMAT) = Dates.Date(s, ccmdateformat)

#helper methods and special cases
parseccm(::Type{<:Any}, v::Missing) = missing
parseccm(::Type{Date}, i::Int) = parseccm(Date, "$i")
parseccm(::Type{Symbol}, s::String) = Symbol(s)
