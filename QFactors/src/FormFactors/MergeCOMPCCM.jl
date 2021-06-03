
#MAIN ENTRY POINT
#final columns for the decompressed file

#merges crsp and compustat
# note refresh merge will be automatically run if any of the other two are run
function loadccm(;
  datapath::String=DATA_PATH,
  ccmname::String = CCM_NAME,
  incsvstream::Function = IN_CSV_STREAM)

  local ccm::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  ccm = incsvstream("$datapath\\$ccmname.csv.lz4") |> CSV.read |> DataFrame

  return ccm

end


function mergecompccm!(comp::DataFrame;
    testoutput::Bool = TEST_OUTPUT)

  #first
  local ccm::DataFrame = loadccm()
  ccm = preprocessccm!(ccm)

  testoutput && CSV.write("output\\ccm.csv", ccm)

  #println("Comp rows before merge: $(size(comp, 1))")
  comp = join(comp, ccm, on=:gvkey, kind=:inner)

  #fmatch the dates
  #println("Comp rows after merge, pre-filter: $(size(comp, 1))")
  comp.keep = trues(size(comp,1))

  comp.keep = (r::DataFrameRow->
    r.adateq ∈ r.linkeffdate:Day(1):r.linkenddate).(eachrow(comp))

  comp = comp[comp.keep,:]
  select!(comp, Not(:keep))
  println("Comp rows after merge, post-filter: $(size(comp, 1))")

  comp = reconcilecompccmintervals!(comp)

  return comp
end


#formats the ccm table in a reasonable way
function preprocessccm!(ccm::DataFrame;
  retainedcolsccm::Vector{Symbol} = RETAINED_COLS_CCM,
  lastdatestring::String = LAST_DATE_STRING)

  #names ot lower case
  names!(ccm, (s::Symbol->Symbol(lowercase(string(s)))).(names(ccm)))

  #fix the end date
  ccm.linkenddt .= (s::MString->
    coalesce(s,"")=="E" ? lastdatestring : s).(ccm.linkenddt)

  filter!(r::DataFrameRow->r.linkprim =='P' || r.linkprim =='C', ccm)
  
  #parse the dates
  ccm.linkeffdate = (s->parseccm(Date, s)).(ccm.linkdt)
  ccm.linkenddate = (s->parseccm(Date, s)).(ccm.linkenddt)

  #only need some of the fields
  ccm = select!(ccm, retainedcolsccm)

  return ccm
end

#this function reconciles the link intervals with the announcement intervals
#see
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
