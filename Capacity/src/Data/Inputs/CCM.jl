

#formats the ccm table in a reasonable way
function prepccm(;
  ccmname::String = PARAM[:ccmname],
  ccmpath::String = PARAM[:ccmpath],
  ccmcolumns::Vector{Symbol} = PARAM[:ccmcolumns],
  lastdatestring::String = "20991231")

  ccm::DataFrame = CSV.File("$ccmpath\\$ccmname.csv") |> DataFrame

  #names ot lower case
  rename!(ccm, (s::String->lowercase(s)).(names(ccm)))
  rename!(ccm, :lpermno=>:permno)

  #fix the end date
  ccm.linkenddt .= (s::MString->
    coalesce(s,"")=="E" ? lastdatestring : s).(ccm.linkenddt)

  #below didn't work for some reason
  #filter!(r::DataFrameRow->((!ismissing(r.gvkey)) && (!ismissing(r.lpermno))), ccm)
  ccm = ccm[completecases(ccm, [:gvkey, :permno]),:] #only want valid links
  ccm.linkprim = (s::String->length(s)==1 ? s[1] : error("invalid linkprim in ccm")).(ccm.linkprim)
  ccm = ccm[(ccm.linkprim .== 'P') .| (ccm.linkprim .== 'C'), :] #avoid secondary issues

  #parse the dates
  ccm.ccmeffdate = (s->parseccm(Date, s)).(ccm.linkdt)
  ccm.ccmenddate = (s->parseccm(Date, s)).(ccm.linkenddt)

  ccm = select!(ccm, ccmcolumns)

  return ccm
end

#WARNING WARNING WARNING- a big difference between here and the leverage changes
#code is the use of inclusive eff and end dates (I think that the begin date was exclusive
#in the leverage changes code)
function mergecompccm(comp::DataFrame;
    maxadateinterval::DatePeriod =  PARAM[:compmaxadateinterval],
    defaultadateinterval::DatePeriod =  PARAM[:compdefaultadateinterval],
    )

  #first
  ccm::DataFrame = prepccm()
  #testoutput && CSV.write("output\\ccm.csv", ccm)

  (isempty(comp[nonunique(comp, [:gvkey, :fdate]), :])) || error("non-unique
  gvkey-fdate records detected, dedup prior to merging ccm")
  @assert comp[nonunique(comp, [:gvkey, :adate]),:] |> isempty

  #sort and preliminary contract
  @assert issorted(comp, [:gvkey, :fdate])
  @assert issorted(comp, [:gvkey, :adate])

  #lead the adate to form the announcement interval, and fix intervals where the end is missing
  #or the interval is too long
  comp.beginadate = comp.adate |> deepcopy
  leadwithin2!(comp, [:beginadate], :gvkey, date=:fdate)
  comp.endadate = incrementdt.(comp.Nbeginadate, Day(-1)) #use inclusive, non-overlapping ranges
  #rename!(comp, :Nbeginadate => :endadate)
  badenddates = view(comp, (comp.endadate .=== missing) .|
    (comp.beginadate .+ maxadateinterval .< comp.endadate), [:beginadate, :endadate])
  badenddates.endadate .=  badenddates.beginadate .+ defaultadateinterval

  #integrity checks
  @assert comp[nonunique(comp, [:gvkey, :endadate]),:] |> isempty
  @assert (comp.endadate .≥ comp.beginadate) |> all
  @assert (comp.endadate .≤ comp.beginadate .+ maxadateinterval) |> all
  @assert issorted(comp, [:gvkey, :endadate])

  comp = innerjoin(comp, ccm, on=:gvkey)
  sort!(comp, [:gvkey, :endadate])
  comp = filter!(comp) do r #drop rows with no interseciton between effective link ranges and adate interals
    (r.ccmeffdate ≤ r.beginadate) &&
    ((r.ccmenddate === missing) || (r.ccmenddate ≥ r.endadate))
  end
  #want the intersection of the link range and the announcement range
  comp.effdate = max.(comp.beginadate, comp.ccmeffdate)
  comp.enddate = min.(comp.endadate, comp.ccmenddate)
  comp.effrange = ((effdate,enddate)->effdate:Day(1):enddate).(comp.effdate, comp.enddate)
  @assert (comp.effdate .≤ comp.enddate) |> all

  #date range checks
  #use the intersection of the link range and the announcement range
  adateranges = ((beginadate, endadate)->
    beginadate:Day(1):endadate).(comp.beginadate,comp.endadate)
  ccmdateranges = ((ccmeffdate, ccmenddate)->
    ccmeffdate:Day(1):ccmenddate).(comp.ccmeffdate,comp.ccmenddate)
  effrangecheck = ((arange, ccmrange)->intersect(arange, ccmrange)).(adateranges, ccmdateranges)
  @assert (comp.effrange .== effrangecheck) |> all

  #integrity checks on the dates
  @assert view(comp, (comp.effdate .=== missing) .| (comp.enddate .=== missing), :) |> isempty
  @assert (comp.enddate .≥ comp.effdate) |> all
  @assert (comp.enddate .≤ comp.effdate .+ maxadateinterval) |> all

  #validates the exclusivity of dates within a grouping scheme
  function validatedateswithingrouping(Fgrp::Symbol)
    @assert issorted(comp, [Fgrp, :effdate])
    @assert issorted(comp, [Fgrp, :enddate])

    #check that the intervals are exclusive
    leadwithin2!(comp, [:effdate], Fgrp, date=:fdate)
    lagwithin2!(comp, [:enddate], Fgrp, date=:fdate)

    #check that the start of the next period is always greater than the end of the preceding
    completeNeffdate = view(comp, completecases(comp, [:Neffdate]), :)
    @assert (completeNeffdate.Neffdate .> completeNeffdate.enddate) |> all

    #same as before but in reverse- check end of preceding is less than the start of the current
    completeLenddate = view(comp, completecases(comp, [:Lenddate]), :)
    @assert (completeLenddate.Lenddate .< completeLenddate.effdate) |> all

    select!(comp, Not([:Neffdate, :Lenddate]))
  end

  #validate with the current gvkey groupings
  validatedateswithingrouping(:gvkey)

  #repeat the validaiton excercise, but this time with :permno
  sort!(comp, [:permno, :fdate])
  validatedateswithingrouping(:permno)


  return comp
end



function parseccm(::Type{T}, s::String;) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)


  return v
end

#the date case
parseccm(::Type{Date}, s::String;
    ccmdateformat::DateFormat = dateformat"yyyymmdd") = Dates.Date(s, ccmdateformat)

#helper methods and special cases
parseccm(::Type{<:Any}, v::Missing) = missing
parseccm(::Type{Date}, i::Int) = parseccm(Date, "$i")
parseccm(::Type{Symbol}, s::String) = Symbol(s)
