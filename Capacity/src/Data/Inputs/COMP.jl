

function prepcomp(;
  compname::String = PARAM[:compname],
  companame::String = PARAM[:companame],
  compqname::String = PARAM[:compqname],
  comppath::String = PARAM[:comppath],
  workingpath::String = PARAM[:workingpath],
  refreshcomp=PARAM[:refreshcomp],
  incsvstream::Function = IN_CSV_STREAM,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  csvextension::String = CSV_EXTENSION,
  binextension::String = BIN_EXTENSION,
  compcolumns::Union{Nothing, Vector{Symbol}} = nothing,
  comprequiredcolumns::Union{Nothing, Vector{Symbol}} = nothing)::DataFrame

  local compa::DataFrame
  local compq::DataFrame
  local comp::DataFrame

  local binpath::String = "$workingpath\\$compname.$binextension"

  throw("ERROR: computstat and IBES data needs to be updated before use")

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcomp
    compq = incsvstream(
      "$comppath\\$compqname.$csvextension") |> CSV.File |> f->DataFrame(f,copycols=false)

    #comp = preprocesscomp(compa, compq)
    comp = preprocesscomp(compq)
    outbinstream(binpath, comp)

  else
    comp = inbinstream("$workingpath\\$compname.$binextension")

  end


  @info "Compustat and CCM linking data loaded and saved into file $compname.$binextension"

  return comp
end




function preprocesscomp(compq::DataFrame;
    columnstokeep::Vector{Symbol}=PARAM[:compcolumnstokeep]
    )

  compq = preprocesscompq(compq)

  #TODO: chekc this!!- is there anything I should use before I delete the below??
  #validatefadates!(comp; Fgroup=:gvkey, Ffdate=:fdate, Fadate=:adate)

  @assert issorted(compq, [:gvkey, :fdate, :adate])
  compq = mergecompccm(compq)

  @info "summary of comp df prior to column pruning"
  println(describe(compq))

  @assert setdiff(columnstokeep, propertynames(compq)) |> isempty
  select!(compq, columnstokeep)
end


#date shaping for compq
function preprocesscompq(compq::DataFrame;
  maxtimetofile::DatePeriod=PARAM[:compmaxtimetofile],
  maxadateinterval::DatePeriod=PARAM[:compmaxadateinterval],
  maxfdateinterval::DatePeriod=PARAM[:compmaxfdateinterval],
  maxannualfdateinterval=PARAM[:compmaxannualfdateinterval],
  defaultadatelag::DatePeriod=PARAM[:compdefaultadatelag],
  bookequityfloor::Float64 = PARAM[:compbookequityfloor],
  yearrange::UnitRange{Int} = PARAM[:compyearrange],
  )

  rename!(compq, (s::String->lowercase(s)).(names(compq)))

  #format dates and rename for clarity
  function formatwrdsdate(p::Pair)
    source, dest = p[1], p[2]
    compq[!, dest] = (s->parsecomp(Date, s)).(compq[!, source])
  end
  formatwrdsdate.([:datadate=>:ddate, :fdateq=>:finaldate, :pdateq=>:preliminarydate, :rdq=>:rdq])
  select!(compq, Not([:datadate, :fdateq, :pdateq]))


  rename!(compq, :indfmt=>:indfmtold)
  compq.indfmt = (s::MString->parsecomp(Symbol, s)).(compq.indfmtold)
  select!(compq, Not(:indfmtold))

  #get rid of some unused columns
  @assert (compq.datafmt .=== "STD") |> all
  @assert (compq.consol .=== "C") |> all
  @assert (compq.popsrc .=== "D") |> all
  select!(compq, Not([:datafmt, :consol, :popsrc,]))

  #usedcompafields::Vector{Symbol} = USED_COMPA_VALUE_FIELDS,


  compq = compq[(!ismissing).(compq.fyearq), :]


  #ddate is the effective record period end date
  #see http://www.stat.rice.edu/~dobelman/courses/482/faq/FAQ00044.htm
  compq.fdate = compq.ddate |> deepcopy

  #dedup the record effective dates
  #use the following logic
  # 1) if the announcement date (rdq, preminary date, final date) is available,
  #   we capture the most recent.
  # 2) for multicurrency records, take :usd
  # 3) assume and assert that the fiscal year changed, so the fiscal year end month fyr
  #   will change. Take the record using the most recent fiscal year accounting.
  compq=transform(groupby(compq, [:gvkey, :fdate]), :fdate=>length=>:recordsingrp)
  compq.curcdqtb = (!).(compq.curcdq .== "USD") #tiebreaker- want these records to appear first in sort
  compq=transform(groupby(compq, [:gvkey, :fdate, :rdq, :preliminarydate, :finaldate,:curcdqtb]),
    :fdate=>length=>:recordsingrpeasy)
  @info "deduping to gvkey-date pairs. $(sum(compq.recordsingrp .> 1)) gvkey-fdate dups remain
    ($(sum(compq.recordsingrpeasy .> 1)) of which cannot be deduped on
    adate-currency criteria)."

  #compute the age in quarrters of each fyr value
  #(possible issue- fyr changes, then changes back. This will mess that up, but it seems like a very
  #remote 3rd order issue)
  sort!(compq, [:gvkey, :fyr, :fdate])
  compq = transform(groupby(compq, [:gvkey, :fyr]), :fdate=>
    (c->collect(1:length(c))) => :qtrssincefyearchange)

  #make sure we still don't have any dups after including fyr iun the grouping
  @assert sum(nonunique(compq, [:gvkey, :fdate, :rdq, :preliminarydate, :finaldate, :fyr])) == 0
  @assert (!ismissing).(compq.qtrssincefyearchange) |> all

  #order the records by desirability
  sort!(compq, [:gvkey, :fdate, :rdq, :preliminarydate, :finaldate, :curcdqtb, :qtrssincefyearchange])
  compq = combine(first, groupby(compq, [:gvkey, :fdate]))

  @assert compq[nonunique(compq, [:gvkey, :fdate]), :] |> isempty
  @info "Deduping complete. Note $(sum(nonunique(compq, [:gvkey, :fyearq, :fqtr]))) duplicate
    gvkey-fyearq-fqtrs remain, though data effective dates (gvkey-fdate) remain unique"

  #clean up temporary fields
  select!(compq, Not([:curcdqtb, :qtrssincefyearchange, :recordsingrp, :recordsingrpeasy]))

  #now splice in adate and Nadate
  ibes = prepibes()
  compq = leftjoin(compq, ibes, on=[:gvkey, :fdate=>:perioddate], validate=(true,true))
  sort!(compq, [:gvkey, :fdate])

  #Note: we do not allow announcement dates after the start of the subsequent fiscal period
  validadate(::Any, ::Any, ::Missing) = missing
  #validadate(fdate::Date, ::Missing, adate::Date) = ifelse(adate≤fdate,missing, adate)
  function validadate(fdate::Date, maxadate::Date, adate::Date)
    (adate ≤ fdate ||
      maxadate < adate) && return missing
    #adate ≤ fdate && return missing
    return adate
  end

  leadwithin2!(compq, [:fdate], :gvkey; date=:fdate)
  lagwithin2!(compq, [:fdate], :gvkey; date=:fdate)
  compq.maxadate = compq.fdate .+ maxtimetofile #limit on the time until the announcement
  for f ∈ [:ibesadate, :rdq, :preliminarydate, :finaldate]
    compq[!,f] .= validadate.(compq.fdate,compq.maxadate,compq[!, f])
  end
  #println(describe(compq))

  #first use RDQ from compustat, as per
  #https://wrds-www.wharton.upenn.edu/pages/support/support-articles/ibes/difference-between-ibes-earnings-announcement-date-and-compustat-announcement-dates/
  #Next try to get the announcement date from ibes
  #follow https://researchfinancial.wordpress.com/2014/02/26/publication-dates-of-annual-reportsearnings/
  #Then pdateq (preliminary data date), then fdateq
  # if none are present, use 60 days
  #=compq = transform(compq,
    [:rdq, :ibesadate, :preliminarydate, :finaldate] =>
      ByRow((rdq, ibesadate, preliminarydate, finaldate)->
        coalesce(rdq, ibesadate, preliminarydate, finaldate, missing)) => :adate)=#
  transform!(compq,
    [:rdq, :ibesadate, :preliminarydate, :finaldate] =>
      ((rdq, ibesadate, preliminarydate, finaldate)->
        coalesce.(rdq, ibesadate, preliminarydate, finaldate, missing)) =>
        :adate)
  select!(compq, Not([:ibesadate, :rdq, :preliminarydate, :finaldate]))


  compq.adateisdefault = falses(size(compq,1))
  #set default adates, but flag them
  badadates = view(compq, compq.adate .=== missing, :)
  badadates.adate .= badadates.fdate .+ defaultadatelag
  #set a default adate with a flag to note that it is a default
  badadates.adateisdefault .= true

  #do not allow overlapping adates NOTE: careful of ordering of these commandes
  leadwithin2!(compq, [:adate], :gvkey; date=:fdate)
  @info "Deleting $(sum(((ad,Nad)->
    (ad !== missing) && (Nad !== missing) && ad≥Nad).(compq.adate, compq.Nadate)))/$(size(compq,1)
    ) bad (overlapping) adates"
  compq.adate[((ad,Nad)->
    (ad !== missing) && (Nad !== missing) && ad≥Nad).(compq.adate, compq.Nadate)] .= missing

  @assert issorted(compq, [:gvkey, :fdate, :adate])
  #sort!(compq, [:gvkey, :fdate, :adate])

  #select nonunique gvkey adate pairs and select the most recent fdate
  validpairs = combine(last, groupby(compq, [:gvkey, :adate]))
  validpairs = validpairs[completecases(validpairs, [:gvkey, :fdate, :adate]),[:gvkey, :fdate, :adate]]
  validpairs.validadate = trues(size(validpairs,1))
  prevrows::Int = size(compq,1)
  #the below join ensures a msising validpairs value for any invalid gvkey-fdate-adate pairing
  compq = leftjoin(compq, validpairs, on=[:gvkey,:fdate, :adate],
    matchmissing=:equal, validate=(false,true))
  sort!(compq, [:gvkey, :fdate, :adate])
  @assert size(compq,1) == prevrows
  @info "Dup adates: $(sum((compq.validadate .=== missing) .& (compq.adate .!== missing)))"
  compq.adate[compq.validadate .=== missing] .= missing


  #drop missing adates
  leadwithin2!(compq, [:adate, :fyearq, :fqtr], :gvkey; date=:fdate) #propogate the deletions forward
  @info "Number of records = $(size(compq,1)), missing adates = $(sum(ismissing.(compq.adate)))"

  #assert adate is within the allowable period
  @assert ((fd, maxad, ad)->(ad===missing) ||
    #(Nfd === missing && ad > fd) ||
    (ad > fd && ad ≤ maxad)).(compq.fdate, compq.maxadate, compq.adate) |> all

  #make sure adates don't overlap
  @assert ((compq.adate .=== missing) .| (compq.Nadate .=== missing) .|
    (compq.adate .< compq.Nadate)) |> all

  #assert next fdate is always greater than past fdate
  @assert ((compq.Nfdate .=== missing) .|  (compq.fdate .< compq.Nfdate)) |> all

  #assert unique gvkey-adate pairs
  @assert view(compq, nonunique(compq, [:gvkey, :adate]) .& (compq.adate .!==missing), :) |> isempty


  #except for missing values, all of the relevant date fields should be in sonsecutive order
  @assert issorted(compq, [:gvkey, :fdate])
  @assert issorted(compq[completecases(compq,[:gvkey, :Nfdate]),:], [:gvkey, :Nfdate])
  @assert issorted(compq[completecases(compq,[:gvkey, :adate]),:], [:gvkey, :adate])
  @assert issorted(compq[completecases(compq,[:gvkey, :Nadate]),:], [:gvkey, :Nadate])
  select!(compq, Not([:Nadate]))

  compa = deepcopy(compq)
  #NOTE: #three alternatives for book equity listed in order of preference
  transform!(compq, [:ceqq, :txditcq, :txdbq] |> AsTable  =>
    ((t)->t.ceqq .+ coalesce.(t.txditcq,t.txdbq, 0.0)) =>
    :bkeceqq)
  transform!(compq, [:seqq, :txditcq, :txdbq, :pstkq, :pstkrq] |> AsTable  =>
    ((t)->t.seqq .+ coalesce.(t.txditcq,t.txdbq, 0.0) .- coalesce.(t.pstkrq,t.pstkq,0.0)) =>
    :bkeseqq)
  transform!(compq, [:atq, :ltq, :txditcq, :txdbq, :pstkq, :pstkrq] |> AsTable  =>
    ((t)->t.atq .- t.ltq .+ coalesce.(t.txditcq,t.txdbq, 0.0) .-
      coalesce.(t.pstkrq,t.pstkq,0.0)) =>
    :bkeaml)
  #take the best metric for bkequity that isnt missing
  assetfloor(x) = max(x,bookequityfloor)
  assetfloor(::Missing) = missing
  transform!(compq, [:bkeceqq, :bkeseqq, :bkeaml] |> AsTable  =>
    ((t) -> coalesce.(t.bkeceqq, t.bkeseqq, t.bkeaml) .|> assetfloor) =>
    :bkequity)
  select!(compq, Not([:bkeceqq, :bkeseqq, :bkeaml]))

  #a couple of alternatives for getting the operating profitability starting
  #at various places on the income statement
  transform!(compq, [:piq, :dpq] |> AsTable =>
    ((t)-> t.piq .+ coalesce.(t.dpq,0.0)) =>
    :oppi)
  transform!(compq, [:ibq, :niq, :dpq, :txtq] |> AsTable =>
    ((t)-> coalesce.(t.ibq,t.niq)  .+ t.txtq .+ coalesce.(t.dpq,0.0)) =>
    :opni)
  transform!(compq, [:revtq, :cogsq, :xsgaq, :xintq] |> AsTable =>
    ((t)-> t.revtq .- t.cogsq .- t.xsgaq .- coalesce.(t.xintq,0.0)) =>
    :oprev)
  transform!(compq, [:oiadpq, :dpq, :xintq] |> AsTable => #oiadpq=ebit
    ((t)->t.oiadpq .+ coalesce.(t.dpq,0.0) .- coalesce.(t.xintq,0.0)) =>
    :opebit)

  #take the best metric for operating profit that isnt missing and divide by bkequity
  #should I average bkequity with Lbkequity?
  transform!(compq, [:oppi, :opni, :oprev, :opebit, :bkequity] |> AsTable => #oiadpq=ebit
    ((t)-> coalesce.(t.oppi, t.opni, t.oprev, t.opebit) ./ t.bkequity) =>
    :op)
  select!(compq, Not([:oppi, :opni, :oprev, :opebit]))
  lagwithin2!(compq, [:atq], :gvkey, date=:fdate)

  #Tcompute the interval of atq, but only if its not too long ago
  transform!(compq, [:atq, :Latq, :fdate, :Lfdate] |> AsTable =>
    ByRow((t)-> (!ismissing(t.Lfdate)) && (t.fdate - maxfdateinterval ≤ t.Lfdate) ?
      log(assetfloor(t.atq)) - log(assetfloor(t.Latq)) : missing) =>
    :inv)

  #annual inv data
  lagwithin2!(compq, [:Latq, :Lfdate], :gvkey, date=:fdate)
  lagwithin2!(compq, [:LLatq, :LLfdate], :gvkey, date=:fdate)
  lagwithin2!(compq, [:LLLatq, :LLLfdate], :gvkey, date=:fdate)
  select!(compq, Not([:LLatq, :LLfdate, :LLLatq, :LLLfdate]))
  rename!(compq, :LLLLatq=>:L4atq, :LLLLfdate=>:L4fdate)
  transform!(compq, [:atq, :L4atq, :fdate, :L4fdate] |> AsTable =>
    ByRow((t)-> (!ismissing(t.L4fdate)) && (t.fdate - maxannualfdateinterval ≤ t.L4fdate) ?
      log(assetfloor(t.atq)) - log(assetfloor(t.L4atq)) : missing) =>
    :inv4)

  #transform!(compq,[:pi, :revtq, :cogsq, :xsgaq, :xrdq, :dpq]

  ####checks on metrics
  mapprox(x,y) = x ≈ y
  mapprox(::Missing, ::Missing) = true
  mapprox(::Any, ::Missing) = false
  mapprox(::Missing, ::Any) = false

  #apply all three defintions of book equity at once
  compa.bkequity = (((ceqq, seqq, atq, ltq, pstkrq, pstkq)->
    coalesce(ceqq, seqq-coalesce(pstkrq, pstkq,0.0), atq-ltq-coalesce(pstkrq, pstkq,0.0))).(
    compa.ceqq,  compa.seqq, compa.atq, compa.ltq, compa.pstkrq, compa.pstkq) .+
    ((txditcq::MFloat64, txdbq::MFloat64)-> coalesce(txditcq, txdbq, 0.0)).(
    compa.txditcq,compa.txdbq)) .|> (x)->max(x,bookequityfloor)
  @assert mapprox.(compq.bkequity, compa.bkequity) |> all

  #apply the 4 defintions of op
  compa.oppi = ((piq, dpq)->(piq + coalesce(dpq,0.0))).(
    compa.piq, compa.dpq)
  compa.opni = ((ibq, niq, txtq, dpq)->(coalesce(ibq, niq) + txtq + coalesce(dpq,0.0))).(
    compa.ibq, compa.niq, compa.txtq, compa.dpq)
  compa.oprev = ((revtq, cogsq, xsgaq, xintq)->(revtq - cogsq - xsgaq - coalesce(xintq,0.0))).(
    compa.revtq, compa.cogsq, compa.xsgaq, compa.xintq)
  compa.opebit = ((oiadpq, xintq, dpq)->(oiadpq  + coalesce(dpq,0.0) - coalesce(xintq,0.0))).(
    compa.oiadpq, compa.xintq, compa.dpq)
  #use compq for bkequity in the below since we already tested bkequity
  compa.op = coalesce.(compa.oppi, compa.opni, compa.oprev, compa.opebit) ./ compq.bkequity
  @assert mapprox.(compq.op, compa.op) |> all

  compa.atqadj = ((x)->max(x,bookequityfloor)).(compa.atq)
  lagwithin2!(compa, [:atqadj], :gvkey, date=:fdate)
  @assert view(compa, (compa.Latqadj .!== missing) .& (compa.Lfdate .=== missing),:) |> isempty
  compa.Latqadj[((fdate, Lfdate, Latqadj)->
    (Latqadj !== missing) && (fdate - maxfdateinterval > Lfdate)
    ).(compa.fdate, compa.Lfdate, compa.Latqadj)] .= missing
  compa.inv = (log.(compa.atqadj) .- log.(compa.Latqadj)) .|> finiteormissing
  #@info "Num different: $(sum((!).(mapprox.(compq.inv, compa.inv))))"
  @assert mapprox.(compq.inv, compa.inv) |> all

  ####end checks
  filter!(r::DataFrameRow->year(r.fdate) ∈ yearrange, compq) #restrict to desired years

  #can't work with records with missing adates
  @info "Deleting $(sum(compq.adate .=== missing)) records with missing adates"
  compq = compq[compq.adate .!== missing, :]
  #drop most shifted dates since they may no longer be valid or acurate
  #keep lagged fdate varient since it timestamps the lagged assets for the investment factor
  #but rename it

  select!(compq, Not([:Nfdate]))
  rename!(compq, :Lfdate=>:Linvfdate, :L4fdate=>:L4invfdate)

  #println(describe(compq))
  #error("stopping here")

  return compq
end


parsecomp(::Type{Date}, s::String, compdateformat::DateFormat=dateformat"yyyymmdd"
  ) = Dates.Date(s, compdateformat)

#helper methods and special cases
parsecomp(::Type{<:Any}, v::Missing) = missing
parsecomp(::Type{Date}, i::Int, compdateformat::DateFormat=dateformat"yyyymmdd") = parsecomp(
  Date, "$i",compdateformat)
parsecomp(::Type{Symbol}, s::String) = Symbol(s)

#parse a value where the type is known: Generic case
function parsecomp(::Type{T}, s::String;
    ##NOTE: uncomment when we have the missings parsedmissings = COMP_PARSED_MISSINGS
    ) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  #NOTE: uncomment when we have the missings
  #(!ismissing(v)) && (v ∈ parsedmissings) && (v=missing)

  return v
end


function parsecomp(::Type{YearQuarter}, s::String)
  local yq::MYearQuarter

  if length(s) ≠ 6
    yq = missing
  else
    y::MInt = parsecomp(Int, "$(s[1:4])")
    q::MInt = parsecomp(Int, "$(s[6])")
    (q==1) || (q==2) || (q==3) || (q==4) || (q=missing)

    yq = (!ismissing(y)) && (!ismissing(q)) ? YearQuarter(y,q) : missing
  end


  return yq
end
