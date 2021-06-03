
prepcomp(;refreshcomp=true) = loadcompdata(refreshcomp=refreshcomp)

function loadcompdata(;companame=COMP_A_NAME, compqname = COMP_Q_NAME, compname = COMP_NAME,
   datapath=DATA_PATH, refreshcomp=true,
   incsvstream::Function = IN_CSV_STREAM,
   injlsstream::Function = IN_JLS_STREAM,
   outjlsstream::Function = OUT_JLS_STREAM)::DataFrame

  local compa::DataFrame
  local compq::DataFrame
  local comp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcomp || (!isfile("$datapath\\$compname.jls.lz4"))
    compa = incsvstream("$datapath\\$companame.csv.lz4") |> CSV.read |> DataFrame
    compq = incsvstream("$datapath\\$compqname.csv.lz4") |> CSV.read |> DataFrame

    comp = preprocesscomp(compa, compq)
    comp = mergecompccm!(comp)

    outjlsstream("$datapath\\$compname.jls.lz4") do s
      serialize(s, comp) #NOTE: DEBUG ONLY
    end

    #NOTE: FINAL VERSION
    #preprocesscomp(compa, compq)
    #serialize("$datapath\\$compname.jls", comp)
  else
    injlsstream("$datapath\\$compname.jls.lz4") do s
      comp = deserialize(s) #NOTE: DEBUG ONLY
    end
  end

  @info "Compustat and CCM linking data loaded and saved into file $compname.jls.lz4"

  return comp
end
#helper function to clean up names in the quarterly file
function properqname(s::Symbol)
  str::String = string(s)

  if str[end] == 'q' || (length(str) ≥ 3 && str[(end-2):end] == "qtr")
    return s
  end

  s = Symbol(s,:q)
  return s
end

#helper function to convert a fiscal year-quarter to a date
@inline function fyearqtr2date(yq::YearQuarter, fyr::Int)
  ym::YearMonth = YearMonth(yq.y, fyr)
  offsetmonths::Int = -(4 - yq.q)*3
  ym += offsetmonths
  return eom(ym)
end


#this processes both the q and the a comp frames, with as much commonly used code as possible
function preprocesscomp(compa::DataFrame, compq::DataFrame;
  compdateformat::DateFormat = COMP_DATE_FORMAT,
  retainedcolscomp::Vector{Symbol} = RETAINED_COLS_COMP,
  testoutput::Bool = TEST_OUTPUT)

  #do operations that are common to both the annual and the quarterly data
  compa = preprocesscomp!(compa)
  compq = preprocesscomp!(compq)

  compa = preprocesscompa!(compa)
  compq = preprocesscompq!(compq)


  comp::DataFrame = join(
    compa, compq, on=[:gvkey=>:gvkeyq, :fyear=>:fyearq], kind=:inner)

  #NOTE: The below is not actually done in the paper, but I think its important
  #(The paper only uses the annual data)
  comp.iaq = ((iaq,ia)->ismissing(iaq) ? ia : iaq).(comp.iaq, comp.ia)

  #form the initial interval candidates
  #these are accurate if the link doesn't changes (adjust for this case later)
  comp = formannouncementintervals!(comp)
  comp = select!(comp,retainedcolscomp)

  #println(describe(comp))
  if testoutput
    d = Date(2015,12,31)
    testoutput && CSV.write("output\\comp_$d.csv", comp[((bd,ed)->
      (d > bd) && (d ≤ ed)).(comp.begindate, comp.enddate), :])
  end


  return comp
end

#operations that are common to both the annual and the quarterly data
function preprocesscomp!(comp::DataFrame)
  #make all column names lower case
  names!(comp, (s::Symbol->Symbol(lowercase(string(s)))).(names(comp)))

  #format dates
  comp.ddate = (s->parsecomp(Date, s)).(comp.datadate)
  comp.dyear = (d::MDate-> ismissing(d) ? missing : Dates.year(d)).(comp.ddate)

  comp = comp[(!ismissing).(comp.dyear), :] #require the calendar date


  return comp
end


function preprocesscompa!(compa::DataFrame;
  usedcompafields::Vector{Symbol} = USED_COMPA_VALUE_FIELDS)
  compa = compa[(!ismissing).(compa.fyear), :]


  #NOTE: REMINDER- convert to first day of month for date math
  #Date math works better using the first day of the month,
  #even though the last day of the month is technically correct
  compa.fdate = ((y::Int,m::Int)->
    lastdayofmonth(Date(y,m,1))).(compa.fyear, compa.fyr)

  #delete rows with missing data
  compa = compa[(!ismissing).(compa.at), :]

  #now process cases where firms have duplicate date-firm key pairs
  #the date field is used to affirm the latest record
  #(for annual records, that is the annual date,
  #for quarterly its the announcement date)
  compa = dedupcomp!(compa, [:gvkey, :fyear], usedcompafields, :indfmt, :ddate)

  #sort the data
  sort!(compa, [:gvkey, :fyear])

  ####now make the I/A amount, both in a quarterly and an annual form
  compa.ia = Vector{MFloat64}(undef, size(compa,1))

  lagwithin!(compa, :at, :gvkey, :fyear, sorted=true)
  compa.ia .= (compa.at .- compa.Lat) ./ compa.Lat
  compa.ia .= (f::MFloat64 -> (!ismissing(f)) && isfinite(f) ? f : missing).(compa.ia)

  return compa
end

function preprocesscompq!(compq::DataFrame;
  usedcompqfields::Vector{Symbol} = USED_COMPQ_VALUE_FIELDS)
  #use a consistent naming convention
  names!(compq, (properqname).(names(compq)))

  compq = compq[((cq, fq)->(!ismissing(cq)) && (!ismissing(fq))).(
    compq.datacqtr, compq.datafqtr),:]
  compq.fyearqtr =  (s->parsecomp(YearQuarter, s)).(compq.datafqtr)
  compq.cyearqtr =  (s->parsecomp(YearQuarter, s)).(compq.datacqtr)
  compq.adateq = (s->parsecomp(Date, s)).(compq.rdq)
  compq = compq[(!ismissing).(compq.adateq),:]


  #NOTE: REMINDER- convert to first day of month for date math
  #Date math works better using the first day of the month,
  #even though the last day of the month is technically correct
  compq.fdateq = ((yq::YearQuarter,fyr::Int)->
    fyearqtr2date(yq,fyr)).(compq.fyearqtr, compq.fyrq)

  #don't let the fiscal date be too far from the data date (uncomment to use data date)
  compq = compq[(r->(r.adateq - Month(6) ≤ r.fdateq)).(eachrow(compq)),:]

  #under no circumstances allow restatements past the subsequent quarter's announcement date
  validatefadates!(compq; Fgroup=:gvkeyq, Ffdate=:fyearqtr,
    Fadate=:adateq)

  #now process cases where firms have duplicate date-firm key pairs
  #the date field is used to affirm the latest record
  #(for annual records, that is the annual date,
  #for quarterly its the announcement date)
  #Require a 1:1 relationship between announcement dates and year-quarters
  #Mainly want to avoid interm information disclosures
  #NOTE: could revisit in robustness checks
  compq = dedupcomp!(compq, [:gvkeyq, :fyearqtr], usedcompqfields, :indfmtq, :adateq)
  compq = dedupcomp!(compq, [:gvkeyq, :adateq], usedcompqfields, :indfmtq, :fyearqtr)

  #sort the data
  sort!(compq, [:gvkeyq, :fyearqtr])

  ####now make the I/A amount, both in a quarterly and an annual form
  compq.iaq = Vector{MFloat64}(undef, size(compq,1))

  lagwithin!(compq, :atq, :gvkeyq, :fyearqtr, sorted=true, lags=4)
  compq.iaq .= (compq.atq .- compq.L4atq) ./ compq.L4atq
  compq.iaq .= (f::MFloat64 -> (!ismissing(f)) && isfinite(f) ? f : missing).(compq.iaq)

  ###derive ROE
  #first get shareholders equity
  compq.seq = (r::DataFrameRow-> coalesce(r.seqq, r.ceqq+r.pstkq, r.atq - r.ltq)).(eachrow(compq))
  #now add in deferred taxes and inv credits to get book equity, and subtract out preferred stock
  compq.beq = ((r::DataFrameRow)->
    r.seq + coalesce(r.txditcq,0.0) - coalesce(r.pstkrq, r.pstkq,0.0)).(eachrow(compq))
  lagwithin!(compq, :beq, :gvkeyq, :fyearqtr, sorted=true)
  compq.roeq = compq.ibq ./ compq.Lbeq
  compq.roeq .= (f::MFloat64->(!ismissing(f)) && isfinite(f) ? f : missing).(compq.roeq)

  return compq
end

  #under no circumstances allow restatements past the subsequent quarter's announcement date
function validatefadates!(comp::DataFrame; Fgroup::NSymbol=nothing, Ffdate::NSymbol=nothing,
  Fadate::NSymbol=nothing, sort::Bool = true)::Nothing

  (isnothing(Fgroup)) && error("Fgroup missing from validatefadates")
  (isnothing(Ffdate)) && error("Ffdate missing from validatefadates")
  (isnothing(Fadate)) && error("Fadate missing from validatefadates")

  sort!(comp, [Fgroup, Fadate])
  T::Type = eltype(comp[!,Ffdate])
  comp.maxfdate_at_adate = Vector{Union{T,Missing}}(undef, size(comp,1))

  compgrp::GroupedDataFrame = groupby(comp, Fgroup)
  @mpar for i ∈ 1:length(compgrp)
    subdf::SubDataFrame = compgrp[i]
    maxfdate::T = subdf[1,Ffdate]
    for r ∈ eachrow(subdf)
      (r[Ffdate] > maxfdate) &&  (maxfdate = r[Ffdate])
      r.maxfdate_at_adate = maxfdate
    end
  end

  filter!(comp) do r::DataFrameRow
    r.maxfdate_at_adate ≤ r[Ffdate]
  end

  return nothing
end


function dedupcomp!(comp::DataFrame,
  groupbyfields::Vector{Symbol}, usedcompfields::Vector{Symbol},
  formatfield::Symbol, datefield::Symbol)

  #keepflags
  comp.keep = trues(size(comp,1))

  #use this for deduping
  comp.valuesmissings = Vector{MFloat64}(undef, size(comp,1))

  ncols::Int = size(comp,2)
  scomps::GroupedDataFrame = groupby(comp, groupbyfields)
  @mpar for i ∈ 1:length(scomps)
    scomp::SubDataFrame = scomps[i]

    if size(scomp,1) > 1

      #keep the most recent version
      maxdate::Union{Int, YearQuarter, Date} = maximum(scomp[!,datefield])
      mindate::Union{Int, YearQuarter, Date} = minimum(scomp[!,datefield])
      #println("got here")

      #why could this happen? not clear to me, but it does happen a bit
      if maxdate ≠ mindate
        scomp.keep .= false
      end

      #special case where we have two records of two different formats (indl and fs)
      #in this case, keep the format with the least missing values
      if sum(scomp.keep) == 2
        sscomp::SubDataFrame = view(scomp, scomp.keep, :)
        if (length(unique(sscomp.indfmt)) == 2)

          sscomp.valuesmissings[1] = sum((c->ismissing(sscomp[1,c])).(usedcompfields))
          sscomp.valuesmissings[2] = sum((c->ismissing(sscomp[2,c])).(usedcompfields))

          sscomp.keep .= (sscomp.valuesmissings .== minimum(sscomp.valuesmissings))

          #tie breaker is to keep the INDL record (industrial since we are looking at non-financials)
          if sum(sscomp.keep) > 1
            sscomp.keep .= (sscomp[!, formatfield] .== :INDL)
          end
        end
      end

      #if we still have a dup, delete all the records
      if sum(scomp.keep) > 1
        @warn "dedup failed group fields:
          $groupbyfields with values $((f->scomp[1,f]).(groupbyfields)) Records removed."
        scomp.keep .= false
      end
    end
  end

  comp = comp[comp.keep, :] #bandaid until filter is fixed


  select!(comp, Not([:keep, :valuesmissings]))

  return comp
end


#NOTE: messy function
#basically to merge in crsp, the crsp day or month must occur
#in the intersection of the link dates, the period six months following
#the announcement date, and the next announcement date
#this handles everything but the intersection with the link dates
function formannouncementintervals!(comp::DataFrame;
    lastdate::Date = LAST_DATE, time2stale::T = TIME_2_STALE) where T<:DatePeriod

  sort!(comp, [:gvkey, :adateq])

  comp.Nadateq = Vector{MDate}(undef, size(comp,1))

  scomps = groupby(comp, :gvkey)
  @mpar for i ∈ 1:length(scomps)
    scomp::SubDataFrame = scomps[i]
    Nsub::Int = size(scomp,1)

    for i ∈ 1:(Nsub-1)
      scomp[i,:Nadateq] = scomp[i+1,:adateq]

      #check some things that shouldn't happen
      (scomp[i,:adateq] ≥ scomp[i,:Nadateq]) && @warn (
        "announcement dates out of order.\n $(scomp[i:(i+1),:])")

      (scomp[i,:fyearqtr] ≥ scomp[i+1,:fyearqtr]) && @warn (
        "fiscal year dates out of order.\n $(scomp[i:(i+1),:])")
    end

    scomp[Nsub, :Nadateq] = lastdate
  end


  comp.fdateqmax = (lastdayofmonth).(
    (firstdayofmonth).(comp.fdateq) .+ time2stale)

  comp.begindate = comp.adateq
  comp.enddate = (min).(comp.Nadateq,comp.fdateqmax)

  return comp
end

parsecomp(::Type{Date}, s::String;
    compdateformat::DateFormat = COMP_DATE_FORMAT) = Dates.Date(s, compdateformat)

#helper methods and special cases
parsecomp(::Type{<:Any}, v::Missing) = missing
parsecomp(::Type{Date}, i::Int) = parsecomp(Date, "$i")

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
