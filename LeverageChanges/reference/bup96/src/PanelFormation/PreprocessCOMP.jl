

function prepcomp(;companame=COMP_A_NAME, compqname = COMP_Q_NAME, compname = COMP_NAME,
   datapath=DATA_PATH, refreshcomp=true,
   incsvstream::Function = IN_CSV_STREAM,
   injlsstream::Function = IN_JLS_STREAM,
   outjlsstream::Function = OUT_JLS_STREAM,
   usequarterlydates::Bool=USE_QUARTERLY_DATES,
   parallel::NBool=nothing,
   retainedcolscomp::Vector{Symbol} = RETAINED_COLS_COMP,
   archive::Bool = true)::DataFrame

   (!isnothing(parallel)) && (PARALLEL[] = parallel)

  local compa::DataFrame
  local compq::DataFrame
  local comp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcomp || (!isfile("$datapath\\$compname.jls.lz4"))
    compa = incsvstream("$datapath\\$companame.csv.lz4") |> CSV.read |> DataFrame
    compq = incsvstream("$datapath\\$compqname.csv.lz4") |> CSV.read |> DataFrame

    comp = preprocesscomp(compa, compq, usequarterlydates=usequarterlydates)
    comp = mergecompccm!(comp)
    comp = select!(comp,retainedcolscomp)
    outjlsstream("$datapath\\$compname.jls.lz4") do s
      serialize(s, comp) #NOTE: DEBUG ONLY
    end

    archive && archivefile("$compname.jls.lz4") #never hurts to have a backup

  else
    try
      injlsstream("$datapath\\$compname.jls.lz4") do s
        comp = deserialize(s) #NOTE: DEBUG ONLY
      end
    catch err
      @warn "Could not load $compname. Error: $err"
      @info "Attempting to restore back-up and load from archive..."
      unarchivefile("$compname.jls.lz4")
      injlsstream("$datapath\\$compname.jls.lz4") do s
        comp = deserialize(s) #NOTE: DEBUG ONLY
      end
      @info "Successfully restored from archive."
    end


  end

  @info "Compustat and CCM linking data loaded and saved into file $compname.jls.lz4"

  return comp
end


#helper function to clean up names in the quarterly file
function propername(s::Symbol)
  str::String = string(s)

  if str[end] == 'q' || (length(str) ≥ 3 && str[(end-2):end] == "qtr")
    return s
  end

  s = Symbol(s,:q)
  return s
end

#this processes both the q and the a comp frames, with as much commonly used code as possible
function preprocesscomp(compa::DataFrame, compq::DataFrame;
  compdateformat::DateFormat = COMP_DATE_FORMAT,
  testoutput::Bool = TEST_OUTPUT,
  usequarterlydates::Bool = true)

  local comp::DataFrame
  #local usedcompfields::Vector{Symbol}


  #do operations that are common to both the annual and the quarterly data
  compa = preprocesscomp!(compa) #methods common to compa and compq
  compa = preprocesscompa!(compa) #methods specific to compa

  if usequarterlydates #need to merge this in if we are using it
    compq = preprocesscomp!(compq)
    compq = preprocesscompq!(compq)

    #(if we were using quarterly data it would be joined in here
    comp = join(compa, compq, on=[(:gvkey, :gvkeyq),
      (:fyear, :fyearq), (:indfmt,:indfmtq)], kind=:left)

    #the below performs a lookup and makes sure each row has an
    #actual quarterly announcement date if available
    gdf::GroupedDataFrame = groupby(comp, [:gvkey, :fyear])
    @mpar for i ∈ 1:length(gdf)
      scomp::SubDataFrame = gdf[i]

      for r ∈ eachrow(scomp)
        local adateq::MDate = r.adateq
        (!ismissing(adateq)) && (r.adate = adateq)
        if (size(scomp,1) > 1) && (!ismissing(adateq))
          scomp.adate .= ((ad::MDate,adq::MDate)->
            (ismissing(adq)) ? adateq : ad).(scomp.adate, scomp.adateq)
        end
      end
    end

    comp.adate = ((ad::MDate,adq::MDate)->
      ismissing(adq) ? ad : adq).(comp.adate, comp.adateq)

      #usedcompfields = [USED_COMPA_VALUE_FIELDS; USED_COMPQ_VALUE_FIELDS]
  else
    comp = compa
    #usedcompfields = USED_COMPA_VALUE_FIELDS
  end

  println("rows: $(size(comp,1))")
  filter!(comp) do r::DataFrameRow
      r[:adate] ≥ eom(r[:fyrmonth])
  end

  println("rows adate≥fdate: $(size(comp,1))")
  validatefadates!(comp; Fgroup=:gvkey, Ffdate=:fyear, Fadate=:adate)

  #form the initial interval candidates
  #these are accurate if the link doesn't changes (adjust for this case later)
  comp = formannouncementintervals!(comp)

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

  rename!(comp, :indfmt=>:indfmtold)
  comp.indfmt = (s::MString->parsecomp(Symbol, s)).(comp.indfmtold)
  select(comp, Not(:indfmtold))

  return comp
end



function preprocesscompa!(compa::DataFrame;
  usedcompafields::Vector{Symbol} = USED_COMPA_VALUE_FIELDS)


  #try to use the more efficient year-month format when possible
  compa = compa[(!ismissing).(compa.fyear), :]
  compa.fyrmonth = ((y::Integer,m::Integer)->
    YearMonth(y,m)).(compa.fyear, compa.fyr)
  compa.fdate = (eom).(compa.fyrmonth)


  #sets the fixed lead of each entry
  @inline function ayrmonth(ym::YearMonth)::Date
    Δ::Int = 4 + (ym.y < 2004)
    return eom(ym + Δ)
  end

  compa.adate = (ayrmonth).(compa.fyrmonth)

  compa = dedupcomp!(compa, [:gvkey, :fyear], usedcompafields, :indfmt, :adate)

  #make the cash field
  compa.cash = ((che::MFloat64, ch::MFloat64)-> coalesce(che, ch)).(compa.che, compa.ch)
  compa = compa[(!ismissing).(compa.cash), :] #this is an essential field
  compa.bknegcash = 1.0 .- (compa.cash ./ compa.at)

  #compute #flev field
  #no negative values
  @inline function maxval0!(df::AbstractDataFrame, s::Symbol)
    df[!,s] .= (f::MFloat64->max(f,0.0)).(compa[!,s])
  end

  @mpar for s ∈ [:dltt, :dlc, :ceq, :txditc, :txdb, :pstk]
    maxval0!(compa, s)
  end

  compa.netdebt = compa.dltt .+ compa.dlc .- compa.cash
  compa.bkequity = compa.ceq .- compa.pstk .- ((txditc::MFloat64, txdb::MFloat64)->
    coalesce(txditc, txdb, 0.0)).(compa.txditc,compa.txdb)
  @mpar for s ∈ [:netdebt, :bkequity]
    maxval0!(compa, s)
  end
  compa.bkflev = compa.netdebt ./ (compa.netdebt .+ compa.bkequity)

  #compute liab
  compa.bkliab = compa.lt ./ compa.at

  #gross profit is proxy for operating profits as in Novy Marx 2013
  compa.op = (r::DataFrameRow->ismissing(r.gp) ? r.revt - r.cogs : r.gp).(eachrow(compa))

  #compa = compa[completecases(compa[:, [:bknegcash, :bkflev, :bkliab]]),:]

  #make finite
  @mpar for s ∈ [:bknegcash, :bkflev, :bkliab]
    compa[!,s] = (f::MFloat64->coalesce(isfinite(f),false) ? f : missing).(compa[!,s])
  end

  #get the market value if available
  compa.mkequity = ((mkvalt::MFloat64, pso::MFloat64)->
    coalesce(mkvalt, pso)).(compa.mkvalt, compa.prcc_f .* compa.csho)

  #NOTE: we won't translate these into our measures until after we do the primary merge

  #sort the data
  sort!(compa, [:gvkey, :adate])


  #println(describe(compa))

  return compa
end


function preprocesscompq!(compq::DataFrame;
  usedcompqfields::Vector{Symbol} = USED_COMPQ_VALUE_FIELDS,
  month2stale::Int = MONTHS_2_STALE_SHORT)
  #use a consistent naming convention
  names!(compq, (propername).(names(compq)))

  #helper function to convert a fiscal year-quarter to a datw
  @inline function fyearqtr2date(yq::YearQuarter, fyr::Int)
    ym::YearMonth = YearMonth(yq.y, fyr)
    offsetmonths::Int = -(4 - yq.q)*3
    ym += offsetmonths
    return eom(ym)
  end

  #pre-proess date info
  compq = compq[((fq)->(!ismissing(fq))).(compq.datafqtr),:]
  compq.fyearqtr =  (s->parsecomp(YearQuarter, s)).(compq.datafqtr)
  #compq.fyearqtr =  (s->parsecomp(Symbol, s)).(compq.indfmtq)
  compq.fdateq = ((yq::YearQuarter,fyr::Int)->
    fyearqtr2date(yq,fyr)).(compq.fyearqtr, compq.fyrq)

  #WARNING: We drop all quarters except Q4 here since we only want annnual data
  compq = compq[(yq::YearQuarter->yq.q==4).(compq.fyearqtr),:]

  compq.adateq = (s->parsecomp(Date, s)).(compq.rdq)
  compq = compq[(!ismissing).(compq.adateq),:]
  compq = compq[((yq::YearQuarter, fyearq::Int)->yq.y == fyearq).(compq.fyearqtr, compq.fyearq),:]


  #NOTE: REMINDER- convert to first day of month for date math
  #Date math works better using the first day of the month,
  #even though the last day of the month is technically correct

  #don't let the fiscal date be too far from the announcement date
  period2stale::Month = Month(month2stale)
  compq = compq[(r->(r.adateq - period2stale ≤ r.fdateq)).(eachrow(compq)),:]



  validatefadates!(compq; Fgroup=:gvkeyq, Ffdate=:fyearqtr, Fadate=:adateq)

  compq = dedupcomp!(compq, [:gvkeyq, :fyearq], usedcompqfields, :indfmt, :adateq)
  compq = dedupcomp!(compq, [:gvkeyq, :adateq], usedcompqfields, :indfmt, :fyearq)
  #now process cases where firms have duplicate date-firm key pairs
  #the date field is used to affirm the latest record
  #(for annual records, that is the annual date,
  #for quarterly its the announcement date)
  #compq = dedupcomp!(compq, [:gvkeyq, :adateq], usedcompqfields, :indfmtq, :fyearq)


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

  select!(comp, Not(:maxfdate_at_adate))

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
      #keep the most recent version
      maxdate::Union{Int, YearQuarter, Date} = maximum(scomp[!,datefield])
      mindate::Union{Int, YearQuarter, Date} = minimum(scomp[!,datefield])
      #println("got here")

      #why could this happen? not clear to me, but it does happen a bit
      if maxdate ≠ mindate
        scomp.keep .= false
      end
      #println("got here")
      #scomp.keep .= (yd->isequal(yd,maxdate)).(scomp[!,datefield])

      #special case where we have two records of two different formats (indl and fs)
      #in this case, keep the format with the least missing values
      if sum(scomp.keep) > 1
        sscomp::SubDataFrame = view(scomp, scomp.keep, :)
        Nsscomp::Int = size(sscomp, 1)

        for k ∈ 1:Nsscomp
          sscomp.valuesmissings[k] = sum((c->ismissing(sscomp[k,c])).(usedcompfields))
        end

        minmissing::Int = minimum(sscomp.valuesmissings)
        sscomp.keep .= (sscomp.valuesmissings .== minmissing)

        #tie breaker is to keep the INDL record (industrial since we are looking at non-financials)
        if sum(sscomp.keep) > 1
          sscomp.keep .= ((keep::MBool, form::MSymbol)->
            (!ismissing(form)) && keep && (form == :INDL)).(sscomp.keep, sscomp[!, formatfield])
        end

      end

      #if we still have a dup, delete all the records
      if sum(scomp.keep) > 1
        @warn "dedup failed group fields:
          $groupbyfields with values $((f->scomp[1,f]).(groupbyfields)) Records removed.\n
          $scomp"
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
    lastdate::Date = LAST_DATE, months2stale::Int =  MONTHS_2_STALE_LONG)

  sort!(comp, [:gvkey, :adate])

  comp.Nadate = Vector{MDate}(undef, size(comp,1))

  scomps = groupby(comp, :gvkey)
  @mpar for i ∈ 1:length(scomps)
    scomp::SubDataFrame = scomps[i]
    Nsub::Int = size(scomp,1)

    for i ∈ 1:(Nsub-1)
      scomp[i,:Nadate] = scomp[i+1,:adate]

      #check some things that shouldn't happen
      (scomp[i,:adate] ≥ scomp[i,:Nadate]) && @warn (
        "announcement dates out of order.\n $(scomp[i:(i+1),:])")

      (scomp[i,:fyear] ≥ scomp[i+1,:fyear]) && @warn (
        "fiscal year dates out of order.\n $(scomp[i:(i+1),:])")
    end

    #final vliad interval extends to the present
    scomp[Nsub, :Nadate] = lastdate
  end


  comp.fdatemax = (eom).(comp.fyrmonth .+ months2stale)

  comp.begindate = comp.adate
  comp.enddate = (min).(comp.Nadate,comp.fdatemax)

  return comp
end


parsecomp(::Type{Date}, s::String;
    compdateformat::DateFormat = COMP_DATE_FORMAT) = Dates.Date(s, compdateformat)

#helper methods and special cases
parsecomp(::Type{<:Any}, v::Missing) = missing
parsecomp(::Type{Date}, i::Int) = parsecomp(Date, "$i")
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
