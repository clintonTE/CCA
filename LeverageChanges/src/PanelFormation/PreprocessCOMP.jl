

function prepcomp(;companame=COMP_A_NAME, compqname = COMP_Q_NAME, compname = COMP_NAME,
   datapath=DATA_PATH, refreshcomp=true,
   incsvstream::Function = IN_CSV_STREAM,
   inbinstream::Function = IN_BIN_STREAM,
   outbinstream::Function = OUT_BIN_STREAM,
   csvextension::String = CSV_EXTENSION,
   binextension::String = BIN_EXTENSION,
   usequarterlydates::Bool=USE_QUARTERLY_DATES,
   parallel::NBool=nothing,
   retainedcolscomp::Vector{Symbol} = RETAINED_COLS_COMP,
   archive::Bool = true)::DataFrame

   (!isnothing(parallel)) && (PARALLEL[] = parallel)

  local compa::DataFrame
  local comp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcomp
    compa = incsvstream("$datapath\\$companame.$csvextension") |> CSV.File |> DataFrame
    usequarterlydates && (compq = incsvstream("$datapath\\$compqname.$csvextension") |> CSV.File |> DataFrame)

    comp = preprocesscomp(compa)
    comp = mergecompccm!(comp)
    comp = select!(comp,retainedcolscomp)
    outbinstream("$datapath\\$compname.$binextension", comp)

  else
    comp = inbinstream("$datapath\\$compname.$binextension")

  end


  @info "Compustat and CCM linking data loaded and saved into file $compname.$binextension"

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
function preprocesscomp(comp::DataFrame;
  compdateformat::DateFormat = COMP_DATE_FORMAT,
  testoutput::Bool = TEST_OUTPUT)

  #do operations that are common to both the annual and the quarterly data
  comp = preprocesscomp!(comp) #methods common to compa and compq
  comp = preprocesscompa!(comp) #methods specific to compa


  println("rows: $(size(comp,1))")
  filter!(comp) do r::DataFrameRow
      r[:adate] ≥ eom(r[:fyrmonth])
  end

  println("rows adate≥fdate: $(size(comp,1))")
  validatefadates!(comp; Fgroup=:gvkey, Ffdate=:fdate, Fadate=:adate)
  sort!(comp, [:gvkey, :fdate])

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
  rename!(comp, (s::String->lowercase(s)).(names(comp)))

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
  usedcompafields::Vector{Symbol} = USED_COMPA_VALUE_FIELDS,
  yearrange::UnitRange{Int} = YEAR_RANGE[],
  minannouncementdelay::Int = MIN_ANNOUNCEMENT_DELAY)

  #println(describe(compa))
  #try to use the more efficient year-month format when possible
  compa = compa[(!ismissing).(compa.fyear), :]

  compa.fic = (s::MString->parsecomp(Symbol,s)).(compa.fic)
  compa = compa[compa.fic.==:USA, :]
  compa.fyrmonth = ((y::Integer,m::Integer)->
    YearMonth(y,m)).(compa.fyear, compa.fyr)
  compa.fdate = (eom).(compa.fyrmonth)

  filter!(r::DataFrameRow->year(r.fdate) ∈ yearrange, compa) #restrict to desired years


  #sets the fixed lead of each entry
  @inline function ayrmonth(ym::YearMonth)::YearMonth
    Δ::Int = minannouncementdelay + (ym.y < 2004)
    return ym + Δ#eom(ym + Δ)
  end

  sort!(compa, [:gvkey, :fdate])
  compa.ayrmonth = (ayrmonth).(compa.fyrmonth)
  compa.adate = (eom).(compa.ayrmonth)
  compa = dedupcomp(compa, [:gvkey, :fyear], :fyrmonth)
  #compa = dedupcomp(compa, [:gvkey, :fyear], :fyrmonth) #more of a validation check

  #make the cash field
  compa.cash = ((che::MFloat64, ch::MFloat64)-> coalesce(che, ch, 0.0)).(compa.che, compa.ch)
  compa.cash .= (floororx).(compa.cash)
  compa = compa[(!ismissing).(compa.cash), :] #this is an essential field #should this line go away?
  compa.cashat = (boundunit1).(compa.cash ./ compa.at)
  compa.bknegcash = (boundunit).(1.0 .- compa.cashat)


  compa.debt = (floororx).(compa.dltt .+ compa.dlc)
  compa.netdebt = compa.debt .- compa.cash
  compa.bkequity = (compa.ceq .+
    ((txditc::MFloat64, txdb::MFloat64)-> coalesce(txditc, txdb, 0.0)).(compa.txditc,compa.txdb) #.-
    #(pstk->coalesce(pstk,0.0)).(compa.pstk) #WARNING- commented out since excluded from ceq
    #this is different than how Ivo has it in the original
    )
  compa.bkequityold = (compa.ceq .+
    ((txditc::MFloat64, txdb::MFloat64)-> coalesce(txditc, txdb, 0.0)).(compa.txditc,compa.txdb) .-
    (pstk->coalesce(pstk,0.0)).(compa.pstk)
    )
  #compa.bkequity .= compa.bkequityold #WARNING- only enable for testing purposes
  @mpar for s ∈ [:netdebt, :bkequity, :bkequityold, :at, :lt] #can thread this a bit
    compa[!,s] .= (floororx).(compa[!,s])
  end

  #put in Ivo's minimum values
  compa = filter(r::DataFrameRow->
    (!ismissing(r.at)) &&
    (!ismissing(r.bkequity)) &&
    (r.at ≥ MIN_AT) &&
    (r.bkequity ≥ MIN_BKEQUITY), compa)

  compa.bkflev = (boundunit).(compa.netdebt ./ (compa.netdebt .+ compa.bkequity))

  #compute liab
  compa.netliab = compa.lt .- compa.cash
  compa.bkliab = (boundunit).(compa.netliab ./ compa.at)

  #gross profit is proxy for operating profits as in Novy Marx 2013
  #compa.op = (r::DataFrameRow->ismissing(r.gp) ? r.revt - r.cogs : r.gp).(eachrow(compa))
  compa.nop = (boundabsunit).(compa.pi ./ compa.bkequity)
  #compa.nop = compa.op ./ compa.bkequity #this is used for the operating profit ratio


  #make finite
  @mpar for s ∈ [:bknegcash, :bkflev, :bkliab, :nop]
    compa[!,s] = (f::MFloat64->coalesce(isfinite(f),false) ? f : missing).(compa[!,s])
  end

  #get the market value if available
  compa.mkequity = ((mkvalt::MFloat64, pso::MFloat64)->
    coalesce(mkvalt, pso)).(compa.mkvalt, compa.prcc_f .* compa.csho)

  #NOTE: we won't translate these into our measures until after we do the primary merge

  #sort the data
  sort!(compa, [:gvkey, :adate])
  issorted(compa, [:gvkey, :fdate]) || error("compa not sorted on fdate")
  lagwithin2!(compa, [:at], :gvkey, date=:fdate)

  compa.inv = compa.at ./ compa.Lat .- 1.

  #Now make the BRS metrics
  compa.Aat = (r::DataFrameRow->coalesce((r.at + r.Lat)/2, r.at, r.Lat)).(eachrow(compa))
  compa.brse = (compa.sstk .- compa.prstkc .- compa.dv) ./ compa.Aat
  compa.brsd = (r::DataFrameRow->(r.dltis - r.dltr - coalesce(r.dlcch,0.))/r.Aat).(eachrow(compa))
  compa.brsx = compa.brse .+ compa.brsd

  (f::Symbol->finiteormissing!(compa, f)).([:brse, :brsd, :brsx])




  #println(describe(compa))

  return compa
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

function dedupcomp(comp::DataFrame,
  groupbyfields::Vector{Symbol}, datefield::Symbol)

  comp = comp[comp.indfmt .== :INDL,:]

  issorted(comp, [groupbyfields; datefield]) || error("comp not sorted by $datefield")
  #keepflags
  #comp.keep = trues(size(comp,1))

  #this is to match ivo's code

  ncols::Int = size(comp,2)
  scomps::GroupedDataFrame = groupby(comp, groupbyfields)
  @mpar for i ∈ 1:length(scomps)
    scomp::SubDataFrame = scomps[i]
    (size(scomp,1) == 1) && continue #the most likely case

    #keep the most recent version
    maxdate::Union{Int, YearQuarter, Date, YearMonth} = maximum(scomp[!,datefield])
    mindate::Union{Int, YearQuarter, Date, YearMonth} = minimum(scomp[!,datefield])
    #println("got here")

    #Happens with change in fiscal year
    #=if maxdate ≠ mindate
      scomp.keep .= false
    end=#
    #println("got here")
    #scomp.keep .= (yd->isequal(yd,maxdate)).(scomp[!,datefield])

    #special case where we have two records of two different formats (indl and fs)
    #in this case, keep the format with the least missing values
    if (size(scomp,1) > 1) && (maxdate==mindate)
      error("duplicate rows after dropping fs records")
    end

    #if we still have a dup, delete all the records
    #=if sum(scomp.keep) > 1
      @warn "dedup failed group fields:
        $groupbyfields with values $((f->scomp[1,f]).(groupbyfields)) Records removed.\n
        $scomp"
      scomp.keep .= false
    end=#
  end

  #comp = comp[comp.keep, :] #bandaid until filter is fixed


  #select!(comp, Not([:keep, :valuesmissings]))

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
      (scomp[i,:adate] ≥ scomp[i,:Nadate]) && error(
        "announcement dates out of order.\n $(scomp[i:(i+1),:])")

      (scomp[i,:fyear] ≥ scomp[i+1,:fyear]) && error(
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
