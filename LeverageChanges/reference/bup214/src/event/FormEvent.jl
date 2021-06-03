


getdayfields(days::Int) = (d->d<0 ? Symbol("DM",abs(d)) : Symbol("D",d)).(collect((-1*days):days))

function formeventwindows(event::DataFrame, days::Int;
  daysextended::Int = days,
  refresheventcrsp::Bool = true,
  Fret::Symbol = error("Fret is required"),
  validpermnos::Vector{Int} = unique(event[!, :permno]),
  Fdate::Symbol = error("Fdate is required"),
  retbound::NFloat64 = nothing,
  prefix::String = error("prefix is required"),
  daystostaleprice::NInt=nothing,
  yearrange::UnitRange = error("yearrange (for crsp) is required")
  )

  days>0 || error("days must be greater than 0")
  Nevent::Int = size(event,1)
  event.eid = collect(1:Nevent)


  #load the crsp df
  local crsp::DataFrame = prepcrsp(
    crspcolumns = [:date, :permno, Fret, :price, :divamt, :distcd, :numtrd, :vol],
    crsptype = :dailyexpanded,
    refreshcrsp=refresheventcrsp,
    validpermnos = validpermnos,
    prefix = "$(prefix)event",
    yearrange=yearrange) do df::DataFrame
    computenetreturns!(df)
    computevolumes!(df)
    sort!(df,[:permno,:date])
    #df.Tprice = Vector{MFloat64}(missing, size(df,1))
    #(!isnothing(daystostaleprice)) && lagprice!(df)
  end

  #println("namescrsp", names(crsp))
  rename!(crsp, [:divamt=>:divamt2, :distcd=>:distcd2])

  #make an index of possible dates
  crspdates::Vector{Date} = sort(unique(crsp.date))
  Ndates::Int = length(crspdates)


  #the following interlude only applies to the distributions dataframe
  if :exdate ∈ names(event)
    event.days2lastdiv = Vector{MInt}(missing, size(event,1))

    #looks up the date against crspdates. If it finds the date,
    #returns the index, otherwise returns missing
    @inline function crspdateormissing(d::Date)::MInt
      idx::Int = searchsortedfirst(crspdates, d)
      ((idx > Ndates) || (crspdates[idx] ≠ d)) && return missing
      return idx
    end

    #get the difference between the dividend dates in terms of trading days
    @mpar for r::DataFrameRow ∈ eachrow(event)
      ismissing(r.Lexdate) && continue

      Lexidx::MInt = crspdateormissing(r.Lexdate)
      ismissing(Lexidx) && continue

      exidx::MInt = crspdateormissing(r.exdate)
      r.days2lastdiv = exidx - Lexidx
    end
  end

  #winsorize to an fixed vlaue if desired
  deleterows!(crsp, (ismissing).(crsp[!,Fret]))
  if !isnothing(retbound)
    (retbound > 0.0) || error("retbound must be > 0")
    crsp.winsorized = (f::MFloat64->ismissing(f) ? missing : abs(f) > retbound).(crsp[!, Fret])

    crsp[!, Fret] .= (f::Float64->(max(min(f,retbound), -retbound))).(crsp[!, Fret])
  end

  sort!(event, [:permno, Fdate])

  #assume a fixed set of dates for each event
  dateranges::Dict =Dict{Date,NTuple{2, Date}}(
    crspdates[i]=> (crspdates[i-days], crspdates[i+days])
    for i ∈ (1+days):(length(crspdates)-days))

  #get the begin and end date of each range
  filter!(r->haskey(dateranges, r[Fdate]), event)
  event.begindate = (d->dateranges[d][1]).(event[!,Fdate])
  event.enddate = (d->dateranges[d][2]).(event[!,Fdate])

  #now check for overlapping dates
  event.oneex=trues(size(event,1))
  sevents::GroupedDataFrame = groupby(event, :permno)
  for sevent ∈ groupby(event, :permno)
    Nsevent::Int = size(sevent,1)
    for (i,r) ∈ enumerate(eachrow(sevent))
      t₀::Date = r[Fdate]
      (i > 1) &&  (r[:oneex] *= t₀>sevent[i-1,:enddate])
      (i < Nsevent) && (r[:oneex] *= t₀<sevent[i+1,:begindate])
    end
  end

  (:exyear ∈ names(event)) && event[event.exyear .== 2015, :] |> CSV.write(
    "output\\oneex-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")

  #now build the panel
  #ensure we capture the declaration day with the dividends
  lowerdays::Int = isnothing(daystostaleprice) ? days : max(daysextended, daystostaleprice)
  upperdays::Int = daysextended
  eventXday::Vector{NTuple{2,Int}} = vec(collect(Iterators.product(event.eid, -lowerdays:upperdays)))
  newevent::DataFrame = DataFrame(eid= (t->t[1]).(eventXday), day = (t->t[2]).(eventXday))
  event = join(event, newevent, on=:eid)

  dateidx::Dict = Dict(d=>i for (i,d) ∈ enumerate(crspdates))
  Ncrspdates::Int = length(crspdates)
  filter!(event) do r::DataFrameRow #only keep data within ranges
    idx = dateidx[r[Fdate]] + r.day
    return (idx ≥ 1) && (idx ≤ Ncrspdates)
  end
  event.date = (r::DataFrameRow->crspdates[dateidx[r[Fdate]] + r.day]).(eachrow(event))

  Nevent = size(event,1)
  event = join(event, crsp, on=[:permno, :date], kind=:left)
  @assert Nevent == size(event,1) #make sure we picked the correct keys

  #qc on the returns and net returns
  (s->finiteormissing!(event, s)).([Fret, :price])

  println("################# eventsize: $(size(event)) crspsize: $(size(crsp))")
  #println(describe(crsp))
  #println(describe(event))
  sort!(event, [:eid, :day])

  return event
end


function formevent(;refreshevent::Bool = true,
    eventpath::String = EVENT_PATH,
    refreshseocrsp::Bool=true,
    refreshdistcrsp::Bool=true,
    eventname::String = EVENT_NAME,
    inbinstream::Function = IN_BIN_STREAM,
    outbinstream::Function = OUT_BIN_STREAM,
    binextension::String = BIN_EXTENSION,
    seodays::Int = SEO_DAYS,
    seoyearrange::UnitRange = SEO_YEAR_RANGE,
    distdays::Int = DIST_DAYS,
    distdaysextended::Int = DIST_DAYS_EXTENDED,
    distyearrange::UnitRange = DIST_YEAR_RANGE
    )::NTuple{2,DataFrame}

  local bb::DataFrame, dist::DataFrame
  eventname = "$eventname-$(REPLICATION_TYPE[])-$(DIST_TYPE[])"

  if   refreshevent
    bb = readbbseo()
    #bb |> CSV.write("output\\bbevent-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")

    bb = formeventwindows(bb,seodays,
      refresheventcrsp=refreshseocrsp,
      validpermnos=unique(bb.permno),
      Fret = SEO_RET,
      Fdate=:date,
      prefix="seo",
      yearrange=seoyearrange
      )



    if DIST_TYPE[] ∈ [:primaryshort, :primaryfull]
      dist = primaryreaddist(prefix="", refreshdist=refreshdistcrsp)
    elseif DIST_TYPE[] ∈ [:oldshort, :oldfull]
      dist = oldreaddist()
    else
      @assert false #shouldn't happen due to previous checks
    end

    dist = preprocessdist(dist)

    (DIST_TYPE[] == :primaryshort || DIST_TYPE[]==:oldshort) && (
      dist = dist[3:5:10000,:])  #WARNING TEST DEBUG ONLY

    dist = formeventwindows(dist,distdays,
      daysextended=distdaysextended,
      refresheventcrsp=refreshdistcrsp,
      validpermnos=unique(dist.permno),
      Fret = DIST_RET,
      Fdate=:exdate,
      retbound = DIST_RET_BOUND,
      prefix="dist-$(DIST_TYPE[])",
      daystostaleprice=DIST_DAYS_TO_STALE_PRICE,
      yearrange=distyearrange
      )
    #println("Size of file: $(length(unique(dist.eid)))")
    #println(describe(dist))
    #error("stop")
    processdist!(dist)

    outbinstream("$eventpath\\$eventname.$binextension", (bb,dist))
  else
    (bb,dist) = inbinstream("$eventpath\\$eventname.$binextension")
  end

  #bb |> CSV.write("output\\cleanbb-$(REPLICATION_TYPE[]).csv")

  #if DIST_TYPE[] = :cleandist
  dist[dist.exyear .== 2015,:] |> CSV.write(
    "output\\cleandist-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-2015.csv")
  #println(describe(dist[dist.exyear .== 2015,:]))
  return (bb,dist)
end
