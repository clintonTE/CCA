


getdayfields(days::Int) = (d->d<0 ? Symbol("DM",abs(d)) : Symbol("D",d)).(collect((-1*days):days))

function formeventwindows(event::DataFrame, days::Int;
  refresheventcrsp::Bool = true,
  Fret::Symbol = error("Fret is required"),
  validpermnos::Vector{Int} = unique(event[!, :permno]),
  Fdate::Symbol = error("Fdate is required"),
  retbound::NFloat64 = nothing,
  prefix = error("prefix is required"),
  daystostaleprice::NInt=nothing
  )

  days>0 || error("days must be greater than 0")
  Nevent::Int = size(event,1)
  event.eid = collect(1:Nevent)


  #load the crsp df
  local crsp::DataFrame = prepcrsp(
    crspcolumns = [:date, :permno, Fret, :price],
    crsptype = REPLICATION_TYPE[],
    refreshcrsp=refresheventcrsp,
    validpermnos = validpermnos,
    prefix = "$(prefix)event") do df::DataFrame
    computenetreturns!(df)
    sort!(df,[:permno,:date])

    #df.Tprice = Vector{MFloat64}(missing, size(df,1))
    #(!isnothing(daystostaleprice)) && lagprice!(df)

  end

  #make the appropriate sorts
  crspdates::Vector{Date} = sort(unique(crsp.date))

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
  event.nooverlap=trues(size(event,1))
  sevents::GroupedDataFrame = groupby(event, :permno)
  for sevent ∈ groupby(event, :permno)
    Nsevent::Int = size(sevent,1)
    for (i,r) ∈ enumerate(eachrow(sevent))
      t₀::Date = r[Fdate]
      (i > 1) &&  (r[:nooverlap] *= t₀>sevent[i-1,:enddate])
      (i < Nsevent) && (r[:nooverlap] *= t₀<sevent[i+1,:begindate])
    end
  end

  #now build the panel
  #ensure we capture the declaration day with the dividends
  lowerdays::Int = isnothing(daystostaleprice) ? days : max(days, daystostaleprice)
  upperdays::Int = days
  eventXday::Vector{NTuple{2,Int}} = vec(collect(Iterators.product(event.eid, -lowerdays:upperdays)))
  newevent::DataFrame = DataFrame(eid= (t->t[1]).(eventXday), day = (t->t[2]).(eventXday))
  event = join(event, newevent, on=:eid)

  dateidx::Dict = Dict(d=>i for (i,d) ∈ enumerate(crspdates))
  event.date = (r::DataFrameRow->crspdates[dateidx[r[Fdate]] + r.day]).(eachrow(event))

  Nevent = size(event,1)
  event = join(event, crsp, on=[:permno, :date], kind=:left)
  @assert Nevent == size(event,1) #make sure we picked the correct keys

  #qc on the returns and net returns
  (s->finiteormissing!(event, s)).([Fret, :price])

  println("################# eventsize: $(size(event)) crspsize: $(size(crsp))")
  println(describe(crsp))
  println(describe(event))

  return event
end



function formevent(;refreshevent::Bool = true,
    eventpath::String = EVENT_PATH,
    refreshseocrsp::Bool=true,
    refreshdistcrsp::Bool=true,
    eventname::String = EVENT_NAME,
    injlsstream::Function = IN_JLS_STREAM,
    outjlsstream::Function = OUT_JLS_STREAM,
    seodays::Int = SEO_DAYS,
    distdays::Int = DIST_DAYS
    )::NTuple{2,DataFrame}

  local bb::DataFrame, dist::DataFrame
  eventname = "$eventname-$(REPLICATION_TYPE[])"

  if   refreshevent
    bb = readbbseo()
    bb |> CSV.write("output\\bbevent-$(REPLICATION_TYPE[]).csv")
    dist = readcrspevent()

    #table8(dist, panela=true)

    #keep track of permnos so that we only refresheventcrsp once
    #validpermnos::Vector{Int} = unique([bb.permno; dist.permno])
    bb = formeventwindows(bb,seodays,
      refresheventcrsp=refreshseocrsp,
      validpermnos=unique(bb.permno),
      Fret = SEO_RET,
      Fdate=:date,
      prefix="seo"
      )

    #dist = dist[1:100,:]  #WARNING WARNING WARNING TEST DEBUG ONLY
    dist = formeventwindows(dist,distdays,
      refresheventcrsp=refreshdistcrsp,
      validpermnos=unique(dist.permno),
      Fret = DIST_RET,
      Fdate=:exdate,
      retbound = DIST_RET_BOUND,
      prefix="dist",
      daystostaleprice=DIST_DAYS_TO_STALE_PRICE
      )
    #println("Size of file: $(length(unique(dist.eid)))")
    #println(describe(dist))
    #error("stop")
    processcrspevent!(dist)

    outjlsstream("$eventpath\\$eventname.jls.lz4") do s
      serialize(s, (bb,dist))
    end
  else
    injlsstream("$eventpath\\$eventname.jls.lz4") do s
      (bb,dist) = deserialize(s) #NOTE: DEBUG ONLY
    end
  end

  bb |> CSV.write("output\\cleanbb-$(REPLICATION_TYPE[]).csv")
  #dist |> CSV.write("output\\cleandist-$(REPLICATION_TYPE[]).csv")
  return (bb,dist)
end
