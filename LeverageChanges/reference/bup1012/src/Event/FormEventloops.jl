


getdayfields(days::Int) = (d->d<0 ? Symbol("DM",abs(d)) : Symbol("D",d)).(collect((-1*days):days))

function formeventwindows(event::DataFrame, days::Int;
  refresheventcrsp::Bool = true,
  Fret::Symbol = error("Fret is required"),
  validpermnos::Vector{Int} = unique(event[!, :permno]),
  Fdate::Symbol = error("Fdate is required"),
  retbound::NFloat64 = nothing,
  prefix = error("prefix is required")
  )

  days>0 || error("days must be greater than 0")
  Nevent::Int = size(event,1)
  event.eid = collect(1:Nevent)

  #load the crsp df
  local crsp::DataFrame = prepcrsp(
    crspcolumns = [:date, :permno, :ret, :net, :price, :Lprice],
    crsptype = REPLICATION_TYPE[],
    refreshcrsp=refresheventcrsp,
    validpermnos = validpermnos,
    prefix = "$(prefix)event") do df::DataFrame
    computenetreturns!(df)
    sort!(df,[:permno,:date])

    #lag the price in the simplest fashion
    df.Lprice = Vector{MFloat64}(missing, size(df,1))
    for sdf ∈ groupby(df, :permno)
      Nsdf::Int = size(sdf,1)
      (Nsdf ≤ 1) && continue

      sdf.Lprice[2:end] .= sdf.price[1:(end-1)]
    end
  end

  #make the appropriate sorts
  crspdates::Vector{Date} = sort(unique(crsp.date))

  #winsorize to an fixed vlaue if desired
  deleterows!(crsp, (ismissing).(crsp[!,Fret]))
  if !isnothing(retbound)
    (retbound > 0.0) || error("retbound must be > 0")
    crsp.net .= (f::Float64->(max(min(f,retbound), -retbound))).(crsp.net)
    crsp.ret .= (f::Float64->(max(min(f,retbound), -retbound))).(crsp.ret)
  end

  sort!(event, [:permno, Fdate])

  #assume a fixed set of date for each event
  dateranges::Dict =Dict{Date,NTuple{2, Date}}(
    crspdates[i]=> (crspdates[i-days], crspdates[i+days])
    for i ∈ (1+days):(length(crspdates)-days))

  #get the begin and end date of each range
  filter!(r->haskey(dateranges, r[Fdate]), event)
  event.begindate = (d->dateranges[d][1]).(event[!,Fdate])
  event.enddate = (d->dateranges[d][2]).(event[!,Fdate])
  event.tokeep=trues(size(event,1))

  sevents::GroupedDataFrame = groupby(event, :permno)

  #now check for overlapping dates
  for sevent ∈ groupby(event, :permno)
    Nsevent::Int = size(sevent,1)
    for (i,r) ∈ enumerate(eachrow(sevent))
      t₀::Date = r[Fdate]
      (i > 1) &&  (r[:tokeep] *= t₀>sevent[i-1,:enddate])
      (i < Nsevent) && (r[:tokeep] *= t₀<sevent[i+1,:begindate])
    end
  end
  event = event[event.tokeep, :]

  #now build the panel
  eventXday::Vector{NTuple{2,Int}} = vec(collect(Iterators.product(event.eid, -days:days)))
  newevent::DataFrame = DataFrame(eid= (t->t[1]).(eventXday), day = (t->t[2]).(eventXday))
  event = join(event, newevent, on=:eid)

  dateidx::Dict = Dict(d=>i for (i,d) ∈ enumerate(crspdates))
  event.date = (r::DataFrameRow->crspdates[dateidx[r[Fdate]] + r.day]).(eachrow(event))

  Nevent = size(event,1)
  #event = join(event, crsp, on=[:permno, :date], kind=:left)
  scrsps
  @assert Nevent == size(event,1) #make sure we picked the correct keys

  #qc on the returns and net returns
  (s->finiteormissing!(event, s)).([:ret, :net, :price, :Lprice])

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

    #keep track of permnos so that we only refresheventcrsp once
    #validpermnos::Vector{Int} = unique([bb.permno; dist.permno])
    bb = formeventwindows(bb,seodays,
      refresheventcrsp=refreshseocrsp,
      validpermnos=unique(bb.permno),
      Fret = SEO_RET,
      Fdate=:date,
      prefix="seo"
      )
    dist = formeventwindows(dist,distdays,
      refresheventcrsp=refreshdistcrsp,
      validpermnos=unique(dist.permno),
      Fret = DIST_RET,
      Fdate=:exdate,
      retbound = DIST_RET_BOUND,
      prefix="dist"
      )
    print("process: ")
    @time processcrspevent!(dist)

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
