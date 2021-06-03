

#const RETAINED_COLS_BB_SEO = []
const BB_DATE_FORMAT = dateformat"mm/dd/yyyy"
const BB_EVENT_NAMES = ["seo1", "seo2"]

function cleanbbnames!(bb::DataFrame; )
  bbnames::Vector{String} = (string).(names(bb))

  bbnames .= (lowercase).(bbnames)
  bbnames .= (s->replace(s, r"[^a-z0-9]"=>"")).(bbnames)
  names!(bb, (Symbol).(bbnames))


end

function filterbbeventdates!(bb::DataFrame)
  #WARNING: STUB -Need to acquire data and put in the filter for 10 days between filing and offer date
  return bb
end

function getidentifiersforbb(bb::DataFrame)
  local ccm::DataFrame = loadccm()
  ccm = preprocessccm!(ccm)

  bb.eid = collect(1:size(bb,1)) #assign an id to each event
  ccm.cusip = (Symbol).(ccm.cusip)
  bb.cusip = (Symbol).(bb.cusip)

  bb = join(bb, ccm, on=:cusip=>:cusip)
  filter!(r::DataFrameRow -> r.date ∈ r.linkeffdate:Day(1):r.linkenddate, bb)
  @assert allunique(bb.eid) #no dups

  rename!(bb, :lpermno=>:permno) #sanitize permno

  return bb
end

function readbbseo(fnames::Vector{String} = BB_EVENT_NAMES;
  path::String = EVENT_PATH,
  bbdateformat = BB_DATE_FORMAT)

  local bb::DataFrame = DataFrame()

  #read in the bloomberg files
  for (i, f) ∈ enumerate(fnames)
    bbin::DataFrame = CSV.read("$path\\$f.csv") |> DataFrame
    cleanbbnames!(bbin)
    bb = vcat(bb, bbin)
  end

  #println("got here")

  filterbbeventdates!(bb)
  bb[!, :date] = (s->Date(s, bbdateformat)).(bb.announceddate)
  bb = getidentifiersforbb(bb)
end

function formeventwindows(event::DataFrame;
  days=7,
  refresheventcrsp::Bool = true)

  Nevent::Int = size(event,1)

  #create the fields that will hold the returns
  local dayfields::Vector{Symbol} = (d->d<0 ? Symbol("DM",abs(d)) : Symbol("D",d)).(collect((-1*days):days))
  for f::Symbol ∈ dayfields
    event[!, f] = Vector{MFloat64}(missing, Nevent)
  end

  validpermnos::Vector{Int} = unique(event[!, :permno])

  #load the crsp df
  local crsp::DataFrame = prepcrsp(
    crspcolumns = [:date, :permno, :ret, :net],
    crsptype = REPLICATION_TYPE[],
    refreshcrsp=refresheventcrsp,
    validpermnos = validpermnos,
    prefix = "event") do df::DataFrame
    computenetreturns!(df)
  end

  #make the appropriate sorts
  crspdates::Vector{Date} = sort(unique(crsp.date))
  sort!(crsp,[:permno,:date])
  sort!(event, [:permno, :date])

  validdateranges::Dict =Dict(crspdates[i]=> (crspdates[i-days], crspdates[i+days])
    for i ∈ (1+days):(length(crspdates)-days))

  #allocate space for the information we will get for the study
  event.begindate = Vector{Union{Missing,Date}}(missing, Nevent)
  event.enddate = Vector{Union{Missing,Date}}(missing, Nevent)
  event.tokeep=trues(Nevent)

  #index the dataframe rows
  eventidx::Dict{Int, SubDataFrame} = Dict{Int, SubDataFrame}(
    sevent.permno[1] => sevent for sevent ∈ groupby(event, :permno))

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  dateidx::Dict = Dict(d=>i for (i,d) ∈ enumerate(crspdates))

  #for each valid permno we have and each event for those permnos
  @mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    for r::DataFrameRow ∈ eachrow(eventidx[scrsp.permno[1]])


      t₀::Date = r.date
      Nscrsp::Int = size(scrsp,1)

      #locate the appropriate entry, and beform first set of validity checks
      i₀::Int = searchsortedfirst(scrsp.date, t₀)
      #println("i₀: $i₀ Nscrsp: $Nscrsp")
      println("i₀: $i₀
        days: $days
        length(scrsp.date): $(length(scrsp.date))
        scrsp.date: $(scrsp.date[i₀])
        t₀: $t₀
      ")
      if ((i₀ - days < 1) || #also takes care of certain "not found" cases
        (i₀ + days > Nscrsp) ||
        (scrsp.date[i₀] ≠ t₀) ||
        (!haskey(validdateranges, t₀)))

        r.tokeep = false
        error("break")
        continue
      end

      #check that the range is complete
      (tₗ::Date,tₕ::Date) = validdateranges[t₀]
      if (tₗ ≠ scrsp.date[i₀-days]) || (tₗ ≠ scrsp.date[i₀+days])
        r.tokeep = false
        continue
      end

      #now copy the return info
      for (i::Int, f::Symbol) ∈ enumerate(dayfields)
        iₜ::Int = i₀ + (i - days)
        r[f] = scrsp.net[iₜ]
      end


    end
  end

  #drop the rows that missed the prior check
  event = event[event.tokeep, :]

  #drop row with incomplete returns
  event = event[completecases(event[:, dayfields]), :]

  return event
end

function formevent(;refresheventcrsp::Bool = true)
  bb::DataFrame = readbbseo()
  bb |> CSV.write("output\\cleanevent.csv")
  bb = formeventwindows(bb)
  println(size(bb))


end
