

const DCLR_EX_DAYS = 9
const DIST_DAYS = 12
const DIST_DAYS_EXTENDED = 40
const DIST_RET = :net
const DIST_RET_BOUND = 0.2

const DIST_MIN_ALLOWED = 1#DIST_DAYS*2
const DIST_DAYS_TO_STALE_PRICE = 22
const DIST_PRIMARY_YEAR_RANGE = YEAR_RANGE[]
const DIST_YEAR_RANGE = 1926:2018



function oldreaddist(;fname::String = CRSP_DISTRIBUTIONS_NAME,
  path::String=EVENT_PATH,
  incsvstream::Function = IN_CSV_STREAM,
  csvextension::String = CSV_EXTENSION,
  crspdateformat::DateFormat = CRSP_DATE_FORMAT)
  local dist::DataFrame


  dist = incsvstream("$path\\$fname.$csvextension") |> CSV.read |> DataFrame
  names!(dist, (s::Symbol->Symbol(lowercase(string(s)))).(names(dist)))
  rename!(dist, [:dclrdt=>:dclrdate, :exdt=>:exdate])

  @warn "Using crsp event file from crsp, but without the checks against main crsp.
    (There are a bunch of checks that were used against the main crsp.
    Check the 11/20 backup for details.)"

  return dist
end

function preprocessdist(dist::DataFrame;
      distdays::Int = DIST_DAYS)::DataFrame
  for f ∈ [:dclrdate, :exdate]
    (Date <: eltype(dist[!,f])) && continue
    dist[!,f] = (i::MInt->ismissing(i) ? missing : Date("$i", CRSP_DATE_FORMAT)).(dist[!,f])
  end


  #drop some bad data
  dist = dist[completecases(dist[!, [:permno, :exdate, :distcd]]), :]

  #get the last dividend for table 10
  sort!(dist, [:permno, :exdate])
  dist.Lexdate = Vector{MDate}(missing, size(dist,1))
  sdist::SubDataFrame = view(dist, (cd::MInt-> (cd ≥ 1000) && (cd < 2000)).(dist.distcd),:)
  ssdists::GroupedDataFrame = groupby(sdist, :permno)
  @mpar for i ∈ 1:length(ssdists)
    ssdist::SubDataFrame = ssdists[i]
    for (i::Int, r::DataFrameRow) ∈ enumerate(eachrow(ssdist))
      (i==1) && continue
      r.Lexdate = ssdist[i-1,:exdate]
    end
  end

  #get the interval since the last dividend
  #dist.days2lastdiv = ((Ldt::MDate, dt::Date)->
  #  (ismissing(Ldt) ? missing : (dt-Ldt).value)).(dist.Lexdate, dist.exdate)

  #drop some more bad data
  dist = dist[completecases(dist[!, [:dclrdate, :divamt]]), :]

  #time between exdate and dclrdate is important
  dist.dclrexdays = (d::Day->d.value).(dist.exdate .- dist.dclrdate)
  dist.exyear = (year).(dist.exdate)


  #filter!(r::DataFrameRow->r.exyear ∈ yearrange, dist) #restrict to desired years

  #ensures integrity
  dist = dist[dist.dclrexdays .≥ 0, :]

  return dist
end

function primaryreaddist(;  prefix = error("prefix is required"),
  path::String=EVENT_PATH,
  incsvstream::Function = IN_CSV_STREAM,
  crspdateformat::DateFormat = CRSP_DATE_FORMAT,
  refreshdist::Bool=true,
  distyearrange::UnitRange = DIST_YEAR_RANGE
  )

  if REPLICATION_TYPE[] == :monthly
    @warn "REPLICATION TYPE set to monthly but monthly is not supported for event.
      Event processing will run using daily data."
  end

  local dist::DataFrame = prepcrsp(
    crspcolumns = [:date, :permno, :divamt, :distcd, :dclrdt],
    crsptype = :dailyexpanded,
    refreshcrsp=refreshdist,
    validpermnos = nothing,
    prefix = "altdist",
    yearrange=distyearrange) do df::DataFrame

    filter!(r::DataFrameRow->(!ismissing(r.distcd)), df)
    #df.Tprice = Vector{MFloat64}(missing, size(df,1))
    #(!isnothing(daystostaleprice)) && lagprice!(df)
  end

  rename!(dist, [:date=>:exdate, :dclrdt=>:dclrdate])

  years::Vector{Int} = unique((year).(dist.exdate))
  if (minimum(distyearrange) ≠ minimum(years)) || (maximum(distyearrange) ≠ maximum(years))
    @error("Distribution year range not found in input dataframe.
      Check DIST_YEAR_RANGE vs the range contained in the dataframe pointed at by
      CRSPD_EXPANDED_NAME.")
  end

  println(describe(dist))

  return dist
end

function processdist!(dist::DataFrame,
    dclrexdays::Int = DCLR_EX_DAYS,
    Fret::Symbol = DIST_RET,
    retbound::Float64=DIST_RET_BOUND,
    minallowed = DIST_MIN_ALLOWED,
    daystostaleprice::Int=DIST_DAYS_TO_STALE_PRICE,
    days::Int = DIST_DAYS,
    daysextended::Int = DIST_DAYS_EXTENDED,
    primaryyearrange::UnitRange = DIST_PRIMARY_YEAR_RANGE,
    distyieldtol::Float64=10^-5)


  local sdists::GroupedDataFrame

  dist.primaryyearrange = (y::Int->y ∈ primaryyearrange).(dist.exyear)

  #now compute the dividend yield from the declaration date
  dist.yield = Vector{MFloat64}(undef, size(dist,1))
  dist.rawyield = Vector{MFloat64}(undef, size(dist,1))
  dist.yield2 = Vector{MFloat64}(undef, size(dist,1))
  dist.winsorizedyield = falses(size(dist,1))
  dist.dclrextradingdays = Vector{MInt}(undef, size(dist,1))
  dist.lowerdays = Vector{MInt}(undef, size(dist,1))
  dist.upperdays = Vector{MInt}(undef, size(dist,1))

  @assert issorted(dist, [:eid, :date])

  #lowerdays::Int = max(daystostaleprice,days,daysextended)
  #upperdays::Int = max(days, daysextended)
  sdists = groupby(dist, :eid)
  @mpar for i ∈ 1:length(sdists)
    sdist::SubDataFrame = sdists[i]
    #some basic integrity checks

    #find the upper and lower limits of the window
    local lowerdays::Int = -minimum(sdist.day)
    local upperdays::Int = maximum(sdist.day)
    sdist.lowerdays .= lowerdays
    sdist.upperdays .= upperdays

    @assert maximum(sdist.divamt) == minimum(sdist.divamt)
    @assert maximum(sdist.dclrexdays) == minimum(sdist.dclrexdays)
    @assert sdist.dclrdate[1] == sdist.date[lowerdays + 1] - Day(sdist.dclrexdays[1])

    #find the row for the declaration date
    #per the paper, we take the later of price at declaration or 22 days prior
    dclridx::Int = searchsortedfirst(sdist.date, sdist.dclrdate[1])
    (dclridx>size(sdist,1)) && continue
    ((sdist.day[dclridx]) < -daystostaleprice) && (dclridx = lowerdays+1-daystostaleprice)

    #compute the yield and winsorize if necessary
    yield::MFloat64 = sdist.divamt[1]/sdist.price[dclridx]
    sdist.rawyield .= ismissing(yield) || (!isfinite(yield)) ? missing : yield
    yield2::MFloat64 = sdist.divamt2[lowerdays+1]/sdist.price[dclridx]

    #no negative yields
    (!ismissing(yield)) && (yield < 0.0) && (yield=missing)
    (!ismissing(yield2)) && (yield2 < 0.0) && (yield2=missing)

    #if we can't get the yield from the event file, try the proper crsp file (yield2)
    ((ismissing(yield)) && (!ismissing(yield2)) && (!ismissing(dsdist.distcd2[lowerdays+1])) &&
      (dsdist.distcd[lowerdays+1] == sdist.distcd2[lowerdays+1]) && (yield=yield2))
    (ismissing(yield)) && continue

    sdist.dclrextradingdays .=  abs(sdist.day[dclridx])
    sdist.winsorizedyield .= yield > retbound
    sdist.yield .= min(yield,retbound)
    (!ismissing(yield2)) && (sdist.yield2 .= min(yield2,retbound))

  end
  #now determine the dividend yields from the exiration

  filter!(r::DataFrameRow->((r.day ≥ -daysextended) && (r.day ≤ daysextended)), dist)


  dist.primary = Vector{Bool}(undef, size(dist,1))
  dist.inrange = Vector{Bool}(undef, size(dist,1))
  dist.match = Vector{Bool}(undef, size(dist,1))
  dist.ordinary = Vector{Bool}(undef, size(dist,1))
  dist.dclrexok = Vector{Bool}(undef, size(dist,1))

  sdists = groupby(dist, :eid)
  @mpar for i ∈ 1:length(sdists)
    sdist = sdists[i]
    r::DataFrameRow = sdist[1, :] #examine fields common to each eid

    sdist.ordinary .= r.distcd == 1232
    sdist.dclrexok .= r.dclrexdays ≥ dclrexdays
    sdist.inrange .= (sdist.day .≥ -days) .* (sdist.day .≤ days)
    notmissing::Int = sum((!ismissing).(sdist[sdist.inrange,Fret]))
    sdist.match .= (notmissing ≥ minallowed) && (!ismissing(r.yield))

  end

  dist.primary .= (dist.ordinary .*
    dist.match .*
    dist.oneex .*
    dist.dclrexok .*
    dist.inrange .*
    dist.primaryyearrange)

  #println("inrange: $(sum(dist.inrange))")

  #some more integrity checks
  yieldmismatches::Int = 0
  avgyieldmismatch::Float64 = 0.0
  distcdmismatches::Int = 0
  dist.mismatches = falses(size(dist,1))
  for sdist ∈ groupby(dist, :eid)
    r::DataFrameRow = sdist[sdist.lowerdays[1]+1,:]
    @assert r.day == 0 #true by the prior filter

    if ((!ismissing(r.yield2)) && (!ismissing(r.yield)) && (r.distcd == r.distcd2) &&
      (abs(r.yield2 - r.yield)>distyieldtol))

      avgyieldmismatch += abs(r.yield2 - r.yield)
      yieldmismatches += 1
      sdist.mismatches .= true
    end

    if   ((!ismissing(r.distcd2)) && (!ismissing(r.distcd)) &&
      (r.distcd2 ≠ r.distcd))
      distcdmismatches += 1
    #  sdist.mismatches .= true
    end
  end
  avgyieldmismatch /= yieldmismatches

  #volume data for figure 4
  dist.dolvol = dist.vol .* dist.price
  dist.lnumtrd = (log).(dist.numtrd)
  dist.lvol = (log).(dist.vol)
  dist.ldolvol = (log).(dist.dolvol)
  (s::Symbol->finiteormissing!(dist, s)).([:lvol, :lnumtrd, :ldolvol])


  println("yieldmismatches: $yieldmismatches, avgyieldmismatch=$avgyieldmismatch
    distcdmismatches: = $distcdmismatches")

  #dist[(max(end-20_000,1)):end,:] |> CSV.write(
  #  "output\\distest-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")

  if occursin("short", string(DIST_TYPE[]))
    dist |> CSV.write(
      "output\\disttest-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")
  end
  dist[dist.mismatches,:] |> CSV.write(
    "output\\distmismatches-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")
end
