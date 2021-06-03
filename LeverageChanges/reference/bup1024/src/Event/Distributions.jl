

const DCLR_EX_DAYS = 9
const DIST_DAYS = 12
const DIST_RET = :net
const DIST_RET_BOUND = 0.2

const DIST_MIN_ALLOWED = 1#DIST_DAYS*2
const DIST_DAYS_TO_STALE_PRICE = 22


function getidentifierscrspevent(dist::DataFrame; Fdate::Symbol = :exdate)
  local ccm::DataFrame = loadccm()
  ccm = preprocessccm!(ccm)

  dist.eid = collect(1:size(dist,1)) #assign an id to each event

  dist = join(dist, ccm, on=:permno=>:lpermno)
  filter!(r::DataFrameRow -> r[Fdate] ∈ r.linkeffdate:Day(1):r.linkenddate, dist)
  @assert allunique(dist.eid) #no dups

  #rename!(dist, :lpermno=>:permno) #sanitize permno

  return dist
end

function readcrspevent(;fname::String = CRSP_DISTRIBUTIONS_NAME,
  path::String=EVENT_PATH,
  incsvstream::Function = IN_CSV_STREAM,
  crspdateformat::DateFormat = CRSP_DATE_FORMAT,
  Fdate::Symbol=:exdate,
  distdays::Int = DIST_DAYS)

  dist::DataFrame = incsvstream("$path\\$fname.csv.lz4") |> CSV.read |> DataFrame

  names!(dist, (s::Symbol->Symbol(lowercase(string(s)))).(names(dist)))

  rename!(dist, [:dclrdt=>:dclrdate, :exdt=>:exdate, :paydt=>:paydate])

  for f ∈ [:dclrdate, :exdate, :paydate]
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

  #ensures integrity
  dist = dist[dist.dclrexdays .≥ 0, :]


  #dist = getidentifierscrspevent(dist, Fdate=Fdate) NOTE: We will need to call this method- just not here

  return dist
end

function processcrspevent!(dist::DataFrame,
    dclrexdays::Int = DCLR_EX_DAYS,
    Fret::Symbol = DIST_RET,
    retbound::Float64=DIST_RET_BOUND,
    minallowed = DIST_MIN_ALLOWED,
    daystostaleprice::Int=DIST_DAYS_TO_STALE_PRICE,
    days::Int = DIST_DAYS,
    distyieldtol::Float64=10^-5)

  local sdists::GroupedDataFrame

  #now compute the dividend yield from the declaration date
  dist.yield = Vector{MFloat64}(undef, size(dist,1))
  dist.yield2 = Vector{MFloat64}(undef, size(dist,1))
  dist.winsorizedyield = falses(size(dist,1))
  dist.dclrextradingdays = Vector{MInt}(undef, size(dist,1))

  lowerdays::Int = max(daystostaleprice,days)
  sdists = groupby(dist, :eid)
  @mpar for i ∈ 1:length(sdists)
    sdist::SubDataFrame = sdists[i]
    #some basic integrity checks
    @assert maximum(sdist.divamt) == minimum(sdist.divamt)
    @assert size(sdist,1) == days + lowerdays + 1
    @assert maximum(sdist.dclrexdays) == minimum(sdist.dclrexdays)
    @assert sdist.dclrdate[1] == sdist.date[lowerdays + 1] - Day(sdist.dclrexdays[1])

    #find the row for the declaration date
    dclridx::Int = searchsortedfirst(sdist.date, sdist.dclrdate[1])
    (dclridx ≤ 0 || dclridx>size(sdist,1) || sdist.date[dclridx] ≠ sdist.dclrdate[1]) && continue

    #compute the yield and winsorize if necessary
    yield::MFloat64 = sdist.divamt[1]/sdist.price[dclridx]
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

  filter!(r::DataFrameRow->(r.day ≥ -days), dist)


  dist.primary = Vector{Bool}(undef, size(dist,1))
  dist.match = Vector{Bool}(undef, size(dist,1))
  dist.ordinary = Vector{Bool}(undef, size(dist,1))
  dist.dclrexok = Vector{Bool}(undef, size(dist,1))

  sdists = groupby(dist, :eid)
  @mpar for i ∈ 1:length(sdists)
    sdist = sdists[i]
    r::DataFrameRow = sdist[1, :] #examine fields common to each eid

    notmissing::Int = sum((!ismissing).(sdist[!,Fret]))
    sdist.match .= (notmissing ≥ minallowed) && (!ismissing(r.yield))
    sdist.ordinary .= r.distcd == 1232
    sdist.dclrexok .= r.dclrexdays ≥ dclrexdays
  end

  dist.primary .= dist.ordinary .* dist.match .* dist.oneex .* dist.dclrexok

  #some more integrity checks
  yieldmismatches::Int = 0
  avgyieldmismatch::Float64 = 0.0
  distcdmismatches::Int = 0
  dist.mismatches = falses(size(dist,1))
  for sdist ∈ groupby(dist, :eid)
    r::DataFrameRow = sdist[days+1,:]
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

  println("yieldmismatches: $yieldmismatches, avgyieldmismatch=$avgyieldmismatch
    distcdmismatches: = $distcdmismatches")

  dist[(max(end-20_000,1)):end,:] |> CSV.write(
    "output\\distest-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")

  dist[dist.mismatches,:] |> CSV.write(
    "output\\distmismatches-$(REPLICATION_TYPE[])-$(DIST_TYPE[]).csv")
end
