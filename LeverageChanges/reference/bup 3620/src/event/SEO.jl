

#const RETAINED_COLS_BB_SEO = []
const BB_DATE_FORMAT = dateformat"mm/dd/yyyy"
const BB_EVENT_NAMES = ["seo-nasdaq", "seo-otc", "seo-ny"]
const BB_EVENT_RETAINED = [:date, :adate, :offertype, :filedate, :permno, :gvkey, :offersizem,
  :primarysharesm, :securitytype, :primaryexchange, :issuerticker]
const SEO_DAYS = 7
const SEO_RET = :ret
const SEO_MIN_ALLOWED = SEO_DAYS*2
const SEO_YEAR_RANGE = YEAR_RANGE[]


function cleanbbnames!(bb::DataFrame; )
  bbnames::Vector{String} = (string).(names(bb))

  bbnames .= (lowercase).(bbnames)
  bbnames .= (s->replace(s, r"[^a-z0-9]"=>"")).(bbnames)
  rename!(bb, ((oldn::Symbol, newn::Symbol)->oldn=>newn).(names(bb), (Symbol).(bbnames)))


end

function filterbbeventdates!(bb::DataFrame, bbdateformat = BB_DATE_FORMAT,
    yearrange::UnitRange{Int} = YEAR_RANGE[])
  #initial filters to ensure complete data
  filter!(r::DataFrameRow->
    (!ismissing(r.offersizem)) && (r.offersizem > 0.) &&
    (!ismissing(r.offertype)) && ("IPO" ≠ r.offertype[1:3]) &&
    (!occursin("Secondary", r.offertype)) &&
    (!ismissing(r.primaryexchange)) && ("OTC" ≠ r.primaryexchange[1:3]) &&
    (length(r.cusip) == 9) &&
    (!ismissing(r.issuerticker)) &&
    (!ismissing(r.primarysharesm)) && (r.primarysharesm > 0.) &&
    (r.offerstage == "Trading" )
    , bb)

  #begin processing the dates
  for f ∈ [:launchdate, :filingdate, :filingtermdate,
      :pricingdate, :tradedateus, :effectivedate, :announceddate]

    bb[!, f] = (s::MString->ismissing(s) ? missing : Date(s, bbdateformat)).(bb[!,f])
  end

  rename!(bb, :announceddate=>:adate)
  bb[!, :filedate] = (r::DataFrameRow->coalesce(r.filingtermdate, r.filingdate, r.launchdate)).(eachrow(bb))
  bb[!, :date] = (r::DataFrameRow->coalesce(r.pricingdate, r.tradedateus, r.effectivedate)).(eachrow(bb))

  #issue and file date must be at least 10 days appart
  filter!(r::DataFrameRow->
    (!ismissing(r.filedate)) &&  (!ismissing(r.date)) && ((r.date - r.filedate).value ≥ 10) &&
    (!ismissing(r.adate)) && (r.date > r.adate)
    , bb)


  filter!(r::DataFrameRow->year(r.adate) ∈ yearrange, bb) #restrict to desired years
  #WARNING: STUB -Need to acquire data and put in the filter for 10 days between filing and offer date
  return bb
end

function getidentifiersbb(bb::DataFrame)
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
  bbretained::Vector{Symbol} = BB_EVENT_RETAINED)

  local bb::DataFrame = DataFrame()

  #read in the bloomberg files
  for (i, f) ∈ enumerate(fnames)
    bbin::DataFrame = CSV.read("$path\\$f.csv") |> DataFrame
    cleanbbnames!(bbin)
    bb = vcat(bb, bbin)
  end

  #println("got here")

  filterbbeventdates!(bb)

  bb = getidentifiersbb(bb)
  select!(bb, bbretained)
end
