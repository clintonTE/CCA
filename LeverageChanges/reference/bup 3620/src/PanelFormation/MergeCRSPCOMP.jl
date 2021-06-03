
#MAIN ENTRY POINT
#final columns for the decompressed file
#TODO: Fix the assertion checks

#merges crsp and compustat
# note refresh merge will be automatically run if any of the other two are run
#pre-op allows for a function to be applied before univ is written
function mergecrsp(crsp::Union{DataFrame,Nothing}=nothing,
  comp::Union{DataFrame,Nothing}=nothing;
  crsptype::Symbol = :monthly,
  refreshcomp::Bool = true,
  refreshmerge::Bool = true,
  datapath::String=DATA_PATH,
  prefix::String = "",
  univname::String = "$(prefix)$(crsptype==:monthly ? "CRSP-M" : "CRSP-D")_$(COMP_NAME)_$UNIV_NAME",
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  binextension::String = BIN_EXTENSION,
  validatemerged::Bool = VALIDATE_MERGED,
  retainedcolsuniv::Vector{Symbol} = [names(comp); names(crsp); :compid],
  cachedata::Bool = true
  )::DataFrame

  local univ::DataFrame

  #as of now, refreshmerged must be true if any of the other flags are true
  refreshmerge = refreshmerge || refreshcrsp || refreshcomp

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshmerge
    univ = mergecrspcomp!(crsp, comp, validatemerged=validatemerged)
    select!(univ, retainedcolsuniv)

    if TEST_OUTPUT
      d = Date(2015,12,31)
      CSV.write("output\\$(univname)_$d.csv", univ[(cd->cd==d).(univ.date), :])
    end

    if cachedata
      outbinstream("$datapath\\$univname.$binextension", univ)
    end

  else
    univ = inbinstream("$datapath\\$univname.$binextension")
  end

  @info "CRSP and Compustat merged in file $(univname).$binextension."

  return univ
end

#this is core code to prepare the merge
function matchcompid!(crsp::DataFrame, comp::DataFrame;
    mergetype::Symbol = :full,
    Fcrspdate=error("Fcrspdate is required"))

  (issorted(crsp, [:permno, Fcrspdate])) || error(
    "crsp comp merge must be sorted by :permno, $Fcrspdate")

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  scomps::Dict = Dict([scomp.lpermno[1]=>scomp for scomp ∈ groupby(comp, :lpermno)])

  #join(crsp, comp, kind)

@mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    local permno::Int = scrsp.permno[1]

    if haskey(scomps, permno)
      local scomp = scomps[permno]
      for rscomp ∈ eachrow(scomp)
        begindate::Date =  rscomp[:begindate]
        enddate::Date = rscomp[:enddate]
        sscrsp::SubDataFrame = view(scrsp, (d->d>begindate && d≤enddate).(scrsp[!, Fcrspdate]),:)
        if size(sscrsp,1) ≥ 1
          (sum((!ismissing).(sscrsp.compid)) > 0) && error("multiple valid compids for one permno")
          sscrsp.compid .= rscomp[:compid]
        end
      end
    end
  end


end


function mergecrspcomp!(crsp::DataFrame, comp::DataFrame;
    validatemerged::Bool=true,
    joinkind::Symbol = :inner,
    Fcrspdate::Symbol = :date,
    allowdatecollisions::Bool=false)::DataFrame

  local univ::DataFrame
  local scomp::SubDataFrame

  local Ncrsp::Int = size(crsp,1)
  #create a record id to form the dictionary
  comp.compid = Vector{MInt}(1:(size(comp,1)))
  comp.validinterval = ((begindate,enddate)->
    (begindate+Day(1)):Day(1):(enddate)).(comp.begindate, comp.enddate)

  crsp.compid = Vector{MInt}(undef, Ncrsp)

  #make the necessary allocations ahead of time
  print("matching...time to match: ")
  @time matchcompid!(crsp,comp, Fcrspdate=Fcrspdate)

  #println("sum crsp compid: $(sum(skipmissing(crsp.compid)))")
  println("Matching complete: $(100*(Ncrsp - sum((ismissing).(crsp.compid)))/Ncrsp)% matched.")
  #println("crsp records: $(Ncrsp)")
  (joinkind==:inner) && (deleterows!(crsp, (ismissing).(crsp.compid)))

  #mNOTE: main join
  univ = join(crsp, comp, on=:compid, kind=joinkind)

  validatemerged && testmerged(univ, Fdate=Fcrspdate, allowdatecollisions=allowdatecollisions)
  println("Total records: $(size(univ))")

  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\univ_mid_$d.csv", univ[(cd->cd==d).(univ[!, Fcrspdate]), :])
  end

  return univ
end

#do a bunch of integrity checks on the merged dataset
function testmerged(univ::DataFrame;
  months2stale::Int = MONTHS_2_STALE_LONG,
  Fdate::Symbol = error("Fdate is required"),
  allowdatecollisions::Bool = false)

  @mpar for r ∈ eachrow(univ)
    ismissing(r.gvkey) && continue
    #check date is valid
    d::Date = r[Fdate]
    (d ≥ r.linkeffdate) || error("failed check: (d ≥ r.linkeffdate)")
    (d ≤ r.linkenddate) || error("failed check: (d ≤ r.linkenddate)")
    (d ∈ r.validinterval) || error("failed check: (d ∈ r.validinterval)")
    (d > r.adate) || error("failed check: (d > r.adate)")
    (d ≤ r.Nadate) || error("failed check: (d ≤ r.Nadate)")
    period2stale = Month(months2stale)

    (!(d ≤ lastdayofmonth(firstdayofmonth(r.adate) + period2stale))) && error(
      "lastdayofmonth(firstdayofmonth(r.adate) + period2stale)))
      $r")

    (r.lpermno == r.permno) || error("failed check: (r.lpermno == r.permno)")
  end

  univ = sort(univ, [:permno, :adate, Fdate])

  sunivs::GroupedDataFrame = groupby(univ, :permno)
  @mpar for i ∈ 1:length(sunivs)
    local suniv::SubDataFrame = sunivs[i]
    maxfiscal::Int= 1900
    maxadate::Date = Date(1900,1,1)

    for r ∈ eachrow(suniv)
      ismissing(r.gvkey) && continue
      (!allowdatecollisions) && ((r[Fdate] > maxadate) || error("failed check: r[Fdate] > maxadate
        permno=$(r[:permno]), r[Fdate]=$(r[Fdate]), maxadate=$maxadate"))

      (allowdatecollisions) && ((r[Fdate] ≥ maxadate) || error("failed check: r[Fdate] ≥ maxadate
        permno=$(r[:permno]), r[Fdate]=$(r[Fdate]), maxadate=$maxadate"))

      maxadate = r[Fdate]

      (maxfiscal ≤ r.fyear) ||  error("failed check: maxfiscal ≤ r.fyear")
      (maxfiscal ≠ r.fyear) && (maxfiscal=r.fyear)
    end

  end
end
