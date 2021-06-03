
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
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  validatemerged::Bool = VALIDATE_MERGED,
  retainedcolsuniv::Vector{Symbol} = [names(comp); names(crsp); :compid],
  cachedata::Bool = true
  )::DataFrame

  local univ::DataFrame

  #as of now, refreshmerged must be true if any of the other flags are true
  refreshmerge = refreshmerge || refreshcrsp || refreshcomp

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshmerge || (!isfile("$datapath\\$univname.jls.lz4"))
    univ = mergecrspcomp!(crsp, comp, validatemerged=validatemerged)
    select!(univ, retainedcolsuniv)

    if TEST_OUTPUT
      d = Date(2015,12,31)
      CSV.write("output\\$(univname)_$d.csv", univ[(cd->cd==d).(univ.date), :])
    end

    if cachedata
      outjlsstream("$datapath\\$univname.jls.lz4") do s
        serialize(s, univ)
      end
    end

  else
    injlsstream("$datapath\\$univname.jls.lz4") do s
      univ = deserialize(s)
    end
  end

  @info "CRSP and Compustat merged in file $(univname).jls.lz4."

  return univ
end

#this is core code to prepare the merge
function matchcompid!(crsp::DataFrame, comp::DataFrame; mergetype::Symbol = :full)

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  scomps::Dict = Dict([scomp.lpermno[1]=>scomp for scomp ∈ groupby(comp, :lpermno)])

  @mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    local permno::Int = scrsp.permno[1]

    if haskey(scomps, permno)
      local scomp = scomps[permno]
      for rscrsp ∈ eachrow(scrsp)
        #ctr::Int = 0
        local compidaddress::NInt
        local rscrspdate::Date = rscrsp.date

        #check if the crsp row has a comp row w/ same permno and within a valid interval
        compidaddress = findfirst(
          interval->rscrspdate ∈ interval, scomp.validinterval)

        (!isnothing(compidaddress)) && (rscrsp.compid = scomp.compid[compidaddress])
      end
    end
  end


end


function mergecrspcomp!(crsp::DataFrame, comp::DataFrame;
    validatemerged::Bool=true)::DataFrame

  local univ::DataFrame
  local scomp::SubDataFrame

  local Ncrsp::Int = size(crsp,1)
  #create a record id to form the dictionary
  comp.compid = Vector{MInt}(1:(size(comp,1)))
  comp.validinterval = ((begindate,enddate)->
    (begindate+Day(1)):Day(1):(enddate)).(comp.begindate, comp.enddate)

  crsp.compid = Vector{MInt}(undef, Ncrsp)

  #make the necessary allocations ahead of time
  @time matchcompid!(crsp,comp)

  #println("sum crsp compid: $(sum(skipmissing(crsp.compid)))")
  println("Matching complete: $((Ncrsp - sum((ismissing).(crsp.compid)))/Ncrsp)% matched.")
  #println("crsp records: $(Ncrsp)")
  deleterows!(crsp, (ismissing).(crsp.compid))

  #mNOTE: main join
  univ = join(crsp, comp, on=:compid, kind=:inner)

  validatemerged && testmerged(univ)
  println("Total records: $(size(univ))")

  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\univ_mid_$d.csv", univ[(cd->cd==d).(univ.date), :])
  end

  return univ
end

#do a bunch of integrity checks on the merged dataset
function testmerged(univ::DataFrame;
  months2stale::Int = MONTHS_2_STALE_LONG)

  @mpar for r ∈ eachrow(univ)
    #check date is valid
    d::Date = r.date
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

  univ = sort(univ, [:permno, :adate, :date])

  sunivs::GroupedDataFrame = groupby(univ, :permno)
  @mpar for i ∈ 1:length(sunivs)
    local suniv::SubDataFrame = sunivs[i]
    maxfiscal::Int= 1900
    maxadate::Date = Date(1900,1,1)

    for r ∈ eachrow(suniv)
      (r.date > maxadate) || error("failed check: r.date > maxadate")
      maxadate = r.date

      (maxfiscal ≤ r.fyear) ||  error("failed check: maxfiscal ≤ r.fyear")
      (maxfiscal ≠ r.fyear) && (maxfiscal=r.fyear)
    end

  end
end
