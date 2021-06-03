
#MAIN ENTRY POINT
#final columns for the decompressed file

#merges crsp and compustat
# note refresh merge will be automatically run if any of the other two are run
#pre-op allows for a function to be applied before univ is written
function makeuniv(prewriteop::Function = (noop(df::DataFrame)=df);
  refreshcomp::Bool = true,
  refreshcrsp::Bool = true,
  refreshmerge::Bool = true,
  datapath::String=DATA_PATH,
  univname::String = UNIV_NAME,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  validatemerged::Bool = true)::DataFrame

  local crsp::DataFrame
  local comp::DataFrame
  local univ::DataFrame

  #as of now, refreshmerged must be true if any of the other flags are true
  refreshmerge = refreshmerge || refreshcrsp || refreshcomp || validatemerged

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshmerge || (!isfile("$datapath\\$univname.jls.lz4"))
    crsp = prepcrsp(refreshcrsp=refreshcrsp)
    comp = prepcomp(refreshcomp=refreshcomp)

    univ = mergecrspcomp!(crsp, comp, validatemerged=validatemerged)
    univ = prewriteop(univ)

    if TEST_OUTPUT
      d = Date(2015,12,31)
      CSV.write("output\\univ_$d.csv", univ[(cd->cd==d).(univ.date), :])
    end

    outjlsstream("$datapath\\$univname.jls.lz4") do s
      serialize(s, univ)
    end

    #NOTE: FINAL VERSION
    #preprocesscomp(compa, compq)
    #serialize("$datapath\\$compname.jls", comp)
  else
    injlsstream("$datapath\\$univname.jls.lz4") do s
      univ = deserialize(s) #NOTE: DEBUG ONLY
    end
  end

  @info "CRSP and Compustat merged in file $(univname).jls.lz4."

  return univ
end

#this is core code to prepare the merge
function matchrid!(crsp::DataFrame,comp::DataFrame)

  #scrsps::GroupedDataFrame = groupby(crsp, :permno)
  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  scomps::Dict = Dict((scomp->scomp.lpermno[1]=>scomp).(groupby(comp, :lpermno)))

  @mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    local permno::Int = scrsp.permno[1]

    if haskey(scomps, permno)
      local scomp = scomps[permno]
      for rscrsp ∈ eachrow(scrsp)
        #ctr::Int = 0
        local ridaddress::NInt
        local rscrspdate::Date = rscrsp.date

        #check if the crsp row has a comp row w/ same permno and within a valid interval
        ridaddress = findfirst(
          interval->rscrspdate ∈ interval, scomp.validinterval)

        (!isnothing(ridaddress)) && (rscrsp.rid = scomp.rid[ridaddress])
      end
    end
  end
end


function mergecrspcomp!(crsp::DataFrame, comp::DataFrame;
    retainedcolsuniv::Vector{Symbol} = RETAINED_COLS_UNIV,
    validatemerged::Bool=true)::DataFrame

  local univ::DataFrame
  local scomp::SubDataFrame

  local Ncrsp::Int = size(crsp,1)
  #create a record id to form the dictionary
  comp.rid = Vector{MInt}(1:(size(comp,1)))
  comp.validinterval = ((begindate,enddate)->
    (begindate+Day(1)):Day(1):(enddate)).(comp.begindate, comp.enddate)

  crsp.rid = Vector{MInt}(undef, Ncrsp)

  #make the necessary allocations ahead of time
  @time matchrid!(crsp, comp)

  println("sum crsp rid: $(sum(skipmissing(crsp.rid)))")
  println("Matching complete: $((Ncrsp - sum((ismissing).(crsp.rid)))/Ncrsp)% matched.")
  #println("crsp records: $(Ncrsp)")

  #mNOTE: main join
  univ = join(crsp, comp, on=:rid, kind=:inner)

  validatemerged && testmerged(univ)
  println("Total records: $(size(univ))")

  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\univ_mid_$d.csv", univ[(cd->cd==d).(univ.date), :])
  end

  select!(univ, retainedcolsuniv)

  println(describe(univ))

  return univ
end

#do a bunch of integrity checks on the merged dataset
function testmerged(univ::DataFrame)

  for r ∈ eachrow(univ)
    #check date is valid
    d::Date = r.date
    @assert d ≥ r.linkeffdate
    @assert d ≤ r.linkenddate
    @assert d ∈ r.validinterval
    @assert d > r.adateq
    @assert d ≤ r.Nadateq
    @assert d ≤ lastdayofmonth(firstdayofmonth(r.adateq) + Month(6))

    @assert r.lpermno == r.permno
  end

  sort!(univ, [:permno, :adateq, :date])

  for suniv ∈ groupby(univ, :permno)
    maxfiscal::YearQuarter = YearQuarter(1900,1)
    maxadate::Date = Date(1900,1,1)

    for r ∈ eachrow(suniv)
      @assert r.date > maxadate
      maxadate = r.date

      @assert maxfiscal ≤ r.fyearqtr
      (maxfiscal ≠ r.fyearqtr) && (maxfiscal=r.fyearqtr)
    end

  end
end
