


function getvalidpermnos(df::DataFrame)::Vector{Int}
  permnosym = :permno ∈ propertynames(df) ? :permno : :lpermno
  validpermnos::Vector{Int} = unique(df[!, permnosym])

  return validpermnos
end

function getvalidpermnos()::Vector{Int}
  local ccm::DataFrame = loadccm()
  ccm = preprocessccm!(ccm)

  return getvalidpermnos(ccm)
end


#do some pre-filtering to reduce the processing time
function prefiltercrsp!(crsp::DataFrame, validpermnos::Vector{Int})::Nothing
  #rename columns to lower case


  crsp.tokeep = trues(size(crsp,1))

  crspgrp::GroupedDataFrame = groupby(crsp, :permno)
  @mpar for i ∈ 1:length(crspgrp)
    scrsp::SubDataFrame = crspgrp[i]

    if scrsp.permno[1] ∉ validpermnos
      scrsp.tokeep .= false
    end
  end

  filter!(r::DataFrameRow->r.tokeep, crsp)
  select!(crsp, Not([:tokeep]))


  return nothing

end

@inline normalizecrspnames!(crsp::DataFrame) = rename!(crsp, (s::String->lowercase(s)).(names(crsp)))

function prepcrsp(op!::Function = (noop(x)=x);
  crsptype = :not_set,
  refreshcrsp=true,
  datapath=DATA_PATH,
  incsvstream::Function = IN_CSV_STREAM,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  csvextension::String = CSV_EXTENSION,
  binextension::String = BIN_EXTENSION,
  validpermnos::Union{Nothing, Vector{Int}} = getvalidpermnos(),
  prefix::String = "",
  crspcolumns::Union{Nothing,Vector{Symbol}} = [:permno, :date, :ret, :lret],
  testoutput::Bool = TEST_OUTPUT,
  yearrange::UnitRange{Int} = YEAR_RANGE[])::DataFrame

  crsptype ∈ [:monthly, :daily, :dailyexpanded] || error("Set crsptype to :monthly or :daily instead of $crsptype")

  local crspname
  local crspparts

  if crsptype == :monthly
    crspname = CRSPM_NAME
    #crspparts = CRSPM_PARTS
  elseif crsptype==:daily
    crspname = CRSPD_NAME
    #crspparts = CRSPD_PARTS
  elseif crsptype==:dailyexpanded
    crspname = CRSPD_EXPANDED_NAME
    #crspparts = CRSPD_EXPANDED_PARTS
  else
    @assert false
  end


  local p::String = "$datapath\\$(prefix)$crspname"
  local crsp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcrsp
    #do the below so we can pre-filter

    pin::String = "$datapath\\$crspname"
    crsp = incsvstream("$pin.$csvextension") |> CSV.File |> DataFrame
    normalizecrspnames!(crsp)
    (!isnothing(validpermnos)) && (prefiltercrsp!(crsp, validpermnos))

    #=if length(crspparts) > 1
      #try to do this is in a space efficient way
      crsptypes::Vector{Type} = (eltype).(eachcol(crsp))
      for crsppart ∈ crspparts[2:end]
        pin = "$datapath\\$crsppart"
        local addedcrsp::DataFrame = CSV.read(incsvstream("$pin.$csvextension"), threaded=CSV_THREADED, types=crsptypes) |> DataFrame
        normalizecrspnames!(addedcrsp)
         (!isnothing(validpermnos)) && (prefiltercrsp!(addedcrsp, validpermnos))
        append!(crsp, addedcrsp)
      end
    end=#

    #filter by year first
    crsp.date = (s->parsecrsp(Date, s)).(crsp.date)
    filter!(r::DataFrameRow->year(r.date) ∈ yearrange, crsp) #restrict to desired years

    crsp = preprocesscrsp(crsp, crsptype)
    op!(crsp)
    select!(crsp, crspcolumns)

    if testoutput
      d = Date(2015,12,31)
      testoutput && CSV.write("output\\$(prefix)_$(crspname)_$d.csv", crsp[(cd->cd==d).(crsp.date), :])
    end


    outbinstream("$p.$binextension", crsp)

  else
    crsp = inbinstream("$p.$binextension")
  end

  @info "CRSP data loaded and/or saved into file $(prefix)$crspname.$binextension"

  return crsp
end

function computecrspreturns!(crsp::DataFrame)
  crsp.ret = (s::MString->parsecrsp(Float64, s)).(crsp.ret)
  crsp.dlret= (s::MString->parsecrsp(Float64, s)).(crsp.dlret)


  rename!(crsp, :ret=>:retlisted)
  crsp.retlisted .= (r->
    !ismissing(r) && (r > -1.0) ? r : missing).(crsp.retlisted)
  crsp.dlret .= (r->
    !ismissing(r) && (r > -1.0) ? r : missing).(crsp.dlret)

  crsp.ret = ((retlisted,dlret)->
    ismissing(retlisted) && ismissing(dlret) ? missing :
      ismissing(dlret) ? retlisted :
      ismissing(retlisted) ? dlret :
      (1+retlisted)*(1+dlret) - 1.
    ).(crsp.retlisted,crsp.dlret)

  #winsorizekeydata!(crsp, fields2winsorize = [:ret])

  #get log returns
  crsp.lret = (log).(1.0 .+ crsp.ret)
  delete!(crsp, (ismissing).(crsp.ret))

  return crsp
end

function filtercrsp!(crsp::DataFrame; filternonfinancials::Bool = false)

    #NOTE: SHRCD==10,11 : US Common Shares
    #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
    delete!(crsp, (!).(completecases(crsp[:,[:shrcd, :exchcd]])))
    local shrcdvalid::Set{Int} = Set{Int}((10,11))
    local exchcdvalid::Set{Int} = Set{Int}((1,2,3))
    filter!(r::DataFrameRow->
      (r.shrcd ∈ shrcdvalid) && (r[:exchcd] ∈ exchcdvalid), crsp)

    #NOTE: r[:SICCD] ∉ 6000:6999 : No financial companies

    if filternonfinancials
        error("Need to get the data again- sic code no longer incldued in crsp file")
        crsp.siccd = (s::MString->parsecrsp(Int, s)).(crsp.siccd)
        local siccdinvalid::Set{Int} = Set{Int}(6000:6999)
        filter!(r::DataFrameRow->
          (!ismissing(r.siccd)) && (r.siccd ∉ siccdinvalid), crsp)
    end

    return crsp
end

function computecrspprices!(crsp::DataFrame, crsptype::Symbol)
  #Now work on the price
  #in this case, use the mean of the two values
  #price cannot be zero
  crsp.prc .= (
    prc->coalesce(iszero(prc),true) ? missing : prc).(crsp.prc)
  crsp.dlprc .= (
    dlprc->coalesce(iszero(dlprc),true) ? missing : dlprc).(crsp.dlprc)

  crsp.price = ((prc,dlprc)->
    ismissing(prc) && ismissing(dlprc) ? missing :
      ismissing(dlprc) ? prc :
      ismissing(prc) ? dlprc :
      (prc+dlprc)/2
  ).(crsp.prc, crsp.dlprc)


  #add in the alternative price if nothing else is available
  if crsptype == :monthly
    crsp.altprc .= (
      altprc->coalesce(iszero(altprc),true) ? missing : altprc).(crsp.altprc)
    crsp.altprc .= (abs).(crsp.altprc)
    crsp.price .= ((price,altprc)->
      ismissing(price) ? altprc : price).(crsp.price,crsp.altprc)
  end


  #negative values computed via bid-ask spread method
  crsp.price .= (abs).(crsp.price)

  return crsp
end

function computecrspmktcap!(crsp::DataFrame)

  #adjust the units and compute the market cap
  crsp.mktcap = Vector{MInt}((round).(crsp.price .* crsp.shrout .* 1000)) #Note- put this in dollars
  crsp.mktcap .= (i::MInt->coalesce(iszero(i),true) ? missing : i).(crsp.mktcap)

  delete!(crsp, (ismissing).(crsp.mktcap))
  return crsp
end

function computenetreturns!(crsp::DataFrame; Findex::Symbol = VWCRSP_INDEX)
  Flindex::Symbol = Symbol(:l, Findex)
  #look up the benchmark against an index
  #crsp doesn't handle irregular intervals very well
  mergevwcrspindex!(crsp)
  crsp.lnet = Vector{MFloat64}(undef, size(crsp,1))

  sort!(crsp, [:permno, :date])

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  print("time to check (net ret formation): ")
  @time @mpar for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]

    scrsp[2:end, :lnet] .= scrsp[2:end, :lret] .- (scrsp[2:end, Flindex] .- scrsp[1:(end-1), Flindex])

  end

  #compute the logged versions
  crsp.net = (exp).(crsp.lnet) .- 1.0


  #=if TEST_OUTPUT
    CSV.write("output\\crspnetseries.csv", crsp[(end-10_000):end, :])
    @info "crspnetseries.csv written."
  end=#

  return nothing

end


function preprocesscrsp(crsp::DataFrame, crsptype::Symbol;
  crspdateformat::DateFormat = CRSP_DATE_FORMAT,
  testoutput::Bool = TEST_OUTPUT,
  usedcrspfields::Vector{Symbol} = USED_CRSP_FIELDS)::DataFrame where T<:DatePeriod

  local scrsps::GroupedDataFrame

  filtercrsp!(crsp)



  #Doing this before the dedup eliminates rows without returns,
  #prices and shares outstanding first. Implicitly requires that
  #these fields exist
  computecrspreturns!(crsp)
  computecrspprices!(crsp, crsptype)

  #now dedup crsp
  dedupctr::Int = 0
  crsp.keeps = trues(size(crsp,1)) |> Vector{Bool}

  scrsps = groupby(crsp, [:permno, :date])
  for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    local partialscrsp::SubDataFrame

    (size(scrsp,1) == 1) && continue

    copyto!(scrsp.keeps, scrsp.distcd .== 1232)
    #scrsp.keeps .= (cd->ifelse(cd == 1232, true, false)).(scrsp.distcd)
    (sum(scrsp.keeps) == 1) && continue

    #if we have too many 1232s, we only will select from these
    partscrsp = (sum(scrsp.keeps) > 1) ? view(scrsp, scrsp.keeps, :) : scrsp

    #try to take the minimum distribution code
    nonmissingdistcds::Vector{Int} = collect(skipmissing(partscrsp.distcd))
    if length(nonmissingdistcds) ≠ 0
      copyto!(partscrsp.keeps, partscrsp.distcd .== minimum(nonmissingdistcds))
    end
    (sum(partscrsp.keeps) == 1) && continue

    #again, see if we have too many  of the minimum distcd, or the distcd was missing
    partscrsp = (sum(partscrsp.keeps) > 1) ? view(partscrsp, partscrsp.keeps, :) : partscrsp

    #finally, check if the rows are identical- if so, keep 1, if not, drop the rows
    partscrsp.keeps .= false
    if eachrowequal(partscrsp)
      partscrsp.keeps[1] = true
    end

  end

  delete!(crsp, (!).(crsp.keeps))

  @assert length(groupby(crsp, [:permno, :date]))==size(crsp,1)

  return crsp
end

#applies a function to a vector of trailing variables
function computetrailing!(F::Function,
  crsp::DataFrame, period::DatePeriod,  Ftarget::Symbol, Fnew::Symbol;
  eomonly::Bool,
  offset::Int=0,
  minpointsperwindow::Int = 0)::Nothing

  local getr₀::Function

  N::Int = size(crsp,1)

  #check and make sure the target is available
  (Fnew ∈ propertynames(crsp)) || (error("$Fnew not found in df"))
  (sum((ismissing).(crsp[!,Fnew])) == size(crsp,1)) ||  (error("$Fnew already contains data"))

  #make sure everything is sorted
  (issorted(crsp, [:permno, :date])) || error("dataframe not sorted by [:permno, :date]")

  #when we are looking at eom values only
  @inline getr₀eomonly(d::Date, dates::AbstractVector{Date}) = searchsortedfirst(
    dates, lastdayofmonth(d) + Day(1) - period)
  @inline getr₀all(d::Date, dates::AbstractVector{Date}) = searchsortedfirst(
    dates, d + Day(1) - period)
  getr₀ = eomonly ? getr₀eomonly : getr₀all

  #split and apply by permno
  completecrsp::SubDataFrame = view(crsp, (!ismissing).(crsp[!,Ftarget]),:)
  scrsps::GroupedDataFrame = groupby(completecrsp, :permno)
  lowestr::Int = minpointsperwindow > 0 ? minpointsperwindow : 1+offset

  #now compute the trailing returns
  @mpar for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    Nscrsp = size(scrsp, 1)
    (lowestr ≤ Nscrsp) || continue #if there are enough  entries to compute the trailing returns
    #isnotmissings::Vector{MInt} = (!ismissing).(scrsp[!,Ftarget])

    #a hack to only compute calcualtions at end of month
    for r ∈ lowestr:Nscrsp
      rₜ::Int = r - offset
      eomonly && (!scrsp[rₜ,:eom]) && continue #if we are only computing end of month values

      r₀::Int = getr₀(scrsp[rₜ, :date], scrsp.date)

      #if the returns are not stale and complete enough, record them
      if rₜ - r₀ + 1 ≥ minpointsperwindow
        scrsp[r, Fnew] = F(scrsp[r₀:rₜ , Ftarget])
      end
    end
  end

  return nothing
end

#computes the trailing returns from the logged returns
function computetrailingreturns!(crsp::DataFrame, period::DatePeriod;
  Flret::Symbol = :lret,
  Fnewret::Symbol = Symbol(:t,:ret,replace("$period", " "=>"")),
  offset::Int = 0,
  minpointsperwindow::Int = 0,
  eomonly::Bool = false
  )::Nothing


  local Flnewret::Symbol = Symbol(:l, Fnewret)
  crsp[!,Flnewret] = Vector{MFloat64}(undef, size(crsp,1))
  computetrailing!(v::AbstractVector{MFloat64}->sum(skipmissing(v)),
    crsp, period, Flret, Flnewret,
    eomonly=eomonly,
    offset=offset,
    minpointsperwindow=minpointsperwindow)

  #create the unlogged version
  crsp[!,Fnewret] = (exp).(crsp[!,Flnewret]) .- 1.0

  return nothing
end

#computes the trailing returns from the logged returns
function computetrailingvol!(crsp::DataFrame, period::DatePeriod, Fret::Symbol;
  Fvol::Symbol = Symbol(:t,:vol,replace("$period", " "=>"")),
  offset::Int = 0,
  minpointsperwindow::Int = 0,
  eomonly::Bool = false
  )::Nothing

  crsp[!,Fvol] = Vector{MFloat64}(undef, size(crsp,1))
  computetrailing!((v::AbstractVector{MFloat64})->std(skipmissing(v)),
    crsp, period, Fret, Fvol,
    eomonly=eomonly,
    offset=offset,
    minpointsperwindow=minpointsperwindow)

  if TEST_OUTPUT
    TEST_OUTPUT && CSV.write(
      "output\\vol_p$(period)_off$(offset).csv",
        crsp[1:20_000, :])
  end


  return nothing
end

function computetrailingvol!(crsp::DataFrame, period::DatePeriod,
    Frets::Vector{Symbol};
  Fvols::Vector{Symbol} = error("Frets is required"),
  offset::Int = 0,
  minpointsperwindow::Int = 0,
  eomonly::Bool = false
  )::Nothing

  for Fvol ∈ Fvols
    crsp[!,Fvol] = Vector{MFloat64}(undef, size(crsp,1))
  end

  #quick parallel routine
  tasks = ((Fret, Fvol)->@spawn(computetrailing!((v::AbstractVector{MFloat64})->std(skipmissing(v)),
    crsp, period, Fret, Fvol,
    eomonly=eomonly,
    offset=offset,
    minpointsperwindow=minpointsperwindow))).(Frets,Fvols)

  (wait).(tasks)


  if TEST_OUTPUT
    TEST_OUTPUT && CSV.write(
      "output\\vol_p$(period)_off$(offset).csv",
        crsp[1:20_000, :])
  end

  return nothing
end

function computevolumes!(crsp::DataFrame)::Nothing
  checkmissing(f::Union{Real,Missing}) = (
    ismissing(f) || (f == 99) || (f < 0) ? missing : f)

  crsp.vol .= (checkmissing).(crsp.vol)

  (String <: eltype(crsp.numtrd)) && (crsp.numtrd = (s::MString->parsecrsp(Int, s)).(crsp.numtrd))
  crsp.numtrd .= (checkmissing).(crsp.numtrd)

  return nothing
end

#make tryparse handle missing values
#Base.tryparse(::Type{<:Any}, v::Missing)::Nothing = nothing

#parse a value where the type is known: Generic case
function parsecrsp(::Type{T}, s::String;
    parsedmissings = CRSP_PARSED_MISSINGS) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  (!ismissing(v)) && (v ∈ parsedmissings) && (v=missing)

  return v
end

#the date case
parsecrsp(::Type{Date}, s::String;
    crspdateformat::DateFormat = CRSP_DATE_FORMAT) = Dates.Date(s, crspdateformat)

#helper methods and special cases
parsecrsp(::Type{<:Any}, v::Missing) = missing
parsecrsp(::Type{Date}, i::Int) = parsecrsp(Date, "$i")
