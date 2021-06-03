


function getvalidpermnos(df::DataFrame)::Vector{Int}
  permnosym = :permno ∈ names(df) ? :permno : :lpermno
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
  names!(crsp, (s::Symbol->Symbol(lowercase(string(s)))).(names(crsp)))

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

function prepcrsp(op!::Function = (noop(x)=x);
  crsptype = :not_set,
  refreshcrsp=true,
  datapath=DATA_PATH,
  incsvstream::Function = IN_CSV_STREAM,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  validpermnos::Union{Nothing, Vector{Int}} = getvalidpermnos(),
  prefix::String = "",
  crspcolumns::Union{Nothing,Vector{Symbol}} = [:permno, :date, :ret, :lret],
  testoutput::Bool = TEST_OUTPUT,
  archive::Bool=true)::DataFrame

  crsptype ∈ [:monthly, :daily] || error("Set crsptype to :monthly or :daily instead of $crsptype")


  local crspname = crsptype == :monthly ? CRSPM_NAME :  CRSPD_NAME
  local crspparts = crsptype == :monthly ? CRSPM_PARTS :  CRSPD_PARTS

  local p::String = "$datapath\\$(prefix)$crspname"
  local crsp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcrsp
    #do the below so we can pre-filter

    pin::String = "$datapath\\$(crspparts[1])"
    crsp = incsvstream("$pin.csv.lz4") |> CSV.read |> DataFrame
    prefiltercrsp!(crsp, validpermnos)

    if length(crspparts) > 1
      #try to do this is in a space efficient way
      for crsppart ∈ crspparts[2:end]
        pin = "$datapath\\$crsppart"
        local addedcrsp::DataFrame = incsvstream("$pin.csv.lz4") |> CSV.read |> DataFrame
         prefiltercrsp!(addedcrsp, validpermnos)
        append!(crsp, addedcrsp)
      end
    end

    crsp = preprocesscrsp(crsp, crsptype)
    op!(crsp)
    select!(crsp, crspcolumns)

    if testoutput
      d = Date(2015,12,31)
      testoutput && CSV.write("output\\$(prefix)_$(crspname)_$d.csv", crsp[(cd->cd==d).(crsp.date), :])
    end

    outjlsstream("$p.jls.lz4") do s
      serialize(s, crsp)
    end

    archive && archivefile("$(prefix)$crspname.jls.lz4") #never hurts to have a backup

  else
    try
      injlsstream("$p.jls.lz4") do s
        crsp = deserialize(s)
      end
    catch err
      @warn "Could not load $(prefix)$crspname. Error: $err"
      @info "Attempting to restore back-up and load from archive..."
      unarchivefile("$(prefix)$crspname.jls.lz4")
      injlsstream("$p.jls.lz4") do s
        crsp = deserialize(s)
      end
      @info "Successfully restored from archive."
    end
  end

  @info "CRSP data loaded and/or saved into file $(prefix)$crspname.jls.lz4."

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
  deleterows!(crsp, (ismissing).(crsp.ret))

  return crsp
end

function filtercrsp!(crsp::DataFrame; filternonfinancials::Bool = false)

    #NOTE: SHRCD==10,11 : US Common Shares
    #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
    deleterows!(crsp, (!).(completecases(crsp[:,[:shrcd, :exchcd]])))
    local shrcdvalid::Set{Int} = Set{Int}((10,11))
    local exchcdvalid::Set{Int} = Set{Int}((1,2,3))
    filter!(r::DataFrameRow->
      (r.shrcd ∈ shrcdvalid) && (r[:exchcd] ∈ exchcdvalid), crsp)

    #NOTE: r[:SICCD] ∉ 6000:6999 : No financial companies
    if filternonfinancials
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
  deleterows!(crsp, (ismissing).(crsp.price))

  return crsp
end

function computecrspmktcap!(crsp::DataFrame)

  #adjust the units and compute the market cap
  crsp.mktcap = Vector{MInt}((round).(crsp.price .* crsp.shrout .* 1000))
  crsp.mktcap .= (i::MInt->coalesce(iszero(i),true) ? missing : i).(crsp.mktcap)

  deleterows!(crsp, (ismissing).(crsp.mktcap))
  return crsp
end

#=function computenetreturnsOLD!(crsp::DataFrame)

  crsp.net = Vector{MFloat64}(undef, size(crsp,1))

  alldates::Vector{Date} = sort(unique(crsp.date))
  Ldates::Dict = Dict((i->alldates[i]=>alldates[i-1]).(2:length(alldates)))
  sort!(crsp, [:permno, :date])

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  print("time to check: ")
  @time @mpar for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]
    for (j,r) ∈ enumerate(eachrow(scrsp))
      (j==1) && continue
      Ldate::Date = Ldates[r.date]
      (Ldate == scrsp.date[j-1]) && (r.net = r.ret - r.ewretd)
    end
  end

  #compute the logged versions
  crsp.net .= (f-> (!ismissing(f)) && (f ≥ -1.0) ? f : missing).(crsp.net)
  crsp.lnet = (log).(1.0 .+ crsp.net)

  return nothing

end=#

function computenetreturns!(crsp::DataFrame)

  #look up the benchmark against an index
  #crsp doesn't handle irregular intervals very well
  mergeewindd!(crsp)
  crsp.lnet = Vector{MFloat64}(undef, size(crsp,1))

  sort!(crsp, [:permno, :date])

  scrsps::GroupedDataFrame = groupby(crsp, :permno)
  print("time to check (net ret formation): ")
  @time @mpar for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]

    scrsp[2:end, :lnet] .= scrsp[2:end, :lret] .- (scrsp[2:end, :lewindd] .- scrsp[1:(end-1), :lewindd])

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
  usedcrspfields::Vector{Symbol} = USED_CRSP_FIELDS,
  months2stale::Int = MONTHS_2_STALE_SHORT)::DataFrame where T<:DatePeriod

  local scrsps::GroupedDataFrame


  crsp.date = (s->parsecrsp(Date, s)).(crsp.date)
  filtercrsp!(crsp)

  #Doing this before the dedup eliminates rows without returns,
  #prices and shares outstanding first. Implicitly requires that
  #these fields exist
  computecrspreturns!(crsp)
  computecrspprices!(crsp, crsptype)

  #now dedup crsp
  dedupctr::Int = 0
  crsp.keep = trues(size(crsp,1))

  scrsps = groupby(crsp, [:permno, :date])
  @mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    if size(scrsp,1) > 1
      #first try to pick records by number missing
      nmissing::Vector{Int} = zeros(Int, size(scrsp,1))
      for f ∈ usedcrspfields
        nmissing .+= (ismissing).(scrsp[!,f])
      end
      minmissing::Int = minimum(nmissing)
      scrsp.keep .= nmissing .== minmissing

      #keep the top record iff it is identical to the others
      if sum(scrsp.keep) > 1
        sscrsp = view(scrsp, scrsp.keep, :)
        nchecks::Int = size(sscrsp,1)
        keep1::Bool = true
        for i ∈ 2:nchecks
          for f ∈ usedcrspfields
            #checks if the given used field is equal to the top row
            (coalesce(sscrsp[i,f],-99.0) ≠ coalesce(sscrsp[1,f], -99.0)) && (keep1 = false)
          end
        end

        sscrsp.keep .= false
        sscrsp.keep[1] = keep1
      end
    end
  end

  deleterows!(crsp, (!).(crsp.keep))

  return crsp
end


#applies a function to a vector of trailing variables
function computetrailing!(F::Function,
  crsp::DataFrame, periods::Int,  Ftarget::Symbol,   Fnew::Symbol;
  offset::Int=0, months2stale::Int = MONTHS_2_STALE_LONG,
  sorted::Bool = false)::Nothing

  periods2stale = Month(months2stale)

  N::Int = size(crsp,1)

  #allocate spcae for the new fields
  crsp[!,Fnew] = Vector{MFloat64}(undef, N)

  #need to sort
  (!sorted) && (sort!(crsp, [:permno, :date]))

  #avoid missing data in the returns for simplicity
  scrspnomissing::SubDataFrame = view(crsp, (!ismissing).(crsp[!,Ftarget]),:)
  scrsps::GroupedDataFrame = groupby(scrspnomissing, :permno)

  #now compute the trailing returns
  @mpar for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    Nscrsp = size(scrsp, 1)
    if periods ≤ Nscrsp #if there are enough  entries to compute the trailing returns
      for r ∈ periods:Nscrsp
        r₀::Int = r-periods+1-offset #this si the first row of the period
        rₜ::Int = r - offset
        staledate::Date = firstdayofmonth(scrsp[r,:date]) - periods2stale
        mindate::Date = scrsp[r₀,:date]

        #if the returns are not stale, record them
        if staledate ≤ mindate
          scrsp[r, Fnew] = F(scrsp[r₀:rₜ , Ftarget])
        end
      end
    end
  end

  return nothing
end

#computes the trailing returns from the logged returns
function computetrailingreturns!(crsp::DataFrame, periods::Int;
  months2stale = MONTHS_2_STALE_LONG,
  Flret::Symbol = :lret,
  Fnewret::Symbol = Symbol(:t,:ret,period),
  offset::Int = 0,
  sorted::Bool = false
  )::Nothing


  local Flnewret::Symbol = Symbol(:l, Fnewret)
  computetrailing!(sum, crsp, periods, Flret, Flnewret,
    months2stale=months2stale, offset=offset, sorted=sorted)

  #create the unlogged version
  crsp[!,Fnewret] = (exp).(crsp[!,Flnewret]) .- 1.0

  return nothing
end

#computes the trailing returns from the logged returns
function computetrailingvol!(crsp::DataFrame, periods::Int;
  months2stale = MONTHS_2_STALE_LONG,
  Fret::Symbol = :ret,
  Fvol::Symbol = Symbol(:t,:vol,periods),
  offset::Int = 0,
  sorted::Bool = false
  )::Nothing

  computetrailing!(std, crsp, periods, Fret, Fvol, months2stale=months2stale, offset=offset)

  if TEST_OUTPUT
    TEST_OUTPUT && CSV.write(
      "output\\vol_p$(periods)_m2s$(months2stale)_off$(offset).csv",
        crsp[(end-20000):end, :])
  end

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
