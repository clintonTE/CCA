
#do some pre-filtering to reduce the processing time
function prefiltercrsp!(crsp::DataFrame, validpermnos::Vector{Int})::Nothing
  #rename columns to lower case


  crsp.tokeep = trues(size(crsp,1))

  crspgrp::GroupedDataFrame = groupby(crsp, :permno)
  Threads.@threads for i ∈ 1:length(crspgrp)
    scrsp::SubDataFrame = crspgrp[i]

    if scrsp.permno[1] ∉ validpermnos
      scrsp[:,:tokeep] .= false
    end
  end

  filter!(r::DataFrameRow->r.tokeep, crsp)
  select!(crsp, Not([:tokeep]))

  return nothing

end

@inline normalizecrspnames!(crsp::DataFrame) = rename!(crsp, (s::String->Symbol(lowercase(s))).(names(crsp)))

#helper to handle default case
#prepcrsp(param::AbstractDict; args...) = prepcrsp(x->x, param; args...)::DataFrame

function prepcrsp(op::Function = (noop(x)=x);
  crsptype=error("crsptype is required"),
  refreshcrsp=true,
  incsvstream::Function = IN_CSV_STREAM,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  csvextension::String = CSV_EXTENSION,
  binextension::String = BIN_EXTENSION,
  crspname::String = string(crsptype),
  crsppath::String = PARAM[:crsppath],
  workingpath::String = PARAM[:workingpath],
  yearrange::Union{AbstractRange, Nothing} = nothing,
  validpermnos::Union{Nothing, Vector{Int}} = nothing,
  binprefix::String = "",
  crspnameout::String = "$binprefix$crspname",
  crspcolumns::Union{Nothing, Vector{Symbol}} = nothing,
  crsprequiredcolumns::Union{Nothing, Vector{Symbol}} = nothing,
  #usestdsharesexchanges::Bool = PARAM[:crspusestdsharesexchanges],
  makebin::Bool = true)::DataFrame


  local binpath::String = "$workingpath\\$crspnameout.$binextension"
  local csvpath::String = "$crsppath\\$crspname.$csvextension"
  local crsp::DataFrame
  local crspdateformat::DateFormat = PARAM[:wrdsdateformat]

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcrsp
    #do the below so we can pre-filter
    crsp = incsvstream(csvpath) |> CSV.File |> f->DataFrame(f,copycols=false)
    normalizecrspnames!(crsp)
    (!isnothing(validpermnos)) && (prefiltercrsp!(crsp, validpermnos))

    #filter by year first
    crsp.date = (s->parsecrsp(Date, s, crspdateformat)).(crsp.date)
    (!isnothing(yearrange)) && filter!(r::DataFrameRow->year(r.date) ∈ yearrange, crsp) #restrict to desired years

    crsp = preprocesscrsp(crsp, crsptype)

    if !isnothing(crsprequiredcolumns)
      crsp = crsp[completecases(crsp, crsprequiredcolumns),:]
    end

    for Fposcol ∈ [:price, :shares, :mc]
      if !(all(crsp[!, Fposcol] .≥ 0.0))
        display(crsp[Fposcol .≤ 0.0,[:permno, :date, :price, :shares, :mc]])
        throw("Found negative values for $Fposcol in crsp")
      end
    end


    crsp = op(crsp)
    #throw("crsp colummns: $(propertynames(crsp))")

    if !isnothing(crspcolumns)
      select!(crsp, crspcolumns)
    end


    makebin && outbinstream(binpath, crsp)

  else
    (!makebin) && (@warn "makebin = false and loading a bin! If this works, beware of stale data....")
    crsp = inbinstream(binpath)
  end

  @info "CRSP data loaded and/or saved into file $crspnameout.$BIN_EXTENSION"

  return crsp
end

#function
#creates crsp returns
function computecrspreturns(crsp::DataFrame)
  crspparsedmissing::Vector{Float64} = PARAM[:crspparsedmissing]
  retlisted::Vector{MFloat64} = (s::MString->parsecrsp(Float64, s, crspparsedmissing)).(crsp.ret)
  dlret::Vector{MFloat64} = (s::MString->parsecrsp(Float64, s, crspparsedmissing)).(crsp.dlret)

  retlisted = (r->
    !ismissing(r) && (r > -1.0) ? r : missing).(retlisted)
  dlret .= (r->
    !ismissing(r) && (r > -1.0) ? r : missing).(dlret)

  ret::Vector{MFloat64} = ((retlisted,dlret)->
    ismissing(retlisted) && ismissing(dlret) ? missing :
      ismissing(dlret) ? retlisted :
      ismissing(retlisted) ? dlret :
      (1+retlisted)*(1+dlret) - 1.
    ).(retlisted,dlret)

  return ret
end

#computes a cumulative index for each stock
function computecrsptotalreturnindex(crsp::AbstractDataFrame; Fret::Symbol=:ret)
  @assert (crsp[!, Fret] .≥ -1.0) |> all
  @assert (crsp[!, Fret] .!== missing) |> all
  @assert issorted(crsp, [:permno, :date])

  rets = crsp[!, [:permno, :date, Fret]] |> deepcopy
  rets.retP1 = rets[!, Fret] .+ 1

  transform!(groupby(rets, :permno), :retP1 => cumprod => :tridx)

  return rets.tridx
end


function computecrspprices(crsp::DataFrame, crsptype::Symbol)
  #Now work on the price
  #in this case, use the mean of the two values
  #price cannot be zero
  crsp.prc .= (
    prc->coalesce(iszero(prc),true) ? missing : prc).(crsp.prc)
  crsp.dlprc .= (
    dlprc->coalesce(iszero(dlprc),true) ? missing : dlprc).(crsp.dlprc)

  @inline ismissingorzero(p::Float64) = iszero(p) ? true : false
  @inline ismissingorzero(::Missing) = true

  price = Vector{MFloat64}(undef, size(crsp,1))
  Threads.@threads for i ∈ 1:length(price)
    if ismissingorzero(crsp.dlprc[i])
      price[i] = crsp.prc[i]
    elseif ismissingorzero(crsp.prc[i])
      price[i] = crsp.dlprc[i]
    else
      price[i] = (crsp.prc[i]+crsp.dlprc[i])/2
    end
  end

  #add in the alternative price if nothing else is available
  #(data available monthly only)
  if crsptype == :monthly
  Threads.@threads for i ∈ 1:length(price)
      if ismissing(price[i]) && (!ismissingorzero(crsp.altprc[i]))
        price[i] = crsp.altprc[i]
      end
    end
  end

  #negative values computed via bid-ask spread method
  price .= (abs).(price)

  return price
end


function filtercrsp!(crsp::DataFrame;)

    #NOTE: SHRCD==10,11 : US Common Shares
    #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
    delete!(crsp, (!).(completecases(crsp[:,[:shrcd, :exchcd]])))
    local shrcdvalid::Set{Int} = Set{Int}((10,11))
    local exchcdvalid::Set{Int} = Set{Int}((1,2,3))
    filter!(r::DataFrameRow->
      (r.shrcd ∈ shrcdvalid) && (r[:exchcd] ∈ exchcdvalid), crsp)


    return crsp
end


computecrspshares(crsp) = crsp.shrout .* 1000.
function computecrspmarketcap(crsp::DataFrame)

  #adjust the units and compute the market cap
  mc::Vector{MInt} = Vector{MInt}((round).(crsp.price .* crsp.shares))
  mc .= (i::MInt->coalesce(iszero(i),true) ? missing : i).(mc)

  return mc
end

function computerangemeasures!(crsp::DataFrame; disttol::Float64=10^-4,
    periods2stale::Day = Day(30))

  issorted(crsp, [:permno, :date]) || error("crsp must be sorted on permno, date")
  posormissing(::Missing) = missing
  posormissing(x::Real) = x>0.0 ? x : missing

  #weird crsp thing that bidlo is a bid under certain circumstances, represented by a neg val
  crsp.lo = (posormissing).(crsp.bidlo)
  crsp.hi = (posormissing).(crsp.askhi)

  lagwithin2!(crsp, [:price], :permno, date=:date, maxnotstale = periods2stale,
    laggedvals = [:Ldistadjclose])

  #adjust the lagged close prices for distributions
  #note a distribuition here is defined as anything that causes the price return
  #to deveiate from the crsp ret, so it can be negative (consider splits)
  priceret::Vector{MFloat64} = crsp.price ./ crsp.Ldistadjclose .- 1.
  dist::Vector{MFloat64} = (crsp.ret .- priceret) .* crsp.Ldistadjclose
  dist[(d->((!ismissing(d)) && (abs(d)<disttol))).(dist)] .= 0.0 #assume such small changes are rounding error
  crsp.Ldistadjclose .-= dist
  crsp.Ldistadjclose .= (posormissing).(crsp.Ldistadjclose)

end



function computedollarvol(crsp::DataFrame)
  #use this to compute dollar volume- take the price midpoint if available
  #println(names(crsp))

  dprice::Vector{MFloat64} = ((open, close)->
    ifelse(open !== missing, (open + close)/2, close)).(crsp.open, crsp.price)

  parsevol(::Missing) = missing
  parsevol(f::Real) = ifelse(f > 0.0, f, missing)

  @assert  (v->(v===missing) || (v≥0)).(crsp.vol) |> all
  vol::Vector{MFloat64} = (parsevol).(crsp.vol)


  return dprice .* crsp.vol
end

#normalizes the dollar volume by the point in time total dollar volume
function computenvol!(crsp::DataFrame, Fdvol::Symbol)

  (:nvol ∈ names(crsp)) && (sum((!ismissing).(crsp.nvol)) > 0) && error(
    ":nvol already contains data")

  crsp.nvol =  Vector{MFloat64}(undef, size(crsp,1))

  #use this to compute dollar volume- take the price midpoint if available
  scrsps::GroupedDataFrame = groupby(crsp, :date)
  for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]
    scrsp.nvol .= scrsp[!, Fdvol]./sum(skipmissing(scrsp[!, Fdvol]))
  end

  return nothing
end



#should probably make this based on time interval
function trailingsigma!(crspd::DataFrame; Fendpoint::Symbol = "Fendpoint is required")


  (:sigma ∉ names(crspd)) || (error(":sigma already in crspd"))
  crspd.sigma = Vector{MFloat64}(undef, size(crspd,1))

  trailingmeasure!(std, crspd, PARAM[:sigmacalendarinterval],PARAM[:sigmaminpoints],
    Ftarget=:ret,
    Fmeasure=PARAM[:sigmaname],
    Fconditional=Fendpoint)::Nothing
end


function preprocesscrsp(crsp::DataFrame, crsptype::Symbol,
  crspdateformat::DateFormat = PARAM[:wrdsdateformat])

  local scrsps::GroupedDataFrame

  #NOTE: SHRCD==10,11 : US Common Shares
  #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
  local shrcdvalid::Set{Int} = Set{Int}((10,11))
  local exchcdvalid::Set{Int} = Set{Int}((1,2,3))

  #=crsp.tokeep = trues(size(crsp,1))
  Threads.@threads for r ∈ eachrow(crsp)
    r.tokeep = ((r.shrcd ∈ shrcdvalid) &&
      (r.exchcd ∈ exchcdvalid) &&
      (!ismissing(r.shrcd)) &&
      (!ismissing(r.exchcd)))
  end
  crsp = crsp[crsp.tokeep,:]=#
  delete!(crsp, (!).(completecases(crsp[:,[:shrcd, :exchcd]])))
  filter!(r::DataFrameRow->(r.shrcd ∈ shrcdvalid) && (r.exchcd ∈ exchcdvalid), crsp)

  #Doing this before the dedup eliminates rows without returns,
  #prices and shares outstanding first. Implicitly requires that
  #these fields exist
  crsp.ret = computecrspreturns(crsp)
  crsp.price = computecrspprices(crsp, crsptype)

  #require valid returns and prices
  crsp = crsp[completecases(crsp, [:ret,:price]),:]

  #following FF UMD, force price and return availability before computing TR
  crsp.tridx = computecrsptotalreturnindex(crsp)

  crsp.lret = (log).(1.0 .+ crsp.ret) #make log returns
  select!(crsp, Not([:prc, :dlprc, :dlret])) #cleanup

  #now dedup crsp
  dedupctr::Int = 0
  crsp.keep = trues(size(crsp,1))

  scrsps = groupby(crsp, [:permno, :date])
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]
    local partialscrsp::SubDataFrame

    (size(scrsp,1) == 1) && continue
    scrsp.keep .= scrsp.distcd .== 1232 #sometimes dups on the distribution type- prefer to keep ord dividends
    (sum(scrsp.keep) == 1) && continue

    #if we have too many 1232s, we only will select from these
    partscrsp = (sum(scrsp.keep) > 1) ? view(scrsp, scrsp.keep, :) : scrsp

    #try to take the minimum distribution code
    nonmissingdistcds::Vector{Int} = collect(skipmissing(partscrsp.distcd))
    if length(nonmissingdistcds) ≠ 0
      partscrsp.keep .= (partscrsp.distcd .== minimum(nonmissingdistcds))
    end
    (sum(partscrsp.keep) == 1) && continue

    #again, see if we have too many  of the minimum distcd, or the distcd was missing
    partscrsp = (sum(partscrsp.keep) > 1) ? view(partscrsp, partscrsp.keep, :) : partscrsp

    #finally, check if the rows are identical- if so, keep 1, if not, drop the rows
    partscrsp.keep .= false
    if eachrowequal(partscrsp)
      partscrsp.keep[1] = true
    end

  end

  delete!(crsp, (!).(crsp.keep))

  #compute market cap here since it is used in selection criteria
  crsp.open = (abs).(crsp.openprc)
  crsp.shares = computecrspshares(crsp)
  crsp.mc = computecrspmarketcap(crsp)


  #=print("time assertion: ")
  @time =#@assert (!any(nonunique(crsp[!,[:permno, :date]])))

  return crsp
end


#parse a value where the type is known: Generic case
function parsecrsp(::Type{T}, s::String, crspparsedmissing::AbstractVector{Float64}) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  (!ismissing(v)) && (v ∈ crspparsedmissing) && (v=missing)

  return v
end


#the date case
function parsecrsp(::Type{Date}, s::String, crspdateformat::DateFormat)
  #try
    Dates.Date(s, crspdateformat)
  #catch err
  #  println("could not parse $s with date format $crspdateformat")
  #  error(err)
  #end
end

#helper methods and special cases
parsecrsp(::Type{<:Any}, ::Missing) = missing
parsecrsp(::Type{<:Any}, ::Missing, ::Any) = missing
parsecrsp(::Type{Date}, i::Int, crspdateformat::DateFormat) = parsecrsp(Date, "$i", crspdateformat)
parsecrsp(::Type{Int}, s::String) = something(tryparse(Int, s), missing)
parsecrsp(::Type{Int}, i::Int) = i
