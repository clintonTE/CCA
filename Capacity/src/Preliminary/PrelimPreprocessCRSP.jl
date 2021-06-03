
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

@inline normalizecrspnames!(crsp::DataFrame) = rename!(crsp, (s::Symbol->Symbol(lowercase(string(s)))).(names(crsp)))

#helper to handle default case
#prepcrsp(param::AbstractDict; args...) = prepcrsp(x->x, param; args...)::DataFrame

function prepcrsp(op!::Function = (noop(x)=x);
  crsptype = :not_set,
  refreshcrsp=true,
  incsvstream::Function = IN_CSV_STREAM,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  csvextension::String = CSV_EXTENSION,
  datapath::String = PARAM[:datapath],
  yearrange::Union{AbstractRange, Nothing} = nothing,
  validpermnos::Union{Nothing, Vector{Int}} = nothing,
  prefix::String = "",
  crspcolumns::Vector{Symbol} = error("crsp columns are required"))::DataFrame

  crsptype ∈ [:monthly, :daily] || error("Set crsptype to :monthly or :daily instead of $crsptype")

  local crspname
  local crspparts

  if crsptype == :monthly
    crspname = PARAM[:crspmname]
    crspparts = PARAM[:crspmparts]
  elseif crsptype==:daily
    crspname = PARAM[:crspdname]
    crspparts = PARAM[:crspdparts]
  else
    @assert false
  end


  local p::String = "$datapath\\$(prefix)$crspname"
  local crsp::DataFrame
  local crspdateformat::DateFormat = PARAM[:wrdsdateformat]

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcrsp
    #do the below so we can pre-filter

    pin::String = "$datapath\\$(crspparts[1])"
    crsp = CSV.read(incsvstream("$pin.$csvextension"), threaded=PARAM[:csvthreaded]) |> DataFrame
    normalizecrspnames!(crsp)
    (!isnothing(validpermnos)) && (prefiltercrsp!(crsp, validpermnos))

    if length(crspparts) > 1
      #try to do this is in a space efficient way
      for crsppart ∈ crspparts[2:end]
        pin = "$datapath\\$crsppart"
        local addedcrsp::DataFrame = CSV.read(
          incsvstream("$pin.$CSV_EXTENSION"), threaded=PARAM[:csvthreaded]) |> DataFrame
        normalizecrspnames!(addedcrsp)
         (!isnothing(validpermnos)) && (prefiltercrsp!(addedcrsp, validpermnos))
        append!(crsp, addedcrsp)
      end
    end

    #filter by year first
    crsp.date = (s->parsecrsp(Date, s, crspdateformat)).(crsp.date)
    (!isnothing(yearrange)) && filter!(r::DataFrameRow->year(r.date) ∈ yearrange, crsp) #restrict to desired years

    crsp = preprocesscrsp(crsp, crsptype)
    op!(crsp)
    select!(crsp, crspcolumns)

    outbinstream("$p.$BIN_EXTENSION", crsp)

  else
    crsp = inbinstream("$p.$BIN_EXTENSION")
  end

  @info "CRSP data loaded and/or saved into file $(prefix)$crspname.$BIN_EXTENSION"

  return crsp
end

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
  @mpar for i ∈ 1:length(price)
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
    @mpar for i ∈ 1:length(price)
      if ismissing(price[i]) && (!ismissingorzero(crsp.altprc[i]))
        price[i] = crsp.altprc[i]
      end
    end
  end

  #negative values computed via bid-ask spread method
  price .= (abs).(price)

  return price
end

function computecrspmktcap(crsp::DataFrame)

  #adjust the units and compute the market cap
  mktcap::Vector{MInt} = Vector{MInt}((round).(crsp.price .* crsp.shrout .* 1000))
  mktcap .= (i::MInt->coalesce(iszero(i),true) ? missing : i).(mktcap)

  return mktcap
end

function preprocesscrsp(crsp::DataFrame, crsptype::Symbol)::DataFrame
  crspdateformat::DateFormat = PARAM[:wrdsdateformat]

  local scrsps::GroupedDataFrame

  #NOTE: SHRCD==10,11 : US Common Shares
  #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
  deleterows!(crsp, (!).(completecases(crsp[:,[:shrcd, :exchcd]])))
  local shrcdvalid::Set{Int} = Set{Int}((10,11))
  local exchcdvalid::Set{Int} = Set{Int}((1,2,3))
  filter!(r::DataFrameRow->
    (r.shrcd ∈ shrcdvalid) && (r[:exchcd] ∈ exchcdvalid), crsp)


  #Doing this before the dedup eliminates rows without returns,
  #prices and shares outstanding first. Implicitly requires that
  #these fields exist
  crsp.ret = computecrspreturns(crsp)
  select!(crsp, Not([:dlret]))
  crsp.lret = (log).(1.0 .+ crsp.ret)

  crsp.price = computecrspprices(crsp, crsptype)
  select!(crsp, Not([:prc, :dlprc]))
  crsp.mktcap = computecrspmktcap(crsp)

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


#parse a value where the type is known: Generic case
function parsecrsp(::Type{T}, s::String, crspparsedmissing::AbstractVector{Float64}) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  (!ismissing(v)) && (v ∈ crspparsedmissing) && (v=missing)

  return v
end

#the date case
parsecrsp(::Type{Date}, s::String, crspdateformat::DateFormat) = Dates.Date(s, crspdateformat)

#helper methods and special cases
parsecrsp(::Type{<:Any}, v::Missing, flotsam) = missing
parsecrsp(::Type{Date}, i::Int, crspdateformat::DateFormat) = parsecrsp(Date, "$i", crspdateformat)
