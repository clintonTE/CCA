
#MAIN ENTRY POINT
#final columns for the decompressed file
prepcrsp(;refreshcrsp=true)::DataFrame = loadcrspdata(refreshcrsp=refreshcrsp)

function loadcrspdata(crspname=CRSP_NAME; datapath=DATA_PATH, refreshcrsp=true,
  incsvstream::Function = IN_CSV_STREAM,
  injlsstream::Function = IN_JLS_STREAM, outjlsstream::Function = OUT_JLS_STREAM)::DataFrame
  local p::String = "$datapath\\$crspname"
  local crsp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcrsp || (!isfile("$p.jls.lz4"))
    crsp = incsvstream("$p.csv.lz4") |> CSV.read |> DataFrame

    #println(describe(crsp))

    #preprocesscrsp!(crsp)

    crsp = preprocesscrsp(crsp)

    outjlsstream("$p.jls.lz4") do s
      serialize(s, crsp)
    end

  else
    injlsstream("$p.jls.lz4") do s
      crsp = deserialize(s)
    end

  end

  @info "CRSP data loaded and saved into file $crspname.jls.lz4."

  return crsp
end

function preprocesscrsp(crsp::DataFrame;
  crspdateformat::DateFormat = CRSP_DATE_FORMAT,
  crspcolumns::Vector{Symbol} = RETAINED_COLS_CRSP,
  crsptype::Symbol = CRSP_TYPE,
  testoutput::Bool = TEST_OUTPUT,
  usedcrspfields::Vector{Symbol} = USED_CRSP_FIELDS,
  time2stale::T = TIME_2_STALE)::DataFrame where T<:DatePeriod

  local scrsps::GroupedDataFrame

  #rename columns to lower case

  names!(crsp, (s::Symbol->Symbol(lowercase(string(s)))).(names(crsp)))

  #parse some columns
  crsp.ret = (s::MString->parsecrsp(Float64, s)).(crsp.ret)
  crsp.dlret= (s::MString->parsecrsp(Float64, s)).(crsp.dlret)
  crsp.siccd = (s::MString->parsecrsp(Int, s)).(crsp.siccd)


  rename!(crsp, :ret=>:retlisted)
  #make sure we capture delisting returns

  #treat -100% returns as missing values
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

  #drop rows with missing data
  #below stopped working
  #deleterows!(crsp, (!).(completecases(crsp[:,[:ret, :SHRCD, :EXCHCD, :SICCD]])))

  crsp = crsp[completecases(crsp[:,[:ret, :shrcd, :exchcd, :siccd]]),:]

  #NOTE: SHRCD==10,11 : US Common Shares
  #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
  #NOTE: r[:SICCD] ∉ 6000:6999 : No financial companies
  #this improves performance
  local shrcdvalid::Set{Int} = Set{Int}((10,11))
  local exchcdvalid::Set{Int} = Set{Int}((1,2,3))
  local siccdinvalid::Set{Int} = Set{Int}(6000:6999)
  filter!(r::DataFrameRow->
    (r[:shrcd] ∈ shrcdvalid) &&
    (r[:exchcd] ∈ exchcdvalid ) &&
    (r[:siccd] ∉ siccdinvalid),
    crsp)

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

  #adjust the units and compute the market cap
  crsp.mktcap = Vector{MInt}((round).(crsp.price .* crsp.shrout .* 1000))
  crsp.mktcap .= (i::MInt->coalesce(iszero(i),true) ? missing : i).(crsp.mktcap)


  #println(describe(crsp))
  crsp.date = (s->parsecrsp(Date, s)).(crsp.date)
  #println(describe(crsp))

  #get log returns
  crsp.lret = (log).(1.0 .+ crsp.ret)

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

  crsp = crsp[crsp.keep,:]



  #need to lag the marketcap
  crsp.Lmktcapdate = Vector{MDate}(undef, size(crsp,1))
  crsp.Lmktcap = Vector{MInt}(undef, size(crsp,1))

  sort!(crsp, [:permno, :date])
  scrsps = groupby(crsp, :permno)
  @mpar for i ∈ 1:length(scrsps)
    local scrsp::SubDataFrame = scrsps[i]

    Nscrsp::Int = size(scrsp,1)
    for i ∈ 2:Nscrsp
      local Lmktcap::MInt
      local Lmktcapdate::MDate

      #get the most recent market cap other than the current
      if !ismissing(scrsp[i-1, :mktcap])
        Lmktcap = scrsp[i-1, :mktcap]
        Lmktcapdate = scrsp[i-1, :date]
      else
        Lmktcap = scrsp[i-1, :Lmktcap]
        Lmktcapdate = scrsp[i-1, :Lmktcapdate]
      end

      #this si the check against stale data
      if (!ismissing(Lmktcapdate)) && (Lmktcapdate + time2stale ≥ scrsp[i,:date])
        scrsp[i, :Lmktcap] = Lmktcap
        scrsp[i, :Lmktcapdate] = Lmktcapdate
      end
    end
  end

  crsp.lLmktcap = (log).(crsp.Lmktcap)
  select!(crsp, crspcolumns)

  if testoutput
    d = Date(2015,12,31)
    testoutput && CSV.write("output\\crsp_$d.csv", crsp[(cd->cd==d).(crsp.date), :])
  end

  return crsp
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
