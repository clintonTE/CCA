
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

    println(describe(crsp))

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
  #println(describe(crsp))

  println("hello15")

  return crsp
end

function preprocesscrsp(crsp::DataFrame;
  crspdateformat::DateFormat = CRSP_DATE_FORMAT,
  crspcolumns::Vector{Symbol} = CRSP_COLUMNS)::DataFrame

  #rename columns to lower case

  names!(crsp, (s::Symbol->Symbol(lowercase(string(s)))).(names(crsp)))

  #parse some columns
  crsp.ret = (s::MString->parsecrsp(Float64, s)).(crsp.ret)
  crsp.dlret= (s::MString->parsecrsp(Float64, s)).(crsp.dlret)
  crsp.siccd = (s::MString->parsecrsp(Int, s)).(crsp.siccd)

  rename!(crsp, :ret=>:retlisted)
  #make sure we capture delisting returns
  crsp.ret = ((retlisted,dlret)->
    ismissing(retlisted) && ismissing(dlret) ? missing :
      ismissing(retlisted) ? dlret :
      ismissing(dlret) ? retlisted :
      (1+retlisted)*(1+dlret) - 1.
    ).(crsp.retlisted,crsp.dlret)

  #drop rows with missing data
  #below stopped working
  #deleterows!(crsp, (!).(completecases(crsp[:,[:ret, :SHRCD, :EXCHCD, :SICCD]])))
  crsp = crsp[completecases(crsp[:,[:ret, :shrcd, :exchcd, :siccd]]),:]

  #NOTE: SHRCD==10,11 : US Common Shares
  #NOTE:  :EXCHCD ∈ [1,2,3] : NYSE, Amex, NASDAQ
  #NOTE: r[:SICCD] ∉ 6000:6999 : No financial companies
  crsp[(r::DataFrameRow->
    (r[:shrcd] ∈ (10,11)) &&
    (r[:exchcd] ∈ (1,2,3)) &&
    (r[:siccd] ∉ 6000:6999)).(eachrow(crsp)),:]

  #println(describe(crsp))
  crsp.date = (s->parsecrsp(Date, s)).(crsp.date)
  select!(crsp, Not(setdiff(names(crsp), crspcolumns)))

  #get log returns
  crsp.lret = (log).(1.0 .+ crsp.ret)

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
