#It would be nice to automate downloading of these files
#also may want some routines to process the metadata and flat files
#does pre-processing of Open990 files

#=
uses constants
const OPEN990_PATH = DATA_PATH * "\\open990"
const OPEN990_FIELD_PATH = OPEN990_PATH * "\\fieldfiles"
const OPEN990_DATE_FORMAT = "yyyymm"
const OPEN990_FIELD_INDEX_NAME = "fieldfileindex"
=#

#loads the Open990 index
function loadfieldfileindex(;path::String = OPEN990_PATH,
  indexname::String = OPEN990_FIELD_INDEX_NAME)

  return CSV.read("$path\\$indexname.csv") |> DataFrame
end

#reads a single file
#takes in the file name, the targeted field name (for grouping)
#and the precedence (higher numbers have greater precedence)
function readfieldfile(filename::String;
  istreamer::Function = I_CSV_STREAMER,
  path::String = OPEN990_FIELD_PATH,
  dateformat::String = OPEN990_DATE_FORMAT,
  #submittedondateformat::String = OPEN990_SUBMITTED__ON_DATE_FORMAT,
  checkunique=true)::DataFrame

  fullpath::String = "$path\\$filename.csv.gz"
  b = istreamer(fullpath) do istream
    b = IOBuffer(read(istream))
    (b)
  end

  try
    df = CSV.read(b, strict=false, header=2,
    types = Dict(:last_updated=>String)) |> DataFrame
  catch err
    println("ERROR: $filename.csv.gz")
    error(err)
  end
  #print(describe(df))

  #now do some light pre-processing

  #strips the month off the fiscal year
  df[:fisyr] = (d::Int->d ÷ 100).(df[:tax_period])
  df[:datefiscal] = (d::Int->Date("$d", dateformat)).(df[:tax_period])

  #properly format the ein
  rename!(df, :ein=>:ein_old)
  df[:ein] = (i::Int->lpad(string(i),9,'0')).(df[:ein_old])

  rename!(df, :submitted_on=>:submitted_on_old)
  df[:submitted_on] = Vector{Date}(df[:submitted_on_old])


  #get unique ropws by taking the newest
  dfmaxsubmitted = by(df, [:ein, :fisyr], submitted_on = :submitted_on => maximum)
  df = join(df, dfmaxsubmitted, on = [:ein, :fisyr, :submitted_on])

  println("Read $(size(df,1)) records from $filename")
  #delete columns which will not be used


  deletecols!(df, setdiff(names(df), [:ein, :value, :fisyr, :submitted_on, :datefiscal]))

  #a spot test of the remaining dups shows ver similar, mostly identical
  #values, so just take the mean
  df = by(df, [:ein, :fisyr, :submitted_on, :datefiscal], value = :value => mean)

  #record the file-level metadata
  df[:filename] = filename


  return df
end

#reads in multiple files
function readfieldfiles(indexdf::DataFrame)::DataFrame

  local futures::Vector{Future} = Vector{Future}()
  local fielddfs::Vector{DataFrame} = Vector{DataFrame}()
  local mappingnames = unique(indexdf[:mappingname])

  for r ∈ eachrow(indexdf)
    push!(futures,
      @spawn readfieldfile(r[:filename]))
  end

  fielddflong::DataFrame = [(fetch).(futures)...;]

  println("Total records read: $(size(fielddflong))")
  fielddf::DataFrame = unstack(fielddflong, :filename, :value)

  #get rid of some dups by taking the newest entries
  dfmax = by(fielddf, [:ein, :fisyr], submitted_on = :submitted_on => maximum)
  fielddf = join(fielddf, dfmax, on = [:ein, :fisyr, :submitted_on])

  #get rid of the remainder by taking the latest fiscal year
  dfmax = by(fielddf, [:ein, :fisyr], datefiscal = :datefiscal=>maximum)
  fielddf = join(fielddf, dfmax, on = [:ein, :fisyr, :datefiscal])

  deletecols!(fielddf, [:datefiscal, :submitted_on])

  #NOTE: I checked and there are a few entries with multiple fisyrs. I elect
  #to select the entries with the newest fiscal year

  @assert size(unique(fielddf[[:ein, :fisyr]]),1) == size(fielddf,1)

  println("Total org-years: $(size(fielddf))")

  #println("uniquesize: $(size(unique(fielddf[[:ein, :fisyr]])))")

  return fielddf
end

function mapfieldfiles!(indexdf::DataFrame, fielddf::DataFrame)::Nothing

  #first build up indices of the mapping

  #println(describe(fielddf))
  mappingfields::Vector{Symbol} = (Symbol).(unique(indexdf[:mappingname]))
  Nmap::Int = length(mappingfields)

  #index that provides teh precedence
  #=precedenceindex::Dict =
    Dict(indexdf[i,:filename]=>indexdf[i, :precedence]
      for i ∈ 1:size(indexdf,1))=#

  #index that gives all of the filenames in order of precedence
  sort!(indexdf, (:mappingname, order(:precedence, rev=true)))
  filenameindex::Dict = Dict{Symbol, Vector{Symbol}}()
  lowprecedencevalues::Vector{Float64} = Vector{Float64}(undef, Nmap)
  for subdf ∈ groupby(indexdf, :mappingname)
    mapfield::Symbol = Symbol(subdf[1,:mappingname])
    filenameindex[mapfield] = (Symbol).(subdf[:filename])

    imapping::Int = findfirst(isequal(mapfield), mappingfields)
    lowprecedencevalues[imapping] = subdf[1,:lowprecedencevalue]
  end

  #allocate the mappingfields. Initially store the highest precedence
  for f ∈ mappingfields
    fielddf[f] = deepcopy(fielddf[filenameindex[f][1]])
  end

  #now do the mapping. Deafult to the first value.
  for r ∈ eachrow(fielddf) #iterate down each row in the new data
    for imapping ∈ 1:Nmap
      mapfield::Symbol = mappingfields[imapping]

      for filename ∈ filenameindex[mapfield][2:end] #iterate in order of precedence

        #if the candidate field is not missing, replace the old field if the old field is missing
        #or its the low precedence value
        if (!ismissing(r[filename])) && (ismissing(r[mapfield]) ||
          (lowprecedencevalues[imapping] == r[mapfield]))

          r[mapfield] = r[filename]
        end
      end
    end
  end

  #CSV.write("$WORKING_PATH\\open990test.csv", fielddf[1:100_000, :])

  #delete the old columns
  for f ∈ mappingfields
    deletecols!(fielddf, filenameindex[f])
  end

  #println(describe(fielddf))
  return nothing
end

function avg(v)::MFloat64
  local f::MFloat64
  if !isempty(skipmissing(v))
    f = mean(skipmissing(v))
    f = isfinite(f) ? f : missing
  else
    f = missing
  end

  return f
end

#does any aggregation functions that are needed and rejoins them to the main dataframe if desired
function aggregatefielddf(fielddf::DataFrame; rejoin::Bool=true,
  functions::Vector = [avg])

  local avgdf::DataFrame = aggregate(
    fielddf[setdiff(names(fielddf), [:fisyr])], :ein, functions)

  if rejoin
    fielddf = join(fielddf, avgdf, on=[:ein], kind=:left)
    return fielddf
  end

  return avgdf
end


#creates the full dataframe if needed
function createfielddf()
  indexdf::DataFrame = loadfieldfileindex()
  fielddf = readfieldfiles(indexdf)
  mapfieldfiles!(indexdf, fielddf)
  fielddf = aggregatefielddf(fielddf)

  return fielddf
end


#spits out the computationally heavy componeny sof the work into a compressed jls
function getfielddf(; path::String = OPEN990_PATH,
  fielddataname::String = OPEN990_FIELD_DATA_NAME,
  refreshopen990data::Bool = true,
  ostreamer::Function = O_JLS_STREAMER,
  istreamer::Function = I_JLS_STREAMER)

  local fielddf::DataFrame

  local fulldatapath::String = "$path\\$fielddataname.jls.gz"

  #recreate the file if none exists or a refresh is called
  if (!isfile(fulldatapath)) || refreshopen990data
    fielddf = createfielddf()

    ostream = ostreamer(fulldatapath)
    serialize(ostream, fielddf)
    close(ostream)
  else
    istream = istreamer(fulldatapath)
    fielddf = deserialize(istream)
    close(istream)
  end

  return fielddf
end

function mergeopen990data!(data::NCCSData; refreshopen990data= true)::Nothing

  fielddf::DataFrame = getfielddf(refreshopen990data=refreshopen990data)
  println(describe(fielddf))
  data.df = join(data.df, fielddf, on=[:ein, :fisyr], kind=:left)

  return nothing
end
