#This file contains a wrapper object and some helper functions
#relating to endowment files



struct NCCSFile
    df::DataFrame
    meta::MetaDictionary

    fileType::Symbol
    fileYear::Int
    fullFile::Bool
end

const NNCCSFile = Union{NCCSFile, Nothing}

#constructs the objects from the metaDictionary
function NCCSFile(meta::MetaDictionary; preProcess::Bool = true,
nccsPath::String = NCCS_PATH, iStreamer::Function = I_CSV_STREAMER)

  inPath::String = "$nccsPath\\$(meta.fileName).csv.gz"

  df::DataFrame = readSingleTypeCSV(inPath)

  #now get the other info we need
  yearMatch::RegexMatch = match(r"\d{4}", meta.fileName)
  fileType::Symbol = Symbol(meta.fileName[(yearMatch.offset+5):(yearMatch.offset+6)])
  fileYear::Int = parse(Int, yearMatch.match)
  fullFile::Bool = contains(meta.fileName, "full")

  data = NCCSFile(df, meta, fileType, fileYear, fullFile)

  preProcess?preProcessNCCSFile!(data):data

  return data
end


  #serializes and compresses an already downloaded NCCS file
function willFail(meta::MetaDictionary; nccsPath::String = NCCS_PATH,
  iStreamer::Function = I_CSV_STREAMER, oStreamer::Function = O_JLS_STREAMER,
  serializeOnly::Bool = true, verbose::Bool = true)::DataFrame

  #these are the column names that we expect
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))

  inPath::String = "$nccsPath\\$(meta.fileName).csv.gz"
  outPath::String = "$nccsPath\\$(meta.fileName).jls.gz"

  if verbose
    println("Serializing $(meta.fileName)...")
  end

  iStream = iStreamer(inPath)
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, rows=2)
  close(iStream)

  #get the column names and make sure each one is in the list of columns we want
  actualCols::Vector{String} = (s::Symbol->lowercase(string(s))).(names(df))

  #form the lookup dictionary for the type. Takes a column name and gets the type, with a default of String
  typeDict::Dict = Dict(
    actualCols[i] => (
    haskey(meta.dictionary, Symbol(actualCols[i])) ? meta.dictionary[Symbol(actualCols[i])]: Union{String, Missing})
    for i::Int ∈ 1:(length(actualCols)))

  #now read in the file
  iStream = iStreamer(inPath)
  #CSV.validate(iStream, header=actualCols, datarow=2, types=typeDict)
  try
    df = CSV.read(iStream, header=actualCols, datarow=2, types=typeDict)
  catch
    error("Try Block worked.")
  end
  close(iStream)

  #remove columns we don't expect or want
  try
    df = df[:,expectedCols]
  catch err
    error("$err\nmeta: $(meta)")
  end

  oStream = oStreamer(outPath)
  serialize(oStream, df)
  close(oStream)

  if verbose
    println("Serialization of $(meta.fileName) complete.")
  end

  if serializeOnly
    df = DataFrame([[nothing]]) #free the memory with an empty dataframe
  end

  return df

end


function changeNCCSColumnType!(data::NCCSFile, col::Symbol, T::Type)::NCCSFile

  if T <: MString #case 1: String or StringUMissing
    return data #do nothing
  end

  rename!(data.df, [col=>:oldCol])
  oldCol::Vector{R where R <: MString} = data.df[:oldCol]
  NValues::Int = length(oldCol)

  if (T <: CategoricalValue{Union{Symbol, Missing}}) #case 3: pool the categorical data
    data.df[:,col] = Vector{Union{Symbol, Missing}}(NValues)
    data.df[:,col] .= ((s::MString->ismissing(s)?missing:Symbol(s)).(oldCol))
    #println("Before: $(typeof(data.df[col]))")
    categorical!(data.df, col)
    #println("After: $(typeof(data.df[col]))")
    data.df[col] = compress(data.df[col])

  elseif T <: Union{MFloat64, MInt, Int, Float64, Bool, MBool} #case 2: concrete numerical type
    #need to parse each entry
    concreteCol::Vector{T} = Vector{T}(NValues)
    if T == Missing #special case where the column is entirely missing
      concreteCol = missings(NValues)
    else
      TSub::Type = typeintersect(T, Union{Float64, Bool, Int}) #get the nonmissable type


      #=@inbounds @simd for i ∈ 1:NValues #parse the column text
        concreteCol[i] = ifelse(ismissing(oldCol[i]),missing, parse("$(oldCol[i])"))
      end=#
      @inbounds @simd for i ∈ 1:NValues
        concreteCol[i] = ismissing(oldCol[i]) ? missing : parse(TSub, oldCol[i])
      end

    end

    data.df[:,col] = concreteCol

  else #Write explicit code for each type of column
    error("Col $col from file $(data.meta.fileName) is of type $T which is not recognized.")
  end

  #delete the old column
  delete!(data.df, :oldCol)

  return data
end



#fixes the type information
function preProcessNCCSFile!(data::NCCSFile)::NCCSFile

  meta::MetaDictionary = data.meta
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))
  actualCols::Vector{Symbol} = (s::Symbol->Symbol(lowercase(string(s)))).(names(data.df))
  names!(data.df, actualCols) #make the names lower case

  #form the lookup dictionary for the type. Takes a column name and gets the type, with a default of String
  typeDict::Dict = Dict(
    actualCols[i] => (
    haskey(meta.dictionary, actualCols[i]) ? meta.dictionary[actualCols[i]]: MString)
    for i::Int ∈ 1:(length(actualCols)))

  #remove columns we don't expect or want
  diff::Vector{Symbol} = setdiff(actualCols, expectedCols)
  if length(diff) > 0
    try
      delete!(data.df, setdiff(actualCols, expectedCols))
    catch err
      error("$err\nmeta: $(meta.fileName)")
    end
  end

  #set the column types. Would have preferred to do this on the read, but got a weird error.
  for col::Symbol ∈ actualCols
    changeNCCSColumnType!(data, col, typeDict[col])
  end

  return data
end


#serializes a NCCSFile object
function serializeNCCSFile(data::NCCSFile; nccsPath::String = NCCS_PATH,
    oStreamer::Function = O_JLS_STREAMER,verbose=true)::NNCCSFile

    outPath::String = "$nccsPath\\$(data.meta.fileName).jls.gz"

    oStream = oStreamer(outPath)
    serialize(oStream, data)
    close(oStream)

    return nothing
end

#wrapper object to serialize a data file object
#For convenience, can return the object if requested
function serializeNCCSFile(meta::MetaDictionary;
  verbose::Bool = true, returnNCCSFile::Bool=true)::NNCCSFile

  #these are the column names that we expect

  if verbose
    println("Serializing $(meta.fileName)...")
  end

  data::NCCSFile = NCCSFile(meta) #make the object
  serializeNCCSFile(data) #serialize the object

  if verbose
    println("Serialization of $(meta.fileName) complete.")
  end

  return returnNCCSFile?data:nothing
end

#deserializes a NCCSFile object
function deserializeNCCSFile(fileName::String; nccsPath::String = NCCS_PATH,
    iStreamer::Function = I_JLS_STREAMER)::NCCSFile

    outPath::String = "$nccsPath\\$(fileName).jls.gz"

    iStream = iStreamer(outPath)
    data::NCCSFile = deserialize(iStream)
    close(iStream)

    return data
end

deserializeNCCSFile(meta::MetaDictionary)::NCCSFile = deserializeNCCSFile(meta.fileName)

#loads an already downloaded NCCS file
function loadNCCSFile(meta::MetaDictionary; nccsPath::String = NCCS_PATH,
  reserializeData::Bool = true, returnNCCSFile::Bool = true)::NNCCSFile

  filePath::String = "$nccsPath\\$(meta.fileName).jls.gz"

  #only re-process the file if necessesary
  if reserializeData || !isfile(filePath)
    data::NNCCSFile = serializeNCCSFile(meta, returnNCCSFile=returnNCCSFile)
  elseif returnNCCSFile
    data = deserializeNCCSFile(meta)
  else #do nothing
    data = nothing #an empty dataframe
  end

  return data
end

#loads and/or serializes a vector of corresponding files
function loadNCCSFiles(metadata::Vector{MetaDictionary}; mode::Symbol = :serial,
  reserializeData::Bool = true,  serializeOnly::Bool = true)::Union{Nothing,Vector{NCCSFile}}

  if mode == :parallel #relaunch the dedicated parallel function
    return loadNCCSFilesPar(metadata::Vector{MetaDictionary}, reserializeData=reserializeData,
      serializeOnly=serializeOnly)
  end

  if serializeOnly
    for meta::MetaDictionary ∈ metadata #do this in a loop format for space conservation
      loadNCCSFile(meta,reserializeData=reserializeData, returnNCCSFile=false)
    end
    return nothing
  else
    return (meta::MetaDictionary->loadNCCSFile(meta,returnNCCSFile=true)).(metadata)
  end
end

#loads many NCCS serial files in parallel
function loadNCCSFilesPar(metadata::Vector{MetaDictionary}; reserializeData::Bool = true,
  serializeOnly::Bool=true)::Union{Nothing, Vector{DataFrame}}

  numFiles::Int = length(metadata)
  pids = workers()
  np::Int = length(pids)

  #divide up the work
  assignments::Vector{Vector{Int}} = ((i::Int)->Vector{Int}()).(1:np)
  for i::Int ∈ 1:numFiles
    push!(assignments[i%np+1], i) #assigns a file to a pid
  end

  #set up a vector for pushing and receiving assignments
  futures::Vector{Future} = Vector{Future}()
  sizehint!(futures, np)

  for i::Int ∈ 1:np
    if length(assignments[i]) > 0
      push!(futures,
        (@spawnat pids[i] loadNCCSFiles(metadata[assignments[i]], reserializeData=reserializeData,
          serializeOnly=serializeOnly)))
    end
  end

  #get the work and return it if desired
  if serializeOnly
    (fetch).(futures)
    return nothing
  else
    return vcat((fetch).(futures)...)
  end
end

#creates a table of summary values
function summarizeNCCSFile(data::NCCSFile)::String
  df::DataFrame = data.df #for convenience since we are just reading, not writing

  fields::Vector{Symbol} = names(df)
  N::Int = size(df,1)
  NFields::Int = length(fields)

  #the names of the columns
  headers::Vector{Symbol} = [:N, :pMissing, :Unique, :Mean, :StdDev, :Skew, :Kurtosis,
    :Min, :p05, :p25, :p50, :p75, :p95, :Max]

  #set up our columns as a dictionary
  columns::OrderedDict = OrderedDict(headers[i] => OrderedDict{String, String}() for i::Int ∈ 1:length(headers))
  (s::Symbol->sizehint!(columns[s], NFields)).(headers)

  fieldStrings::Vector{String} = (s::Symbol->join(split(string(s),"_"),"\\_")).(fields)

  f(x::T where T<:Real) = num2Str(x,2,Ints=true, scaleHurdle = 1000.)

  for i ∈ 1:length(fields)
    #fieldString::String = string(fields[i])

    #first get some type info
    T::Type = data.meta.dictionary[fields[i]]
    TSub::Type = typeintersect(T, Union{Float64, Bool, Int, String, CategoricalValue}) #get the nonmissable type
    TSub = CategoricalValue <: TSub || TSub <: CategoricalValue? Symbol : TSub

    #compute if the type is numeric, then get the column type
    numeric::Bool = (T <: MFloat64) || (T <: MInt)
    VecT = (T <: CategoricalValue) || (CategoricalValue <: T) ? Vector{Union{Symbol,Missing}} : Vector{T}

    #form the rows
    col::VecT = VecT(N)
    col .= VecT(df[fields[i]])
    sansMissing::Vector{TSub} = collect(skipmissing(col))
    if numeric
      push!(columns[:N], fieldStrings[i] => f(length(sansMissing)))
      push!(columns[:pMissing], fieldStrings[i] => f(1.-length(sansMissing)/length(col)))
      push!(columns[:Unique], fieldStrings[i] => f(length(unique(col))))
      push!(columns[:Mean], fieldStrings[i] => f(mean(sansMissing)))
      push!(columns[:StdDev], fieldStrings[i] => f(std(sansMissing)))
      push!(columns[:Skew], fieldStrings[i] => f(skewness(sansMissing)))
      push!(columns[:Kurtosis], fieldStrings[i] => f(kurtosis(sansMissing)))
      push!(columns[:Min], fieldStrings[i] => f(minimum(sansMissing)))
      push!(columns[:p05], fieldStrings[i] => f(percentile(sansMissing,.05)))
      push!(columns[:p25], fieldStrings[i] => f(percentile(sansMissing,.25)))
      push!(columns[:p50], fieldStrings[i] => f(percentile(sansMissing,.50)))
      push!(columns[:p75], fieldStrings[i] => f(percentile(sansMissing,.75)))
      push!(columns[:p95], fieldStrings[i] => f(percentile(sansMissing,.95)))
      push!(columns[:Max], fieldStrings[i] => f(maximum(sansMissing)))
    else
      push!(columns[:N], fieldStrings[i] => f(length(sansMissing)))
      push!(columns[:pMissing], fieldStrings[i] => f(1.-length(sansMissing)/length(col)))
      push!(columns[:Unique], fieldStrings[i] => f(length(unique(col))))
    end
  end

  #convert the dictionaries to table columns
  tableColumns::Vector{TableCol} = (s::Symbol -> TableCol(string(s), columns[s])).(headers)

  table::String = texTable(hcat(tableColumns...), title=data.meta.fileName, titleDescription="Summary Table")

  return table
end

#convenience method since loading many files is memory intensive
summarizeNCCSFile(meta::MetaDictionary) =
  summarizeNCCSFile(loadNCCSFile(meta, reserializeData=false, returnNCCSFile=true))::String

#vectorized method for creating summary tables
function summarizeNCCSFiles(metadata::Vector{MetaDictionary}; mode=:serial)::Vector{String}

  if mode == :parallel #relaunch the dedicated parallel function
    return summarizeNCCSFilesPar(metadata::Vector{MetaDictionary})
  end

  #get the summary tables
  NData::Int = length(metadata)
  summaryTables::Vector{String} = (summarizeNCCSFile).(metadata)

  return summaryTables
end

#helper function to parallelize the formation of the metadata
function summarizeNCCSFilesPar(metadata::Vector{MetaDictionary})::Vector{String}

  numFiles::Int = length(metadata)
  pids = workers()
  np::Int = length(pids)

  #divide up the work
  assignments::Vector{Vector{Int}} = ((i::Int)->Vector{Int}()).(1:np)
  for i::Int ∈ 1:numFiles
    push!(assignments[i%np+1], i) #assigns a file to a pid
  end

  #set up a vector for pushing and receiving assignments
  futures::Vector{Future} = Vector{Future}()
  sizehint!(futures, np)

  for i::Int ∈ 1:np
    if length(assignments[i]) > 0
      push!(futures,
        (@spawnat pids[i] summarizeNCCSFiles(metadata[assignments[i]], mode=:serial)))
    end
  end

  return vcat((fetch).(futures)...) #concatenate the reulsts
end

function writeSummaryNCCSTables(metadata::Vector{MetaDictionary};
  outputPath::String = OUTPUT_PATH, headerName::String = HEADER_NAME, footerName::String = FOOTER_NAME)::Nothing

  tables::Vector{String} = summarizeNCCSFiles(metadata, mode=:parallel)
  writeTables2File(tables, "header.tex", "footer.tex", path=outputPath, outName = "summaryTables.tex")

  return nothing

end
