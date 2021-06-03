#This file contains a wrapper object and some helper functions
#relating to endowment files



mutable struct NCCSFile
    df::DataFrame
    meta::MetaDictionary
end

const NNCCSFile = Union{NCCSFile, Nothing}

#constructs the objects from the metaDictionary
function NCCSFile(meta::MetaDictionary; nccsPath::String = NCCS_PATH,
  iStreamer::Function = I_CSV_STREAMER)

  inPath::String = "$nccsPath\\$(meta.fileName).csv.gz"

  df::DataFrame = readSingleTypeCSV(inPath)

  #check and make sure we have the right file
  yearMatch::RegexMatch = match(r"\d{4}", meta.fileName)
  nonprofitType::Symbol = Symbol(meta.fileName[(yearMatch.offset+5):(yearMatch.offset+6)])
  fileYear::Int = parse(Int, yearMatch.match)
  fullFile::Bool = occursin("full", meta.fileName)

  if nonprofitType ≠ meta.nonprofitType || fileYear ≠ meta.fileYear || fullFile ≠ meta.fullFile
    @warn ("File type information does not match for $(meta.fileName).\n
      nonprofitType: $(nonprofitType) | meta.nonprofitType: $(meta.nonprofitType)\n
      fileYear: $(fileYear)  | meta.fileYear: $(fileYear)\n
      fullFile: $(fullFile) | meta.fullFile: $(fullFile)")
  end

  data::NCCSFile = NCCSFile(df, meta)
  preProcessNCCSFile!(data)

  return data
end



function changeNCCSColumnType!(data::NCCSFile, col::Symbol, T::Type)::NCCSFile

  if T <: MString #case 1: String or StringUMissing
    return data #do nothing
  end

  rename!(data.df, [col=>:oldCol])
  oldCol::Vector{R where R <: MString} = (v->ismissing(v) ? missing : string(v)).(data.df[:oldCol])
  NValues::Int = length(oldCol)

  if (T <: CategoricalValue{Union{Symbol, Missing}}) #case 3: pool the categorical data
    data.df[col] = Vector{Union{Symbol, Missing}}(undef, NValues)
    data.df[col] .= ((s::MString->ismissing(s) ? missing : Symbol(s)).(oldCol))
    categorical!(data.df, col)
    #data.df[col] = compress(data.df[col]) NOTE: Consider using this in the future, currently crashes

  elseif T <: Union{MFloat64, MInt, Int, Float64, Bool, MBool} #case 2: concrete numerical type
    #need to parse each entry
    concreteCol::Vector{T} = Vector{T}(undef, NValues)
    if T == Missing #special case where the column is entirely missing
      concreteCol = missings(NValues)
    else
      TSub::Type = typeintersect(T, Union{Float64, Bool, Int}) #get the nonmissable type


      #=@inbounds @simd for i ∈ 1:NValues #parse the column text
        concreteCol[i] = ifelse(ismissing(oldCol[i]),missing, parse("$(oldCol[i])"))
      end=#
      @inbounds @simd for i ∈ 1:NValues
        concreteCol[i] = ismissing(oldCol[i]) ? missing : something(tryparse(TSub, oldCol[i]), missing)
      end

    end

    data.df[col] = concreteCol

  else #Write explicit code for each type of column
    println("Col $col from file $(data.meta.fileName) is of type $T which is not recognized.")
    error("Col $col from file $(data.meta.fileName) is of type $T which is not recognized.")

  end

  #delete the old column
  deletecols!(data.df, :oldCol)

  return data
end

#fixes the type information
function preProcessNCCSFile!(data::NCCSFile)::NCCSFile

  meta::MetaDictionary = data.meta
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))
  actualCols::Vector{Symbol} = (s::Symbol->Symbol(lowercase(string(s)))).(names(data.df))
  names!(data.df, actualCols) #make the names lower case

  #println("$(meta.fileName) presize:", size(data.df))
  #form the lookup dictionary for the type. Takes a column name and gets the type, with a default of String
  typeDict::Dict = Dict(
    actualCols[i] => (
    haskey(meta.dictionary, actualCols[i]) ? meta.dictionary[actualCols[i]] : MString)
    for i::Int ∈ 1:(length(actualCols)))

  #remove columns we don't expect or want
  diff::Vector{Symbol} = setdiff(actualCols, expectedCols)
  if length(diff) > 0
    try
      deletecols!(data.df, setdiff(actualCols, expectedCols))
    catch err
      error("$err\nmeta: $(meta.fileName)")
    end
  end

  #set the column types. Would have preferred to do this on the read, but got a weird error.
  for col::Symbol ∈ actualCols
    changeNCCSColumnType!(data, col, typeDict[col])
  end

  #println("$(meta.fileName) postsize:", size(data.df))
  fields = names(data.df)
  #println("Names $(meta.fileName): $fields")
  for f ∈ fields
    if f ≠ data.meta.metaTuples[data.meta.index[f]].fieldName
      @warn "Dictionary dows not match (preProcessNCCSFile in NCCSFiles)"
      println("Name: Expected=$f Actual= $(data.meta.metaTuples[data.meta.index[f]].fieldName)")
      println("Type: Expected=$(eltype(data.df[f]))
        Actual= $(data.meta.dictionary[f])")

    end
  end

  data = preFilterNCCSFile!(data)

  return data
end

#removes rows from the dataframe in NCCSData
#function receives a dataframe row and outputs a boolean
function Base.filter!(F::Function, data::NCCSFile)

  data.df = data.df[((r::DataFrameRow)->F(r)).(eachrow(data.df)), :]
  #NOTE: Replace with sinple filter
end

#holds a limited amount of pre-filtering logic to improve later stage data processing performance
#NOTE: drops full files and co files before 1997
function preFilterNCCSFile!(data::NCCSFile; revenueLowerBound::Float64 = REVENUE_LOWER_BOUND)::NCCSFile
  nonprofitType::Symbol = data.meta.nonprofitType

  #=for i ∈ 1:size(data.df,2)
    println("eltype $i: ", eltype(data.df[i]))
  end=#

  N::Int = size(data.df, 1)

  #println("size:", size(data.df))
  if data.meta.fullFile
    filter!(r::DataFrameRow->false, data) #ignore these for now
  elseif nonprofitType == :pc
    #data.df = data.df[((ismissing).(data.df[:fundbal]) .== false) .& (data.df[:fundbal] .> assetLowerBound),:]
    filter!(r::DataFrameRow->ismissing(r[:totrev]) == false && r[:totrev] > revenueLowerBound, data)
  elseif nonprofitType == :co #the below filter doesn't really work due to dq issues
    if data.meta.fileYear ≥ 2000
      filter!(r::DataFrameRow->ismissing(r[:totrev2]) == false && r[:totrev2] > revenueLowerBound, data)
    else
      filter!(r::DataFrameRow->ismissing(r[:grossrec]) == false && r[:grossrec] > revenueLowerBound, data)
    end
  elseif nonprofitType == :pf
    filter!(r::DataFrameRow->ismissing(r[:p1totrev]) == false && r[:p1totrev] > revenueLowerBound, data)
  end
  #=if data.meta.fullFile
    filter!(r::DataFrameRow->false, data.df) #ignore these for now
  elseif nonprofitType == :pc
    #data.df = data.df[((ismissing).(data.df[:fundbal]) .== false) .& (data.df[:fundbal] .> assetLowerBound),:]
    filter!(r::DataFrameRow->ismissing(r[:totrev]) == false && r[:totrev] > revenueLowerBound, data.df)
  elseif nonprofitType == :co #the below filter doesn't really work due to dq issues
    if data.meta.fileYear ≥ 2000
      filter!(r::DataFrameRow->ismissing(r[:totrev2]) == false && r[:totrev2] > revenueLowerBound, data.df)
    else
      filter!(r::DataFrameRow->ismissing(r[:grossrec]) == false && r[:grossrec] > revenueLowerBound, data.df)
    end
  elseif nonprofitType == :pf
    filter!(r::DataFrameRow->ismissing(r[:p1totrev]) == false && r[:p1totrev] > revenueLowerBound, data.df)
  end=#

  return data
end


#serializes a NCCSFile object
function serializeNCCSFile(data::NCCSFile; nccsPath::String = NCCS_PATH,
    oStreamer::Function = O_JLS_STREAMER,verbose::Bool=true)::NNCCSFile

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

  return returnNCCSFile ? data : nothing
end

#deserializes a NCCSFile object
function deserializeNCCSFile(fileName::String; nccsPath::String = NCCS_PATH,
    iStreamer::Function = I_JLS_STREAMER)::NCCSFile

    inPath::String = "$nccsPath\\$(fileName).jls.gz"

    #=iStream = iStreamer(inPath)
    data::NCCSFile = deserialize(iStream)
    close(iStream)=#

    b = iStreamer(inPath) do iStream
      b = IOBuffer(read(iStream))
      (b)
    end

    data::NCCSFile = deserialize(b)

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
  serializeOnly::Bool=true, maxSerializationWorkers::Int = MAX_SERIALIZATION_WORKERS)::Union{Nothing, Vector{DataFrame}}

  numFiles::Int = length(metadata)
  pids = workers()
  np::Int = min(length(pids), maxSerializationWorkers)

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

#define this type. We may want to add more later.
mutable struct NCCSData
  df::DataFrame
end

#deals with removing fields
function removeFields!(data::NCCSData, fields::Vector{Symbol})::NCCSData
  deletecols!(data.df, fields)
  #NOTE: If we add processed metadata, deal with that here
  return data
end
