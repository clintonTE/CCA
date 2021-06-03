#this file holds functions related to data processing




#reads in a CSV type
function readSingleTypeCSV(path::String; colType=Union{Missing, String},
    streamer::Function = (p::String)->GzipDecompressorStream(open(p)))::DataFrame

  #count the number of cols
  iStream = streamer(path)
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, rows=2)
  close(iStream)
  nCols::Int = size(df,2)

  #set the type
  types::Dict = Dict(i=>Union{String,Missing} for i ∈ 1:nCols)
  iStream = streamer(path)
  df = CSV.read(iStream, strings=:raw, categorical=false, types=types)
  close(iStream)

  return df
end

#same as the read version only returns nothing and just checks for errors
function validateSingleTypeCSV(path::String; colType=Union{Missing, String},
  streamer::Function = (p::String)->GzipDecompressorStream(open(p)))::Nothing

  #count the number of cols

  iStream = streamer(path)
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, rows=2)
  close(iStream)
  nCols::Int = size(df,2)

  #set the type

  types::Dict = Dict(i=>Union{String,Missing} for i ∈ 1:nCols)
  iStream = streamer(path)
  CSV.validate(iStream, strings=:raw, categorical=false, types=types)
  close(iStream)

end



#reads in the metadata and creates type tuples
function readMetaIndex(;indexPath::String = DOCUMENTATION_PATH,
  fileName::String = META_INDEX_NAME)::DataFrame

  iStream::IOStream = open("$indexPath\\$fileName.csv")
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false)
  close(iStream)

  return df

end

#froms the vector of Metadata Dictionaries
function formMetadata(indexPath::String = DOCUMENTATION_PATH,
  indexName::String = META_INDEX_NAME)::Vector{MetaDictionary}

  #first get a list of the paths
  indexDF::DataFrame = readMetaIndex()
  names::Vector{String} = indexDF[:name]
  dataNames::Vector{String} = indexDF[:dataName]
  NFiles::Int = length(names)

  metadata::Vector{MetaDictionary} = (MetaDictionary).(dataNames, names)

  return metadata
end

#checks for and removes metadata fields not actually in the file
function validateMetaDictionary!(meta::MetaDictionary, keeperList::Vector{Symbol}; nccsPath::String = NCCS_PATH,
  streamer::Function = I_CSV_STREAMER)::MetaDictionary

  #these are the column names that we expect
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))

  path::String = "$nccsPath\\$(meta.dataFileName).csv.gz"
  iStream = streamer(path)
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, rows=2)
  close(iStream)

  #check first that the column names are correct
  actualCols::Vector{Symbol} = (s::Symbol->Symbol(lowercase(string(s)))).(names(df))

  #get the list of files we don't want
  expectedNotFound::Vector{Symbol} = setdiff(expectedCols, actualCols)
  expectedNotWanted::Vector{Symbol} = setdiff(expectedCols, keeperList)
  fieldsToRemove::Vector{Symbol} = unique([expectedNotFound; expectedNotWanted])

  #remove the offending fields
  if length(fieldsToRemove) > 0
    #=if contains(meta.dataFileName,"2009") || contains(meta.dataFileName,"2008")
      println("$(meta.dataFileName) pre: $(intersect(keys(meta.dictionary), fieldsToRemove))")
    end=#
    removeFields!(meta, fieldsToRemove)
    println("Validated $(meta.dataFileName). Removed $(length(fieldsToRemove)) fields from metadata.")
  else
    println("Validated $(meta.dataFileName).")
  end

  #=if contains(meta.dataFileName,"2009") || contains(meta.dataFileName,"2008")
    println("$(meta.dataFileName) post: $(intersect(keys(meta.dictionary), fieldsToRemove))")
  end=#

  return meta

end

function validateMetadata(metadata::Vector{MetaDictionary}; indexPath::String = DOCUMENTATION_PATH,
    keeperName::String=KEEPER_NAME)::Vector{MetaDictionary}

  metadata = deepcopy(metadata) #its hard to do this without some modifications AND copies
                                #hence err on the side of caution here

  metadata = removeDictionary(metadata, "LEGACY")
  NMeta::Int = length(metadata)

  #get list of kept columns
  iStream::IOStream = open("$indexPath\\$keeperName.csv")
  keeperDF::DataFrame = CSV.read(iStream, strings=:raw, categorical=false)
  close(iStream)
  keeperList::Vector{Symbol} = (Symbol).(keeperDF[:keeper])

  #this can take a while so multi-thread it
  println("Beginning file validation")

  futures::Vector{Future} = Vector{Future}()

  #use this to check orders
  metaNames::Vector{String} = (m::MetaDictionary->m.metaFileName).(metadata)

  #=metaDict::Dict = meta2Index(metadata)
  i2008::Int = metaDict["Core 2008 PC"]

  println("Core 2008 PC pre [$(length(metadata[i2008].metaTuples))]:
    $(keys(metadata[i2008].dictionary))")=#
  #assign and retrieve the work
  for i::Int ∈ 1:NMeta
    push!(futures, (@spawn validateMetaDictionary!(metadata[i], keeperList)))
  end
  metadata .= (fetch).(futures)

  #=
  println("Core 2008 PC post1 [$(length(metadata[i2008].metaTuples))]:
    $(keys(metadata[i2008].dictionary))")=#

  return metadata
end

#this function either forms or loads the metadata
function getMetadata(;refreshMetadata::Bool = true, metadataPath::String = DOCUMENTATION_PATH,
  metadataName::String = METADATA_NAME)

  metadataFilePath::String = "$metadataPath\\$metadataName.jls"

  #load the metadata if we don't need to redo the validation
  if refreshMetadata || !(isfile(metadataFilePath))
    metadata::Vector{MetaDictionary}= formMetadata()
    metadata = validateMetadata(metadata)

    #=  metaDict::Dict = meta2Index(metadata)
    i2008::Int = metaDict["Core 2008 PC"]

    println("Core 2008 PC post2 [$(length(metadata[i2008].metaTuples))]:
      $(keys(metadata[i2008].dictionary))")
    warn("breakpoint")=#

    oStream::IOStream = open(metadataFilePath, "w")
    serialize(oStream, metadata)
    close(oStream)
  else
    iStream::IOStream = open(metadataFilePath)
    metadata = deserialize(iStream)
    close(iStream)
  end

  return metadata
end



function preProcessData(; update::Bool=false, validateDownloads::Bool = true,
  nccsPath::String = NCCS_PATH,  refreshMetadata::Bool = true, refreshData::Bool = true)::Nothing

  #println("got here 1")
  getNCCSFiles(NCCS_PATH, update=update, #=delay=0.001,=# validate=validateDownloads)
  metadata::Vector{MetaDictionary} = getMetadata(refreshMetadata = refreshMetadata)
  loadNCCSFiles(metadata, refresh=refreshData, serializeOnly = true, mode=:parallel)


  return nothing
end
