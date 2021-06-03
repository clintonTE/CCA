#contains functions related to building, loading, and saving metadata
#currently metadata is a vector of MetaDictionary, but it could eventually
#become its own object


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
  fileNames::Vector{String} = indexDF[:name]
  NFiles::Int = length(fileNames)

  metadata::Vector{MetaDictionary} = (MetaDictionary).(fileNames)

  return metadata
end

#checks for and removes metadata fields not actually in the file
function validateMetaDictionary!(meta::MetaDictionary, keeperList::Vector{Symbol};
  nccsPath::String = NCCS_PATH,  streamer::Function = I_CSV_STREAMER)::MetaDictionary

  #these are the column names that we expect
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))

  path::String = "$nccsPath\\$(meta.fileName).csv.gz"
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

    removeFields!(meta, fieldsToRemove)
    println("Validated $(meta.fileName). Removed $(length(fieldsToRemove)) fields from metadata.")
  else
    println("Validated $(meta.fileName).")
  end


  return meta

end

function validateMetadata(metadata::Vector{MetaDictionary};
    indexPath::String = DOCUMENTATION_PATH, keeperName::String=KEEPER_NAME)::Vector{MetaDictionary}

  metadata = deepcopy(metadata) #its hard to do this without some modifications AND copies
                                #hence err on the side of caution here

  #metadata = removeDictionary(metadata, "nccs_core_2013_co_legacy")
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
  metaNames::Vector{String} = (m::MetaDictionary->m.fileName).(metadata)

  #assign and retrieve the work
  for i::Int ∈ 1:NMeta
    push!(futures, (@spawn validateMetaDictionary!(metadata[i], keeperList)))
  end
  metadata .= (fetch).(futures)


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
