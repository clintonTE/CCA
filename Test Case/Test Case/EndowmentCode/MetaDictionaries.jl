#contains methods and associated helper functions for forming a MetaDictionary object
#each object contains the complete metadata of a file

mutable struct MetaDictionary
  fileName::String
  metaTuples::Vector{MetaTuple}
  dictionary::Dict
  index::Dict
end

#return a metadictionary object. The actual dictionary can be passed as a type lookup file
function MetaDictionary(fileName::String, metaTuples::Vector{MetaTuple})::MetaDictionary
  NTuples = length(metaTuples)
  dictionary::Dict = Dict(metaTuples[i].fieldName => metaTuples[i].columnType for i::Int ∈ 1:NTuples)
  index::Dict = Dict(metaTuples[i].fieldName => i for i::Int ∈ 1:NTuples)

  return MetaDictionary(fileName, metaTuples, dictionary, index)
end

#createas a MetaDictionary from the metadata file
function MetaDictionary(fileName::String; metaPath::String = META_PATH,
  subTablePath::String = SUB_TABLES_PATH)::MetaDictionary

  #first open the file
  iStream = open("$metaPath\\$fileName.csv")
  metaDF::DataFrame = CSV.read(iStream, strings=:raw, categorical=false)
  close(iStream)

  NFields::Int = size(metaDF,1)

  #this is where we will put the results
  metadata::Vector{MetaTuple} = Vector{MetaTuple}(NFields)

  #these are the data
  names::Vector{Symbol} = (Symbol).(metaDF[:name])
  types::Vector{Symbol} = (Symbol).(metaDF[:type])
  subTables::Vector{Bool} = (s::String->s=="True"?true:false).(metaDF[:subTable])
  shortDescriptions::Vector{MString} = metaDF[:shortDesc]
  longDescriptions::Vector{MString} = metaDF[:longDesc]


  for i::Int ∈ 1:NFields
    metadata[i] = MetaTuple(fileName, names[i], types[i],
      subTables[i], shortDescriptions[i], longDescriptions[i])
  end

  return MetaDictionary(fileName, metadata)
end

#removes fields from a meta dictionary
function removeField!(meta::MetaDictionary, field::Symbol)::MetaDictionary
  target::Int = meta.index[field]
  meta.metaTuples = meta.metaTuples[(i::Int->i≠target).(1:end)]

  NTuples = length(meta.metaTuples)
  delete!(meta.dictionary, field)
  meta.index = Dict(meta.metaTuples[i].fieldName => i for i::Int ∈ 1:NTuples)

  return meta
end

#removes fields from a meta dictionary
function removeFields!(meta::MetaDictionary, fields::Vector{Symbol})::MetaDictionary
  targets::Vector{Int} = ((s::Symbol)->meta.index[s]).(fields) #get the indices
  meta.metaTuples = meta.metaTuples[(i::Int->i ∉ targets).(1:end)] #delete the data

  NTuples = length(meta.metaTuples) #update the dictionaries
  #(delete!).(meta.dictionary, fields)
  meta.dictionary = Dict(
    meta.metaTuples[i].fieldName => meta.metaTuples[i].columnType for i::Int ∈ 1:NTuples)
  meta.index = Dict(
    meta.metaTuples[i].fieldName => i for i::Int ∈ 1:NTuples)

  return meta
end

function removeDictionary(metadata::Vector{MetaDictionary}, fileName::String)
  targets::Vector{Int} = Vector{Int}()
  NMeta::Int = length(metadata)

  for i::Int ∈ 1:NMeta
    if metadata[i].fileName == fileName
      push!(targets, i)
    end
  end

  if length(targets) < 1
    error("Attempted to delete $fileName but $fileName not found.")
  end

  metadata = metadata[((i::Int)->i ∉ targets).(1:NMeta)]

  return metadata
end



#returns a dictionary for converting data names to meta names
meta2Index(metadata::Vector{MetaDictionary})::Dict =
  Dict(metadata[i].fileName => i for i::Int ∈ 1:length(metadata))
