#contains methods and associated helper functions for forming a MetaDictionary object
#each object contains the complete metadata of a file

mutable struct MetaDictionary
  fileName::String
  metaTuples::Vector{MetaTuple}

  nonprofitType::Symbol
  fileYear::Int
  fullFile::Bool
  fileDescription::String

  dictionary::Dict
  index::Dict
end


#return a metadictionary object. The actual dictionary can be passed as a type lookup file
function MetaDictionary(metaTuples::Vector{MetaTuple}, indexRow::DataFrameRow)::MetaDictionary
  NTuples = length(metaTuples)

  fileName::String = indexRow[:name]
  nonprofitType::Symbol = Symbol(FILE_TYPES_NAME_TO_CODE[indexRow[:type]])
  fileYear::Int = indexRow[:year]
  fullFile::Bool = uppercase(indexRow[:full]) == "TRUE"
  fileDescription::String = indexRow[:description]

  dictionary::Dict = Dict(metaTuples[i].fieldName => metaTuples[i].columnType for i::Int ∈ 1:NTuples)
  index::Dict = Dict(metaTuples[i].fieldName => i for i::Int ∈ 1:NTuples)

  #extract the information from the row and create the object
  return MetaDictionary(fileName, metaTuples, nonprofitType,
    fileYear, fullFile, fileDescription, dictionary, index)
end

#createas a MetaDictionary from the metadata file
function MetaDictionary(indexRow::DataFrameRow;
   metaPath::String = META_PATH,
  subTablePath::String = SUB_TABLES_PATH)::MetaDictionary

  fileName::String = indexRow[:name]
  #first open the file
  metaDF::DataFrame = CSV.read("$metaPath\\$(fileName).csv", strings=:raw, categorical=false)
  #metaDF::DataFrame = load("$metaPath\\$(fileName).csv") |> DataFrame

  NFields::Int = size(metaDF,1)

  #this is where we will put the results
  metaTuples::Vector{MetaTuple} = Vector{MetaTuple}(undef, NFields)

  #these are the data
  names::Vector{Symbol} = (Symbol).(metaDF[:name])
  types::Vector{Symbol} = (Symbol).(metaDF[:type])
  subTables::Vector{Bool} = (s::String->s=="True" ? true : false).(metaDF[:subTable])
  shortDescriptions::Vector{MString} = metaDF[:shortDesc]
  longDescriptions::Vector{MString} = metaDF[:longDesc]


  for i::Int ∈ 1:NFields
    metaTuples[i] = MetaTuple(fileName, names[i], types[i],
      subTables[i], shortDescriptions[i], longDescriptions[i])
  end

  return MetaDictionary(metaTuples, indexRow)
end

#rebuilds the metadata table from the dictionary
function metaTable(meta::MetaDictionary)::DataFrame
  N::Int = length(meta.metaTuples)

  #pre-allocate
  fieldNames::Vector{Symbol} = Vector{Symbol}(undef, N)
  types::Vector{Type} = Vector{Type}(undef, N)
  typeNames::Vector{String} = Vector{String}(undef, N)
  shortDescriptions::Vector{MString} = Vector{MString}(undef, N)
  longDescriptions::Vector{MString} = Vector{MString}(undef, N)
  numeric::Vector{MBool} = Vector{MBool}(undef, N)

  @inbounds @simd for i ∈ 1:N
    fieldNames[i] = meta.metaTuples[i].fieldName
    types[i] = meta.metaTuples[i].columnType
    shortDescriptions[i] = meta.metaTuples[i].shortDescription
    longDescriptions[i] = meta.metaTuples[i].longDescription

    #need to figure out the types
    typeNames[i] =
      ifelse(types[i] <: MFloat64, "Float64",
        ifelse(types[i] <: MInt, "Int",
          ifelse(types[i] <: MBool, "Bool",
            ifelse(types[i] <: MString, "String",
              ifelse((types[i] <: CategoricalValue) || (CategoricalValue <: types[i]), "Conformed", "UNK")))))
    numeric[i] = ifelse(types[i] <: MFloat64 || types[i] <: MInt, true, false)
  end

  #now build the dataframe
  df::DataFrame = DataFrame(fieldNames = fieldNames, shortDescriptions=shortDescriptions,
    longDescriptions=longDescriptions, types=typeNames, numeric=numeric)

  #now write the file-level metadata
  df[:fileNames] = Vector{String}(undef, N)
  df[:nonprofitTypes] = Vector{Symbol}(undef, N)
  df[:fileYears] = Vector{Int}(undef, N)
  df[:fullFiles] = Vector{Bool}(undef, N)
  df[:fileDescriptions] = Vector{MString}(undef, N)


  df[:fileNames] .= meta.fileName
  df[:nonprofitTypes] .= meta.nonprofitType
  df[:fileYears] .= meta.fileYear
  df[:fullFiles] .= meta.fullFile
  df[:fileDescriptions] .= meta.fileDescription

  return df
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
