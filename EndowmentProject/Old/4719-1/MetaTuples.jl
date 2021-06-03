# holds all functions related to the type MetaTuple



#all information about each field
struct MetaTuple{T<:Union{String, Bool, Conformed, Missing, Int, Float64}}
  fieldName::Symbol
  columnType::Type
  columnDict::S where S<:Union{Dict, Nothing}
  shortDescription::U where U<:MString
  longDescription::V where V<:MString
end

#for atomic type fields
function MetaTuple{T}(fieldName::Symbol, shortDescription::U where U<:MString,
  longDescription::V where V<:MString, ::Type{T})::MetaTuple{T}

  if !(T <: Union{Bool, String, Missing, Float64, Int})
    error("Invalid atomic metatuple type $T")
  end
  return MetaTuple{T}(fieldName, T, nothing, shortDescription, longDescription)
end

#for conformed fields
function MetaTuple{T}(fieldName::Symbol, valueDict::Dict, shortDescription::MString,
  longDescription::MString, ::Type{T})::MetaTuple{T}

  if !(T<:Union{Conformed, Missing})
    error("Invalid conformed metatuple type $T")
  end

  return MetaTuple{T}(fieldName::Symbol, CategoricalValue{T}, valueDict, shortDescription, longDescription)
end

#convenience method which creates the dictionary of valid values from the dataframe
function MetaTuple{T}(fieldName::Symbol, columnDF::DataFrame,
    shortDescription::MString, longDescription::MString, ::Type{T})::MetaTuple{T}

  NValues::Int = size(columnDF, 1)

  #convert the value column to symbols
  values::Vector{Union{Symbol, Missing}}=
    ((s::U where U<:Union{Int, Float64, String, Missing})->ismissing(s)?missing:Symbol(s)).(columnDF[:value])

  #convert the description columns to nothing

  columnDict = Dict(values[i] => columnDF[i,:description] for i::Int ∈ 1:NValues)

  return MetaTuple(fieldName, columnDict, shortDescription, longDescription, T)

end

#creates a MetaTuple purely using the fields from the metadata index
function MetaTuple(fileName::String, fieldName::Symbol, fieldType::Symbol, subTableInd::Bool,
    shortDescription::MString, longDescription::MString; subTablesPath::String = SUB_TABLES_PATH)::MetaTuple

    if !subTableInd
      fieldType::Type = Union{VAR_TYPES[fieldType], Missing}
      return MetaTuple(fieldName, shortDescription, longDescription, fieldType)
    else
      fieldType = Union{Conformed, Missing}

      #get the subtable info
      iStream::IOStream = open("$subTablesPath\\$fileName\\$fileName - $fieldName.csv")
      columnDF::DataFrame = CSV.read(iStream, strings=:raw, categorical=false)
      close(iStream)

      return MetaTuple(fieldName, columnDF, shortDescription, longDescription, fieldType)
    end

end
