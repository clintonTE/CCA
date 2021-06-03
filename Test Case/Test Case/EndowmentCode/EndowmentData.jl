#this file holds functions related to data processing




#reads in a CSV type
function readSingleTypeCSV(path::String; colType=Union{Missing, String},
    streamer::Function = I_CSV_STREAMER)::DataFrame

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
  streamer::Function = I_CSV_STREAMER)::Nothing

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


#this fucntion is almsot a script which analyzes data
function preProcessData(; refreshData::Bool=true, validateDownloads::Bool = true,
  nccsPath::String = NCCS_PATH,  refreshMetadata::Bool = true, reserializeData::Bool = true)::Nothing

  #downloadNCCSCSVs(NCCS_PATH, refreshData=refreshData, #=delay=0.001,=# validate=validateDownloads)
  metadata::Vector{MetaDictionary} = getMetadata(refreshMetadata = refreshMetadata)
  #=loadNCCSFiles(metadata, reserializeData=reserializeData, serializeOnly = true, mode=:parallel)
  summarizeNCCSFile(testFile)=#
  #writeSummaryNCCSTables(metadata[25:27])

  #indexDF::DataFrame = readMetaIndex()
  #metaDict::MetaDictionary = MetaDictionary(indexDF[1,:name])
  #println("fileName: $(metaDict.fileName)")
  future = @spawnat workers()[1] willFail(metadata[1])

  fetch(future)

  #=println("fileName: $(metadata[1].fileName)")
  willFail(metadata[1])=#

  return nothing
end
