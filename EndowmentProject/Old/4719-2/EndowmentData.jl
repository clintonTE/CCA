#this file holds functions related to data processing




#reads in a CSV type
function readSingleTypeCSV(path::String; colType=Union{Missing, String},
    streamer::Function = I_CSV_STREAMER)::DataFrame

  #count the number of cols
  #iStream = streamer(path)
  #df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, silencewarnings=true, rows=2)
  #close(iStream)

  b::IOBuffer = gzopen(path) do f
    s::String = "$(readline(f))\n$(readline(f))\n$(readline(f))"
    println
    b = IOBuffer(s)
    (b)
  end
  df::DataFrame = CSV.read(b, strings=:raw,
    categorical=false, silencewarnings=true, strict=false)


  nCols::Int = size(df,2)
  #println("nCols: $nCols")

  #set the type
  types::Dict = Dict(i=>Union{String,Missing} for i ∈ 1:nCols)
  iStream = streamer(path)
  #df = load(iStream,  colparsers=types) |> DataFrame
  b = IOBuffer(read(iStream))
  df = CSV.read(#=iStream=#b, strings=:raw, categorical=false, types=types,
     silencewarnings=true, strict=false)
  #df = uCSV.read(iStream, coltypes=types, allowmissing=true, escape='\\', quotes='\"')
  close(iStream)

  #println("$(df[1:3, :NAME])")

  return df
end


#same as the read version only returns nothing and just checks for errors
function validateSingleTypeCSV(path::String; colType=Union{Missing, String},
  streamer::Function = I_CSV_STREAMER,
  missingStrings::Vector{String}=MISSING_STRINGS,
  verifyFields::Vector{Symbol} = VERIFY_FIELDS)::Nothing

  #count the number of cols

  #iStream = streamer(path)
  b::IOBuffer = gzopen(path) do f
    b = IOBuffer("$(readline(f))\n$(readline(f))")
    (b)
  end
  df::DataFrame = CSV.read(b, strings=:raw, categorical=false,
    silencewarnings=true, strict=false)

  dfNames::Vector{Symbol} = (s::Symbol->Symbol(lowercase(string(s)))).(names(df))
  for f::Symbol ∈ verifyFields
    if f ∉ dfNames
      @error "Field $f not found in file $path, names found: $dfNames"
    end
  end


  nCols::Int = size(df,2)

  #set the type
  types::Dict = Dict(i=>Union{String,Missing} for i ∈ 1:nCols)
  iStream = streamer(path)
  CSV.read(iStream, strings=:raw, categorical=false, types=types,
    missingstrings=missingStrings, strict=false, silencewarnings=true)
  #df = load(iStream,  colparsers=types) |> DataFrame
  close(iStream)

  return nothing

end



#this fucntion is almsot a script which analyzes data
function preProcessData(; refreshData::Bool=true, validateDownloads::Bool = true,
  nccsPath::String = NCCS_PATH,  refreshMetadata::Bool = true, reserializeData::Bool = true,
  refreshSummary = true)::Vector{MetaDictionary}

  downloadNCCSCSVs(NCCS_PATH, refreshData=refreshData, #=delay=0.001,=#
    validate=validateDownloads, mode=:parallel)
  metadata::Vector{MetaDictionary} = getMetadata(refreshMetadata = refreshMetadata)
  loadNCCSFiles(metadata, reserializeData=reserializeData, serializeOnly = true, mode=:parallel)
  consolidateSummary(metadata, refreshSummary=refreshSummary)

  #=println("fileName: $(metadata[1].fileName)")
  willFail(metadata[1])=#

  return metadata
end
