#this file holds functions specific to acquiring the NCCS data



function getNCCSFile(url::String, fileName::String, outPath::String;
  outName::String = fileName, validate::Bool=true, update::Bool=false)::Nothing

  if update && isfile("$outPath\\$outName.gz")# if the file exists
    if validate #validate if desired
      try #validate the file
        validateSingleTypeCSV("$outPath\\$outName.gz")
        println("Download Validated: $outName")
        return nothing
      end
    else
      return nothing
    end
  end

  #download the file
  #oStream::GZipStream = GZip.open("$outPath\\$outName.gz", "w")
  oStream::GZipStream = GZip.open("$outPath\\$outName.gz", "w")
  #weird error if using the transcode codec

  #acquire the data
  HTTP.open(io::IO->
    while !eof(io)
      write(oStream, readavailable(io))
    end,
    "GET", "$url/$fileName", retries=4, connection_limit=16,
    pipeline_limit=32, readtimeout=60)
  GZip.close(oStream)

  #check if its a valid file
  if validate
    try
      validateSingleTypeCSV("$outPath\\$outName.gz")
      println("Download Validated: $outName")
    catch err
      warn("DL Error with $outName. \nErrorMessage:\n$err")
    end
  end



  return nothing
end

#vectorized version (useful for parallelization)
function getNCCSFiles(urls::Vector{String}, fileNames::Vector{String},
  outPath::String, mode::Symbol = :serial;
  outNames::Vector{String} = fileNames, delay::Float64 = 0.0,
  validate::Bool = true, update::Bool=false)::Nothing



  #short-circuit for parallel operation
  if mode == :parallel
    return getNCCSFilesPar(urls, fileNames, outPath, outNames=outNames,
    delay=delay, validate=validate, update=update)
  end

  for i::Int ∈ 1:length(urls)

    if delay ≥ .001 #allows for a delay to avoid slamming the site
      sleep(delay)
    end

    getNCCSFile(urls[i], fileNames[i], outPath, validate=validate, update=update)
  end

  return nothing
end

function getNCCSFilesPar(urls::Vector{String}, fileNames::Vector{String},
  outPath::String; outNames::Vector{String} = fileNames,
  delay::Float64 = 0.0, validate::Bool = true, update::Bool=false)::Nothing

  numFiles::Int = length(urls)
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
        (@spawnat pids[i] getNCCSFiles(urls[assignments[i]],
          fileNames[assignments[i]], outPath,
          outNames=outNames, delay=delay, validate=validate, update=update)))
    end
  end

  #get the work
  (fetch).(futures)

  return nothing
end

function getNCCSFilesIndex(;inPath::String = DATA_PATH, filesIndexName::String = FILES_INDEX_NAME)::DataFrame
  return CSV.read("$inPath\\$filesIndexName.csv", strings=:raw,
    categorical=false)
end

function getNCCSFiles(outPath::String; mode::Symbol=:parallel,
  delay::Float64=0.0, validate::Bool = true, update::Bool=false)

  #read in the metadata
  fileDF::DataFrame = getNCCSFilesIndex()
  numFiles::Int = size(fileDF,1)

  #build the urls
  urls::Vector{String} = ((urlStem::String, urlYear::Int)->
    "$urlStem/$urlYear").(fileDF[:urlStems], fileDF[:urlYears])
  fileNames::Vector{String} =  fileDF[:fileNames]

  return getNCCSFiles(urls, fileNames, outPath, mode, delay=delay,
    validate=validate, update=update)

end

#serializes and compresses an already downloaded NCCS file
function serializeNCCSFile(meta::MetaDictionary; nccsPath::String = NCCS_PATH,
  iStreamer::Function = I_CSV_STREAMER, oStreamer::Function = O_JLS_STREAMER,
  serializeOnly::Bool = true, verbose::Bool = true)::DataFrame

  #these are the column names that we expect
  expectedCols::Vector{Symbol} = collect(keys(meta.dictionary))

  inPath::String = "$nccsPath\\$(meta.dataFileName).csv.gz"
  outPath::String = "$nccsPath\\$(meta.dataFileName).jls.gz"

  if verbose
    println("Serializing $(meta.dataFileName)...")
  end

  iStream = iStreamer(inPath)
  df::DataFrame = CSV.read(iStream, strings=:raw, categorical=false, rows=2)
  close(iStream)

  #get the column names and make sure each one is in the list of columns we want
  actualCols::Vector{String} = (s::Symbol->lowercase(string(s))).(names(df))

  #form the lookup dictionary for the type. Takes a column name and gets the type, with a default of String
  typeDict::Dict = Dict(
    actualCols[i] => (
    haskey(meta.dictionary, actualCols[i]) ? meta.dictionary[Symbol(actualCols[i])]: Union{String, Missing})
    for i::Int ∈ 1:(length(actualCols)))

  #now read in the file
  iStream = iStreamer(inPath)
  #CSV.validate(iStream, header=actualCols, datarow=2, types=typeDict)
  df = CSV.read(iStream, header=actualCols, datarow=2, types=typeDict)
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
    println("Serialization of $(meta.dataFileName) complete.")
  end

  if serializeOnly
    df = DataFrame([[nothing]]) #free the memory with an empty dataframe
  end

  return df

end

#loads an already downloaded NCCS file
function loadNCCSFile(meta::MetaDictionary; nccsPath::String = NCCS_PATH,
  refresh::Bool = true, iStreamer::Function = I_JLS_STREAMER, serializeOnly::Bool = true)

  filePath::String = "$nccsPath\\$(meta.dataFileName).jls.gz"

  #only re-process the file if necessesary
  if refresh || !isfile(filePath)
    df::DataFrame = serializeNCCSFile(meta, serializeOnly=serializeOnly)
  elseif !serializeOnly
    iStream = iStreamer(filePath)
    df = deserialize(iStream)
    close(iStream)
  else #do nothing
    df = DataFrame([[nothing]]) #an empty dataframe
  end

  return df
end

#loads and/or serializes a vector of corresponding files
function loadNCCSFiles(metadata::Vector{MetaDictionary}; mode::Symbol = :serial,
  refresh::Bool = true, returnFiles::Bool = false, serializeOnly::Bool = true)::Vector{DataFrame}

  if mode == :parallel #relaunch the dedicated parallel function
    return loadNCCSFilesPar(metadata::Vector{MetaDictionary}, refresh=refresh, serializeOnly=serializeOnly)
  end

  if serializeOnly
    (meta::MetaDictionary->loadNCCSFile(meta,refresh=refresh, serializeOnly=serializeOnly)).(metadata)
    return Vector{DataFrame}()
  else
    return (meta::MetaDictionary->loadNCCSFile(meta,refresh=refresh)).(metadata)
  end
end

#loads many NCCS serial files in parallel
function loadNCCSFilesPar(metadata::Vector{MetaDictionary}; refresh::Bool = true,
  serializeOnly::Bool=true)::Vector{DataFrame}

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
        (@spawnat pids[i] loadNCCSFiles(metadata[assignments[i]], refresh=refresh, serializeOnly=serializeOnly)))
    end
  end

  #get the work
  return vcat((fetch).(futures)...)
end
