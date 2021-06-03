#this file holds functions specific to acquiring the NCCS data



function downloadNCCSCSV(url::String, fileName::String, outPath::String,
  outName::String = fileName; validate::Bool=true,
  refreshData::Bool=true, maxAttempts::Int = 12)::Nothing




  if !refreshData && isfile("$outPath\\$(outName).csv.gz")# if the file exists
    if validate #validate if desired
    try #validate the file
        validateSingleTypeCSV("$outPath\\$outName.csv.gz")
        println("Download Validated: $outName")
        return nothing
      catch err
        #do nothing and do not return anything
      end
    else
      return nothing
    end
  end

  attempt::Int = 0
  retry::Bool = true
  while attempt < maxAttempts && retry #loops and tries to download all the files
    attempt+=1
    retry = false
    #download the file
    oStream::GZipStream = GZip.open("$outPath\\$outName.csv.gz", "w")
    #b::IOBuffer = IOBuffer()
    #weird error if using the transcode codec

    #acquire the data
    HTTP.open(io::IO->
      while !eof(io)
        write(oStream, readavailable(io))
      end,
      "GET", "$url/$fileName", retries=4)#, connection_limit=16,
      #pipeline_limit=32, readtimeout=120)


    #bString::String = String(take!(b))
    #TODO: Might need to uncomment #bString = replace(bString, "\0"=>"")
    #write(oStream, bString)
    GZip.close(oStream)

    #check if its a valid file
    try
      validateSingleTypeCSV("$outPath\\$outName.csv.gz")
      println("Download Validated: $outName")
    catch err
      @warn "DL Error with $outName. $(
        attempt<maxAttempts ? "Retrying ($attempt/$maxAttempts)" : "Max attempts Reached for file $outName.csv.")
        \nErrorMessage:\n$err\nurl:$url/$fileName"
      retry = true
    end

  end

  return nothing
end

#vectorized version (useful for parallelization)
function downloadNCCSCSVs(urls::Vector{String}, fileNames::Vector{String},
  outPath::String, outNames::Vector{String}, mode::Symbol = :serial,
  ;delay::Float64 = 0.0,
  validate::Bool = true, refreshData::Bool=true)::Nothing



  #short-circuit for parallel operation
  if mode == :parallel
    return downloadNCCSCSVsPar(urls, fileNames, outPath, outNames,
    delay=delay, validate=validate, refreshData=refreshData)
  end

  for i::Int ∈ 1:length(urls)

    if delay ≥ .001 #allows for a delay to avoid slamming the site
      sleep(delay)
    end

    downloadNCCSCSV(urls[i], fileNames[i], outPath, outNames[i],
      validate=validate, refreshData=refreshData)
  end

  return nothing
end

#useful for parallelization
function downloadNCCSCSVsPar(urls::Vector{String}, fileNames::Vector{String},
  outPath::String, outNames::Vector{String};
  delay::Float64 = 0.0, validate::Bool = true, refreshData::Bool=true)::Nothing

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
        (@spawnat pids[i] downloadNCCSCSVs(urls[assignments[i]],
          fileNames[assignments[i]], outPath,
          outNames[assignments[i]], delay=delay, validate=validate,
          refreshData=refreshData)))
    end
  end

  #get the work
  (fetch).(futures)

  return nothing
end

function readNCCSCSVsIndex(;inPath::String = DATA_PATH, filesIndexName::String = FILES_INDEX_NAME)::DataFrame
  return CSV.read("$inPath\\$filesIndexName.csv", categorical=false)
  #return load("$inPath\\$filesIndexName.csv") |> DataFrame
end

function downloadNCCSCSVs(outPath::String; mode::Symbol=:parallel,
  delay::Float64=0.0, validate::Bool = true, refreshData::Bool=true)

  #read in the metadata
  fileDF::DataFrame = readNCCSCSVsIndex()
  numFiles::Int = size(fileDF,1)

  #build the urls
  urls::Vector{String} = ((urlStem::String, urlYear::Int)->
    "$urlStem/$urlYear").(fileDF[:urlStems], fileDF[:urlYears])
  fileNames::Vector{String} =  fileDF[:fileNames]
  outNames::Vector{String} =  fileDF[:outNames]

  return downloadNCCSCSVs(urls, fileNames, outPath, outNames, mode, delay=delay,
    validate=validate, refreshData=refreshData)

end
