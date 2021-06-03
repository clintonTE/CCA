#this simple script gets IRS pdf data



function downloadIRSPDF(url::String, fileName::String, outPath::String,
  outName::String = fileName;
  refreshData::Bool=true, maxAttempts::Int = 12)::Nothing


  if !refreshData && isfile("$outPath\\$(outName).pdf")# if the file exists
    return nothing
  end

  attempt::Int = 0
  retry::Bool = true
  while attempt < maxAttempts && retry #loops and tries to download all the files
    attempt+=1
    retry = false
    #download the file
    oStream::IOStream = open("$outPath\\$outName.pdf", "w")
    #weird error if using the transcode codec

    #acquire the data
    try
      HTTP.open(io::IO->
        while !eof(io)
          write(oStream, readavailable(io))
        end,
        "GET", "$url/$fileName", retries=4, connection_limit=16,
        pipeline_limit=32, readtimeout=120)
      close(oStream)
    catch err
      @warn "DL Error with $outName. $(
        attempt<maxAttempts ? "Retrying ($attempt/$maxAttempts)" : "Max attempts Reached for file $outName.pdf.")
        \nErrorMessage:\n$err"
      retry = true
    end

    println("Downloaded IRS file $outName.")
  end

  return nothing
end

#vectorized version (useful for parallelization)
function downloadIRSPDFs(urls::Vector{String}, fileNames::Vector{String},
  outPath::String, outNames::Vector{String}, mode::Symbol = :serial,
  ;delay::Float64 = 0.0, refreshData::Bool=true)::Nothing



  #short-circuit for parallel operation
  if mode == :parallel
    return downloadIRSPDFs(urls, fileNames, outPath, outNames,
    delay=delay, refreshData=refreshData)
  end

  for i::Int ∈ 1:length(urls)

    if delay ≥ .001 #allows for a delay to avoid slamming the site
      sleep(delay)
    end

    downloadIRSPDF(urls[i], fileNames[i], outPath, outNames[i],
      refreshData=refreshData)
  end

  return nothing
end

#useful for parallelization
function downloadIRSPDFsPar(urls::Vector{String}, fileNames::Vector{String},
  outPath::String, outNames::Vector{String};
  delay::Float64 = 0.0, refreshData::Bool=true)::Nothing

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
        (@spawnat pids[i] downloadIRSPDFs(urls[assignments[i]],
          fileNames[assignments[i]], outPath,
          outNames[assignments[i]], delay=delay,
          refreshData=refreshData)))
    end
  end

  #get the work
  (fetch).(futures)

  return nothing
end

function readIRSPDFsIndex(;inPath::String = DOCUMENTATION_PATH,
  filesIndexName::String = IRS_INDEX_NAME)::DataFrame
  return CSV.read("$inPath\\$filesIndexName.csv",# strings=:raw,
    categorical=false)
end

function downloadIRSPDFs(;outPath::String = IRS_PATH, mode::Symbol=:parallel,
  delay::Float64=0.0, refreshData::Bool=true)

  #read in the metadata
  fileDF::DataFrame = readIRSPDFsIndex()
  numFiles::Int = size(fileDF,1)

  #build the file names and urls
  fileNames::Vector{String} =  fileDF[:filePrefix]
  fileNames .*= ((string).(fileDF[:urlYears]))
  fileNames .*= ".pdf"
  urls::Vector{String} = fileDF[:urlStems]

  #now build the outnames
  outNames::Vector{String} = fileDF[:outPrefix]
  outNames .*= "_" .* ((string).(fileDF[:urlYears]))
  outNames .*= "_" .* fileDF[:outSuffix]

  return downloadIRSPDFs(urls, fileNames, outPath, outNames, mode, delay=delay,
    refreshData=refreshData)

end
