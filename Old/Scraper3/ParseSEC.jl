
using  DataFrames, Distributions, StatsBase, ZipFile

const WORKERS = 9 #seems like logical cores * (3/4) is the magic value
const CLEAN_WORKERS = true

if nworkers() < WORKERS || (nworkers()==1 && WORKERS == 1)
  addprocs(nworkers()==1?WORKERS-nworkers()+1:WORKERS-nworkers())
end

@everywhere using DataFrames, ZipFile

@everywhere if !isdefined(:constDef)
  const constDef = true
  #other constants here

  const debugHFMod = true #switch for performance impairing checks
  const dataFolder = "data"
  const dataPath = pwd() * "\\$dataFolder"
  const numColumns = 9
  const sep = '+' #readtable is finicy and so a character no one uses is important
  const extractURL = "https://www.sec.gov/files/dera/data/financial-statement-data-sets/"
  const forceLoad = false
end




#parallelizes the downloading and processing
@everywhere function _parTasks(y::Int, q::Int)
  address::String = extractURL * "$y" * "q" * "$q.zip"
  fileName::String = "$y" * "q$q"
  numFileName::String = "$y" * "q$q" * "num"

  if !isfile("$dataPath\\$numFileName.jls") || forceLoad == true
    println("JLS binary not found. Attempting to load data directly.")
    #attempt to download the data files for the appropriate range
    try
      download(address,"$dataPath\\$fileName.zip")
    catch err
      println("Error: Failed to download ", address, "\n Message: ",err)
    end

    #unzip the appropriate content
    zHandle = ZipFile.Reader("$dataPath\\$fileName.zip")
    oStream::IOStream = open("$dataPath\\$numFileName.txt", "w+")
    for f in zHandle.files
      if f.name == "num.txt"
        b::Vector{UInt8} = Vector{UInt8}(f.uncompressedsize)
        readbytes!(f,b)
        write(oStream,b)
      end
    end

    close(oStream)

    #write the binary
    oStream = open("$dataPath\\$numFileName.dat","w+")
    iStream::IOStream = open("$dataPath\\$numFileName.txt")

    lineMisses = 0
    #Read in and pre-process the data
    for l in eachline(iStream)
      lineArr::Vector{String} = split(l,'\t')

      #makes this a formal NA statement
      while length(lineArr) < numColumns
        push!(lineArr, "NA")
      end
      if length(lineArr) == numColumns
        write(oStream,join(lineArr[1:numColumns],sep))
        write(oStream, '\n')
      else
        lineMisses += 1
      end
    end
    close(iStream)
    close(oStream)

    #writes the data frame object as a binary
    oStream = open("$dataPath\\$numFileName.jls", "w+")
    serialize(oStream, readtable("$dataPath\\$numFileName.dat", separator = sep))
    close(oStream)

    println("Pre-processed $numFileName.\n", lineMisses, " lines unprocessed\n")
  end

  #load the data from the binary
  iStream = open("$dataPath\\$numFileName.jls")
  finDat::DataFrame = deserialize(iStream)
  close(iStream)

  ##do parallel tasks with the data here
end

function loadFile(yearBegin::Int=2016, yearEnd::Int=2016; noParallel=true)
  try mkdir(dataFolder) end

  #for each row and year
  yearByQuarter::Matrix{Int} = Matrix{Int}((yearEnd-yearBegin + 1)*4,2)
  for y::Int in yearBegin:yearEnd, q::Int in 1:4
    yearByQuarter[(y-yearBegin)*4+q, 1] = y
    yearByQuarter[(y-yearBegin)*4+q, 2] = q
  end

  #preprocess and execute any parallel tasks
  if !noParallel
    pmap(_parTasks, yearByQuarter[:,1], yearByQuarter[:,2])
  end


  for y::Int in yearBegin:yearEnd, q::Int in 1:4
    #form the file names
    address::String = extractURL * "$y" * "q" * "$q.zip"
    fileName::String = "$y" * "q$q"
    numFileName::String = "$y" * "q$q" * "num"

    #load the data from the binary
    iStream = open("$dataPath\\$numFileName.jls")
    finDat::DataFrame = deserialize(iStream)
    close(iStream)

    ####do serial tasks with the data here

    println(finDat[1:3,:])

  end

end

function cleanworkers()
  myProcs = procs()
  if CLEAN_WORKERS
    for i in 2:length(myProcs)
      rmprocs(myProcs[i])
    end
  end
end

@time loadFile(noParallel = false)
cleanworkers()
