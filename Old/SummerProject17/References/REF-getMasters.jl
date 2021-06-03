
using Requests, FTPClient, GZip, DataFrames

const WORKERS = 9 #seems like logical cores * (3/4) is the magic value
const CLEAN_WORKERS = true

if nworkers() < WORKERS || (nworkers()==1 && WORKERS == 1)
  addprocs(nworkers()==1?WORKERS-nworkers()+1:WORKERS-nworkers())
end

@everywhere using DataFrames

@everywhere if !isdefined(:constDef)
  const constDef = true
  const indexUrl = """http://www.sec.gov/Archives/edgar/full-index/"""
  const rootUrl = """http://www.sec.gov/Archives/"""
  const dataDir = "SEC10QK\\"
  const typeListToGet = ["10-Q","10-K"]
  const numMasterCols = 5
end

function getMaster(yearBegin::Int=2016, yearEnd::Int=2016)
  try mkdir(dataDir)  end

  #get the master file for each quarter
  for y::Int in yearBegin:yearEnd, q::Int in 1:4
    address::String = indexUrl * "$y" * "/QTR" * "$q" * "/master.gz"
    fileName = pwd() * "\\" * dataDir * "\\$y" * "QTR$q" * "master"
    println(fileName)
    try
      download(address, fileName * ".gz")
    catch err
      println(err)
    end

    #open and extract
    gHandle = GZip.open(fileName * ".gz")
    write(fileName, readlines(gHandle,chomp=false))

    close(gHandle)

  end


end


function parseMasterAsCSV2(yearBegin::Int=2016, yearEnd::Int=2016)

  #open and parse the master file for each quarter
  for y::Int in yearBegin:yearEnd, q::Int in 1:4

    fileName = pwd() * "\\" * dataDir * "\\$y" * "QTR$q" * "master"

    inFile::IOStream = open(fileName)
    outFile::IOStream = open(fileName*".txt", "w+")


    #get the number of lines to skip. #Hardcodes thes tructure
    for l in eachline(inFile)
      if ismatch(Regex(".txt"), l)
        lineArr::Vector{String} = split(l,'|')
        if (isnumber(lineArr[1]) == true && length(lineArr) == numMasterCols
          && sum(map(x::String->length(strip(x))>0?1:0, lineArr))==numMasterCols)
          write(outFile,join(lineArr,'\t'))
          write(outFile, '\n')
        end
      end
    end


    close(inFile)
    close(outFile)

  end
end

#helper function for pmap. Requires the year and quarter
@everywhere function _downloadFromMaster(y::Int, q::Int)
  mstrName::String = pwd() * "\\" * dataDir * "$y" * "QTR$q" * "master.txt"
  logFileName::String = "$(Dates.format(now(),"yymmdd_THMS"))_PID$(myid())_Log.txt"
  logFile::IOStream = open(logFileName, "w+")
  typeList::Vector{String} = deepcopy(typeListToGet) #make a local copy
  numFilesDLed = 0
  numErrors = 0
  numFiles = 0

  println(logFile, "Begin Processing ",mstrName)

  master::DataFrame = readtable(mstrName, header=false,separator='\t')
  println("\nSpun up process ", myid(),
    "\nDownloading from ", mstrName,
    "\nLog name: ", logFileName)
  for i in 1:(size(master,1)) #assumes first column is the CIK code
    if sum(map(x::String -> x==master[i,3]?1:0,typeList)) > 0
      outName::String = pwd() * "\\" * dataDir * "\\$y" * "QTR$q" *"_" * "$(master[i,1])" * ".txt"
      try
        download(rootUrl * "$(master[i,end])", outName)
        numFilesDLed += 1
      catch err
        println(logFile, "ERROR: ", err, "\nFailed to download ", outName)
        numErrors += 1
      end
      numFiles += 1
    end
  end
  println(logFile, "Successfully downloaded ", mstrName,
    "\n Success Count: ", numFilesDLed, "\nFailure Count:", numErrors,
    "\nTotal Files:", numFiles)
  println("Successfully downloaded ", mstrName,
    "\n Success Count: ", numFilesDLed, "\nFailure Count:", numErrors,
    "\nTotal Files:", numFiles)
end

function getSECFiles(yearBegin::Int=2016, yearEnd::Int=2016)


  #for each row and year
  yearByQuarter::Matrix{Int} = Matrix{Int}((yearEnd-yearBegin + 1)*4,2)
  for y::Int in yearBegin:yearEnd, q::Int in 1:4
    yearByQuarter[(y-yearBegin)*4+q, 1] = y
    yearByQuarter[(y-yearBegin)*4+q, 2] = q
  end

  pmap(_downloadFromMaster , yearByQuarter[:,1], yearByQuarter[:,2])
  println("All downloads attempted")

  cleanworkers()

end

function cleanworkers()
  myProcs = procs()
  if CLEAN_WORKERS
    for i in 2:length(myProcs)
      rmprocs(myProcs[i])
    end
  end
end


#getMaster()
#parseMasterAsCSV2()
getSECFiles()
