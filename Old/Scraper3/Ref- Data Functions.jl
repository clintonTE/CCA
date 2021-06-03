
module HFMod
using  DataFrames, Distributions, StatsBase, GZip, JLD#, ArrayFire (Needs to be updated to 0.6)

if !isdefined(:constDef)
  const constDef = true
  #other constants here

  const debugHFMod = true #switch for performance impairing checks
  const dataPath = pwd() * "\\data"
  const dataFileName = "HFPerformance_2015-2017"

end

function PreProcessHFData()


  #this makes sure we have a binary of the data (Improves load times 3-4x)
  if !isfile("$dataPath\\$dataFileName.jls")
    #extract from the zip file
    gHandle::GZipStream = GZip.open("$dataPath\\$dataFileName.gz")
    write("$dataPath\\$dataFileName.csv", readlines(gHandle,chomp=false))
    close(gHandle)

    #write the binary
    oStream::IOStream = open("$dataPath\\$dataFileName.jls","w+")
    serialize(oStream, readtable("$dataPath\\$dataFileName.csv"))
    close(oStream)
  end

  #load the binary
  iStream::IOStream = open("$dataPath\\$dataFileName.jls")
  HFData::DataFrame = deserialize(iStream)
  close(iStream)

  #type conversions here
  #=for i = 1:size(HFData,1)
    HFData[i,[:]:date] = Date("$(HFData[i,:date])",Dates.DateFormat("yyyymmdd"))
  end=#
  rename!(HFData,:date,:date_old)
  HFData[:,:date] = map((x::Int)->Date("$x",Dates.DateFormat("yyyymmdd")),HFData[:,:date_old])

  println(HFData[1:3,:])

  #save the binary
  oStream = open("$dataPath\\$(dataFileName)_clean.jls","w+")
  serialize(oStream, HFData)
  close(oStream)

  #println(HFData[])


end

@time PreProcessHFData()

end
