#this file contains a list of file specific functions which each produce a dataframe
#The objective here is to minimize lines of code specific to each file
#Each function takes in the input dataframe only. A linker function selects
#the correct function from a symbol and a constant lookup table

#the function naming convention is [filetype]_[year(s)], where all refers to all for the type
#examples:
# all_all, pf_all, pf_2003-2009_2011
# refers to all files, all private foundation files,
# and the 2003-2009 and 2011 private foundation files
# the first function is always all_all, which is the only one with no
# dataframe argument

#holds the funciton linkings necessary to process the data
const PROC_FUNCTIONS = Dict{String, Vector{Function}}(
  "nccs_core_2015_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015, pf_2011_2015],
  "nccs_core_2014_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015, pf_2011_2015],
  "nccs_core_2013_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015, pf_2011_2015],
  "nccs_core_2012_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015, pf_2011_2015],
  "nccs_core_2011_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015, pf_2011_2015],
  "nccs_core_2010_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2009_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2008_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2007_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2006_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2005_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2004_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2003_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2002_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2001_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_2000_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_1999_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_1998_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015, pf_1998_2015],
  "nccs_core_1997_pf" => [all_all, pf_all, pf_1992_2015, pf_1997_2015],
  "nccs_core_1996_pf" => [all_all, pf_all, pf_1992_2015],
  "nccs_core_1995_pf" => [all_all, pf_all, pf_1992_2015],
  "nccs_core_1994_pf" => [all_all, pf_all, pf_1992_2015],
  "nccs_core_1992_pf" => [all_all, pf_all, pf_1992_2015],
  "nccs_core_1991_pf" => [all_all, pf_all],
  "nccs_core_1990_pf" => [all_all, pf_all],
  "nccs_core_1989_pf" => [all_all, pf_all],
  "nccs_core_2015_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2014_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2013_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2012_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2011_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2010_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2009_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015, pc_2009_2015],
  "nccs_core_2008_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2007_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2006_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2005_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2004_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2003_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2002_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2001_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_2000_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_1999_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_1998_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_1997_pc" => [all_all, pc_all, pc_1991_2015, pc_1997_2015],
  "nccs_core_1996_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1995_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1994_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1993_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1992_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1991_pc" => [all_all, pc_all, pc_1991_2015],
  "nccs_core_1990_pc" => [all_all, pc_all],
  "nccs_core_1989_pc" => [all_all, pc_all],
  "nccs_core_2015_co" => [all_all, co_all, co_2000_2015, co_2011_2015],
  "nccs_core_2014_co" => [all_all, co_all, co_2000_2015, co_2011_2015],
  "nccs_core_2013_co" => [all_all, co_all, co_2000_2015, co_2011_2015],
  "nccs_core_2012_co" => [all_all, co_all, co_2000_2015, co_2011_2015],
  "nccs_core_2011_co" => [all_all, co_all, co_2000_2015, co_2011_2015],
  "nccs_core_2010_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2009_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2008_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2007_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2006_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2005_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2004_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2003_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2002_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2001_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_2000_co" => [all_all, co_all, co_2000_2015],
  "nccs_core_1999_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1998_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1997_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1996_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1995_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1994_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1993_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1992_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1991_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1990_co" => [all_all, co_all, co_1989_1999],
  "nccs_core_1989_co" => [all_all, co_all, co_1989_1999]
)


#holds the destination column types
#const DATA_TYPES = DICT{Symbol, Type}(:ein=>Symbol, :fisyr=>Int)



function processNCCS(data::NCCSFile)::DataFrame
  processFuncs::Vector{Function} = PROC_FUNCTIONS[data.meta.fileName]

  #println(data.df[1:3, :name])

  #run the initial function
  outDF::DataFrame = processFuncs[1](data)

  for f::Function ∈ processFuncs[2:end]
    outDF = f(data, outDF)
  end

  return outDF
end

#convenience method since loading many files is memory intensive
processNCCS(meta::MetaDictionary) =
  processNCCS(loadNCCSFile(meta, reserializeData=false, returnNCCSFile=true))::DataFrame


#this function coordinates the processing of the data files
function processNCCS(metadata::Vector{MetaDictionary}, mode::Symbol)::DataFrame

  local outDF::DataFrame #contains the merged dataframe for output

  if mode == :serial #terminal branch, summarize the data

    #apply the transformation rules and concatenate
    fileOutDFs::Vector{DataFrame} = (processNCCS).(metadata)
    outDF = vcat(fileOutDFs...)

  elseif mode == :parallel #relaunch the dedicated parallel function
    numFiles::Int = length(metadata)
    pids = workers()
    np::Int = min(length(pids), MAX_MAPPING_WORKERS)

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
        #gets the metadata index for all files to be processed for a given pid
        push!(futures,
          (@spawnat pids[i] processNCCS(metadata[assignments[i]], :serial)))
      end
    end
    outDF = vcat((fetch).(futures)...) #concatenate the results
    return outDF
  else
    @error "Incorrect mode for function processNCCS"
  end


  return outDF

end

#serializes a NCCSData object
function serializeNCCSData(data::NCCSData;
    dataPath::String = DATA_PATH,
    dataName::String = PROCESSED_DATA_NAME,
    oStreamer::Function = O_JLS_STREAMER)::Nothing

    outPath::String = "$dataPath\\$(dataName).jls.gz"

    oStream = oStreamer(outPath)
    serialize(oStream, data)
    close(oStream)

    return nothing
end

#serializes a NCCSData object
function deserializeNCCSData(;
    dataPath::String = DATA_PATH,
    dataName::String = PROCESSED_DATA_NAME,
    iStreamer::Function = I_JLS_STREAMER)::NCCSData

    inPath::String = "$dataPath\\$(dataName).jls.gz"

    iStream = iStreamer(inPath)
    data::NCCSData = deserialize(iStream)
    close(iStream)

    return data
end

#serializes a NCCSData object
#for now, assume the data are uncompressed
function writeNCCSDataCSV(data::NCCSData;
    processedDataName::String = PROCESSED_DATA_NAME,
    dataPath::String = DATA_PATH)::Nothing

    outPath::String = "$dataPath\\$(processedDataName).csv"

    CSV.write(outPath, data.df)
    return nothing
end

#this is the intial funciton
#main purpose is to filter the metadictionaries to only files
#we are equipped to process
function NCCSData(metadata::Vector{MetaDictionary};
  mode::Symbol=:parallel, refreshProcessed::Bool=true, dataPath::String = DATA_PATH,
  processedDataName::String = PROCESSED_DATA_NAME,
  refreshNCCSDataCSV::Bool = false)::NCCSData

  if refreshProcessed || !isfile("$dataPath\\$processedDataName.jls.gz")

    #get the index of files
    metaIndex::Dict = meta2Index(metadata)

    #get the list of files where we have the transformation rules
    filesToProcess::Vector{String} = collect(keys(PROC_FUNCTIONS))
    metadata = metadata[(s::String->metaIndex[s]).(filesToProcess)]

    #now write out the data
    outDF::DataFrame = processNCCS(metadata, mode)
    data::NCCSData = NCCSData(outDF)
    serializeNCCSData(data)
  else
    data = deserializeNCCSData()
  end

  #since this is the main datafile, we may want to write it out
  if refreshNCCSDataCSV
    writeNCCSDataCSV(data)
  end



  return data
end
