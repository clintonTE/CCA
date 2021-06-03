using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end



using Endowment, DataFrames, Distributed, Finometrics, CodecZlib, CSV, GZip#, CSVFiles#, CSV

const WORKERS=1

const CLEAN_WORKERS=false

#set up for parallel operations
initializePar(WORKERS, reset=CLEAN_WORKERS)

#a list of high-level functions to run at the script level
function endowmentTasks()::Nothing
  #downloadIRSPDFs()

  #metadata::Vector{MetaDictionary} = preProcessData(refreshData=false,
  #  validateDownloads=false,  refreshMetadata = false, reserializeData = false,
  #  refreshSummary = true)

  #NCCSData(metadata, refreshProcessed=true, refreshNCCSDataCSV=false, mode=:parallel)
  #NOTE: Persistance plots= return persistance plots
  makeDescriptive(refreshFilter = false, refreshDescriptiveFigures = false,
    refreshDescriptiveTables = false, doOneOffTasks = false, refreshβ=false,
    refreshSP500=false, refreshPersistancePlots=true)

  return nothing
end

@time endowmentTasks()


#cleanup the extra processes if desired
if CLEAN_WORKERS
  cleanup()
end
