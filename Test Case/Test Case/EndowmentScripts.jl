module EndowmentScripts

const Nothing = Void #NOTE: 0.7 compatabiltiy hack, delete on ugrade

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


using Endowment, CSV, DataFrames

const WORKERS=0

const CLEAN_WORKERS=false

#set up for parallel operations
initializePar(WORKERS, reset=CLEAN_WORKERS)


#a list of high-level functions to run at the script level
function endowmentTasks()::Nothing
  preProcessData(refreshData=false, validateDownloads=false,
    refreshMetadata = true, reserializeData = true)

  return nothing
end

@time endowmentTasks()


#cleanup the extra processes if desired
if CLEAN_WORKERS
  cleanup()
end
end
