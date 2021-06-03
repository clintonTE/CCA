module EndowmentScripts

const Nothing = Void #NOTE: 0.7 compatabiltiy hack, delete on ugrade

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


using Endowment, CSV, CodecZlib, DataFrames

const WORKERS=2
const CLEAN_WORKERS=false

#set up for parallel operations
initializePar(WORKERS, reset=CLEAN_WORKERS)


#a list of high-level functions to run at the script level
function endowmentTasks()::Nothing
  preProcessData(update=true, validateDownloads=false, refreshMetadata = false, refreshData = false)

  return nothing
end

@time endowmentTasks()

#cleanup the extra processes if desired
if CLEAN_WORKERS
  cleanup()
end
end
