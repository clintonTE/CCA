using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end



using Endowment, DataFrames, Distributed, Finometrics, CodecZlib, CSV, GZip#, CSVFiles#, CSV

const WORKERS=5

const CLEAN_WORKERS=false

#set up for parallel operations
initializePar(WORKERS, reset=CLEAN_WORKERS)

#a list of high-level functions to run at the script level
function endowmentTasks()::Nothing
  #downloadIRSPDFs()
  #NOTE:
  #function in endowment data
  #RefreshData: Refreshes the data from NCCS. Takes a long time. Main func in downloadnccsdata
  #Validatedownloads: checks for valid data files for NCCS. Takes a long time.
  #  Checks for a few common fields. Main func in endowmentdata
  #refreshSummary: runs function consolidateSummary in NCCSSummary. Refreshes metadata files.
  #reserialize data: rewrites all teh data files as serial files. Takes a long time.

  #=metadata::Vector{MetaDictionary} = preProcessData(refreshData=false,
    validateDownloads=false,  refreshMetadata = true, reserializeData = true,
    refreshSummary = false)

  #NOTE:
  #refreshProcessed=refreshes the specific field mapping and some basic transformation rules
  #refreshNCCSDataCSV= Writes the data to a CSV. Careful- CSV is uncompressed and may be large
  NCCSData(metadata, refreshProcessed=true, refreshNCCSDataCSV=false, mode=:parallel)=#


  #NOTE:
  #refreshFilter= redo the processing of return, wealth, and similar computations
  #refreshDescriptiveFigures= draw basic return and wealth graphs
  #refreshDescriptiveTables=write summary and return tables
  #doOneOffTasks= executes onoff function in descriptiveFigures file
  #refreshβ= compute regressions of returns on S&P, contributions on S&P, retunrs on contribtuions etc
  #refreshSP500=re-read and merge S&P500 returns
  #refreshpersistence= return persistance plots
  #filternoncentral= filters out organizations classified as noncentral
  #  (values in the code dictionary file)
  #refreshtopinstitutions= prints the top institutions by wealth and worst perfoming institutions
  #refreshopen990data= refreshes the open990 data from the downloads (~60s parallel)
  reportandanalyze(refreshFilter = false,
    refreshopen990data = false,
    refreshDescriptiveFigures = false,
    refreshDescriptiveTables = false,
    doOneOffTasks = false,
    refreshβ=false,
    refreshoutsidedata=false,
    refreshpersistence=false,
    refreshtopinstitutions=false,
    refreshsensitivity=true)

  return nothing
end

@time endowmentTasks()

#@time getfielddf(refreshopen990data=true)


#cleanup the extra processes if desired
if CLEAN_WORKERS
  cleanup()
end
