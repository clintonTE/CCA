using Revise, Pkg

#TODO: Check daily data, debug whatever is wrong w/ the table
#if this is hard maybe do the summary table and walk from there
#other idea is to try valuewieghting

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using LeverageChanges

#LeverageChanges.recompressfiles()
#LeverageChanges.recompressfile("CRSP-D")
#LeverageChanges.recompressfile("CRSP-M")

#refreshcomp (PreprocessCOMP): refreshes the compustat data from the csv.lz4 file
#refreshcrsp (PreprocessCRSP): refreshes the crsp data from the csv.lz4 file
#refreshmerge[1] (MergeCRSPCOMP): refreshes the compustat-crsp merge.
#refreshfactors (FormFactors): refreshes the compustat-crsp merge.
#[1]True if refreshcomp, refreshcrsp, or validatemerged are true
#[2]True if refreshcomp, refreshcrsp, refreshmerge, or validatemerged are true
#hello
#usequarterlydates

function leveragechangesscript()
  #LeverageChanges.recompressfile("CCM-A")
  #LeverageChanges.recompressfile("CRSP-D-I")
  #LeverageChanges.recompressfile("CRSP-D-II")
  #LeverageChanges.recompressfile("CRSP-M")
  #LeverageChanges.recompressfile("CCM")
  #LeverageChanges.recompressfile("COMP-A")

  LeverageChanges.replicate(
    refreshcomp = false,
    refreshtvolmcapcrsp = false,
    refreshtretcrsp = false,
    refreshdataseries=false,
    refreshportfolios=false,
    parallel = false,
    trialrunonly = false,
    table1 = false,
    table2 = false,
    table3 = false,
    table4 = false,
    table5 = false,
    table6 = false,
    refreshevent=true,
    refresheventcrsp=true)

  #LeverageChanges.formevent()
end

#LeverageChanges.archivefile("lev-daily.jls.lz4")
#LeverageChanges.unarchivefile("lev-daily.jls.lz4")
#LeverageChanges.recompressfile("CCM")
@time leveragechangesscript()
#@time LeverageChanges.formevent(refresheventcrsp = true)

#winsorize vol252d? Maybe mess with months2stale short, or otherwise
#adjust for gaps in the returns/benchmark
