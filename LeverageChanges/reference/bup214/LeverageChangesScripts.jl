using Revise, Pkg
#TODO continue reconciliation
#TODO valide new merging algorithm- make sure tables don't change


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
#


function leveragechangesscript()
  #LeverageChanges.recompressfile("CCM-A")
  #LeverageChanges.recompressfile("CRSP-D-I")
  #LeverageChanges.recompressfile("CRSP-D-II")
  #LeverageChanges.recompressfile("CRSP-M")
  #LeverageChanges.recompressfile("CCM")
  #LeverageChanges.recompressfile("COMP-A")

  SOUP_2_NUTS::Bool = false

  LeverageChanges.replicate(
    refreshcomp = (false || SOUP_2_NUTS),
    refreshtvolmcapcrsp = (false || SOUP_2_NUTS),
    refreshtretcrsp = (false || SOUP_2_NUTS),
    refreshdataseries=(false || SOUP_2_NUTS),
    refreshportfolios=(true || SOUP_2_NUTS),
    parallel = (true || SOUP_2_NUTS),
    trialrunonly = false, #WARNING MANUAL set to false on soup-2-nuts run

    #full set analysis tables
    table1 = (true || SOUP_2_NUTS),
    table2 = (true || SOUP_2_NUTS),
    table3 = (true || SOUP_2_NUTS),
    table4 = (true || SOUP_2_NUTS),
    table5 = (true || SOUP_2_NUTS),
    table6 = (true || SOUP_2_NUTS),

    #event analysis
    disttype=:primaryfull, #valid=[:primaryfull, :primaryshort, :oldfull, :oldshort]
    refreshseocrsp=(true  || SOUP_2_NUTS),
    refreshdistcrsp=(true  || SOUP_2_NUTS),
    refreshevent=(true  || SOUP_2_NUTS),
    refreshcumex=(true|| SOUP_2_NUTS),

    table7 = (true || SOUP_2_NUTS),
    table8 = (true || SOUP_2_NUTS),
    figure1 = (true  || SOUP_2_NUTS),
    figure2 = (true  || SOUP_2_NUTS),
    table9 = (true || SOUP_2_NUTS),
    figure3 = (true  || SOUP_2_NUTS),
    table10 = (true  || SOUP_2_NUTS),
    figure4 = (true  || SOUP_2_NUTS),

    tablea1 = (true || SOUP_2_NUTS),
    tablea2 = (true  || SOUP_2_NUTS),

    yearrange = 1961:2019,
    outsuffix = "all",
    writepaneltocsv=false)

    println(LeverageChanges.YEAR_RANGE[])


end
#@time leveragechangesscript()


#=LeverageChanges.recompressfile("CCM")
LeverageChanges.recompressfile("event\\crspdistributions")
LeverageChanges.recompressfile("CRSP-D")
LeverageChanges.recompressfile("COMP-A")
LeverageChanges.recompressfile("COMP-Q")=#
#LeverageChanges.recompressfile("CRSP-M")

function leveragechangestestscript()

  panel = LeverageChanges.formdataseries(
    refreshcomp = false,
    refreshtvolmcapcrsp = false,
    refreshtretcrsp = false,
    refreshdataseries=true,
    trialrunonly = false)
end

@time leveragechangestestscript()
