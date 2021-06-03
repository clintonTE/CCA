using Revise
using Pkg
#TODO 1Run soup2nuts 2) continue reconciliation

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using LeverageChanges

@warn("I'm not sure the fyear always corresponds to the actual fiscal year end-
  if so, I should investigate using ddate instead. Seemed to work well in capacity.
    Also check out the ccm reconciliation- make sure the endpoints work")
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
    table1 = (false || SOUP_2_NUTS),
    table2 = (false || SOUP_2_NUTS),
    table3 = (false || SOUP_2_NUTS),
    table4 = (false || SOUP_2_NUTS),
    table5 = (false || SOUP_2_NUTS),
    table6 = (false || SOUP_2_NUTS),

    #event analysis
    disttype=:primaryfull, #valid=[:primaryfull, :primaryshort, :oldfull, :oldshort]
    refreshseocrsp=(false  || SOUP_2_NUTS),
    refreshdistcrsp=(false  || SOUP_2_NUTS),
    refreshevent=(false  || SOUP_2_NUTS),
    refreshcumex=(false|| SOUP_2_NUTS),

    table7 = (false || SOUP_2_NUTS),
    table8 = (false|| SOUP_2_NUTS),
    figure1 = (false  || SOUP_2_NUTS),
    figure2 = (false  || SOUP_2_NUTS),
    table9 = (false || SOUP_2_NUTS),
    figure3 = (false  || SOUP_2_NUTS),
    table10 = (false  || SOUP_2_NUTS),
    figure4 = (false  || SOUP_2_NUTS),

    tablea1 = (false || SOUP_2_NUTS),
    tablea2 = (false  || SOUP_2_NUTS),

    yearrange = 1961:2013,
    outsuffix = "all",
    writepaneltocsv=false)

    println(LeverageChanges.YEAR_RANGE[])


end

#LeverageChanges.recompressfile("CCM")
@time leveragechangesscript()

#using Dates
#function gendatedict(yearrange::AbstractRange=1961:2018, period::DatePeriod)
#  mindate = Date()
