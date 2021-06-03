using Revise, Pkg

#NOTE: Consider getting both indstl and fs records in ccm table

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

Pkg.activate(pwd())
using QFactors

#QFactors.recompressqfactorfiles()
#QFactors.recompressqfactorfile("CRSP-D")
#QFactors.recompressqfactorfile("CRSP-M")

#refreshcomp (PreprocessCOMP): refreshes the compustat data from the csv.lz4 file
#refreshcrsp (PreprocessCRSP): refreshes the crsp data from the csv.lz4 file
#refreshmerge[1] (MergeCRSPCOMP): refreshes the compustat-crsp merge.
#refreshfactors (FormFactors): refreshes the compustat-crsp merge.
#[1]True if refreshcomp, refreshcrsp, or validatemerged are true
#[2]True if refreshcomp, refreshcrsp, refreshmerge, or validatemerged are true

#=full flow (updated 7/27)
Key: function (FileName)
... means substantial code and/or subordinate (minor) function calls

formfactors (FormFactors)
  makeuniv (MergeCRSPCOMP)
    loadcrspdata (PreprocessCRSP)
      preprocesscrsp
        ...
    loadcompdata (PreprocessCOMP)
      preprocesscomp(2 arg)
        preprocesscomp(1 arg)
        ...
        dedupcomp!
        ...
        formannouncementintervals!
      mergecompccm (MergeCOMPCCM)
        loadccm
        preprocessccm!
        ...
        reconcilecompccmintervals!
    mergecrspcomp!
      testmerged
  formportfolios!


=#

function qfactorscript()
  QFactors.formfactors(
    refreshcomp = true,
    refreshcrsp = false,
    refreshmerge = true,
    validatemerged = true,
    refreshfactors = true,
    parallel = true)

  QFactors.analyzeqfactors()
end

@time qfactorscript()
