
function replicate(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  refreshdataseries::Bool = true,
  refreshportfolios::Bool = true,
  parallel::NBool = nothing,
  trialrunonly::Bool = false,
  table1=true)::Nothing

  local panel::DataFrame

  (!isnothing(parallel)) && (PARALLEL[] = parallel)
  REPLICATION_TYPE[] = trialrunonly ? :monthly : :daily

  #not strictly necessary, but helps stop doubling up on the multi-threading


  #this encompasses all pre-processing activities
  if sum([refreshcomp, refreshtretcrsp, refreshtvolmcapcrsp,
    refreshdataseries, refreshportfolios]) > 0

    panel = constructportfolios(
      refreshcomp = refreshcomp,
      refreshtvolmcapcrsp = refreshtvolmcapcrsp,
      refreshtretcrsp = refreshtretcrsp,
      refreshdataseries=refreshdataseries,
      refreshportfolios=refreshportfolios,
      trialrunonly = trialrunonly)
  end

  table1 && maketable1(trialrunonly=trialrunonly)

  return nothing
end
