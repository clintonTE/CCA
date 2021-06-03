
function replicate(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  refreshdataseries::Bool = true,
  refreshportfolios::Bool = true,
  parallel::NBool = nothing,
  trialrunonly::Bool = false,
  table1::Bool=true,
  table2::Bool=true,
  table3::Bool=true,
  table4::Bool=true,
  table5::Bool=true,
  table6::Bool=true,
  refresheventcrsp::Bool=true,
  refreshevent::Bool=true)::Nothing

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
  elseif table1 + table2 + table3 + table4 + table5 + table6 > 0
    panel = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  end


  table1 && maketable1(panel, trialrunonly=trialrunonly)
  table2 && maketable2(panel, trialrunonly=trialrunonly)
  table3 && maketable3(panel, trialrunonly=trialrunonly)
  table4 && maketable4(panel, trialrunonly=trialrunonly)
  table5 && maketable5(panel, trialrunonly=trialrunonly)
  table6 && maketable6(panel, trialrunonly=trialrunonly)

  refreshevent && formevent(refresheventcrsp = refresheventcrsp)

  return nothing
end
