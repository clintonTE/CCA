
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

  table7::Bool=true,
  table8::Bool=true,
  refreshseocrsp::Bool=true,
  refreshdistcrsp::Bool=true,
  refreshevent::Bool=true,
  figure1::Bool = true)::Nothing

  local panel::DataFrame
  local bb::DataFrame
  local dist::DataFrame

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

  #######################panel study
  table1 && maketable1(panel, trialrunonly=trialrunonly)
  table2 && maketable2(panel, trialrunonly=trialrunonly)
  table3 && maketable3(panel, trialrunonly=trialrunonly)
  table4 && maketable4(panel, trialrunonly=trialrunonly)
  table5 && maketable5(panel, trialrunonly=trialrunonly)
  table6 && maketable6(panel, trialrunonly=trialrunonly)


  ##################################event study
  println("got here")
  (bb,dist) = formevent(refreshevent=refreshevent,
  refreshseocrsp=refreshseocrsp,
  refreshdistcrsp=refreshdistcrsp)
  table7 && maketable7(bb)
  table8 && maketable8(dist)
  figure1 && makefigure1(dist)


  return nothing
end
