
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
  yearrange::Union{UnitRange{Int},Nothing} = nothing,
  outsuffix::Union{Nothing,String} = nothing,

  disttype::Symbol=:full,
  table7::Bool=true,
  table8::Bool=true,
  refreshseocrsp::Bool=true,
  refreshdistcrsp::Bool=true,
  refreshevent::Bool=true,
  refreshcumex::Bool= true,
  figure1::Bool = true,
  figure2::Bool = true,
  figure3::Bool = true,
  figure4::Bool = true,
  table9::Bool=true,
  table10::Bool= true,

  tablea1::Bool = true,
  tablea2::Bool = true)::Nothing


  local panel::DataFrame
  local bb::DataFrame
  local dist::DataFrame
  local cumex::DataFrame

  (!isnothing(parallel)) && (PARALLEL[] = parallel)
  REPLICATION_TYPE[] = trialrunonly ? :monthly : :daily
  (!isnothing(yearrange)) && (YEAR_RANGE[] = yearrange)

  (disttype ∈ [:primaryfull, :primaryshort, :oldfull, :oldshort]) || error("disttype must be ∈ [:primaryfull, :primaryshort, :oldfull, :oldshort]")
  DIST_TYPE[] = disttype

  (!isnothing(outsuffix)) && (OUT_SUFFIX[] = outsuffix)

  #not strictly necessary, but helps stop doubling up on the multi-threading
  @info "Beginning replication"


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
  elseif table1 + table2 + table3 + table4 + table5 + table6 + tablea1 > 0
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
  println("Beginning Event Analysis")
  (bb,dist) = formevent(refreshevent=refreshevent,
    refreshseocrsp=refreshseocrsp,
    refreshdistcrsp=refreshdistcrsp)
  cumex = formcumex(dist, refreshcumex = refreshcumex)
  table7 && maketable7(bb)
  table8 && maketable8(dist)
  figure1 && makefigure1(dist)
  figure2 && makefigure2(dist)
  table9 && maketable9(cumex)
  figure3 && makefigure3(cumex)
  table10 && maketable10(cumex)
  figure4 && makefigure4(dist)

  #####################################appendix
  @info "Begining appendix tables"
  tablea1 && maketablea1(panel, trialrunonly=trialrunonly)
  tablea2 && maketablea2(dist)

  return nothing
end
