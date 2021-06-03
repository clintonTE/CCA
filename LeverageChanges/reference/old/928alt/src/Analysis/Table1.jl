

function maketable1(;trialrunonly::Bool=trialrunonly)::Nothing
  #load the data we will work with
  panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  table1panelD(panel)

  return nothing
end


#NOTE: We require
function table1panelD(panel::DataFrame)
  local cqports::CQuartileConstructions = leverageconstructions()
  Fgrps = Fgroups(cqports)

  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(T1_PANEL_D_ROWNAME) * length(Fgrps))

  #make the regression specifications
  @time for rname ∈ T1_PANEL_D_ROWNAME
    for grp ∈ Fgrps
      push!(fss, FactorSpecification(panel, T1_PANEL_D_RET[rname], :date, :permno,
        Ffactors=nothing,
        Fcharacteristics = [grp],
        Fcontrols = T1_PANEL_D_CONTROLS[rname],
        name=Symbol("Row$(rname)_Col$(grp)"),
        skiptimeseriesdfcheck=true)) #only doing crossectional regressions
    end
  end

  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]
  frs::Vector{FactorRegression} = (fs->characteristicFM(fs, parallel=false)).(fss)

end
