


function maketable2(panel::DataFrame; trialrunonly::Bool=trialrunonly,
    runpanela::Bool=true,
    runpanelb::Bool=true,
    runpanelc::Bool=true,
    runpaneld::Bool=true
    )::Nothing
  #load the data we will work with
  #panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  sort!(panel, [:permno, :date])
  apanel::DataFrame = panel[panel.annflag,:]

  ################common data
  cqports::PortfolioConstructions = leverageconstructions()
  T2_FNAMES::Vector{Symbol} = CQUART_FNAMES_D
  T2_FVALS::Vector{Symbol} = CQUART_FVALS_D
  T2_FGROUPS::Vector{Symbol} = CQUART_FGROUPS_D
  T2_FNAMES_ALT::Vector{Symbol} = CQUART_FNAMES_ALT_D
  T2_FGROUPS_ALT::Vector{Symbol} = CQUART_FGROUPS_ALT_D

  ################Table A
  T2_PANEL_A_ROWKEYS::Vector{Symbol} = [:mean, :std, :stdxt]

  T2_PANEL_A_ROWLABEL::Dict = Dict(
    :mean=>"Mean",
    :std=>"SD",
    :stdxt=>"SD(X,T)")


  T2_PANEL_A_IDX::Dict = Dict(
    :mean=>1,
    :std=>2,
    :stdxt=>3)
  T2_PANEL_A_NAME::String = "t2panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"


  runpanela && tablexpanelA(apanel,
    fvals = T2_FVALS,
    rowkeys=T2_PANEL_A_ROWKEYS,
    rowlabel=T2_PANEL_A_ROWLABEL,
    idx=T2_PANEL_A_IDX,
    tablename=T2_PANEL_A_NAME)


  #####################TABLE B
  T2_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :leverageQmeans,
    :quartilenum,
    :LLsd,
    :sd,
    :Δsd,
    :LLret12m,
    :ret12m,
    :Δret12m]

  T2_PANEL_B_ROWLABEL::Dict = Dict(
    :leverageQmeans=>"Leverage Q Means",
    :quartilenum=>"Quartile \\#",
    :LLsd=>"\\smallskip SD Net Lagged\$^2\$",
    :sd=>"SD Net Lead",
    :Δsd=>"SD Net 2Y Delta",
    :LLret12m=>"\\smallskip Compound Raw Lagged\$^2\$",
    :ret12m=>"Compound Raw Lead",
    :Δret12m=>"Compound Raw 2Y Delta")


  T2_PANEL_B_IDX::Dict = Dict(
    T2_PANEL_B_ROWKEYS[i]=>i for i ∈ 1:(length(T2_PANEL_B_ROWKEYS)))

  #the below is to help with future abstraction if necessary
  T2_PANEL_B_VALCOLIDX::Dict= Dict(f=>f for f ∈ T2_FVALS)

  #these are pairs of lagged and unlagged fields
  T2_PANEL_B_VOLFIELDS::NTuple{2,Symbol} = (:LLvol252dnet, :vol252dnet)
  T2_PANEL_B_RETFIELDS::NTuple{2,Symbol} = (:LLret12m, :ret12m)

  T2_PANEL_B_NAME::String = "t2panelb-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  runpanelb && tablexpanelB(apanel,
    fvals = T2_FVALS,
    rowkeys=T2_PANEL_B_ROWKEYS,
    rowlabel=T2_PANEL_B_ROWLABEL,
    idx=T2_PANEL_B_IDX,
    valcolidx = T2_PANEL_B_VALCOLIDX,
    tablename=T2_PANEL_B_NAME,
    volfields=T2_PANEL_B_VOLFIELDS,
    retfields=T2_PANEL_B_RETFIELDS)

  ###################PANEL C############


  T2_PANEL_C_ROWKEYS::Vector{Symbol} = [:constant, :plusxmkt, :plussmbhml, :plusrmwcma, :plusumd]

  T2_PANEL_C_ROWLABEL::Dict = Dict(
    :constant=>"Constant",
    :plusxmkt=>"+XMKT",
    :plussmbhml=>"+SMB+HML",
    :plusrmwcma=>"+RMW+CMA",
    :plusumd=>"+UMD")

  T2_PANEL_C_CONTROLS::Dict = Dict(
    :constant=>CONTROL0,
    :plusxmkt=>[:mkt],
    :plussmbhml=>[:mkt, :P_ff3m_smb,	:P_ff3m_hml],
    :plusrmwcma=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma],
    :plusumd=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd])

  T2_PANEL_C_FOCALS::Dict = Dict(
    :constant=>T2_FNAMES,
    :plusxmkt=>T2_FNAMES,
    :plussmbhml=>T2_FNAMES,
    :plusrmwcma=>T2_FNAMES,
    :plusumd=>T2_FNAMES)

  T2_PANEL_C_NAME::String = "t2panelc-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T2_PANEL_C_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpanelc && tablexpanelC(panel,
    fvals=T2_FVALS,
    fgroups=T2_FGROUPS,
    Fnames = T2_FNAMES,
    rowkeys=T2_PANEL_C_ROWKEYS,
    rowlabel=T2_PANEL_C_ROWLABEL,
    controls=T2_PANEL_C_CONTROLS,
    focals=T2_PANEL_C_FOCALS,
    tablename=T2_PANEL_C_NAME,
    minrowspergroup = T2_PANEL_C_MIN_ROWS
  )

  #########################PANEL D#########################
  T2_PANEL_D_CONTROLS::Dict = Dict(
    :quartiles=>CONTROL0,
    :quartiles4ctrl=>T1_CONTROL4,
    :quartiles10ctrl=>T1_CONTROL10,
    :values=>CONTROL0,
    :values4ctrl=>T1_CONTROL4,
    :values10ctrl=>T1_CONTROL10,
    :quartiles_ac=>CONTROL0,
    :quartiles_lret=>T1_CONTROL10)

  T2_PANEL_D_RET::Dict = Dict(
    :quartiles=>:ret,
    :quartiles4ctrl=>:ret,
    :quartiles10ctrl=>:ret,
    :values=>:ret,
    :values4ctrl=>:ret,
    :values10ctrl=>:ret,
    :quartiles_ac=>:ret,
    :quartiles_lret=>:lret)

  T2_PANEL_D_ROWKEYS::Vector{Symbol} = [:quartiles, :quartiles4ctrl, :quartiles10ctrl, :values,
    :values4ctrl, :values10ctrl, :quartiles_ac, :quartiles_lret]

  T2_PANEL_D_ROWLABEL::Dict = Dict(
    :quartiles=>"Leverage Quartiles",
    :quartiles4ctrl=>"... +4 vars",
    :quartiles10ctrl=>"... +10 vars",
    :values=>"Leverage Values",
    :values4ctrl=>"... +4 vars",
    :values10ctrl=>"... +10 vars",
    :quartiles_ac=>"\\midrule Asset-Controlled Quartiles",
    :quartiles_lret=>"Log Returns +10 vars")

  T2_PANEL_D_FOCALS::Dict = Dict(
    :quartiles=>T2_FGROUPS,
    :quartiles4ctrl=>T2_FGROUPS,
    :quartiles10ctrl=>T2_FGROUPS,
    :values=>T2_FVALS,
    :values4ctrl=>T2_FVALS,
    :values10ctrl=>T2_FVALS,
    :quartiles_ac=>T2_FGROUPS_ALT,
    :quartiles_lret=>T2_FGROUPS)

  T2_PANEL_D_NAME::String = "t2paneld-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T2_PANEL_D_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpaneld && tablexpanelD(panel,
      fvals=T2_FVALS,
      fgroups=T2_FGROUPS,
      ret = T2_PANEL_D_RET,
      rowkeys=T2_PANEL_D_ROWKEYS,
      rowlabel=T2_PANEL_D_ROWLABEL,
      controls=T2_PANEL_D_CONTROLS,
      focals=T2_PANEL_D_FOCALS,
      tablename=T2_PANEL_D_NAME,
      minrowspergroup = T2_PANEL_D_MIN_ROWS
    )

  return nothing
end
