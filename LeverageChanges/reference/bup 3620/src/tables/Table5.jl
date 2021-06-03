

function maketable5(panel::DataFrame; trialrunonly::Bool=trialrunonly,
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
  T5_FNAMES::Vector{Symbol} = CQUART_FNAMES_RMD
  T5_FVALS::Vector{Symbol} = CQUART_FVALS_RMD
  T5_FGROUPS::Vector{Symbol} = CQUART_FGROUPS_RMD
  T5_FNAMES_ALT::Vector{Symbol} = CQUART_FNAMES_ALT_RMD
  T5_FGROUPS_ALT::Vector{Symbol} = CQUART_FGROUPS_ALT_RMD

  ################Table A
  T5_PANEL_A_ROWKEYS::Vector{Symbol} = [:mean, :std, :stdxt]

  T5_PANEL_A_ROWLABEL::Dict = Dict(
    :mean=>"Mean",
    :std=>"SD",
    :stdxt=>"SD(X,T)")


  T5_PANEL_A_IDX::Dict = Dict(
    :mean=>1,
    :std=>2,
    :stdxt=>3)
  T5_PANEL_A_NAME::String = "t5panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"


  runpanela && tablexpanelA(apanel,
    fvals = T5_FVALS,
    rowkeys=T5_PANEL_A_ROWKEYS,
    rowlabel=T5_PANEL_A_ROWLABEL,
    idx=T5_PANEL_A_IDX,
    tablename=T5_PANEL_A_NAME)


  #####################TABLE B
  T5_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :leverageQmeans,
    :quartilenum,
    :LLsd,
    :sd,
    :Δsd,
    :LLret12m,
    :ret12m,
    :Δret12m]

  T5_PANEL_B_ROWLABEL::Dict = Dict(
    :leverageQmeans=>"Leverage Q Means",
    :quartilenum=>"Quartile \\#",
    :LLsd=>"\\smallskip SD Net Lagged\$^2\$",
    :sd=>"SD Net Lead",
    :Δsd=>"SD Net 2Y Delta",
    :LLret12m=>"\\smallskip Compound Raw Lagged\$^2\$",
    :ret12m=>"Compound Raw Lead",
    :Δret12m=>"Compound Raw 2Y Delta")


  T5_PANEL_B_IDX::Dict = Dict(
    T5_PANEL_B_ROWKEYS[i]=>i for i ∈ 1:(length(T5_PANEL_B_ROWKEYS)))

  #the below is to help with future abstraction if necessary
  T5_PANEL_B_VALCOLIDX::Dict= Dict(f=>f for f ∈ T5_FVALS)

  #these are pairs of lagged and unlagged fields
  T5_PANEL_B_VOLFIELDS::NTuple{2,Symbol} = (:LLvol252dnet, :vol252dnet)
  T5_PANEL_B_RETFIELDS::NTuple{2,Symbol} = (:LLret12m, :ret12m)

  T5_PANEL_B_NAME::String = "t5panelb-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  runpanelb && tablexpanelB(apanel,
    fvals = T5_FVALS,
    rowkeys=T5_PANEL_B_ROWKEYS,
    rowlabel=T5_PANEL_B_ROWLABEL,
    idx=T5_PANEL_B_IDX,
    valcolidx = T5_PANEL_B_VALCOLIDX,
    tablename=T5_PANEL_B_NAME,
    volfields=T5_PANEL_B_VOLFIELDS,
    retfields=T5_PANEL_B_RETFIELDS)

  ###################PANEL C############


  T5_PANEL_C_ROWKEYS::Vector{Symbol} = [:constant, :plusxmkt, :plussmbhml, :plusrmwcma, :plusumd]

  T5_PANEL_C_ROWLABEL::Dict = Dict(
    :constant=>"Constant",
    :plusxmkt=>"+XMKT",
    :plussmbhml=>"+SMB+HML",
    :plusrmwcma=>"+RMW+CMA",
    :plusumd=>"+UMD")

  T5_PANEL_C_CONTROLS::Dict = Dict(
    :constant=>CONTROL0,
    :plusxmkt=>[:mkt],
    :plussmbhml=>[:mkt, :P_ff3m_smb,	:P_ff3m_hml],
    :plusrmwcma=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma],
    :plusumd=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd])

  T5_PANEL_C_FOCALS::Dict = Dict(
    :constant=>T5_FNAMES,
    :plusxmkt=>T5_FNAMES,
    :plussmbhml=>T5_FNAMES,
    :plusrmwcma=>T5_FNAMES,
    :plusumd=>T5_FNAMES)

  T5_PANEL_C_NAME::String = "t5panelc-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T5_PANEL_C_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpanelc && tablexpanelC(panel,
    fvals=T5_FVALS,
    fgroups=T5_FGROUPS,
    Fnames = T5_FNAMES,
    rowkeys=T5_PANEL_C_ROWKEYS,
    rowlabel=T5_PANEL_C_ROWLABEL,
    controls=T5_PANEL_C_CONTROLS,
    focals=T5_PANEL_C_FOCALS,
    tablename=T5_PANEL_C_NAME,
    minrowspergroup = T5_PANEL_C_MIN_ROWS
  )

  #########################PANEL D#########################
  T5_PANEL_D_CONTROLS::Dict = Dict(
    :quartiles=>CONTROL0,
    :quartiles4ctrl=>T1_CONTROL4,
    :quartiles10ctrl=>T1_CONTROL10,
    :values=>CONTROL0,
    :values4ctrl=>T1_CONTROL4,
    :values10ctrl=>T1_CONTROL10,
    :quartiles_ac=>CONTROL0,
    :quartiles_lret=>T1_CONTROL10)

  T5_PANEL_D_RET::Dict = Dict(
    :quartiles=>:ret,
    :quartiles4ctrl=>:ret,
    :quartiles10ctrl=>:ret,
    :values=>:ret,
    :values4ctrl=>:ret,
    :values10ctrl=>:ret,
    :quartiles_ac=>:ret,
    :quartiles_lret=>:lret)

  T5_PANEL_D_ROWKEYS::Vector{Symbol} = [:quartiles, :quartiles4ctrl, :quartiles10ctrl, :values,
    :values4ctrl, :values10ctrl, :quartiles_ac, :quartiles_lret]

  T5_PANEL_D_ROWLABEL::Dict = Dict(
    :quartiles=>"Leverage Quartiles",
    :quartiles4ctrl=>"... +4 vars",
    :quartiles10ctrl=>"... +10 vars",
    :values=>"Leverage Values",
    :values4ctrl=>"... +4 vars",
    :values10ctrl=>"... +10 vars",
    :quartiles_ac=>"\\midrule Asset-Controlled Quartiles",
    :quartiles_lret=>"Log Returns +10 vars")

  T5_PANEL_D_FOCALS::Dict = Dict(
    :quartiles=>T5_FGROUPS,
    :quartiles4ctrl=>T5_FGROUPS,
    :quartiles10ctrl=>T5_FGROUPS,
    :values=>T5_FVALS,
    :values4ctrl=>T5_FVALS,
    :values10ctrl=>T5_FVALS,
    :quartiles_ac=>T5_FGROUPS_ALT,
    :quartiles_lret=>T5_FGROUPS)

  T5_PANEL_D_NAME::String = "t5paneld-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  #T5_PANEL_D_MIN_ROWS::Int = 24
  T5_PANEL_D_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpaneld && tablexpanelD(panel,
      fvals=T5_FVALS,
      fgroups=T5_FGROUPS,
      ret = T5_PANEL_D_RET,
      rowkeys=T5_PANEL_D_ROWKEYS,
      rowlabel=T5_PANEL_D_ROWLABEL,
      controls=T5_PANEL_D_CONTROLS,
      focals=T5_PANEL_D_FOCALS,
      tablename=T5_PANEL_D_NAME,
      minrowspergroup = T5_PANEL_D_MIN_ROWS
    )

  return nothing
end
