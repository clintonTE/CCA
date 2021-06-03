


function maketablea1(panel::DataFrame; trialrunonly::Bool=trialrunonly,
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
  TA1_FNAMES::Vector{Symbol} = CQUART_FNAMES_D
  TA1_FVALS::Vector{Symbol} = CQUART_FVALS_D
  TA1_FGROUPS::Vector{Symbol} = CQUART_FGROUPS_D
  TA1_FNAMES_ALT::Vector{Symbol} = CQUART_FNAMES_ALT_D
  TA1_FGROUPS_ALT::Vector{Symbol} = CQUART_FGROUPS_ALT_D



  ###################PANEL C############
  TA1_PANEL_A_ROWKEYS::Vector{Symbol} = [
    :constantbrsd, :plusxmktbrsd, :plusumdbrsd,
    :constantbrse, :plusxmktbrse, :plusumdbrse,
    :constantbrsx, :plusxmktbrsx, :plusumdbrsx,
  ]

  TA1_PANEL_A_ROWLABEL::Dict = Dict(
    :constantbrsd=>"BRS \$\\Delta\$Debt",
    :plusxmktbrsd=>"BRS+XMKT",
    :plusumdbrsd=>"BRS+FF5+UMD",
    :constantbrse=>"\\midrule BRS \$\\Delta\$Equity",
    :plusxmktbrse=>"BRS+XMKT",
    :plusumdbrse=>"BRS+FF5+UMD",
    :constantbrsx=>"\\midrule BRS \$\\Delta\$External",
    :plusxmktbrsx=>"BRS+XMKT",
    :plusumdbrsx=>"BRS+FF5+UMD",
    )

  TA1_PANEL_A_CONTROLS::Dict = Dict(
    :constantbrsd=>[:P_brsd],
    :plusxmktbrsd=>[:P_brsd, :mkt],
    :plusumdbrsd=>[:P_brsd, :mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd],
    :constantbrse=>[:P_brse],
    :plusxmktbrse=>[:P_brse, :mkt],
    :plusumdbrse=>[:P_brse, :mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd],
    :constantbrsx=>[:P_brsx],
    :plusxmktbrsx=>[:P_brsx, :mkt],
    :plusumdbrsx=>[:P_brsx, :mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd],
    )

  TA1_PANEL_A_FOCALS::Dict = Dict(
  :constantbrsd=>TA1_FNAMES,
  :plusxmktbrsd=>TA1_FNAMES,
  :plusumdbrsd=>TA1_FNAMES,
  :constantbrse=>TA1_FNAMES,
  :plusxmktbrse=>TA1_FNAMES,
  :plusumdbrse=>TA1_FNAMES,
  :constantbrsx=>TA1_FNAMES,
  :plusxmktbrsx=>TA1_FNAMES,
  :plusumdbrsx=>TA1_FNAMES)

  TA1_PANEL_A_NAME::String = "ta1panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  TA1_PANEL_A_MIN_ROWS::Int = 24

  runpanela && tablexpanelC(panel,
    fvals=TA1_FVALS,
    fgroups=TA1_FGROUPS,
    Fnames = TA1_FNAMES,
    rowkeys=TA1_PANEL_A_ROWKEYS,
    rowlabel=TA1_PANEL_A_ROWLABEL,
    controls=TA1_PANEL_A_CONTROLS,
    focals=TA1_PANEL_A_FOCALS,
    tablename=TA1_PANEL_A_NAME,
    minrowspergroup = TA1_PANEL_A_MIN_ROWS
  )

  #########################PANEL B#########################
  TA1_PANEL_B_CONTROLS::Dict = Dict(
    :brsdquartiles=>[:brsd],
    :brsdquartiles10ctrl=>[T1_CONTROL10; :brsd],
    :brsequartiles=>[:brse],
    :brsequartiles10ctrl=>[T1_CONTROL10; :brse],
    :brsxquartiles=>[:brsx],
    :brsxquartiles10ctrl=>[T1_CONTROL10; :brsx],
    :brsdvalues=>[:brsd],
    :brsdvalues10ctrl=>[T1_CONTROL10; :brsd],
    :brsevalues=>[:brse],
    :brsevalues10ctrl=>[T1_CONTROL10; :brse],
    :brsxvalues=>[:brsx],
    :brsxvalues10ctrl=>[T1_CONTROL10; :brsx],)

  TA1_PANEL_B_RET::Dict = Dict(
    :brsdquartiles=>:ret,
    :brsdquartiles10ctrl=>:ret,
    :brsequartiles=>:ret,
    :brsequartiles10ctrl=>:ret,
    :brsxquartiles=>:ret,
    :brsxquartiles10ctrl=>:ret,
    :brsdvalues=>:ret,
    :brsdvalues10ctrl=>:ret,
    :brsevalues=>:ret,
    :brsevalues10ctrl=>:ret,
    :brsxvalues=>:ret,
    :brsxvalues10ctrl=>:ret)

  TA1_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :brsdquartiles,
    :brsdquartiles10ctrl,
    :brsequartiles,
    :brsequartiles10ctrl,
    :brsxquartiles,
    :brsxquartiles10ctrl,
    :brsdvalues,
    :brsdvalues10ctrl,
    :brsevalues,
    :brsevalues10ctrl,
    :brsxvalues,
    :brsxvalues10ctrl,]

  TA1_PANEL_B_ROWLABEL::Dict = Dict(
    :brsdquartiles=>"Leverage Quartiles+\$\\Delta\$Debt",
    :brsdquartiles10ctrl=>"... +10 vars",
    :brsequartiles=>"Leverage Quartiles+\$\\Delta\$Equity",
    :brsequartiles10ctrl=>"... +10 vars",
    :brsxquartiles=>"Leverage Quartiles+\$\\Delta\$External",
    :brsxquartiles10ctrl=>"... +10 vars",
    :brsdvalues=>"\\midrule Leverage Values+\$\\Delta\$Debt",
    :brsdvalues10ctrl=>"... +10 vars",
    :brsevalues=>"Leverage Values+\$\\Delta\$Equity",
    :brsevalues10ctrl=>"... +10 vars",
    :brsxvalues=>"Leverage Values+\$\\Delta\$External",
    :brsxvalues10ctrl=>"... +10 vars")

  TA1_PANEL_B_FOCALS::Dict = Dict(
    :brsdquartiles=>TA1_FGROUPS,
    :brsdquartiles10ctrl=>TA1_FGROUPS,
    :brsequartiles=>TA1_FGROUPS,
    :brsequartiles10ctrl=>TA1_FGROUPS,
    :brsxquartiles=>TA1_FGROUPS,
    :brsxquartiles10ctrl=>TA1_FGROUPS,
    :brsdvalues=>TA1_FVALS,
    :brsdvalues10ctrl=>TA1_FVALS,
    :brsevalues=>TA1_FVALS,
    :brsevalues10ctrl=>TA1_FVALS,
    :brsxvalues=>TA1_FVALS,
    :brsxvalues10ctrl=>TA1_FVALS,  )

  TA1_PANEL_B_NAME::String = "ta1panelb-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  TA1_PANEL_B_MIN_ROWS::Int = 24

  runpanelb && tablexpanelD(panel,
      fvals=TA1_FVALS,
      fgroups=TA1_FGROUPS,
      ret = TA1_PANEL_B_RET,
      rowkeys=TA1_PANEL_B_ROWKEYS,
      rowlabel=TA1_PANEL_B_ROWLABEL,
      controls=TA1_PANEL_B_CONTROLS,
      focals=TA1_PANEL_B_FOCALS,
      tablename=TA1_PANEL_B_NAME,
      minrowspergroup = TA1_PANEL_B_MIN_ROWS
    )

  return nothing
end
