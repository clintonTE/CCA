
function maketable4(panel::DataFrame; trialrunonly::Bool=trialrunonly,
    runpanela::Bool=true,
    )::Nothing
  #load the data we will work with
  #panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  sort!(panel, [:permno, :date])

  ################common data
  cqports::PortfolioConstructions = leverageconstructions()
  T4_FNAMES::Vector{Symbol} = CQUART_FNAMES_FD
  T4_FVALS::Vector{Symbol} = CQUART_FVALS_D
  T4_FGROUPS::Vector{Symbol} = CQUART_FGROUPS_FD




  #########################PANEL A#########################
  T4_PANEL_A_ROWKEYS::Vector{Symbol} = [:constant, :plusxmkt, :plussmbhml, :plusrmwcma, :plusumd]

  T4_PANEL_A_ROWLABEL::Dict = Dict(
    :constant=>"Constant",
    :plusxmkt=>"+XMKT",
    :plussmbhml=>"+SMB+HML",
    :plusrmwcma=>"+RMW+CMA",
    :plusumd=>"+UMD")

  T4_PANEL_A_CONTROLS::Dict = Dict(
    :constant=>CONTROL0,
    :plusxmkt=>[:mkt],
    :plussmbhml=>[:mkt, :P_ff3m_smb,	:P_ff3m_hml],
    :plusrmwcma=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma],
    :plusumd=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd])

  T4_PANEL_A_FOCALS::Dict = Dict(
    :constant=>T4_FNAMES,
    :plusxmkt=>T4_FNAMES,
    :plussmbhml=>T4_FNAMES,
    :plusrmwcma=>T4_FNAMES,
    :plusumd=>T4_FNAMES)

  T4_PANEL_A_NAME::String = "t4panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T4_PANEL_A_MIN_ROWS::Int = 24

  runpanela && tablexpanelC(panel,
    fvals=T4_FVALS,
    fgroups=T4_FGROUPS,
    Fnames = T4_FNAMES,
    rowkeys=T4_PANEL_A_ROWKEYS,
    rowlabel=T4_PANEL_A_ROWLABEL,
    controls=T4_PANEL_A_CONTROLS,
    focals=T4_PANEL_A_FOCALS,
    tablename=T4_PANEL_A_NAME,
    minrowspergroup = T4_PANEL_A_MIN_ROWS
  )
  return nothing
end
