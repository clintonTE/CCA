

function maketable6(panel::DataFrame; trialrunonly::Bool=trialrunonly,
    )::Nothing
  #load the data we will work with
  #panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  sort!(panel, [:permno, :date])

  ################common data


  ###################PANEL A############
  T6_PANEL_A_ROWKEYS::Vector{Symbol} = [:baseline,
    :bknegcash,
    :mknegcash,
    :Dbkflev,
    :Dmkflev,
    :Dbkliab,
    :Dmkliab,
    :Dbknegcash,
    :Dmknegcash]
  T6_PANEL_A_ROWGROUPNAME::Dict = Dict(
    :baseline=>:baseline,
    :bknegcash => :P_bknegcash,
    :mknegcash => :P_mknegcash,
    :Dbkflev => :P_Dbkflev,
    :Dmkflev => :P_Dmkflev,
    :Dbkliab => :P_Dbkliab,
    :Dmkliab => :P_Dmkliab,
    :Dbknegcash => :P_Dbknegcash,
    :Dmknegcash => :P_Dmknegcash,)

  T6_PANEL_A_ROWLABEL::Dict = Dict(
    :baseline=>"none",
    :bknegcash=>"\\midrule bknegcash",
    :mknegcash=>"mknegcash",
    :Dbkflev=>"\\midrule\$\\Delta\$bkflev",
    :Dmkflev=>"\$\\Delta\$mkflev",
    :Dbkliab=>"\$\\Delta\$bkliab",
    :Dmkliab=>"\$\\Delta\$mkliab",
    :Dbknegcash=>"\$\\Delta\$bknegcash",
    :Dmknegcash=>"\$\\Delta\$mknegcash",)

  ALL_ROW_FOCALS::Vector{Symbol} = [:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd]
  T6_PANEL_A_TEST_FOCALS::Vector{Symbol} = [:P_bknegcash, :P_mknegcash, :P_Dmknegcash,
    :P_Dbknegcash, :P_Dmkflev, :P_Dbkflev, :P_Dmkliab, :P_Dbkliab]


  T6_PANEL_A_NAME::String = "t6panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T6_PANEL_A_MIN_ROWS::Int = 24
  T6_PANEL_A_COLKEYS = [:stat, :ret, :mkt, :P_ff5_hml, :P_ff5_smb,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd]
  T6_PANEL_A_COLLABEL = Dict(
    :stat=>"stat",
    :ret=>"see left",
    :mkt=>"xmkt",
    :P_ff5_hml=>"hml",
    :P_ff5_smb=>"smb",
    :P_ff5_rmw=>"rmw",
    :P_ff5_cma=>"cma",
    :P_ff3m_umd=>"umd")

  table6panela(panel,
    rowgroupname = T6_PANEL_A_ROWGROUPNAME,
    rowkeys=T6_PANEL_A_ROWKEYS,
    rowlabel=T6_PANEL_A_ROWLABEL,
    collabel=T6_PANEL_A_COLLABEL,
    testfocals=T6_PANEL_A_TEST_FOCALS,
    tablename=T6_PANEL_A_NAME,
    minrowspergroup = T6_PANEL_A_MIN_ROWS,
    allrowfocals = ALL_ROW_FOCALS,
    colkeys = T6_PANEL_A_COLKEYS
  )


  return nothing
end

function table6panela(panel::DataFrame;
    rowgroupname::Dict= error("portname is required"),
    rowkeys::Vector{Symbol} = error("rowkeys is required"),
    rowlabel::Dict = error("rowlabel is required"),
    colkeys::Vector{Symbol} = error("colkeys is required"),
    collabel::Dict = error("collabel is required"),
    testfocals::Vector{Symbol} = error("testfocals is required"),
    tablename::String = error("table name is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS,
    minrowspergroup::Int = error("minrowspergroup is required"),
    allrowfocals::Vector{Symbol} = error("allrowfocals is required"))

  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)


  #################################
  #build a small  panel of portfolios
  allcols::Vector{Symbol} = [:date; allrowfocals; testfocals]
  minipanel::DataFrame = unique(panel[!, allcols])
  idcols::Vector{Symbol} = setdiff(allcols, testfocals)
  minipanel = stack(minipanel, testfocals, idcols, variable_name = :portname, value_name = :ret)
  sort!(minipanel, [:portname, :date])



  #make the regression specifications
  fs::FactorSpecification = FactorSpecification(minipanel, :ret, :date, :portname,
    Ffactors = allrowfocals,
    Fcharacteristics=nothing,
    Fcontrols = Vector{Symbol}(),
    name=Symbol("t6"),
    skipcrosssectiondfcheck=true,
    minrowspergroup=minrowspergroup)


  frindex::Dict = BJSFF(fs)
  #println(keys(frindex))

  #now build the table
  Ncols::Int = length(colkeys)
  Nrows::Int = length(rowkeys) * 3
  colnames::Vector{Vector{String}} = [(colkey->collabel[colkey]).(colkeys)]
  colidx::Dict = Dict(colkey=>i for (i,colkey) ∈ enumerate(colkeys))

  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)
  descrownames::Vector{String} = Vector{String}(undef, Nrows)

  #now populate the the table
  for (i,rowkey) ∈ enumerate(rowkeys)
    local r::Int = (i-1)*3 + 1 #gets the row index
    descrownames[r] = rowlabel[rowkey]
    descrownames[r+1] = ""
    descrownames[r+2] = ""

    desccontent[r][colidx[:stat]] = "\\text{alpha}"
    desccontent[r+1][colidx[:stat]] = "\\text{T(alpha)}"
    desccontent[r+2][colidx[:stat]] = "\\text{\$R^2\$ \\%}"

    #now get the core content
    for (c,colkey) ∈ enumerate(colkeys)
      #a couple of special cases where the fr object doesn't provide any content
      (colkey == :stat) && continue
      if (rowkey == :baseline) && (colkey == :ret)
        desccontent[r][colidx[colkey]] = ""
        desccontent[r+1][colidx[colkey]] = ""
        desccontent[r+2][colidx[colkey]] = ""
        continue
      end

      frkey::Symbol = Symbol("$(rowgroupname[rowkey])_col$(colkey)")
      fr::FactorRegression = frindex[frkey]
      αidx::Int = fr.coefidx[:intercept]

      #pull the content
      α::Float64 = fr.λᵢ[αidx]
      t::Float64 = fr.λᵢ[αidx] / fr.σᵢ[αidx]
      αR²::Float64 = fr.otherstats[:R²]

      desccontent[r][colidx[colkey]] = brnum(α*100, decimals)
      desccontent[r+1][colidx[colkey]] = brnum(t, decimals)
      desccontent[r+2][colidx[colkey]] = brnum(αR²*100, 0)
    end
  end
  alignmentstring::String = "ll rrrrrrr"

  #display(desccontent)
  #construct the tex table
  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    alignmentstring=alignmentstring)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end
