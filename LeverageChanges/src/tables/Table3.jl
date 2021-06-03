
function maketable3(panel::DataFrame; trialrunonly::Bool=trialrunonly,
    runpanela::Bool=true,
    runpanelb::Bool=true,
    )::Nothing
  #load the data we will work with
  #panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  sort!(panel, [:permno, :date])
  apanel::DataFrame = panel[panel.annflag,:]
  #apanel = apanel[apanel.fyear .≤ 2016,:] #WARNING DELETE THIS!!!!
  apanel.year = (year).(apanel.date)
  ################common data
  cqports::PortfolioConstructions = leverageconstructions()
  T3_FNAMES::Vector{Symbol} = CQUART_FNAMES_D
  T3_FVALS::Vector{Symbol} = CQUART_FVALS_D
  T3_FGROUPS::Vector{Symbol} = CQUART_FGROUPS_D
  T3_FNAMES_ALT::Vector{Symbol} = CQUART_FNAMES_ALT_D
  T3_FGROUPS_ALT::Vector{Symbol} = CQUART_FGROUPS_ALT_D




  #########################PANEL A#########################

  T3_PANEL_A_ROWKEYS::Vector{Symbol} = [:mknegcash, :bknegcash, :mkflev,
    :bkflev, :mkliab, :bkliab]

  T3_PANEL_A_ROWLABEL::Dict = Dict(
    :mknegcash=>"\$\\Delta\$mknegcash",
    :bknegcash=>"\$\\Delta\$bknegcash",
    :mkflev=>"\\\\ \$\\Delta\$mkflev",
    :bkflev=>"\$\\Delta\$bkflev",
    :mkliab=>"\\\\ \$\\Delta\$mkliab",
    :bkliab=>"\$\\Delta\$bkliab")

  T3_PANEL_A_FOCALS::Dict = Dict(
    :mknegcash=>[:Lvol1yearnet, :LLvol1yearnet, :LDmknegcash],
    :bknegcash=>[:Lvol1yearnet, :LLvol1yearnet, :LDbknegcash],
    :mkflev=>[:Lvol1yearnet, :LLvol1yearnet, :LDmkflev],
    :bkflev=>[:Lvol1yearnet, :LLvol1yearnet, :LDbkflev],
    :mkliab=>[:Lvol1yearnet, :LLvol1yearnet, :LDmkliab],
    :bkliab=>[:Lvol1yearnet, :LLvol1yearnet, :LDbkliab])

  T3_PANEL_A_CONTROLS::Dict = Dict(
    :mknegcash=>CONTROL0,
    :bknegcash=>CONTROL0,
    :mkflev=>CONTROL0,
    :bkflev=>CONTROL0,
    :mkliab=>CONTROL0,
    :bkliab=>CONTROL0)

  T3_PANEL_A_RET::Dict = Dict(
    :mknegcash=>:vol1yearnet,
    :bknegcash=>:vol1yearnet,
    :mkflev=>:vol1yearnet,
    :bkflev=>:vol1yearnet,
    :mkliab=>:vol1yearnet,
    :bkliab=>:vol1yearnet)

  T3_PANEL_A_NAME::String = "t3panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T3_PANEL_A_MIN_ROWS = 5

  runpanela && table3panelA(apanel,
      fvals=T3_FVALS,
      fgroups=T3_FGROUPS,
      ret = T3_PANEL_A_RET,
      rowkeys=T3_PANEL_A_ROWKEYS,
      rowlabel=T3_PANEL_A_ROWLABEL,
      controls=T3_PANEL_A_CONTROLS,
      focals=T3_PANEL_A_FOCALS,
      tablename=T3_PANEL_A_NAME,
      minrowspergroup = T3_PANEL_A_MIN_ROWS
    )

    #########################PANEL B#########################

    T3_PANEL_B_ROWKEYS::Vector{Symbol} = [:mknegcash, :bknegcash, :mkflev,
      :bkflev, :mkliab, :bkliab]

    T3_PANEL_B_ROWLABEL::Dict = Dict(
      :mknegcash=>"\$\\Delta\$mknegcash",
      :bknegcash=>"\$\\Delta\$bknegcash",
      :mkflev=>"\\\\ \$\\Delta\$mkflev",
      :bkflev=>"\$\\Delta\$bkflev",
      :mkliab=>"\\\\ \$\\Delta\$mkliab",
      :bkliab=>"\$\\Delta\$bkliab")

    T3_PANEL_B_FOCALS::Dict = Dict(
      :mknegcash=>[:Lvol1yearnet, :LLvol1yearnet, :LDmknegcash],
      :bknegcash=>[:Lvol1yearnet, :LLvol1yearnet, :LDbknegcash],
      :mkflev=>[:Lvol1yearnet, :LLvol1yearnet, :LDmkflev],
      :bkflev=>[:Lvol1yearnet, :LLvol1yearnet, :LDbkflev],
      :mkliab=>[:Lvol1yearnet, :LLvol1yearnet, :LDmkliab],
      :bkliab=>[:Lvol1yearnet, :LLvol1yearnet, :LDbkliab])

    T3_PANEL_B_CONTROLS::Dict = Dict(
      :mknegcash=>CONTROL0,
      :bknegcash=>CONTROL0,
      :mkflev=>CONTROL0,
      :bkflev=>CONTROL0,
      :mkliab=>CONTROL0,
      :bkliab=>CONTROL0)

    T3_PANEL_B_RET::Dict = Dict(
      :mknegcash=>:vol1yearnet,
      :bknegcash=>:vol1yearnet,
      :mkflev=>:vol1yearnet,
      :bkflev=>:vol1yearnet,
      :mkliab=>:vol1yearnet,
      :bkliab=>:vol1yearnet)

    T3_PANEL_B_CLUSTERS::Vector{Symbol} = [:permno, :fyear] #QUESTION: maybe fyrmonth instead?

    T3_PANEL_B_Σ::Vector{Function} = [homoskedasticΣ!,
      lin::FMLM->clusteredΣ!(lin, clusters=[lin.clusters[1]]),
      lin::FMLM->clusteredΣ!(lin, clusters=[lin.clusters[2]]),
      modifiedwhiteΣ!,
      lin::FMLM->clusteredΣ!(lin, clusters=lin.clusters#=, testequivelance=true=#)]

    T3_PANEL_B_NAME::String = "t3panelb-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

    T3_PANEL_B_MIN_ROWS = 5

    runpanelb && table3panelB(apanel,
        fvals=T3_FVALS,
        fgroups=T3_FGROUPS,
        ret = T3_PANEL_B_RET,
        rowkeys=T3_PANEL_B_ROWKEYS,
        rowlabel=T3_PANEL_B_ROWLABEL,
        controls=T3_PANEL_B_CONTROLS,
        focals=T3_PANEL_B_FOCALS,
        tablename=T3_PANEL_B_NAME,
        minrowspergroup = T3_PANEL_B_MIN_ROWS,
        fclusters = T3_PANEL_B_CLUSTERS,
        errorfunctions = T3_PANEL_B_Σ
      )
  return nothing
end


function table3panelA(panel::DataFrame;
    fvals::Vector{Symbol} = error("fvals is required"),
    ret::Dict = error("ret is required"),
    rowkeys::Vector{Symbol} = error("rowkeys is required"),
    rowlabel::Dict = error("rowlabel is required"),
    focals::Dict = error("focals is required"),
    controls::Dict = error("controls is required"),
    tablename::String = error("table name is required"),
    fgroups::Vector{Symbol} = error("fgroups is required"),
    tablepath::String = TABLE_PATH,
    colnames::Vector{Vector{String}} = [["Method",
      "\$\\alpha\$", "\$\\rho_1\$", "\$\\rho_2\$", "\$\\gamma\$",
      "\\% Years Positive"]],
    decimals::Int = DEFAULT_DECIMALS,
    minrowspergroup::Int = error("minrowspergroup is required"))

  Nrows::Int = length(rowkeys)*2
  Ncols::Int = length(colnames[1])
  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)

  #################################
  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(rowkeys) * length(fgroups))

  #make the regression specifications
  for rname ∈ rowkeys
    fs::FactorSpecification = FactorSpecification(
      panel, ret[rname], :fyear, :permno,
      skiptimeseriesdfcheck = true,
      minrowspergroup = minrowspergroup,
      Ffactors=nothing,
      Fcharacteristics = focals[rname],
      Fcontrols = controls[rname],
      name=Symbol("row$(rname)"))
    push!(fss, fs)

    #=if rname == :mkliab
      treg = characteristicFM(fs)

      println("all size: $(size(fs.df,1))")
      for sdf ∈ fs.eachdate

        println(sdf[1,fs.Fdate], ": ", size(sdf,1))
      end
      error("stop")

    end=#
  end

  #Run the regressions
  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]

  frs::Vector{FactorRegression} = Vector{FactorRegression}(undef, length(fss))
  for i ∈ 1:length(fss)
    frs[i] = characteristicFM(fss[i], parallel=parallel)
  end
  frsidx::Dict = Dict(fr.name=>fr for fr ∈ frs)
  #println(collect(keys(frsidx)))
  #make the table
  descrownames::Vector{String} = [(s->[rowlabel[s], ""]).(rowkeys)...;]
  desccontent::Vector{Vector{String}} = (i->
    Vector{String}(undef, Ncols)).(1:Nrows)

  for (i,rname) ∈ enumerate(rowkeys)
    fr::FactorRegression = frsidx[Symbol("row$(rname)")]
    r::Int = 2*(i-1)+1
    desccontent[r][1] = "\\text{Coefficients}"
    desccontent[r+1][1] = "\\text{(NW(2) T-stat)}"

    #coefficients
    desccontent[r][2] = "$(num2str(fr.λᵢ[fr.coefidx[:intercept]]))"
    desccontent[r][3] = "$(num2str(fr.λᵢ[fr.coefidx[focals[rname][1]]]))"
    desccontent[r][4] = "$(num2str(fr.λᵢ[fr.coefidx[focals[rname][2]]]))"
    desccontent[r][5] = "$(num2str(fr.λᵢ[fr.coefidx[focals[rname][3]]]))"

    #SEs
    desccontent[r+1][2] = tstatstr(fr, decimals=decimals, i=fr.coefidx[:intercept])
    desccontent[r+1][3] = tstatstr(fr, decimals=decimals, i=fr.coefidx[focals[rname][1]])
    desccontent[r+1][4] = tstatstr(fr, decimals=decimals, i=fr.coefidx[focals[rname][2]])
    desccontent[r+1][5] = tstatstr(fr, decimals=decimals, i=fr.coefidx[focals[rname][3]])

    perλ::Int = Int(round(sum(fr.λ[:,fr.coefidx[focals[rname][3]]] .≥ 0)/size(fr.λ,1)*100))
    desccontent[r][6] = num2str(perλ, 0)
    desccontent[r+1][6] = ""
  end

  alignmentstring::String = "ll rrrrr"

  tt::String = textable(colnames=colnames,  descrownames=descrownames,
    desccontent=desccontent, alignmentstring=alignmentstring)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function table3panelB(panel::DataFrame;
    fvals::Vector{Symbol} = error("fvals is required"),
    ret::Dict = error("ret is required"),
    rowkeys::Vector{Symbol} = error("rowkeys is required"),
    rowlabel::Dict = error("rowlabel is required"),
    focals::Dict = error("focals is required"),
    controls::Dict = error("controls is required"),
    tablename::String = error("table name is required"),
    fgroups::Vector{Symbol} = error("fgroups is required"),
    tablepath::String = TABLE_PATH,
    colnames::Vector{Vector{String}} = [[
      "\$\\rho_1\$", "\$\\rho_2\$", "\$\\gamma\$",
      "OLS", "CLR-i", "CLR-t", "WHT", "CLR-i\\&t"]],
    decimals::Int = DEFAULT_DECIMALS,
    minrowspergroup::Int = error("minrowspergroup is required"),
    errorfunctions::Vector{Function} = error("error functions are required"),
    fclusters = error("fclusters are required")
    )

  Nrows::Int = length(rowkeys)
  Ncols::Int = length(colnames[1])
  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)

  #################################
  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(rowkeys))

  #make the regression specifications
  for rname ∈ rowkeys
    fs::FactorSpecification = FactorSpecification(
      skiptimeseriesdfcheck = true,
      minrowspergroup = minrowspergroup,
      panel, ret[rname], :fyear, :permno,
      Ffactors=nothing,
      Fcharacteristics = focals[rname],
      Fcontrols = controls[rname],
      name=Symbol("row$(rname)"))
    push!(fss, fs)
  end


  #Run the regressions
  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]

  vfrs::Vector{Vector{FactorRegression}} = []
  sizehint!(vfrs, length(fss))
  for i ∈ 1:length(fss)
    frs::Vector{FactorRegression} = panelFM(fss[i],
      errorfunctions=errorfunctions,
      clustersyms=fclusters)
    push!(vfrs, frs)
  end

  #gets the row
  frsidx::Dict = Dict(rname=>vfrs[j] for (j,rname) ∈ enumerate(rowkeys))
  colidx::Dict = Dict(:ρ₁=>1, :ρ₂=>2, :γ=>3, :ols=>4, :cross=>5, :ts=>6, :wht=>7, :bth=>8)

  #println(collect(keys(frsidx)))
  #make the table
  descrownames::Vector{String} = [(s->rowlabel[s]).(rowkeys);]
  desccontent::Vector{Vector{String}} = (i->
    Vector{String}(undef, Ncols)).(1:Nrows)

  for (r,rname) ∈ enumerate(rowkeys)
    frs::Vector{FactorRegression} = frsidx[rname]

    #coefficients
    desccontent[r][colidx[:ρ₁]] = "$(num2str(frs[1].λᵢ[frs[1].coefidx[focals[rname][1]]]))"
    desccontent[r][colidx[:ρ₂]] = "$(num2str(frs[1].λᵢ[frs[1].coefidx[focals[rname][2]]]))"
    desccontent[r][colidx[:γ]] = "$(num2str(frs[1].λᵢ[frs[1].coefidx[focals[rname][3]]]))"

    #SEs
    desccontent[r][colidx[:ols]] = tstatstr(frs[1], decimals=decimals, i=frs[1].coefidx[focals[rname][3]])
    desccontent[r][colidx[:cross]] = tstatstr(frs[2], decimals=decimals, i=frs[2].coefidx[focals[rname][3]])
    desccontent[r][colidx[:ts]] = tstatstr(frs[3], decimals=decimals, i=frs[3].coefidx[focals[rname][3]])
    desccontent[r][colidx[:wht]] = tstatstr(frs[4], decimals=decimals, i=frs[4].coefidx[focals[rname][3]])
    desccontent[r][colidx[:bth]] = tstatstr(frs[5], decimals=decimals, i=frs[5].coefidx[focals[rname][3]])

  end

  #alignmentstring::String = "l rrrrrrrr"

  tt::String = textable(colnames=colnames,  descrownames=descrownames,
    desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function tstatstr(fr::FactorRegression; decimals::Int=DEFAULT_DECIMALS, i::Int=2,
  )::String
  tstat::Float64 = fr.λᵢ[i]/fr.σᵢ[i]
  return "($(num2str(tstat, decimals)))"
end
