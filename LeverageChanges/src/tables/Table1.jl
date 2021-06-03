

const DISABLE_FM_PARALLEL = true
FM_M = Matrix{Float64}
FM_V = Vector{Float64}
FM_QR = Matrix{Float64}

DEF_MIN_ROWS_T1 = 12 #24

function maketable1(panel::DataFrame; trialrunonly::Bool=trialrunonly,
    runpanela::Bool=true,
    runpanelb::Bool=true,
    runpanelc::Bool=true,
    runpaneld::Bool=true,
    )::Nothing
  #load the data we will work with

  sort!(panel, [:permno, :date])
  apanel::DataFrame = panel[panel.annflag,:]

  ################common data
  cqports::PortfolioConstructions = leverageconstructions()
  T1_FNAMES::Vector{Symbol} = CQUART_FNAMES
  T1_FVALS::Vector{Symbol} = CQUART_FVALS
  T1_FGROUPS::Vector{Symbol} = CQUART_FGROUPS
  T1_FNAMES_ALT::Vector{Symbol} = CQUART_FNAMES_ALT
  T1_FGROUPS_ALT::Vector{Symbol} = CQUART_FGROUPS_ALT

  ################Table A
  T1_PANEL_A_ROWKEYS::Vector{Symbol} = [:mean, :std, :stdxt]

  T1_PANEL_A_ROWLABEL::Dict = Dict(
    :mean=>"Mean",
    :std=>"SD",
    :stdxt=>"SD(X,T)")


  T1_PANEL_A_IDX::Dict = Dict(
    :mean=>1,
    :std=>2,
    :stdxt=>3)
  T1_PANEL_A_NAME::String = "t1panela-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"


  runpanela && tablexpanelA(apanel,
    fvals = T1_FVALS,
    rowkeys=T1_PANEL_A_ROWKEYS,
    rowlabel=T1_PANEL_A_ROWLABEL,
    idx=T1_PANEL_A_IDX,
    tablename=T1_PANEL_A_NAME)


  #####################TABLE B
  T1_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :leverageQmeans,
    :quartilenum,
    :LLsd,
    :sd,
    :Δsd,
    :LLret1year,
    :ret1year,
    :Δret1year]

  T1_PANEL_B_ROWLABEL::Dict = Dict(
    :leverageQmeans=>"Leverage Q Means",
    :quartilenum=>"Quartile \\#",
    :LLsd=>"\\\\ SD Net Lagged\$^2\$",
    :sd=>"SD Net Lead",
    :Δsd=>"SD Net 2Y Delta",
    :LLret1year=>"\\\\ Compound Raw Lagged\$^2\$",
    :ret1year=>"Compound Raw Lead",
    :Δret1year=>"Compound Raw 2Y Delta")


  T1_PANEL_B_IDX::Dict = Dict(
    T1_PANEL_B_ROWKEYS[i]=>i for i ∈ 1:(length(T1_PANEL_B_ROWKEYS)))

  #the below is to help with future abstraction if necessary
  T1_PANEL_B_VALCOLIDX::Dict= Dict(f=>f for f ∈ T1_FVALS)

  #these are pairs of lagged and unlagged fields
  T1_PANEL_B_VOLFIELDS::NTuple{2,Symbol} = (:LLvol1yearnet, :vol1yearnet)
  T1_PANEL_B_RETFIELDS::NTuple{2,Symbol} = (:LLret1year, :ret1year)

  T1_PANEL_B_NAME::String = "t1panelb-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  runpanelb && tablexpanelB(apanel,
    fvals = T1_FVALS,
    rowkeys=T1_PANEL_B_ROWKEYS,
    rowlabel=T1_PANEL_B_ROWLABEL,
    idx=T1_PANEL_B_IDX,
    valcolidx = T1_PANEL_B_VALCOLIDX,
    tablename=T1_PANEL_B_NAME,
    volfields=T1_PANEL_B_VOLFIELDS,
    retfields=T1_PANEL_B_RETFIELDS)

  ###################PANEL C############


  T1_PANEL_C_ROWKEYS::Vector{Symbol} = [:constant, :plusxmkt, :plussmbhml, :plusrmwcma, :plusumd]

  T1_PANEL_C_ROWLABEL::Dict = Dict(
    :constant=>"Constant",
    :plusxmkt=>"+XMKT",
    :plussmbhml=>"+SMB+HML",
    :plusrmwcma=>"+RMW+CMA",
    :plusumd=>"+UMD")

  T1_PANEL_C_CONTROLS::Dict = Dict(
    :constant=>CONTROL0,
    :plusxmkt=>[:mkt],
    :plussmbhml=>[:mkt, :P_ff3m_smb,	:P_ff3m_hml],
    :plusrmwcma=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma],
    :plusumd=>[:mkt, :P_ff5_smb,	:P_ff5_hml,	:P_ff5_rmw,	:P_ff5_cma, :P_ff3m_umd])

  T1_PANEL_C_FOCALS::Dict = Dict(
    :constant=>T1_FNAMES,
    :plusxmkt=>T1_FNAMES,
    :plussmbhml=>T1_FNAMES,
    :plusrmwcma=>T1_FNAMES,
    :plusumd=>T1_FNAMES)

  T1_PANEL_C_NAME::String = "t1panelc-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"

  T1_PANEL_C_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpanelc && tablexpanelC(panel,
    fvals=T1_FVALS,
    fgroups=T1_FGROUPS,
    Fnames = T1_FNAMES,
    rowkeys=T1_PANEL_C_ROWKEYS,
    rowlabel=T1_PANEL_C_ROWLABEL,
    controls=T1_PANEL_C_CONTROLS,
    focals=T1_PANEL_C_FOCALS,
    tablename=T1_PANEL_C_NAME,
    minrowspergroup = T1_PANEL_C_MIN_ROWS
  )

  #########################PANEL D#########################
  T1_PANEL_D_CONTROLS::Dict = Dict(
    :quartiles=>CONTROL0,
    :quartiles4ctrl=>T1_CONTROL4,
    :quartiles10ctrl=>T1_CONTROL10,
    :values=>CONTROL0,
    :values4ctrl=>T1_CONTROL4,
    :values10ctrl=>T1_CONTROL10,
    :quartiles_ac=>CONTROL0,
    :quartiles_lret=>T1_CONTROL10)

  T1_PANEL_D_RET::Dict = Dict(
    :quartiles=>:ret,
    :quartiles4ctrl=>:ret,
    :quartiles10ctrl=>:ret,
    :values=>:ret,
    :values4ctrl=>:ret,
    :values10ctrl=>:ret,
    :quartiles_ac=>:ret,
    :quartiles_lret=>:lret)

  T1_PANEL_D_ROWKEYS::Vector{Symbol} = [:quartiles, :quartiles4ctrl, :quartiles10ctrl, :values,
    :values4ctrl, :values10ctrl, :quartiles_ac, :quartiles_lret]

  T1_PANEL_D_ROWLABEL::Dict = Dict(
    :quartiles=>"Leverage Quartiles",
    :quartiles4ctrl=>"... +4 vars",
    :quartiles10ctrl=>"... +10 vars",
    :values=>"Leverage Values",
    :values4ctrl=>"... +4 vars",
    :values10ctrl=>"... +10 vars",
    :quartiles_ac=>"\\midrule Asset-Controlled Quartiles",
    :quartiles_lret=>"Log Returns +10 vars")

  T1_PANEL_D_FOCALS::Dict = Dict(
    :quartiles=>T1_FGROUPS,
    :quartiles4ctrl=>T1_FGROUPS,
    :quartiles10ctrl=>T1_FGROUPS,
    :values=>T1_FVALS,
    :values4ctrl=>T1_FVALS,
    :values10ctrl=>T1_FVALS,
    :quartiles_ac=>T1_FGROUPS_ALT,
    :quartiles_lret=>T1_FGROUPS)

  T1_PANEL_D_NAME::String = "t1paneld-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[])"
  T1_PANEL_D_MIN_ROWS::Int = DEF_MIN_ROWS_T1

  runpaneld && tablexpanelD(panel,
      fvals=T1_FVALS,
      fgroups=T1_FGROUPS,
      ret = T1_PANEL_D_RET,
      rowkeys=T1_PANEL_D_ROWKEYS,
      rowlabel=T1_PANEL_D_ROWLABEL,
      controls=T1_PANEL_D_CONTROLS,
      focals=T1_PANEL_D_FOCALS,
      tablename=T1_PANEL_D_NAME,
      minrowspergroup = T1_PANEL_D_MIN_ROWS
    )

  return nothing
end

function tablexpanelA(panel::DataFrame;
    fvals::Vector{Symbol} = error("PANELA: fvals is required"),
    rowkeys::Vector{Symbol} = error("PANELA: rowkeys is required"),
    rowlabel::Dict = error("PANELA: rowlabel is required"),
    idx::Dict = error("PANELA: idx is required"),
    tablename::String = error("PANELA: table name is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS)
  #local cqports::CQuartileConstructions = leverageconstructions()


  #NOTE: These are of the form mknegcash ~ α + β*mkt ...



  #################################

  #make the table
  colnames::Vector{Vector{String}} = [(string).(fvals)]
  descrownames::Vector{String} = (s->rowlabel[s]).(rowkeys)

  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, length(fvals))).(1:length(rowkeys))
  for c ∈ 1:length(fvals)
    local sspanels::GroupedDataFrame
    local ngroups::Int

    f::Symbol = fvals[c]
    allfields::Vector{Symbol} = [:date, :permno, f]
    spanel::SubDataFrame = view(panel, completecases(panel[!,allfields]), allfields)

    cmean::Float64 = mean(spanel[!, f])
    cstd::Float64 = std(spanel[!, f])

    cstdx::Float64 = 0.0
    ngroups = 0
    sspanels = groupby(spanel, :date)
    for sspanel ∈ sspanels
      cstdxi::Float64 = std(sspanel[!,f])
      if isfinite(cstdxi)
        cstdx += cstdxi
        ngroups+=1
      end
    end
    cstdx /= ngroups

    cstdt::Float64 = 0.0
    ngroups = 0
    sspanels = groupby(spanel, :permno)
    for sspanel ∈ groupby(spanel, :permno)
      cstdti::Float64 = std(sspanel[!,f])
      if isfinite(cstdti)
        cstdt += cstdti
        ngroups+=1
      end
    end
    cstdt /= ngroups

    desccontent[idx[:mean]][c] = brnum(cmean, decimals)
    desccontent[idx[:std]][c] = "($(num2str(cstd, decimals)))"
    desccontent[idx[:stdxt]][c] = "($(num2str(cstdx, decimals)), $(num2str(cstdt, decimals)))"
  end

  tt::String = textable(colnames=colnames, descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function tablexpanelB(panel::DataFrame;
    fvals::Vector{Symbol} = error("PANELB: fvals is required"),
    #fgroups::Vector{Symbol} = error("PANELB: fgroups is required"),
    rowkeys::Vector{Symbol} = error("PANELB: rowkeys is required"),
    rowlabel::Dict = error("PANELB: rowlabel is required"),
    idx::Dict = error("PANELB: idx is required"),
    valcolidx::Dict = error("PANELB: valcolidx is required"),
    tablename::String = error("PANELB: table name is required"),
    volfields::NTuple{2, Symbol}=error("PANELB: volfields is required"),
    retfields::NTuple{2, Symbol}=error("PANELB: retfields is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = 1)
  local cqports::CQuartileConstructions = leverageconstructions()

  #println(describe(panel))

  #################################

  #will be using this a lot
  @inline skipmean(v::AbstractVector{<:MFloat64})::Float64 = mean(skipmissing(v))
  @inline n2s(x::T where T<:MFloat64, d::Real=decimals) = num2str(x, d, scalefactor=100.)
  @inline br(x::T where T<:MFloat64, d::Real=decimals) = brnum(x, d, scalefactor=100.)

  #=println("Test neg cash: (")
  sdf::SubDataFrame = view(panel, (g->(!ismissing(g)) && g==1).(panel[!, cqports[1].Fgroup]),:)
  print("$(skipmean(sdf[!, T1_PANEL_B_RETFIELDS[1]])), ")
  sdf = view(panel, (g->(!ismissing(g)) && g==4).(panel[!, cqports[1].Fgroup]),:)
  print("$(skipmean(sdf[!, T1_PANEL_B_RETFIELDS[1]])))")=#

  #make the table
  colnames::Vector{Vector{String}} = [[(s->string(s)).(fvals)...;]]
  widthcolnames::Vector{Vector{Int}} = [[(s->2).(fvals)...;]]
  descrownames::Vector{String} = (s->rowlabel[s]).(rowkeys)
  alignmentstring::String = string(" l", join(["r@{ }r" for i ∈ 1:length(fvals)]))

  #pre-allocate for content
  numcols::Int = length(fvals)*2
  desccontent::Vector{Vector{String}} = (
    i->Vector{String}(undef, numcols)).(1:length(rowkeys))



  volscale::Float64 = REPLICATION_TYPE[] == :daily ? 252^0.5 : 12^0.5

  for i ∈ 1:length(fvals)
    c::Int = 1+(i-1)*2

    local sspanelq1::SubDataFrame
    local sspanelq4::SubDataFrame
    local allfields::Vector{Symbol}

    #acquire the fields we will use
    f::Symbol = fvals[i]
    name::Symbol = Symbol("P_", f)
    Fgroup::Symbol = cqports[name].Fgroup
    #println("f: $f, name: $name, Fgroup: $Fgroup")
    basicfields::Vector{Symbol} = [:date; :permno; Fgroup]
    local spanel::SubDataFrame = view(panel, completecases(panel[!,basicfields]),:)

    ###Value Rows
    Fval::Symbol = valcolidx[f]
    #allfields = [basicfields; Fval]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1val::Float64 = mean(skipmissing(sspanelq1[!, Fval]))
    q4val::Float64 = mean(skipmissing(sspanelq4[!, Fval]))
    desccontent[idx[:leverageQmeans]][c] = (
      "$(brnum(q1val, 2))")
    desccontent[idx[:leverageQmeans]][c+1] = ("$(brnum(q4val, 2))")

    desccontent[idx[:quartilenum]][c] = ("\\text{Q1}")
    desccontent[idx[:quartilenum]][c+1] = ("\\text{Q4}")

    ###SD Rows
    #allfields = [basicfields; T1_PANEL_B_VOLFIELDS]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1LLvol::Float64 = skipmean(sspanelq1[!, volfields[1]]) * volscale
    q4LLvol::Float64 = skipmean(sspanelq4[!, volfields[1]]) * volscale
    desccontent[idx[:LLsd]][c] = ("$(n2s(q1LLvol))")
    desccontent[idx[:LLsd]][c+1] = ("$(n2s(q4LLvol))")

    q1vol::Float64 = skipmean(sspanelq1[!, volfields[2]]) * volscale
    q4vol::Float64 = skipmean(sspanelq4[!, volfields[2]]) * volscale
    desccontent[idx[:sd]][c] = ("$(n2s(q1vol))")
    desccontent[idx[:sd]][c+1] = ("$(n2s(q4vol))")

    Δq1vol::Float64 = skipmean(sspanelq1[!, volfields[2]] .-
      sspanelq1[!, volfields[1]]) * volscale
    Δq4vol::Float64 = skipmean(sspanelq4[!, volfields[2]] .-
      sspanelq4[!, volfields[1]]) * volscale
    desccontent[idx[:Δsd]][c] = ("$(br(Δq1vol))")
    desccontent[idx[:Δsd]][c+1] = ("$(br(Δq4vol))")

    ###Raw Ret Rows
    #allfields = [basicfields; T1_PANEL_B_RETFIELDS]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1LLret::Float64 = skipmean(sspanelq1[!, retfields[1]])
    q4LLret::Float64 = skipmean(sspanelq4[!, retfields[1]])
    desccontent[idx[:LLret1year]][c] = ("$(n2s(q1LLret))")
    desccontent[idx[:LLret1year]][c+1] = ("$(n2s(q4LLret))")

    q1ret::Float64 = skipmean(sspanelq1[!, retfields[2]])
    q4ret::Float64 = skipmean(sspanelq4[!, retfields[2]])
    desccontent[idx[:ret1year]][c] = ("$(n2s(q1ret))")
    desccontent[idx[:ret1year]][c+1] = ("$(n2s(q4ret))")

    Δq1ret::Float64 = skipmean(sspanelq1[!, retfields[2]] .-
      sspanelq1[!, retfields[1]])
    Δq4ret::Float64 = skipmean(sspanelq4[!, retfields[2]] .-
      sspanelq4[!, retfields[1]])
    desccontent[idx[:Δret1year]][c] = ("$(br(Δq1ret))")
    desccontent[idx[:Δret1year]][c+1] = ("$(br(Δq4ret))")


  end

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    widthcolnames=widthcolnames,
    alignmentstring=alignmentstring)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function tablexpanelC(panel::DataFrame;
    fvals::Vector{Symbol} = error("PANELC: fvals is required"),
    Fnames::Vector{Symbol} = error("PANELC: Fnames is required"),
    rowkeys::Vector{Symbol} = error("PANELC: rowkeys is required"),
    rowlabel::Dict = error("PANELC: rowlabel is required"),
    focals::Dict = error("PANELC: focals is required"),
    controls::Dict = error("PANELC: controls is required"),
    tablename::String = error("PANELC: table name is required"),
    fgroups::Vector{Symbol} = error("PANELC: fgroups is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS,
    minrowspergroup::Int = error("minrowspergroup is required"))

  #local cqports::CQuartileConstructions = leverageconstructions()

  #println(Fnames)

  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)


  #################################
  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(rowkeys) * length(fgroups))

  #need to re-shape for this exercise
  idcols::Vector{Symbol} = unique([:date; values(controls)...;])
  measurecols::Vector{Symbol} = Fnames
  allcols::Vector{Symbol} = unique([:date; values(controls)...; values(focals)...;])
  minipanel::DataFrame = unique(panel[!,allcols])
  minipanel = stack(minipanel, measurecols, idcols, variable_name = :portname, value_name = :ret)

  #@inline parseportname(s::String) = parse(Int, s)
  #minipanel.portname = (parseportname).(Vector{String}(minipanel.portname)) #keep this as a int for consistency
  minipanel.portname = (Symbol).(Vector{String}(minipanel.portname))
  sort!(minipanel, [:portname, :date])

  #println("allcols: $allcols")
  #println("minipanelun: ", unique(minipanel[!, :portname]))
  #println(describe(minipanel))
  #println("Size minipanel: $(size(minipanel))")

  minipanel |> CSV.write("output\\minipanelc-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[]).csv")

  #make the regression specifications
  for rname ∈ rowkeys
    sallcols::Vector{Symbol} = [:ret; :date; :portname; controls[rname]]
    sminipanel = view(minipanel, completecases(minipanel[!,sallcols]), sallcols)
    #println("sminipanelun: ", describe(sminipanel))
    fs = FactorSpecification(sminipanel, :ret, :date, :portname,
      Ffactors = nothing,
      Fcharacteristics=nothing,
      Fcontrols = controls[rname],
      name=Symbol("row$(rname)"),
      skipcrosssectiondfcheck=true,
      skiptimeseriesdfcheck=false,
      minrowspergroup=minrowspergroup,
      )
    #lock(spinlock)
    push!(fss, fs) #only doing time series regressions
    #unlock(spinlock)
  end

  #run the regressions
  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]
  frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #generate and push all the regressions
  #each call to manyFF pulls a row of regressions
  for i ∈ 1:length(fss)
    sfrs::Vector{FactorRegression} = manyFF(fss[i], FM_M, FM_V,
      qrtype=FM_QR, parallel=parallel)

    push!(frs, sfrs...)
  end

  #println((fr->fr.name).(frs))

  #make the table
  colnames::Vector{Vector{String}} = [
    (string).(fvals), (i->"Alpha (T-stat)").(1:length(fvals))]
  descrownames::Vector{String} = (s->rowlabel[s]).(rowkeys)


  desccontent::Vector{Vector{String}} = factordesccontent(
    fr::FactorRegression->coeftstat(fr, i=1), #fss=fss,
    frs=frs, rowkeys = rowkeys, colindex=focals)

  tt::String = textable(colnames=colnames, descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function tablexpanelD(panel::DataFrame;
    fvals::Vector{Symbol} = error("fvals is required"),
    ret::Dict = error("ret is required"),
    rowkeys::Vector{Symbol} = error("rowkeys is required"),
    rowlabel::Dict = error("rowlabel is required"),
    focals::Dict = error("focals is required"),
    controls::Dict = error("controls is required"),
    tablename::String = error("table name is required"),
    fgroups::Vector{Symbol} = error("fgroups is required"),
    tablepath::String = TABLE_PATH,
    minrowspergroup::Int = error("minrowspergroup is required"))


  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)
  #@warn "TEST RUN ONLY- fix this before production!!!"
  #panel = panel[(d->year(d)>2010).(panel.date),:]

  #################################
  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(rowkeys) * length(fgroups))

  #make the regression specifications
  for rname ∈ rowkeys
    for focal ∈ focals[rname]
      fs::FactorSpecification = FactorSpecification(panel, ret[rname], :date, :permno,
        Ffactors=nothing,
        Fcharacteristics = [focal],
        Fcontrols = controls[rname],
        minrowspergroup=minrowspergroup,
        skiptimeseriesdfcheck=true,
        name=Symbol("row$(rname)_col$(focal)"))


      #=if occursin("mkliab", string(focal))
        fsreg = characteristicFM(fs, FM_M, FM_V,
          parallel=parallel, qrtype=FM_QR)
        print(size(fsreg.λ), ": ")
        #println(fsreg.λ)
        error("stop")

      end=#

      #lock(spinlock)
      push!(fss, fs) #only doing crossectional regressions
      #unlock(spinlock)
    end
  end

  #Run the regressions

  #for testing purposes
  if length(unique((d->year(d)).(panel.date))) < 10
    csvcols::Vector{Symbol} = [:date; :permno; ret[rowkeys[1]];
      focals[rowkeys[1]][1]; controls[rowkeys[1]]]
    panel[completecases(panel[!,csvcols]), csvcols] |> CSV.write(
      "output\\paneld-$(REPLICATION_TYPE[])-$(OUT_SUFFIX[]).csv")
    @info "test file written - warning- may need to remove restriction"
  end

  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]
  frs::Vector{FactorRegression} = Vector{FactorRegression}(undef, length(fss))
  for i ∈ 1:length(fss)
    frs[i] = characteristicFM(fss[i], FM_M, FM_V,
      parallel=parallel, qrtype=FM_QR)
  end

  #make the table
  colnames::Vector{Vector{String}} = [
    (string).(fvals), (i->"Coef (T-stat)").(1:length(fvals))]
  descrownames::Vector{String} = (s->rowlabel[s]).(rowkeys)
  desccontent::Vector{Vector{String}} = factordesccontent(#fss=fss,
    frs=frs,  rowkeys = rowkeys, colindex=focals)

  tt::String = textable(colnames=colnames,
    descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function coeftstat(fr::FactorRegression; decimals::Int=DEFAULT_DECIMALS, i::Int=2,
  )::String
  λᵢ::Float64 = fr.λᵢ[i] * 100
  tstat::Float64 = fr.λᵢ[i]/fr.σᵢ[i]

  outstring = "$(num2str(λᵢ, decimals))\\;($(num2str(tstat, decimals)))"

  return brnum(λᵢ, numstring=outstring)
end

#this is a kernel function to generate the content matrix under a number of different scenarios
function factordesccontent(f::Function = coeftstat;
    #fss::Vector{FactorSpecification} = error("fss is required for factordesccontent"),
    frs::Vector{FactorRegression} = error("frs is required for factordesccontent"),
    rowkeys::Vector{Symbol} = error("Rowkeys is required"),
    colindex::Dict = error("colindex is required for factordesccontent"))

    Nrows::Int = length(rowkeys)
    Ncols::Int = length(collect(values(colindex))[1])

    #pre-allocate
    desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)

    #convenience dictionaries
    frsidx = Dict(fr.name=>fr for fr ∈ frs)

    rowcolnameidx = Dict(
      (rname::Symbol, focal::Symbol)=>Symbol("row$(rname)_col$(focal)")
      for rname::Symbol ∈ rowkeys for focal::Symbol ∈ colindex[rname])

    #write the content matrix
    for (r,rname) ∈ enumerate(rowkeys)
      for (c, focal) ∈ enumerate(colindex[rname])
        fr::FactorRegression = frsidx[rowcolnameidx[(rname, focal)]]
        desccontent[r][c] = f(fr)
      end
    end

    return desccontent
  end

#creates the blue/red numbers used in much fo the tables
@inline function brnum(f::Real,
    decimals::Int = DEFAULT_DECIMALS;
    Ints::Bool = false,
    scalefactor=1.,
    numstring::MString = num2str(f, decimals, Ints=Ints, scalefactor=scalefactor))

  local outstring::String

  if f≥0
    outstring = "\\textcolor{blue}{$numstring}"
  else
    outstring = "\\textcolor{red}{$numstring}"
  end

  return outstring
end

brnum(::Missing, args...; keyargs...)::MString = missing
