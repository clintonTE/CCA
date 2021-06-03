

const CONTROL0 = Vector{Symbol}()
const T1_CONTROL4 = [:Lbm; :Lmkequity; :LDat; :Lop]
const T1_CONTROL10 = [T1_CONTROL4; :L12ret12m; :L24ret12m;
  :fyret12m; :L12vol252d; :L24vol252d; :fyvol252d]


function maketable1(;trialrunonly::Bool=trialrunonly)::Nothing
  #load the data we will work with
  panel::DataFrame = constructportfolios(refreshportfolios=false, trialrunonly=trialrunonly)
  #table1panelA(panel)
  table1panelB(panel)
  #table1panelC(panel)
  #table1panelD(panel)

  return nothing
end


function table1panelA(panel::DataFrame;
    tablepath::String = TABLE_PATH,
    tablename::String = "t1panela-$(REPLICATION_TYPE[])",
    decimals::Int = DEFAULT_DECIMALS)
  local cqports::CQuartileConstructions = leverageconstructions()
  Fgrps::Vector{Symbol} = Fgroups(cqports)

  #NOTE: These are of the form mknegcash ~ α + β*mkt ...

  ################INPUTS############
  T1_PANEL_A_ROWKEYS::Vector{Symbol} = [:mean, :std, :stdxt]

  T1_PANEL_A_ROWLABEL::Dict = Dict(
    :mean=>"Mean",
    :std=>"SD",
    :stdxt=>"SD(X,T)")


  T1_PANEL_A_IDX::Dict = Dict(
    :mean=>1,
    :std=>2,
    :stdxt=>3)

  #################################

  #make the table
  colnames::Vector{Vector{String}} = [(string).(CQUART_FVALS)]
  descrownames::Vector{String} = (s->T1_PANEL_A_ROWLABEL[s]).(T1_PANEL_A_ROWKEYS)

  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, length(CQUART_FVALS))).(1:length(T1_PANEL_A_ROWKEYS))
  for c ∈ 1:length(CQUART_FVALS)
    local sspanels::GroupedDataFrame
    local ngroups::Int

    f::Symbol = CQUART_FVALS[c]
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

    desccontent[T1_PANEL_A_IDX[:mean]][c] = num2str(cmean, decimals)
    desccontent[T1_PANEL_A_IDX[:std]][c] = "($(num2str(cstd, decimals)))"
    desccontent[T1_PANEL_A_IDX[:stdxt]][c] = "($(num2str(cstdx, decimals)), $(num2str(cstdt, decimals)))"
  end

  tt::String = textable(colnames=colnames, descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function table1panelB(panel::DataFrame;
    tablepath::String = TABLE_PATH,
    tablename::String = "t1panelb-$(REPLICATION_TYPE[])",
    decimals::Int = 1)
  local cqports::CQuartileConstructions = leverageconstructions()

  ################INPUTS############
  T1_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :leverageQmeans,
    :quartilenum,
    :LLsd,
    :sd,
    :Δsd,
    :LLret12m,
    :ret12m,
    :Δret12m]

  T1_PANEL_B_ROWLABEL::Dict = Dict(
    :leverageQmeans=>"Leverage Q Means",
    :quartilenum=>"Quartile \\#",
    :LLsd=>"\\smallskip SD Net Lagged\$^2\$",
    :sd=>"SD Net Lead",
    :Δsd=>"SD Net 2Y Lead",
    :LLret12m=>"\\smallskip Compound Raw Lagged\$^2\$",
    :ret12m=>"Compound Raw Lead",
    :Δret12m=>"Compound Raw 2Y Delta")


  T1_PANEL_B_IDX::Dict = Dict(T1_PANEL_B_ROWKEYS[i]=>i for i ∈ 1:(length(T1_PANEL_B_ROWKEYS)))

  #the below is to help with future abstraction if necessary
  T1_PANEL_B_VALCOLIDX::Dict= Dict(f=>f for f ∈ CQUART_FVALS)

  #these are pairs of lagged and unlagged fields
  T1_PANEL_B_VOLFIELDS::NTuple{2,Symbol} = (:LLvol252d, :vol252d)
  T1_PANEL_B_RETFIELDS::NTuple{2,Symbol} = (:LLret12m, :ret12m)

  #################################

  #make the table
  colnames::Vector{Vector{String}} = [[(s->[string(s), ""]).(CQUART_FVALS)...;]]
  pop!(colnames[1])
  widthcolnames::Vector{Vector{Int}} = [[(s->[2,1]).(CQUART_FVALS)...;]]
  pop!(widthcolnames[1])
  descrownames::Vector{String} = (s->T1_PANEL_B_ROWLABEL[s]).(T1_PANEL_B_ROWKEYS)

  #pre-allocate for content
  numcols::Int = length(CQUART_FVALS)*3-1
  desccontent::Vector{Vector{String}} = (
    i->Vector{String}(undef, numcols)).(1:length(T1_PANEL_B_ROWKEYS))

  #will be using this a lot
  @inline skipmean(v::AbstractVector{T where T<:MFloat64})::Float64 = mean(skipmissing(v))
  @inline n2s(x::T where T<:MFloat64, d::Real=decimals) = num2str(x, d, scalefactor=100.)

  volscale::Float64 = REPLICATION_TYPE[] == :daily ? 252^0.5 : 12^0.5

  for i ∈ 1:length(CQUART_FVALS)
    c::Int = 1+(i-1)*3

    local sspanelq1::SubDataFrame
    local sspanelq4::SubDataFrame
    local allfields::Vector{Symbol}

    #acquire the fields we will use
    f::Symbol = CQUART_FVALS[i]
    name::Symbol = Symbol("P_", f)
    Fgroup::Symbol = cqports[name].Fgroup
    basicfields::Vector{Symbol} = [:date; :permno; Fgroup]
    local spanel::SubDataFrame = view(panel, completecases(panel[!,basicfields]),:)

    ###Value Rows
    Fval::Symbol = T1_PANEL_B_VALCOLIDX[f]
    #allfields = [basicfields; Fval]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1val::Float64 = mean(skipmissing(sspanelq1[!, Fval]))
    q4val::Float64 = mean(skipmissing(sspanelq4[!, Fval]))
    desccontent[T1_PANEL_B_IDX[:leverageQmeans]][c] = (
      "$(n2s(q1val, 1))")
    desccontent[T1_PANEL_B_IDX[:leverageQmeans]][c+1] = ("$(n2s(q4val, 1))")

    desccontent[T1_PANEL_B_IDX[:quartilenum]][c] = ("\\text{Q1}")
    desccontent[T1_PANEL_B_IDX[:quartilenum]][c+1] = ("\\text{Q4}")

    ###SD Rows
    #allfields = [basicfields; T1_PANEL_B_VOLFIELDS]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1LLvol::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_VOLFIELDS[1]]) * volscale
    q4LLvol::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_VOLFIELDS[1]]) * volscale
    desccontent[T1_PANEL_B_IDX[:LLsd]][c] = ("$(n2s(q1LLvol))")
    desccontent[T1_PANEL_B_IDX[:LLsd]][c+1] = ("$(n2s(q4LLvol))")

    q1vol::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_VOLFIELDS[2]]) * volscale
    q4vol::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_VOLFIELDS[2]]) * volscale
    desccontent[T1_PANEL_B_IDX[:sd]][c] = ("$(n2s(q1vol))")
    desccontent[T1_PANEL_B_IDX[:sd]][c+1] = ("$(n2s(q4vol))")

    Δq1vol::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_VOLFIELDS[2]] .-
      sspanelq1[!, T1_PANEL_B_VOLFIELDS[1]]) * volscale
    Δq4vol::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_VOLFIELDS[2]] .-
      sspanelq4[!, T1_PANEL_B_VOLFIELDS[1]]) * volscale
    desccontent[T1_PANEL_B_IDX[:Δsd]][c] = ("$(n2s(Δq1vol))")
    desccontent[T1_PANEL_B_IDX[:Δsd]][c+1] = ("$(n2s(Δq4vol))")

    ###Raw Ret Rows
    #allfields = [basicfields; T1_PANEL_B_RETFIELDS]
    sspanelq1 = view(spanel, spanel[!,Fgroup] .== 1, :)
    sspanelq4 = view(spanel, spanel[!,Fgroup] .== 4, :)

    q1LLret::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_RETFIELDS[1]])
    q4LLret::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_RETFIELDS[1]])
    desccontent[T1_PANEL_B_IDX[:LLret12m]][c] = ("$(n2s(q1LLret))")
    desccontent[T1_PANEL_B_IDX[:LLret12m]][c+1] = ("$(n2s(q4LLret))")

    q1ret::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_RETFIELDS[2]])
    q4ret::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_RETFIELDS[2]])
    desccontent[T1_PANEL_B_IDX[:ret12m]][c] = ("$(n2s(q1ret))")
    desccontent[T1_PANEL_B_IDX[:ret12m]][c+1] = ("$(n2s(q4ret))")

    Δq1ret::Float64 = skipmean(sspanelq1[!, T1_PANEL_B_RETFIELDS[2]] .-
      sspanelq1[!, T1_PANEL_B_RETFIELDS[1]])
    Δq4ret::Float64 = skipmean(sspanelq4[!, T1_PANEL_B_RETFIELDS[2]] .-
      sspanelq4[!, T1_PANEL_B_RETFIELDS[1]])
    desccontent[T1_PANEL_B_IDX[:Δret12m]][c] = ("$(n2s(Δq1ret))")
    desccontent[T1_PANEL_B_IDX[:Δret12m]][c+1] = ("$(n2s(Δq4ret))")

    if c+1≠numcols #hacks together a spacer column
      desccontent[T1_PANEL_B_IDX[:leverageQmeans]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:quartilenum]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:LLsd]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:sd]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:Δsd]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:LLret12m]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:ret12m]][c+2] = ("")
      desccontent[T1_PANEL_B_IDX[:Δret12m]][c+2] = ("")
    end

  end

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function table1panelC(panel::DataFrame;
    tablepath::String = TABLE_PATH,
    tablename::String = "t1panelc-$(REPLICATION_TYPE[])",
    decimals::Int = DEFAULT_DECIMALS)
  local cqports::CQuartileConstructions = leverageconstructions()
  Fgrps::Vector{Symbol} = Fgroups(cqports)

  #NOTE: These are of the form mknegcash ~ α + β*mkt ...

  ################INPUTS############
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

  T1_PANEL_C_FOCAL::Dict = Dict(
    :constant=>CQUART_FNAMES,
    :plusxmkt=>CQUART_FNAMES,
    :plussmbhml=>CQUART_FNAMES,
    :plusrmwcma=>CQUART_FNAMES,
    :plusumd=>CQUART_FNAMES)

  #################################
  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(T1_PANEL_C_ROWKEYS) * length(Fgrps))

  #need to re-shape for this exercise
  idcols = unique([:date; values(T1_PANEL_C_CONTROLS)...;])
  measurecols = CQUART_FNAMES
  allcols = unique([:date; values(T1_PANEL_C_CONTROLS)...; values(T1_PANEL_C_FOCAL)...;])
  minipanel::DataFrame = unique(panel[!,allcols])
  minipanel = melt(minipanel, idcols, measurecols, variable_name = :portname, value_name = :ret)
  #println(describe(minipanel))
  #println("Size minipanel: $(size(minipanel))")

  #make the regression specifications
  for rname ∈ T1_PANEL_C_ROWKEYS
    sallcols::Vector{Symbol} = [:ret; :date; :portname; T1_PANEL_C_CONTROLS[rname]]
    sminipanel = view(minipanel, completecases(minipanel[!,sallcols]), sallcols)
    push!(fss, FactorSpecification(sminipanel, :ret, :date, :portname,
      Ffactors = nothing,
      Fcharacteristics=nothing,
      Fcontrols = T1_PANEL_C_CONTROLS[rname],
      name=Symbol("row$(rname)"),
      skipcrosssectiondfcheck=true)) #only doing time series regressions
  end

  #run the regressions
  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]
  frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #generate and push all the regressions
  #each call to manyFF pulls a row of regressions
  (fs->push!(frs, manyFF(fs, parallel=false)...)).(fss)

  #make the table
  colnames::Vector{Vector{String}} = [
    (string).(CQUART_FVALS), (i->"Alpha (T-stat)").(1:length(CQUART_FVALS))]
  descrownames::Vector{String} = (s->T1_PANEL_C_ROWLABEL[s]).(T1_PANEL_C_ROWKEYS)

  desccontent::Vector{Vector{String}} = factordesccontent(
    fr::FactorRegression->coeftstat(fr, i=1), fss=fss, frs=frs,
    rowkeys = T1_PANEL_C_ROWKEYS, colindex=T1_PANEL_C_FOCAL)

  tt::String = textable(colnames=colnames, descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function table1panelD(panel::DataFrame;
    tablepath::String = TABLE_PATH,
    tablename::String = "t1paneld-$(REPLICATION_TYPE[])")
  local cqports::CQuartileConstructions = leverageconstructions()
  Fgrps::Vector{Symbol} = Fgroups(cqports)


  ################INPUTS############
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

  T1_PANEL_D_FOCAL::Dict = Dict(
    :quartiles=>CQUART_FGROUPS,
    :quartiles4ctrl=>CQUART_FGROUPS,
    :quartiles10ctrl=>CQUART_FGROUPS,
    :values=>CQUART_FVALS,
    :values4ctrl=>CQUART_FVALS,
    :values10ctrl=>CQUART_FVALS,
    :quartiles_ac=>CQUART_FGROUPS_ALT,
    :quartiles_lret=>CQUART_FGROUPS)
  #################################

  fss::Vector{FactorSpecification} = Vector{FactorSpecification}()
  sizehint!(fss, length(T1_PANEL_D_ROWKEYS) * length(Fgrps))

  #make the regression specifications
  @time for rname ∈ T1_PANEL_D_ROWKEYS
    for focal ∈ T1_PANEL_D_FOCAL[rname]
      push!(fss, FactorSpecification(panel, T1_PANEL_D_RET[rname], :date, :permno,
        Ffactors=nothing,
        Fcharacteristics = [focal],
        Fcontrols = T1_PANEL_D_CONTROLS[rname],
        name=Symbol("row$(rname)_col$(focal)"),
        skiptimeseriesdfcheck=true)) #only doing crossectional regressions
    end
  end

  #Run the regressions
  #WARNING: When BLAS paralllelization issues are solved replace with parallel=PARALLEL[]
  frs::Vector{FactorRegression} = (fs->characteristicFM(fs, parallel=false)).(fss)

  #make the table
  colnames::Vector{Vector{String}} = [
    (string).(CQUART_FVALS), (i->"Coef (T-stat)").(1:length(CQUART_FVALS))]
  descrownames::Vector{String} = (s->T1_PANEL_D_ROWLABEL[s]).(T1_PANEL_D_ROWKEYS)
  desccontent::Vector{Vector{String}} = factordesccontent(fss=fss, frs=frs,
    rowkeys = T1_PANEL_D_ROWKEYS, colindex=T1_PANEL_D_FOCAL)

  tt::String = textable(colnames=colnames,
    descrownames=descrownames, desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end


function coeftstat(fr::FactorRegression; decimals::Int=DEFAULT_DECIMALS, i::Int=2,
  )::String
  λᵢ::Float64 = fr.λᵢ[i] * 100
  tstat::Float64 = fr.λᵢ[i]/fr.σᵢ[i]
  return "$(num2str(λᵢ, decimals))\\;($(num2str(tstat, decimals)))"
end

#this is a kernel function to generate the content matrix under a number of different scenarios
function factordesccontent(f::Function = coeftstat;
    fss::Vector{FactorSpecification} = error("fss is required for factordesccontent"),
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
