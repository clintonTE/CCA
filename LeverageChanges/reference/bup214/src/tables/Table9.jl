#WARNING- look, think about market caps and scaling v log
#treat thsi table as a one-off since it has no common elements w/ other tables
function maketable9(cumex::AbstractDataFrame;
    panela::Bool=true,
    panelb::Bool=true,
    panelc::Bool=true)


  scumex::SubDataFrame = view(cumex, cumex.primary, :)

  ###################PANEL A############

  T9_PANEL_A_NAME::String = "t9panela-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  T9_PANEL_A_ROWKEYS::Vector{Symbol} = [
    :all,
    :lowyield,
    :midyield,
    :highyield]

  T9_PANEL_A_ROWLABEL::Dict = Dict(
    :all => "All",
    :lowyield=> "\\midrule Low Div Yield \$(\\delta < 0.75\\%)\$",
    :midyield=> "Mid Div Yield \$(0.75\\% < \\delta < 1.5\\%)\$",
    :highyield=> "High Div Yield \$(1.5\\%) < \\delta\$",)

  T9_PANEL_A_COLKEYS::Vector{Symbol} = [
    :yield,
    :riskcum,
    :riskex,
    :riskdiff,
    :rewardcum,
    :rewardex,
    :rewarddiff,
    :N]

  T9_PANEL_A_COLLABEL::Dict = Dict(
    :yield=>"DivYield",
    :riskcum=>"Cum",
    :riskex=>"Ex",
    :riskdiff=>"Diff",
    :rewardcum=>"Cum",
    :rewardex=>"Ex",
    :rewarddiff=>"Diff",
    :N=>"N")

  T9_PANEL_A_COLNAMES::Vector{Vector{String}} =[
    ["", "``Risk'': \$\\sum{|r|}\$", "``Reward'': \$\\sum{r}\$", ""],
    (s::Symbol->T9_PANEL_A_COLLABEL[s]).(T9_PANEL_A_COLKEYS)]

  T9_PANEL_A_WIDTHCOLNAMES::Vector{Vector{Int}} = [[1,3,3,1], ones(Int,length(T9_PANEL_A_COLKEYS))]

  panela && table9panela(scumex,
      tablename = T9_PANEL_A_NAME,
      rowkeys = T9_PANEL_A_ROWKEYS,
      rowlabel = T9_PANEL_A_ROWLABEL,
      colkeys = T9_PANEL_A_COLKEYS,
      colnames = T9_PANEL_A_COLNAMES,
      widthcolnames = T9_PANEL_A_WIDTHCOLNAMES
      )


  ###################PANEL B/C############

  T9_PANEL_BC_ROWKEYS::Vector{Symbol} = [
    :interceptδ,
    :yieldδ,
    :interceptδ1δ,
    :yieldδ1δ,
    :exdayret]

  T9_PANEL_BC_ROWLABEL::Dict = Dict(
    :interceptδ => "Intercept",
    :yieldδ => "Div Yield \$\\delta\$",
    :interceptδ1δ => "\\midrule Intercept",
    :yieldδ1δ => "\$\\delta/(1-\\delta)\$",
    :exdayret => "Ex-Day Return"
    )

  T9_PANEL_BC_COLKEYS::Vector{Symbol} = [
    :coef,
    :se,
    :tstat,
    :stdcoef,
    :mwse,
    :mwtstat]

  T9_PANEL_BC_COLLABEL::Dict = Dict(
    :coef=>"Coef",
    :se=>"StdErr",
    :tstat=>"T-stat",
    :stdcoef=>"StdCoef",
    :mwse=>"StdErr",
    :mwtstat=>"T-stat")

  T9_PANEL_BC_COLNAMES::Vector{Vector{String}} =[
    ["", "Heteroskedastic"],
    (s::Symbol->T9_PANEL_BC_COLLABEL[s]).(T9_PANEL_BC_COLKEYS)]

  T9_PANEL_BC_WIDTHCOLNAMES::Vector{Vector{Int}} = [[4,2], ones(Int,length(T9_PANEL_BC_COLKEYS))]

  T9_PANEL_BC_DECIMALS = 3

  #Now do panel b
  T9_PANEL_B_YSYM = :Dabs
  T9_PANEL_B_YSYM_UNIT = :UDabs
  T9_PANEL_B_NAME::String = "t9panelb-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  panelb && table9panelbc(scumex,
    Ysym = T9_PANEL_B_YSYM,
    Ysymunit = T9_PANEL_B_YSYM_UNIT,
    tablename = T9_PANEL_B_NAME,
    rowkeys = T9_PANEL_BC_ROWKEYS,
    rowlabel = T9_PANEL_BC_ROWLABEL,
    colkeys = T9_PANEL_BC_COLKEYS,
    colnames = T9_PANEL_BC_COLNAMES,
    widthcolnames = T9_PANEL_BC_WIDTHCOLNAMES,
    decimals = T9_PANEL_BC_DECIMALS
    )

  #Now do panel c
  T9_PANEL_C_YSYM = :Dret
  T9_PANEL_C_YSYM_UNIT = :UDret
  T9_PANEL_C_NAME::String = "t9panelc-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"
  panelc && table9panelbc(scumex,
    Ysym = T9_PANEL_C_YSYM,
    Ysymunit = T9_PANEL_C_YSYM_UNIT,
    tablename = T9_PANEL_C_NAME,
    rowkeys = T9_PANEL_BC_ROWKEYS,
    rowlabel = T9_PANEL_BC_ROWLABEL,
    colkeys = T9_PANEL_BC_COLKEYS,
    colnames = T9_PANEL_BC_COLNAMES,
    widthcolnames = T9_PANEL_BC_WIDTHCOLNAMES,
    decimals = T9_PANEL_BC_DECIMALS
    )
end



function table9panelbc(cumex::AbstractDataFrame;
  Ysym::Symbol = error("Ysym is required"),
  Ysymunit::Symbol = error("Ysymunit is required"),
  tablename::String = error("tablename is required"),
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Symbol} = error("colkeys is required"),
  colnames::Vector{Vector{String}} = error("colnames is required"),
  widthcolnames::Vector{Vector{Int}} = error("widthcolnames is required"),
  tablepath::String = TABLE_PATH,
  decimals::Int = DEFAULT_DECIMALS,
  )

  #pre-declare some regresion variables taht will be assigned multiple times
  local Xexpr::FMExpr
  local Xnames::Vector{Symbol}
  local reg::FMLM
  local se::Matrix{Float64}
  local semw::Matrix{Float64}
  local scumex::SubDataFrame

  #housekeeping and pre-allocation
  colidx::Dict = Dict(colkey=>i for (i,colkey) ∈ enumerate(colkeys))
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, sum(widthcolnames[1]))).(1:length(rowkeys))
  rowidx::Dict = Dict(rowkey=>desccontent[i] for (i,rowkey) ∈ enumerate(rowkeys))
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)

  scumex = view(cumex, :, [
    :yield, :yieldd1d, :exdayret, :Uyield, :Uyieldd1d, :Uexdayret, Ysym, Ysymunit])
  scumex = view(scumex, completecases(scumex), :)

  #constructs a row for the table from the regression object and the name of the coefficient
  function regrow!(reg::FMLM,  Fcoef::Symbol, descrow::Vector{String},
    Σ::Matrix{Float64}, Σmw::Matrix{Float64})

    ind::Int = findfirst(isequal(Fcoef), reg.Xnames)
    @assert Fcoef == reg.Xnames[ind]

    #display(Σ)

    coef::Float64 = reg.β[ind]
    sigma::Float64 = Σ[ind,ind] ^0.5
    tstat::Float64 = coef/sigma
    mwsigma::Float64 = Σmw[ind,ind]^0.5
    mwtstat::Float64 = coef/mwsigma

    descrow[colidx[:coef]] = brnum(coef, decimals, scalefactor=100.)
    descrow[colidx[:se]] = num2str(sigma, decimals, scalefactor=100.)
    descrow[colidx[:tstat]] = brnum(tstat, decimals, scalefactor=1.)
    descrow[colidx[:mwse]] = num2str(mwsigma, decimals, scalefactor=100.)
    descrow[colidx[:mwtstat]] = brnum(mwtstat, decimals, scalefactor=1.)
  end


  #do the bivariate regressions
  Xexpr = :yield
  Xnames= [:intercept, :yield]
  reg = FMLM(scumex,  Xexpr, Ysym, Xnames=Xnames, Yname=Ysym, containsmissings=false)
  se = homoskedasticΣ!(reg)
  semw = modifiedwhiteΣ!(reg)
  regrow!(reg, :intercept, rowidx[:interceptδ], se, semw)
  regrow!(reg, :yield, rowidx[:yieldδ], se, semw)

  #unit-standardized bivariate regression
  Xexpr = :Uyield
  reg = FMLM(scumex,  Xexpr, Ysymunit, Xnames=Xnames, Yname=Ysym, containsmissings=false)
  @assert (reg.Xnames[2] == :yield) #name check
  rowidx[:yieldδ][colidx[:stdcoef]] = brnum(reg.β[2], decimals, scalefactor=100.)
  rowidx[:interceptδ][colidx[:stdcoef]] = ""

  #multi-variate regression
  Xexpr = Meta.parse("yieldd1d + exdayret")
  Xnames= [:intercept, :yieldd1d, :exdayret]
  reg = FMLM(scumex,  Xexpr, Ysym, Xnames=Xnames, Yname=Ysym, containsmissings=false)
  se = homoskedasticΣ!(reg)
  semw = modifiedwhiteΣ!(reg)
  regrow!(reg, :intercept, rowidx[:interceptδ1δ], se, semw)
  regrow!(reg, :yieldd1d, rowidx[:yieldδ1δ], se, semw)
  regrow!(reg, :exdayret, rowidx[:exdayret], se, semw)

  #unit-standardized multivariate regression
  Xexpr = Meta.parse("Uyieldd1d + Uexdayret")
  reg = FMLM(scumex,  Xexpr, Ysymunit, Xnames=Xnames, Yname=Ysym, containsmissings=false)
  @assert (reg.Xnames[2] == :yieldd1d) && (reg.Xnames[3] == :exdayret)
  rowidx[:yieldδ1δ][colidx[:stdcoef]] = brnum(reg.β[2], decimals, scalefactor=100.)
  rowidx[:exdayret][colidx[:stdcoef]] = brnum(reg.β[3], decimals, scalefactor=100.)
  rowidx[:interceptδ1δ][colidx[:stdcoef]] = ""

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    #alignmentstring=alignmentstring,
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function table9panela(cumex::AbstractDataFrame;
  tablename::String = error("tablename is required"),
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Symbol} = error("colkeys is required"),
  colnames::Vector{Vector{String}} = error("colnames is required"),
  widthcolnames::Vector{Vector{Int}} = error("widthcolnames is required"),
  tablepath::String = TABLE_PATH,
  decimals::Int = DEFAULT_DECIMALS,
  )

  colidx::Dict = Dict(colkey=>i for (i,colkey) ∈ enumerate(colkeys))
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, sum(widthcolnames[1]))).(1:length(rowkeys))
  rowidx::Dict = Dict(rowkey=>desccontent[i] for (i,rowkey) ∈ enumerate(rowkeys))
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)



  @inline function tablerow!(df::AbstractDataFrame, descrow::Vector{String})
    #first get our datasets


    #get the values
    yield::Float64 = mean(df.yield)
    cumrisk::Float64 = mean(df.cumabs)
    exrisk::Float64 = mean(df.exabs)
    Δrisk::Float64 = cumrisk - exrisk

    cumret::Float64 = mean(df.cumret)
    exret::Float64 = mean(df.exret)
    Δret::Float64 = cumret - exret
    N::Int = size(df,1)

    #write the row
    descrow[colidx[:yield]] = num2str(yield, decimals, scalefactor=100.0)
    descrow[colidx[:riskcum]] = brnum(cumrisk, decimals, scalefactor=100.0)
    descrow[colidx[:riskex]] = brnum(exrisk, decimals, scalefactor=100.0)
    descrow[colidx[:riskdiff]] = brnum(Δrisk, decimals, scalefactor=100.0)
    descrow[colidx[:rewardcum]] = brnum(cumret, decimals, scalefactor=100.0)
    descrow[colidx[:rewardex]] = brnum(exret, decimals, scalefactor=100.0)
    descrow[colidx[:rewarddiff]] = brnum(Δret, decimals, scalefactor=100.0)
    descrow[colidx[:N]] = num2str(N, decimals, Ints=true)
  end

  #write the row using the appropriate data
  tablerow!(cumex, rowidx[:all])

  scumex = view(cumex,cumex.Gyield .== :lowyield, :)
  tablerow!(scumex, rowidx[:lowyield])

  scumex = view(cumex,cumex.Gyield .== :medyield, :)
  tablerow!(scumex, rowidx[:midyield])

  scumex = view(cumex,cumex.Gyield .== :highyield, :)
  tablerow!(scumex, rowidx[:highyield])

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    #alignmentstring=alignmentstring,
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function aggregateoverdays(F::Function, dist::AbstractDataFrame, dayrange::AbstractRange;
  Fgroup::Symbol = :eid, selectcols::Vector{Symbol} = names(dist), suffix::String = "",
  colnames::Vector{Symbol} = names(dist), Fcriteria::Symbol = :day)

  sdist::SubDataFrame = view(dist, (g::Union{MInt, Symbol}->
    (!ismissing(g)) && g ∈ dayrange).(dist[!,Fcriteria]), selectcols)


  agg::DataFrame = aggregate(groupby(sdist, Fgroup), F)
  rename!(agg, colnames)



  return agg
end
