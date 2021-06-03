#WARNING- look, think about market caps and scaling v log
#treat thsi table as a one-off since it has no common elements w/ other tables
function maketable10(cumex::AbstractDataFrame)

  T10_ROWKEYS::Vector{Symbol} = [
    :all1232,
    :lowyield,
    :medyield,
    :highyield,
    :taxdiv,
    :nontaxdiv,
    :monthlydiv,
    :quarterlydiv,
    :semiandannualdiv,
    :specialdiv,
    :nonspecialdiv,
    :lastdiv40,
    :lastdiv4080,
    :lastdiv80,
    :compustat,
    :noncompustat,
    :bkliab13,
    :bkliab1323,
    :bkliab23,
    :assets100,
    :assets1003000,
    :assets3000]

  T10_ROWLABEL::Dict = Dict(
    :all1232 => "All (Code1232)",
    :lowyield => "\\midrule Low Yield \$(\\delta < 0.75\\%)\$",
    :medyield => "Mid Yield \$(0.75\\% < \\delta < 1.5\\%)\$",
    :highyield => "High Yield \$(1.5\\%) < \\delta\$",
    :taxdiv => "\\midrule Taxable Dividend",
    :nontaxdiv => "Non-Taxable",
    :monthlydiv => "\\midrule Monthly Dividend",
    :quarterlydiv => "Quarterly",
    :semiandannualdiv => "Semi- and Annual",
    :specialdiv => "\\midrule Special Dividend",
    :nonspecialdiv => "Not Special",
    :lastdiv40 => "\\midrule Last Dividend <40days",
    :lastdiv4080 => "40-80days",
    :lastdiv80 => ">80days",
    :compustat => "\\midrule w/ Compustat (t-1)",
    :noncompustat => "w/o Compustat",
    :bkliab13 => "\\midrule bkliab < 1/3",
    :bkliab1323 => "bkliab 1/3 - 2/3",
    :bkliab23 => "bkliab > 2/3",
    :assets100 => "\\midrule Small (Assets < \\\$100m)",
    :assets1003000 => "Medium (Assets \\\$100m-\\\$3,000m)",
    :assets3000 => "Large (Assets < \\\$3000m)"
    )

  T10_COLKEYS::Vector{Symbol} = [
    :in,
    :volcoef,
    :voltstat,
    :retcoef,
    :rettstat,
    :N]

  T10_COLLABEL::Dict = Dict(
    :in=>"in \\%",
    :volcoef=>"Coef",
    :voltstat=>"T-stat",
    :retcoef=>"Coef",
    :rettstat=>"T-stat",
    :N=>"N")

  T10_COLNAMES::Vector{Vector{String}} =[
    ["", "Cum-Ex Change from Dividend", ""],
    ["Avg \$\\delta\$", "Volatility", "Avg Return", ""],
    (s::Symbol->T10_COLLABEL[s]).(T10_COLKEYS)]

  T10_WIDTHCOLNAMES::Vector{Vector{Int}} = [
    [1,4,1],
    [1,2,2,1],
    ones(Int,length(T10_COLKEYS))]

  T10_ALIGNMENTCOLNAMES =   [
    ["c", "c", "c"],
    ["c", "c", "c", "c"],
    ["r" for i in 1:length(T10_COLKEYS)]]

  T10_DECIMALS = 1
  T10_COEF_DECIMALS = 2

  T10_ABS_YSYM = :Dabs
  T10_RET_YSYM = :Dret
  T10_NAME::String = "t10-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  table10(cumex,
    Yabssym = T10_ABS_YSYM,
    Yretsym = T10_RET_YSYM,
    tablename = T10_NAME,
    rowkeys = T10_ROWKEYS,
    rowlabel = T10_ROWLABEL,
    colkeys = T10_COLKEYS,
    colnames = T10_COLNAMES,
    widthcolnames = T10_WIDTHCOLNAMES,
    decimals = T10_DECIMALS,
    coefdecimals = T10_COEF_DECIMALS,
    alignmentcolnames=T10_ALIGNMENTCOLNAMES
    )

end



function table10(cumex::AbstractDataFrame;
  Yabssym::Symbol = error("Yvolsym is required"),
  Yretsym::Symbol = error("Yretsym is required"),
  tablename::String = error("tablename is required"),
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Symbol} = error("colkeys is required"),
  colnames::Vector{Vector{String}} = error("colnames is required"),
  widthcolnames::Vector{Vector{Int}} = error("widthcolnames is required"),
  tablepath::String = TABLE_PATH,
  decimals::Int = error("decimals is required"),
  coefdecimals::Int = error("coefdecimals is required"),
  alignmentcolnames::Vector{Vector{String}} = error("alignmentcolnames is required")
  )

  #some of the regression constants and parameters should overlap
  local Xexpr = :yield
  local Xnames= [:intercept, :yield]

  #housekeeping and pre-allocation
  colidx::Dict = Dict(colkey=>i for (i,colkey) ∈ enumerate(colkeys))
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, sum(widthcolnames[1]))).(1:length(rowkeys))
  rowidx::Dict = Dict(rowkey=>desccontent[i] for (i,rowkey) ∈ enumerate(rowkeys))
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)

  #apply the conditions that will apply to all rows
  #NOTE: Already filtered in formation process
  #scumex = view(cumex, dist.match .* dist.oneex .* dist.dclrexok, :)
  scumex = view(cumex, completecases(cumex[!, [:yield, Yretsym, Yabssym]]), :)

  #constructs a row for the table from the regression object and the name of the coefficient
  @inline function regrow!(rdf::AbstractDataFrame,  descrow::Vector{String})
    #start by running the regressions
    retreg = FMLM(rdf,  Xexpr, Yretsym, Xnames=Xnames, Yname=Yretsym, containsmissings=false)
    retΣ::Matrix{Float64} = homoskedasticΣ!(retreg)

    absreg = FMLM(rdf,  Xexpr, Yabssym, Xnames=Xnames, Yname=Yabssym, containsmissings=false)
    absΣ::Matrix{Float64} = homoskedasticΣ!(absreg)

    #make sure the names align
    @assert retreg.Xnames[2] == :yield
    @assert absreg.Xnames[2] == :yield

    #get the values for the table row
    in::Float64 = mean(rdf.yield)
    volcoef::Float64 = absreg.β[2]
    voltstat::Float64 = absreg.β[2] / absΣ[2,2]^0.5
    retcoef::Float64  = retreg.β[2]
    rettstat::Float64 = retreg.β[2] / retΣ[2,2]^0.5
    N::Float64 = size(rdf,1)

    #format and write the output
    descrow[colidx[:in]] = num2str(in, decimals, scalefactor=100.)
    descrow[colidx[:volcoef]] = brnum(volcoef, coefdecimals, scalefactor=100.)
    descrow[colidx[:voltstat]] = brnum(voltstat, decimals, scalefactor=1.)
    descrow[colidx[:retcoef]] = brnum(retcoef, coefdecimals, scalefactor=100.)
    descrow[colidx[:rettstat]] = brnum(rettstat, decimals, scalefactor=1.)
    descrow[colidx[:N]] = num2str(N, 0, Ints=true)
  end

  #first row is all rows (same result as before)
  scumex1232::SubDataFrame = view(scumex, scumex.ordinary, :)
  regrow!(scumex1232, rowidx[:all1232])

  #performs the regressions within each group and builds the corresponding table rows
  @inline function reggroupby!(df::AbstractDataFrame, rowidx::Dict, Fgroup::Symbol)
    local sdf::SubDataFrame = view(df, (!ismissing).(df[!,Fgroup]), :)

    (size(sdf,1) > 10) || error(
      "sdf below minimum size for effective groupings of column $Fgroup")

    #for each group classification
    for ssdf::SubDataFrame ∈ groupby(sdf, Fgroup)
      local grp::Symbol = ssdf[1,Fgroup]
      haskey(rowidx, grp) || error(
        "Group classification $grp from column $Fgroup not found as key in rowidx")

      regrow!(ssdf, rowidx[grp])
    end
  end

  #build the table rows
  reggroupby!(scumex1232, rowidx, :Gyield)
  reggroupby!(scumex, rowidx, :Gtaxdiv)
  reggroupby!(scumex, rowidx, :Gfreqdiv)
  reggroupby!(scumex, rowidx, :Gspecialdiv)
  reggroupby!(scumex, rowidx, :Glastdiv)
  reggroupby!(scumex1232, rowidx, :Gcompustat)
  reggroupby!(scumex1232, rowidx, :Gbkliab)
  reggroupby!(scumex1232, rowidx, :Gassets)

  #println(desccontent)

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    alignmentcolnames=alignmentcolnames,
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end
