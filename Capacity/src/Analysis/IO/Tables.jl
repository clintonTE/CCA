
verificationstats=[:cor, :betamb, :spear, #=:kend, =# :N]
verificationlabels=Dict{Symbol, Any}(
  :cor=>"Pearson Cor",
    :cort=>"t-stat",
    :corhac=>"HAC z",
    :cormbbz=>"MBB z",
    :cormbbp=>"MBB p",
    :corar1h0z=>"AR(1) \$H_0\$ z",
    :corar1h0p=>"AR(1) \$H_0\$ p",
  :spear=>"\\midrule Spearman Cor",
    :speart=>"t-stat",
    :spearmbbz=>"MBB z",
    :spearmbbp=>"MBB p",
    :spearar1h0z=>"AR(1) \$H_0\$ z",
    :spearar1h0p=>"AR(1) \$H_0\$ p",
  :kend=>"\\midrule Kendall's Tau",
    :kendz=>"z-stat",
    :kendmbbz=>"MBB z",
    :kendmbbp=>"MBB p",
    :kendar1h0z=>"AR(1) \$H_0\$ z",
    :kendar1h0p=>"AR(1) \$H_0\$ p",
  :betabm=>"\\midrule Beta", #bm=benchmark v measure
    :betabmmw=>"HC z",
    :betabmhac=>"HAC z",
    :betabmmbbz=>"MBB z",
    :betabmmbbp=>"MBB p",
    :betabmar1h0z=>"AR(1) \$H_0\$ z",
    :betabmar1h0p=>"AR(1) \$H_0\$ p",
  :betamb=>"\\midrule Beta", #bm=benchmark v measure
    :betambmw=>"HC z",
    :betambhac=>"HAC z",
    :betambmbbz=>"MBB z",
    :betambmbbp=>"MBB p",
    :betambar1h0z=>"AR(1) \$H_0\$ z",
    :betambar1h0p=>"AR(1) \$H_0\$ p",
  :N=>"\\midrule Months",
  :comom=>Dict(
    :lG_WXplmomentum_LP366d=>"Mom", #industry adjusted
    :lGVP_WXplmomentum_LP366d=>"Mom Participation",
    :lGP_WXplmomentum_LP366d=>"Mom Pt-Weighted",
    :lG_WXplmomentum121s_LP366d=>"Mom", #non-industry adjusted
    :lGVP_WXplmomentum121s_LP366d=>"Mom Participation",
    :lGP_WXplmomentum121s_LP366d=>"Mom Pt-Weighted",
    :M_WXplmomentum121s_LP366d=>"\$Mom/|ME|\$",
    :GlG_WXplmomentum121s_LP366d=>"\$\\Delta\$ Mom",
    :GlGP_WXplmomentum121s_LP366d=>"\$\\Delta\$ Mom Pt-Weighted",
    :GM_WXplmomentum121s_LP366d=>"\$\\Delta Mom/|ME|\$",
  ),
  :measureadj=>Dict(
#    :lG_A_aumequityXw_WXmomentum121_LP366d=>"\$w=\\beta\$",
    :lG_A_aumequityXabsw_WXmomentum121s_LP366d=>"\$w=|\\beta_i|\$",
    :lG_A_aumequityXabscorw_WXmomentum121s_LP366d=>"\$w=|\\beta_i|_{>0.1}\$",
  ),
  :umdadj=>Dict(
  #  :lG_A_aumequityXw_F_ff3m_umd_LP366d=>"\$w=\\beta\$",
    :lG_A_aumequityXabsw_F_ff3m_umd_LP366d=>"\$w=|\\beta_i|\$",
    :lG_A_aumequityXabscorw_F_ff3m_umd_LP366d=>"\$w=|\\beta_i|_{>0.1}\$",
  ),
)

verificationsubstats=Dict{Symbol, Vector{Symbol}}(
  :cor=>[#=:cort, =#:corhac, :cormbbz, :cormbbp, :corar1h0z, :corar1h0p],
  :spear=>[#=:speart,=# :spearmbbz, :spearmbbp, :spearar1h0z, :spearar1h0p],
  :kend=>[#=:kendz,=# :kendmbbz, :kendmbbp, :kendar1h0z, :kendar1h0p],
  :betabm=>[#=:betamw,=# :betabmhac, :betabmmbbz, :betabmmbbp, :betabmar1h0z, :betabmar1h0p],
  :betamb=>[#=:betamw,=# :betambhac, :betambmbbz, :betambmbbp, :betambar1h0z, :betambar1h0p],
  :N=>Symbol[])

verificationformats = Dict(
  :cor=>(;scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:"),
  :cort=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :corhac=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :cormbbz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :cormbbp=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),
  :corar1h0z=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :corar1h0p=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),

  :spear=>(;scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:"),
  :speart=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :spearmbbz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :spearmbbp=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}") ,
  :spearar1h0z=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :spearar1h0p=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),

  :kend=>(;scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:"),
  :kendz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :kendmbbz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :kendmbbp=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),
  :kendar1h0z=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :kendar1h0p=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),

  :betabm=>(;scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:"),
  :betabmmw=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betabmhac=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betabmmbbz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betabmmbbp=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),
  :betabmar1h0z=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betabmar1h0p=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),

  :betamb=>(;scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:"),
  :betambmw=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betambhac=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betambmbbz=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betambmbbp=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),
  :betambar1h0z=>(;scalefactor=1.0, decimals=2, prefix="\\mathit{", suffix="}\\:\\:"),
  :betambar1h0p=>(;scalefactor=1.0, decimals=3, prefix="\\mathit{", suffix="}"),
  :N=>(;scalefactor=1.0, decimals=0, prefix="", suffix=""),
)

verificationrows = [(s->[s;verificationsubstats[s];]).(verificationstats)...;]

#=mfbenchmarkcols = [
  :lG_A_aumequityXw_WXmomentum121_LP366d
  :lG_A_aumequityXw_F_ff3m_umd_LP366d
  :lG_A_aumequityXabsw_WXmomentum121_LP366d,
  :lG_A_aumequityXabscorw_WXmomentum121_LP366d,
]=#


#=comombenchmarkcols = [ #NOTE- these are the industry adjusted versions
  :lG_WXplmomentum_LP366d,
  :lGVP_WXplmomentum_LP366d,
  :lGP_WXplmomentum_LP366d,]=#
stdcomombenchmarkcols = [
  :lG_WXplmomentum121_LP366d,
  :lGP_WXplmomentum121_LP366d,
  :M_WXplmomentum121_LP366d,
  ]
Gcomombenchmarkcols = [
  :GlG_WXplmomentum121_LP366d,
  :GlGP_WXplmomentum121_LP366d,
  :GM_WXplmomentum121_LP366d,
  ]

fundbenchmarkcols = Dict(
  :umdadj=>[
    #:lG_A_aumequityXw_F_ff3m_umd_LP366d,
    :lG_A_aumequityXabsw_F_ff3m_umd_LP366d,
    :lG_A_aumequityXabscorw_F_ff3m_umd_LP366d,],
  :measureadj => [
    #:lG_A_aumequityXw_WXmomentum121_LP366d,
    :lG_A_aumequityXabsw_WXmomentum121s_LP366d,
    :lG_A_aumequityXabscorw_WXmomentum121s_LP366d
    ])

fundmeasurecol = :lG_WXmomentum121s_LP366d

#create a table of the comomomentum results
function comomentumtables()
  comomentumtable(
    comomtablenamesuffix="",
    comombenchmarkcols=stdcomombenchmarkcols,
    benchmark=:comom)
  comomentumtable(
    comomtablenamesuffix="G",
    comombenchmarkcols=Gcomombenchmarkcols,
    benchmark=:Gcomom)
end
function comomentumtable(;
    tablepath = PARAM[:tablepath],
    measurefilename = PARAM[:tabcomommeasurefilename],
    rid = PARAM[:tabcomomrid],
    resultsname = PARAM[:tabcomomresultsname],
    analysispath = PARAM[:analysispath],
    comomtablename = PARAM[:tabcomomtablename],
    comomtablenamesuffix,
    comombenchmarkcols,
    benchmark)

  comom = CSV.File("$analysispath\\$measurefilename\\$rid\\$(rid)_$resultsname.csv") |> DataFrame
  comom.measure = comom.measure .|> Symbol
  comom.benchmark = comom.benchmark .|> Symbol
  gcomom = groupby(comom, [:measure, :benchmark]) #this effectivelly keys the df

  #get the headers
  colnames::Vector{Vector{String}} = [comombenchmarkcols .|> c->verificationlabels[:comom][c]]
  descrownames::Vector{String} = verificationrows .|> c->verificationlabels[c]

  Ncols = maximum(length.(colnames))
  Nrows = length(descrownames)
  desccontent = [Vector{String}(undef, Ncols) for i ∈ 1:Nrows]

  #extracts the values for a column of the output table from a row of the input
  function verrow2str(r::DataFrameRow)
    outcol = Vector{String}()
    for s ∈ verificationrows
      f = verificationformats[s]
      valstring = num2str(r[s], f.decimals, scalefactor=f.scalefactor)
      push!(outcol, "$(f.prefix)$(valstring)$(f.suffix)")
    end
    return outcol
  end
  #  (s->num2str(r[s], verificationdecimals[s])).(verificationrows)

  function verrow2str(r::AbstractDataFrame)
    @assert nrow(r) == 1
    verrow2str(r[1,:])
  end



  #fill in the content, double checking each row and column entry
  for (c,col) ∈ enumerate(comombenchmarkcols)
    @assert colnames[1][c] === verificationlabels[:comom][col]
    colcontent = verrow2str(gcomom[(;measure=col, benchmark=benchmark)])
    for (r, valstring) ∈ enumerate(colcontent)
      @assert descrownames[r] === verificationlabels[verificationrows[r]]
      desccontent[r][c] = valstring
    end
  end

  #now compile the content
  t = textable(;colnames, descrownames, desccontent)

  writetextable(t, path=tablepath, outname="$(comomtablenamesuffix)$(comomtablename).tex")

  return nothing
end

function fundtables()

  fundtable(fundadjtype=:umdadj)
  fundtable(fundadjtype=:measureadj)
end

#this creates a table of verification results
#I seperate by the weighting adjustment, either French's UMD or
  #returns of the actual momentum measure
function fundtable(;
    tablepath = PARAM[:tablepath],
    measurefilename = PARAM[:tabfundmeasurefilename],
    rid = PARAM[:tabfundrid],
    analysispath = PARAM[:analysispath],
    fundmeasurecol::Symbol=fundmeasurecol,
    resultsname::String=PARAM[:tabfundresultsname],
    fundadjtype::Symbol,
    fundtablename::String = PARAM[:tabfundtablename],
    fundbenchmarkcols::Vector{Symbol} = fundbenchmarkcols[fundadjtype],
    )

  #load the results
  fund = CSV.File("$analysispath\\$measurefilename\\$rid\\$(rid)_$resultsname.csv") |> DataFrame
  fund.label = fund.label .|> Symbol
  fund.measure = fund.measure .|> Symbol
  fund.benchmark = fund.benchmark .|> Symbol
  gfund = groupby(fund, [:label, :measure, :benchmark]) #this effectivelly keys the df

  #get the headers
  colnames::Vector{Vector{String}} = [["\$\\beta\$-Weighted MF AUM Growth",
    "\$\\beta\$-Weighted HF AUM Growth"],
    [fundbenchmarkcols .|> c->verificationlabels[fundadjtype][c];
      fundbenchmarkcols .|> c->verificationlabels[fundadjtype][c]]]
  widthcolnames = [[length(fundbenchmarkcols), length(fundbenchmarkcols)],
    ones(Int, length(fundbenchmarkcols) * 2)]
  alignmentcolnames = [["c" for c ∈ 1:length(colnames[1])],
    ["r" for c ∈ 1:length(colnames[2])]]
  descrownames::Vector{String} = verificationrows .|> c->verificationlabels[c]
  alignmentstring= "l | rr | rr"

  Ncols = maximum(length.(colnames))
  Nfundcols = length(fundbenchmarkcols)
  Nrows = length(descrownames)
  desccontent = [Vector{String}(undef, Ncols) for i ∈ 1:Nrows]

  #extracts the values for a column of the output table from a row of the input
  function verrow2str(r::DataFrameRow)
    outcol = Vector{String}()
    for s ∈ verificationrows
      f = verificationformats[s]
      valstring = num2str(r[s], f.decimals, scalefactor=f.scalefactor)
      push!(outcol, "$(f.prefix)$(valstring)$(f.suffix)")
    end
    return outcol
  end
  #  (s->num2str(r[s], verificationdecimals[s])).(verificationrows)

  function verrow2str(r::AbstractDataFrame)
    @assert nrow(r) == 1
    verrow2str(r[1,:])
  end

  #fill in the content, double checking each row and column entry
  #for this table, we are populating both mutual funds and hedge funds
  #but each fund type uses the same benchmark column adjustment
  for (cmf,col) ∈ enumerate(fundbenchmarkcols)
    chf = cmf + Nfundcols
    @assert colnames[2][cmf] === verificationlabels[fundadjtype][col]
    @assert colnames[2][chf] === verificationlabels[fundadjtype][col]
    colmfcontent = verrow2str(gfund[(;label=:mf, measure=fundmeasurecol, benchmark=col)])
    colhfcontent = verrow2str(gfund[(;label=:hf, measure=fundmeasurecol, benchmark=col)])

    #now iterate over the columns of content
    for (r, (mfvalstring, hfvalstring)) ∈ enumerate(zip(colmfcontent, colhfcontent))
      @assert descrownames[r] === verificationlabels[verificationrows[r]]
      desccontent[r][cmf] = mfvalstring
      desccontent[r][chf] = hfvalstring
    end
  end

  #now compile the content
  t = textable(;colnames, descrownames, desccontent,
    widthcolnames, alignmentstring, alignmentcolnames)

  writetextable(t, path=tablepath, outname="$fundtablename-$fundadjtype.tex")
end

#function univariatelltable(;measurecols, returncols)
