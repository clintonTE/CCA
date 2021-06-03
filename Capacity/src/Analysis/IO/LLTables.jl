

################helper functions for the tables and regressions

function regmbb(X::Vector{TXcol},y; bssamples, bswidth, β=nothing) where TXcol<:AbstractVector

  #draw the coefficients for the bootstrap
  function bsstatfunc(sim)
    y = sim[1]
    X = reduce(hcat, sim[2:end])

    #try
      b = cholesky!(Symmetric(X'*X))\(X'*y)
    #=catch
      println("X[1:5,:]: $(X[1:5,:])
        y[1:5]: $(y[1:5])")
      throw("cholesky factorization failed")

      #return bsstatfunc(sim)
    end=#
    return Dict{Symbol, Vector{Float64}}(
      :beta => b,
      #:betamr => ρ * std(m)/std(r)
      )
  end

  function bsaggfunc(bs)
    bsmat = reduce(hcat,bs)
    aggfunc = (;se=(v->std(skipmissing(v))).(eachrow(bsmat)), p=nonparametricp.(eachrow(bsmat)))
  end

  return mbbootstrap([y, X...,], bsstatfunc, bsaggfunc, Vector{Float64}, B=bssamples, w=bswidth)
end

#helper funciton to handle the likely scenario of an X matrix instead of a vector
regmbb(X::Matrix, y; bssamples, bswidth) = regmbb(X |> eachcol |> collect, y; bssamples, bswidth)

function regnullar1(X::AbstractMatrix,y; β, nullar1simulations, integrated)

  function nullar1statfunc(Xy)
    Xs = @view Xy[:, 1:(end-1)]
    ys = @view Xy[:, end]
    return Dict(:beta => cholesky!(Symmetric(Xs'*Xs))\(Xs'*ys),)
  end


  stats::Dict{Symbol, Vector{Float64}} = Dict(:beta=>β)
  nullar1aggfunc(M,statname) = (;
    se=std(M, dims=1),
    p=nonparametricnullp(M, stats[statname]))

  integrated && return comparearima11undernull(X,y, nullar1statfunc,
    nullar1aggfunc, Nsimulations=nullar1simulations)
  return comparear1undernull(X,y, nullar1statfunc,
    nullar1aggfunc, Nsimulations=nullar1simulations)
end



#####################

#the below function and info is for formatting numbers
function llformat(β;scalefactor, decimals, prefix, suffix)
  valstring = num2str(β, decimals; scalefactor)
  return "$(prefix)$(valstring)$(suffix)"
end

llbetaformat =(; scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:")
llpercentageformat =(; scalefactor=100.0, decimals=2, prefix="", suffix="\\:\\:")
llseformat = (;scalefactor=1.0, decimals=2, prefix="\\mathit{(", suffix=")}")
llNformat = (;scalefactor=1.0, decimals=0, prefix="", suffix="")

#put all of the field information here
llfields = Dict{Symbol, Any}(
  :measure=>Dict(
    :L1g_WXmomentum121s=>(;
      id=:L1g_WXmomentum121s,
      label="\$\\%gMOM_{t-1}\$",
      format=llpercentageformat,
      se=(;
        id=:L1g_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1p_WXmomentum121s=>(;
      id=:L1p_WXmomentum121s,
      label="\$\\%Par_{t-1}\$",
      format=llpercentageformat,
      se=(;
        id=:L1p_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1ag_WXmomentum121s=>(;
      id=:L1ag_WXmomentum121s,
      label="\\%gMOM (active)",
      format=llpercentageformat,
      se=(;
        id=:L1ag_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1cgg_WXmomentum121s=>(;
      id=:L1cgg_WXmomentum121s,
      label="\\%gMOM (returns)",
      format=llpercentageformat,
      se=(;
        id=:L1cgg_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1GMT_WXmomentum121s=>(;
      id=:L1GMT_WXmomentum121s,
      label="\$GMT_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:L1GMT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1nGMT_WXmomentum121s=>(;
      id=:L1nGMT_WXmomentum121s,
      label="\$GMT_{t-1}^{net}\$",
      format=llbetaformat,
      se=(;
        id=:L1nGMT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :dtXL1lMT_WXmomentum121s=>(;
      id=:dtXL1lMT_WXmomentum121s,
      label="dtXL1lMT",
      format=llbetaformat,
      se=(;
        id=:dtXL1lMT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1GAN_WXmomentum121s=>(;
      id=:L1GAN_WXmomentum121s,
      label="GAN",
      format=llbetaformat,
      se=(;
        id=:L1GAN_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1GEON_WXmomentum121s=>(;
      id=:L1GEON_WXmomentum121s,
      label="GEON",
      format=llbetaformat,
      se=(;
        id=:L1GEON_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1GDEON_WXmomentum121s=>(;
      id=:L1GDEON_WXmomentum121s,
      label="GDEON",
      format=llbetaformat,
      se=(;
        id=:L1GDEON_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1M_WXmomentum121s=>(;
      id=:L1M_WXmomentum121s,
      label="\$M^{L+S}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:L1M_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1MT_WXmomentum121s=>(;
      id=:L1MT_WXmomentum121s,
      label="\$M^{tot}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:L1MT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1AN_WXmomentum121s=>(;
      id=:L1AN_WXmomentum121s,
      label="\$MF^{all}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:L1AN_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L1DEON_WXmomentum121s=>(;
      id=:L1DEON_WXmomentum121s,
      label="\$MF^{dom. eq.}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:L1DEON_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :dtXL1M_WXmomentum121s=>(;
      id=:dtXL1M_WXmomentum121s,
      label="\$dtM^{L+S}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:dtXL1M_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :dtXL1MT_WXmomentum121s=>(;
      id=:dtXL1MT_WXmomentum121s,
      label="\$dtM^{tot}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:dtXL1MT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :dtXL1AN_WXmomentum121s=>(;
      id=:dtXL1AN_WXmomentum121s,
      label="\$dtMF^{all}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:dtXL1AN_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :dtXL1DEON_WXmomentum121s=>(;
      id=:dtXL1DEON_WXmomentum121s,
      label="\$dtMF^{dom. eq.}_{t-1}\$",
      format=llbetaformat,
      se=(;
        id=:dtXL1DEON_WXmomentum121s_se,
        label="",
        format=llseformat)),
      :L12M_WXmomentum121s=>(;
        id=:L12M_WXmomentum121s,
        label="\$M^{L+S}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:L12M_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :L12MT_WXmomentum121s=>(;
        id=:L12MT_WXmomentum121s,
        label="\$M^{tot}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:L12MT_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :L12AN_WXmomentum121s=>(;
        id=:L12AN_WXmomentum121s,
        label="\$MF^{all}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:L12AN_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :L12DEON_WXmomentum121s=>(;
        id=:L12DEON_WXmomentum121s,
        label="\$MF^{dom. eq.}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:L12DEON_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :dtXL12M_WXmomentum121s=>(;
        id=:dtXL12M_WXmomentum121s,
        label="\$dtM^{L+S}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:dtXL12M_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :dtXL12MT_WXmomentum121s=>(;
        id=:dtXL12MT_WXmomentum121s,
        label="\$dtM^{tot}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:dtXL12MT_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :dtXL2AN_WXmomentum121s=>(;
        id=:dtXL12AN_WXmomentum121s,
        label="\$dtMF^{all}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:dtXL12AN_WXmomentum121s_se,
          label="",
          format=llseformat)),
      :dtXL12DEON_WXmomentum121s=>(;
        id=:dtXL12DEON_WXmomentum121s,
        label="\$dtMF^{dom. eq.}_{t-12}\$",
        format=llbetaformat,
        se=(;
          id=:dtXL12DEON_WXmomentum121s_se,
          label="",
          format=llseformat)),
    :L1F_ff5_mktrf=>(;
      id=:L1F_ff5_mktrf,
      label="\$\\%Mkt_{t-1}\$",
      format=llpercentageformat,
      se=(;
        id=:L1F_ff5_mktrf_se,
        label="",
        format=llseformat)),
    :L1ng_WXmomentum121s=>(;
      id=:L1ng_WXmomentum121s,
      label="\$\\%gMOM_{t-1}-Mkt_{t-1}\$",
      format=llpercentageformat,
      se=(;
        id=:L1ng_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2GMT_WXmomentum121s=>(;
      id=:L2GMT_WXmomentum121s,
      label="\$GMT_{t-2}\$",
      format=llbetaformat,
      se=(;
        id=:L2GMT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2nGMT_WXmomentum121s=>(;
      id=:L2nGMT_WXmomentum121s,
      label="\$GMT_{t-2}^{net}\$",
      format=llbetaformat,
      se=(;
        id=:L2nGMT_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2p_WXmomentum121s=>(;
      id=:L2p_WXmomentum121s,
      label="\$\\%Par_{t-2}\$",
      format=llpercentageformat,
      se=(;
        id=:L2p_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2g_WXmomentum121s=>(;
      id=:L2g_WXmomentum121s,
      label="\$\\%gMOM_{t-2}\$",
      format=llpercentageformat,
      se=(;
        id=:L2g_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2ng_WXmomentum121s=>(;
      id=:L2ng_WXmomentum121s,
      label="\$\\%gMOM_{t-2}-Mkt_{t-2}\$",
      format=llpercentageformat,
      se=(;
        id=:L2ng_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2GEON_WXmomentum121s=>(;
      id=:L2GEON_WXmomentum121s,
      label="\$\\frac{A_{t-1}-A_{t-2}}{AUM^{eq}_{t-2}}\$",
      format=llbetaformat,
      se=(;
        id=:L2GEON_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L2F_ff5_mktrf=>(;
      id=:L2F_ff5_mktrf,
      label="\$\\%Mkt_{t-2}\$",
      format=llpercentageformat,
      se=(;
        id=:L2F_ff5_mktrf_se,
        label="",
        format=llseformat)),
    :L3g_WXmomentum121s=>(;
      id=:L3g_WXmomentum121s,
      label="\$\\%gMOM_{t-3}\$",
      format=llpercentageformat,
      se=(;
        id=:L3g_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :L3F_ff5_mktrf=>(;
      id=:L3F_ff5_mktrf,
      label="\$Mkt_{t-3}\$",
      format=llpercentageformat,
      se=(;
        id=:L3F_ff5_mktrf_se,
        label="",
        format=llseformat)),

    ######for contemporaneous return relationships
    :R_WXmomentum121s=>(;
      id=:R_WXmomentum121s,
      label="\$mom_{t}\$",
      format=llbetaformat,
      se=(;
        id=:R_WXmomentum121s_se,
        label="",
        format=llseformat)),
    :WXmomentum121s_retlevgross=>(;
      id=:WXmomentum121s_retlevgross,
      label="\$mom_{t}^{L+S}\$",
      format=llbetaformat,
      se=(;
        id=:WXmomentum121s_retlevgross_se,
        label="",
        format=llseformat)),
    :WXmomentum121s_long=>(;
      id=:WXmomentum121s_long,
      label="\$mom_{t}^{L}\$",
      format=llbetaformat,
      se=(;
        id=:WXmomentum121s_long_se,
        label="",
        format=llseformat)),
    :WXmomentum121s_short=>(;
      id=:WXmomentum121s_short,
      label="\$mom_{t}^{S}\$",
      format=llbetaformat,
      se=(;
        id=:WXmomentum121s_short_se,
        label="",
        format=llseformat)),
    :F_ff3m_umd=>(;
      id=:F_ff3m_umd,
      label="\$umd_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff3m_umd_se,
        label="",
        format=llseformat)),
    :F_ff3m_mktrf=>(;
      id=:F_ff3m_mktrf,
      label="\$mktrf_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff3m_mktrf_se,
        label="",
        format=llseformat)),
    :F_ff3m_smb=>(;
      id=:F_ff3m_smb,
      label="\$smb_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff3m_smb_se,
        label="",
        format=llseformat)),
    :F_ff3m_hml=>(;
      id=:F_ff3m_hml,
      label="\$hml_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff3m_hml_se,
        label="",
        format=llseformat)),
    :F_ff5_mktrf=>(;
      id=:F_ff5_mktrf,
      label="\$MktRf_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff5_mktrf_se,
        label="",
        format=llseformat)),
    :F_ff5_smb=>(;
      id=:F_ff5_smb,
      label="\$smb_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff5_smb_se,
        label="",
        format=llseformat)),
    :F_ff5_hml=>(;
      id=:F_ff5_hml,
      label="\$hml_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff5_hml_se,
        label="",
        format=llseformat)),
    :F_ff5_rmw=>(;
      id=:F_ff5_rmw,
      label="\$rmw_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff5_rmw_se,
        label="",
        format=llseformat)),
    :F_ff5_cma=>(;
      id=:F_ff5_cma,
      label="\$cma_{t}\$",
      format=llbetaformat,
      se=(;
        id=:F_ff5_cma_se,
        label="",
        format=llseformat)),
  ),
  :description=>Dict(
    :Nmonths=>(;
      id=:Nmonths,
      label="\\midrule Months",
      format=llNformat),
    :Nweeks=>(;
      id=:Nweeks,
      label="\\midrule Weeks",
      format=llNformat),
  ),
)

#start of a table spec
tablelt112_3672spec = (;
  tablename="lllt112_3672",
  Fmeasures=[
    [:L1MT_WXmomentum121s],
    [:L12MT_WXmomentum121s],
    [:L1MT_WXmomentum121s],
    [:L12MT_WXmomentum121s],
    [:L1MT_WXmomentum121s],
    [:L12MT_WXmomentum121s],
    ],
    #[:L1g_WXmomentum121s, :L1F_ff5_mktrf]],
  colnames=[["12-month", "36-month", "72-month"], 1:6 .|> i->"($i)",],
  widthcolnames= [[2,2,2],ones(Int,6),],
  Freturns=[[:FW12R_WXmomentum121s for i ∈ 1:2]...;[:FW36R_WXmomentum121s for i ∈ 1:2]...;[:FW72R_WXmomentum121s for i ∈ 1:2]...;],
  FN=:Nmonths,
  Fse=:se,
  alignmentstring="l | rr | rr | rr"
)

#start of a table spec
tableltdt112_3672spec = (;
  tablename="llltdt112_3672",
  Fmeasures=[
    [:dtXL1MT_WXmomentum121s],
    [:dtXL12MT_WXmomentum121s],
    [:dtXL1MT_WXmomentum121s],
    [:dtXL12MT_WXmomentum121s],
    [:dtXL1MT_WXmomentum121s],
    [:dtXL12MT_WXmomentum121s],
    ],
    #[:L1g_WXmomentum121s, :L1F_ff5_mktrf]],
  colnames=[["12-month", "36-month", "72-month"], 1:6 .|> i->"($i)",],
  widthcolnames= [[2,2,2],ones(Int,6),],
  Freturns=[[:dtXFW12R_WXmomentum121s for i ∈ 1:2]...;[:dtXFW36R_WXmomentum121s for i ∈ 1:2]...;[:dtXFW72R_WXmomentum121s for i ∈ 1:2]...;],
  FN=:Nmonths,
  Fse=:se,
  alignmentstring="l | rr | rr | rr"
)


#start of a table spec
tablest11spec = (;
  tablename="llst11",
  Fmeasures=[
    [:L1F_ff5_mktrf],
    [:L1g_WXmomentum121s],
    [:L1GMT_WXmomentum121s],
    [:L1nGMT_WXmomentum121s],
    [:L1ng_WXmomentum121s],
    [:L1GMT_WXmomentum121s, :L1F_ff5_mktrf]],
    #[:L1g_WXmomentum121s, :L1F_ff5_mktrf]],
  Freturns=[[:FW1R_WXmomentum121s for i ∈ 1:6]...;],
  FN=:Nweeks,
  Fse=:se,
  alignmentstring="l | rrr | rrr"
)

tablest21spec = (;
  tablename="llst21",
  Fmeasures=[
    [:L2F_ff5_mktrf],
    [:L2GMT_WXmomentum121s],
    [:L2g_WXmomentum121s],
    [:L2nGMT_WXmomentum121s],
    [:L2ng_WXmomentum121s],
    [:L2GMT_WXmomentum121s, :L2F_ff5_mktrf]],
    #[:L2g_WXmomentum121s, :L2F_ff5_mktrf]],
  Freturns=[[:FW1R_WXmomentum121s for i ∈ 1:6]...;],
  #=[[:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];
    [:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];],=#
  #[[:FW1R_WXmomentum121s for i ∈ 1:6]...;],
  FN=:Nweeks,
  Fse=:se,
  alignmentstring="l | rrr | rrr"
)

tablestpspec = (;
  tablename="llstp",
  Fmeasures=[
    [:L1p_WXmomentum121s],
    [:L1p_WXmomentum121s, :L1F_ff5_mktrf],
    [:L1p_WXmomentum121s],
    [:L1p_WXmomentum121s, :L1F_ff5_mktrf],
    [:L2p_WXmomentum121s],
    [:L2p_WXmomentum121s, :L2F_ff5_mktrf],
    [:L2p_WXmomentum121s],
    [:L2p_WXmomentum121s, :L2F_ff5_mktrf]],
  colnames=[["1-week return","2-week return", "1-week return","2-week return"], 1:8 .|> i->"($i)",],
  widthcolnames= [[2,2,2,2],ones(Int,8),],
    #[:L2g_WXmomentum121s, :L2F_ff5_mktrf]],
  Freturns=[[:FW1R_WXmomentum121s for i ∈ 1:2]...;[:FW2R_WXmomentum121s for i ∈ 1:2]...;
    [:FW1R_WXmomentum121s for i ∈ 1:2]...;[:FW2R_WXmomentum121s for i ∈ 1:2]...;],
  #=[[:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];
    [:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];],=#
  #[[:FW1R_WXmomentum121s for i ∈ 1:6]...;],
  FN=:Nweeks,
  Fse=:se,
  alignmentstring="l | rr | rr | rr | rr"
)


table00spec = (;
  tablename="ll00",
  Fmeasures=[
    [:R_WXmomentum121s],
    [:WXmomentum121s_long],
    [:WXmomentum121s_short],
    #[:WXmomentum121s_retlevgross], #maybe consider breaking this all up into two tables
    [:F_ff3m_umd],
    [:F_ff3m_mktrf],
    [:F_ff3m_smb],
    [:F_ff3m_hml],],
    #[:F_ff3m_umd, :F_ff3m_mktrf, :F_ff3m_smb, :F_ff3m_hml, ]],
  Freturns=[[:g_WXmomentum121s for i ∈ 1:7]...;],
  #=[[:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];
    [:FW1R_WXmomentum121s; :FW2R_WXmomentum121s; :FW3R_WXmomentum121s];],=#
  #[[:FW1R_WXmomentum121s for i ∈ 1:6]...;],
  FN=:Nweeks,
  Fse=:se,
  alignmentstring="l | rrr | rrrr"
)


function lltables(;analysispath = PARAM[:analysispath],)

  bssamples = PARAM[:tabllbsar1samples]
  nullar1simulations = PARAM[:tabllbsar1samples]

  measure00 = CSV.File(
    "$analysispath\\$(PARAM[:tabll00measurefilename])" *
    "\\$(PARAM[:tabll00fundrid])\\$(PARAM[:tabll00fundrid])_leadlag_series.csv"
    ) |> DataFrame
  ll00bswidth = nrow(measure00)^(1/3) |> ceil |> Int
  lltable(measure00, (X,y; β)->regmbb(X,y, bswidth=ll00bswidth; bssamples).beta.se;
    table00spec...)

  measurelt112_3672 = CSV.File(
    "$analysispath\\$(PARAM[:tablllt112_3672measurefilename])" *
    "\\$(PARAM[:tablllt112_3672fundrid])\\$(PARAM[:tablllt112_3672fundrid])_leadlag_series.csv"
    ) |> DataFrame
  lt112_3672bswidth = nrow(measurelt112_3672)^(1/4) |> ceil |> Int
  lltable(measurelt112_3672, (X,y; β)->regnullar1(X,y;β,nullar1simulations,integrated=true).beta.se;
    tablelt112_3672spec...)

  measureltdt112_3672 = CSV.File(
    "$analysispath\\$(PARAM[:tabllltdt112_3672measurefilename])" *
    "\\$(PARAM[:tabllltdt112_3672fundrid])\\$(PARAM[:tabllltdt112_3672fundrid])_leadlag_series.csv"
    ) |> DataFrame
  ltdt112_3672bswidth = nrow(measureltdt112_3672)^(1/4) |> ceil |> Int
  lltable(measureltdt112_3672, (X,y; β)->regnullar1(X,y;β,nullar1simulations,integrated=false).beta.se;
    tableltdt112_3672spec...)


  measurest11 = CSV.File(
    "$analysispath\\$(PARAM[:tabllst11measurefilename])" *
    "\\$(PARAM[:tabllst11fundrid])\\$(PARAM[:tabllst11fundrid])_leadlag_series.csv"
    ) |> DataFrame
  st11bswidth = nrow(measurest11)^(1/3) |> ceil |> Int
  lltable(measurest11, (X,y; β)->regmbb(X,y, bswidth=st11bswidth; bssamples).beta.se;
    tablest11spec...)

  measurest21 = CSV.File(
    "$analysispath\\$(PARAM[:tabllst21measurefilename])" *
    "\\$(PARAM[:tabllst21fundrid])\\$(PARAM[:tabllst21fundrid])_leadlag_series.csv"
    ) |> DataFrame
  st21bswidth = nrow(measurest21)^(1/3) |> ceil |> Int
  lltable(measurest21, (X,y; β)->regmbb(X,y, bswidth=st21bswidth; bssamples).beta.se;
    tablest21spec...)

  measurestp = CSV.File(
    "$analysispath\\$(PARAM[:tabllstpmeasurefilename])" *
    "\\$(PARAM[:tabllstpfundrid])\\$(PARAM[:tabllstpfundrid])_leadlag_series.csv"
    ) |> DataFrame
  stpbswidth = nrow(measurestp)^(1/3) |> ceil |> Int
  lltable(measurestp, (X,y; β)->regmbb(X,y, bswidth=stpbswidth; bssamples).beta.se;
    tablestpspec...)


end

#note the se function can be anything that takes X and y, so it cna return multpile stats
#this is not meant to be high performance
function llregression(measure,se::Tse; FX,Fy) where Tse<:Function
  #construct the info for each row

  cleanmeasure = view(measure, completecases(measure, [FX; Fy]),:)
  #form X from the measure DF + an intercept column
  X = hcat(cleanmeasure[:, FX] |> Matrix{Float64}, ones(nrow(cleanmeasure)))
  y = cleanmeasure[:, Fy] |> Vector{Float64}

  #compute the estimate and standard error
  β = svd(X)\y
  @assert (FMLM(X,y).β .≈ β) |> all

  se = se(X, y; β)

  N = nrow(cleanmeasure)

  t = β ./ se

  return (; FX, β, t, se, N)
end

#create a table for the lead/lag
function lltable(measure, se::Tse;
    tablename,
    Fmeasures,
    Freturns::TFreturns,
    FN,
    Fse,
    tablepath = PARAM[:tablepath],
    colnames=[1:length(Fmeasures) .|> i->"($i)",],
    widthcolnames= [ones(Int,length(Fmeasures)),],
    alignmentcolnames = [
      (rc-> "c").(colnamerow) for (r,colnamerow) ∈ enumerate(colnames)],
    alignmentstring = string("l",["r" for i ∈ 1:length(Fmeasures)]...),

    #construct the info for each row
    #first extract lightly constructed field objects IN ORDER
    datafields = reduce(vcat, (reduce(vcat, Fmeasures) |> unique) .|> (
      (s)->[llfields[:measure][s], llfields[:measure][s][Fse]])),
    descriptionfields = [
      llfields[:description][FN]
    ],

    ) where {TFreturns, Tse<:Function}


  #get the row ids IN ORDER
  Frowids = (f->f.id).([datafields; descriptionfields])

  #create an index of row information. order is not preserved.
  Frowindex = Dict(f.id=>f for f ∈ [datafields; descriptionfields])

  #double check alignment
  @assert (id->Frowindex[id].id≡id).(Frowids) |> all

  #extract the row labels, using the order preserved in Frowid
  descrownames = (id->Frowindex[id].label).(Frowids)

  #do this in case I want to use a different parsing routine etc
  RHS = Fmeasures
  lhs::Vector{Symbol} = TFreturns <: Symbol ? [Freturns for r ∈ 1:size(RHS,1)] : Freturns

  @assert length(lhs) ≡ size(RHS,1)


  #NOTE: each ROW in the dataframe corresponds to a COLUMN in the table
  #cols preserve order, though I will try not to rely on it
  regressioncols = DataFrame([id=>String[] for id ∈ Frowids])
  for (FX, Fy) ∈ zip(RHS, lhs)
    reg = llregression(measure, se; FX, Fy)

    regressioncol = Dict{Symbol, Any}()
    for (i,Fx) ∈ enumerate(FX)

      #pull the field information from the GLOBAL VARIABLE llfields
      estimatefield=llfields[:measure][Fx]
      @assert estimatefield.id ≡ Fx #double check the field is correct
      regressioncol[estimatefield.id] = llformat(reg.β[i]; estimatefield.format...)

      sefield = llfields[:measure][Fx][Fse]
      @assert sefield.id ≡ Symbol(Fx, :_, Fse)
      regressioncol[sefield.id] = llformat(reg.t[i]; sefield.format...)
    end

    #compute N, and eventually any other descriptive fields
    Nfield= llfields[:description][FN]
    @assert Nfield.id ≡ FN
    regressioncol[Nfield.id] = llformat(reg.N; Nfield.format...)

    #make sure we captured the expected fields
    #note- below is UNORDERED
    @assert issetequal(keys(regressioncol)|>collect,
      [FX; FX .|> (s)->llfields[:measure][s][Fse].id; Nfield.id]) "
        Unexpected and/or missing fields detected in lltable!:
        keys(regressioncol)|>collect: $(keys(regressioncol)|>collect)
        [...]: $([FX; FX .|> (s)->llfields[:measure][s][Fse].id; Nfield])
        "

    #all other columns and rows are blank
    blankfields = Dict(f => "" for f ∈ setdiff(propertynames(regressioncols), keys(regressioncol)))

    #form the regression col
    push!(regressioncols, merge(regressioncol, blankfields))
  end

  Ncols = length(Fmeasures)
  #desccontent = [Vector{String}(undef, Ncols) for i ∈ 1:Nrows]

  @assert (propertynames(regressioncols) .≡ Frowids) |> all
  desccontent = eachcol(regressioncols) |> collect .|> Vector

  t = textable(;colnames, descrownames, desccontent,
    widthcolnames, alignmentstring, alignmentcolnames)

  writetextable(t, path=tablepath, outname="$tablename.tex")
end
