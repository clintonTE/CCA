


#runs the basic panel form regressions
function olscontributionreg(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Lpcontributions"),
    xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
  xspec=Meta.parse("Lpprogramexpense"),
  xnames = [:intercept, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  push!(specs,
  xspec=Meta.parse("Lpprogramexpense"),
  xnames = [:intercept, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end

#runs the basic panel form regressions
function olscontributionreghist(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("pcontributions"),
    xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn3yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpense"),
    xnames = [:intercept, :Exp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn3yr)

  push!(specs,
  xspec=Meta.parse("pcontributions"),
  xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
  xspec=Meta.parse("pprogramexpense"),
  xnames = [:intercept, :Exp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
  xspec=Meta.parse("pcontributions"),
  xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  push!(specs,
  xspec=Meta.parse("pprogramexpense"),
  xnames = [:intercept, :Exp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end

#runs the basic panel form regressions
function olscontributionregcomb(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("pcontributions"),
    xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn3yr)

  push!(specs,
    xspec=Meta.parse("Lpcontributions"),
    xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
  xspec=Meta.parse("pcontributions"),
  xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
  xspec=Meta.parse("pcontributions"),
  xnames = [:intercept, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end


#runs the basic panel form regressions
function olssensitivitypanelregPE(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end



#runs the basic panel form regressions
function olssensitivitypanelregPEhist(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXpcont + pcontributions"),
    xnames = [:intercept, :ContribAvgExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn3yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpenseXpcont + pcontributions + pprogramexpense"),
    xnames = [:intercept, :ContribExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn3yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXpcont + pcontributions"),
    xnames = [:intercept, :ContribAvgExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpenseXpcont + pcontributions + pprogramexpense"),
    xnames = [:intercept, :ContribExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5003yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXpcont + pcontributions"),
    xnames = [:intercept, :ContribAvgExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  push!(specs,
    xspec=Meta.parse("pprogramexpenseXpcont + pcontributions + pprogramexpense"),
    xnames = [:intercept, :ContribExp, :Contrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff33yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end

#runs the basic panel form regressions
function olssensitivitypanelregPS(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  #NOTE: delete the l before public support on next pass
  push!(specs,
    xspec=Meta.parse("publicsupport_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    withinspec = :ein,
    clusterspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLpcont + Lpcontributions + Lpublicsupport"),
    xnames = [:intercept, :LagContribSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("publicsupport_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLpcont + Lpcontributions + Lpublicsupport"),
    xnames = [:intercept, :LagContribSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("publicsupport_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLpcont + Lpcontributions + Lpublicsupport"),
    xnames = [:intercept, :LagContribSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end

#runs the basic panel form regressions
function olssensitivity2pass(df::DataFrame;
  rhs2ndpassfields_avg::Vector{Symbol} = [:pprogramexpense_avg],#, :publicsupport_avg=#],
  rhs2ndpassfields::Vector{Symbol} = [:Lpprogramexpense#=, :Lpublicsupport_avg=#],
  winsorizeprop::Float64 = WINSORIZE_PROP)

  specs::FMSpecs = FMSpecs(6)

  local ??df::DataFrame
  local ??dfcross::DataFrame
  local focalfields::Vector{Symbol}
  #run the first pass regressions to get the betas
  #println("q expense complete: ",
  #  sum(completecases(df[[:lreturn, :Llcontributions, :lprogramexpense_avg]])))

  regressbynonprofit!(df,
    YField = :lreturn,
    XFields = [:Lpcontributions],
    ??Names = [:b0temp1, :b1lRetLCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmsp5001yr,
    XFields = [:Lpcontributions],
    ??Names = [:b0temp2, :b1lRetNetSP500LCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmff31yr,
    XFields = [:Lpcontributions],
    ??Names = [:b0temp3, :b1lRetNetFF3LCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  deletecols!(df, [:b0temp1, :b0temp2, :b0temp3])

  ??fields::Vector{Symbol} = [:b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

  (f->winsorize!(df[f], prop=winsorizeprop)).(??fields)
  focalfields = [:ein; ??fields]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  ??df = unique(df[[focalfields; :categoryfiltered]], focalfields)
  for f ??? rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    ??df = join(??df, unique(subdf), on=:ein, kind=:left)
  end

  plot??tab!(??df, tabfield = :b1lRetLCon)
  #CSV.write("$WORKING_PATH\\betatest.csv", ??df)

  #make sure the joins went as planned
  @assert size(unique(??df),1) == size(??df,1)
  focalfields = [:ein; :fisyr; ??fields; :categoryfiltered]

  ??dfcross = df[focalfields]

  for f ??? rhs2ndpassfields
    subdf::DataFrame = view(df, completecases(df[:,[:ein,:fisyr,f]]), [:ein,:fisyr,f])
    ??dfcross = join(??dfcross, unique(subdf), on=[:ein, :fisyr], kind=:left)
  end

  #make sure all of the betas within each entity are unique
  @assert size(unique(??dfcross),1) == size(??dfcross,1)

  #make the fisyrs into symbols
  rename!(??dfcross, [:fisyr=>:fisyrold])
  ??dfcross[:fisyr] = (y::Int->Symbol("y", y)).(??dfcross[:fisyrold])
  deletecols!(??dfcross, :fisyrold)

  #Regress the betas on the focal variables
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    withinspec = :categoryfiltered,
    yspec=:b1lRetLCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :categoryfiltered,
    yspec=:b1lRetLCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetSP500LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetSP500LCon)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    #withinspec = :fisyr,
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  println(describe(??df))
  println(describe(??dfcross))

  #run the regressions
  ??dfs::Vector{DataFrame} = [??df, ??dfcross, ??df, ??dfcross, ??df, ??dfcross]
  computeFMLMresults!(??dfs, specs)

  return specs
end

#runs the basic panel form regressions
function olssensitivity2passhist(df::DataFrame;
  rhs2ndpassfields_avg::Vector{Symbol} = [:pprogramexpense_avg],#, :publicsupport_avg=#],
  rhs2ndpassfields::Vector{Symbol} = [:pprogramexpense#=, :Lpublicsupport_avg=#],
  winsorizeprop::Float64 = WINSORIZE_PROP)

  specs::FMSpecs = FMSpecs(6)
  df = prefilterforolssensitivity(deepcopy(df), returnfield = :lreturn3yr)

  local ??df::DataFrame
  local ??dfcross::DataFrame
  local focalfields::Vector{Symbol}
  #run the first pass regressions to get the betas
  #println("q expense complete: ",
  #  sum(completecases(df[[:lreturn, :Llcontributions, :lprogramexpense_avg]])))

  regressbynonprofit!(df,
    YField = :lreturn3yr,
    XFields = [:pcontributions],
    ??Names = [:b0temp1, :b1lRet3yrCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmsp5003yr,
    XFields = [:pcontributions],
    ??Names = [:b0temp2, :b1lRet3yrNetSP500Con],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmff33yr,
    XFields = [:pcontributions],
    ??Names = [:b0temp3, :b1lRet3yrNetFF3Con],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  deletecols!(df, [:b0temp1, :b0temp2, :b0temp3])

  ??fields::Vector{Symbol} = [:b1lRet3yrCon, :b1lRet3yrNetSP500Con,
    :b1lRet3yrNetFF3Con]

  (f->winsorize!(df[f], prop=winsorizeprop)).(??fields)
  focalfields = [:ein; ??fields]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  ??df = unique(df[focalfields])
  for f ??? rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    ??df = join(??df, unique(subdf), on=:ein, kind=:left)
  end

  #CSV.write("$WORKING_PATH\\betatest.csv", ??df)

  #make sure the joins went as planned
  @assert size(unique(??df),1) == size(??df,1)
  focalfields = [:ein; :fisyr; ??fields]

  ??dfcross = df[focalfields]

  for f ??? rhs2ndpassfields
    subdf::DataFrame = view(df, completecases(df[:,[:ein,:fisyr,f]]), [:ein,:fisyr,f])
    ??dfcross = join(??dfcross, unique(subdf), on=[:ein, :fisyr], kind=:left)
  end

  #make sure all of the betas within each entity are unique
  @assert size(unique(??dfcross),1) == size(??dfcross,1)

  #make the fisyrs into symbols
  rename!(??dfcross, [:fisyr=>:fisyrold])
  ??dfcross[:fisyr] = (y::Int->Symbol("y", y)).(??dfcross[:fisyrold])
  deletecols!(??dfcross, :fisyrold)

  plot??tab!(??df)

  #Regress the betas on the focal variables
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRet3yrCon)

  push!(specs,
    xspec=Meta.parse("pprogramexpense"),
    xnames = [:intercept, :Exp],
    withinspec = :fisyr,
    yspec=:b1lRet3yrCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRet3yrNetSP500Con)

  push!(specs,
    xspec=Meta.parse("pprogramexpense"),
    xnames = [:intercept, :Exp],
    withinspec = :fisyr,
    yspec=:b1lRet3yrNetSP500Con)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRet3yrNetFF3Con)

  push!(specs,
    xspec=Meta.parse("pprogramexpense"),
    xnames = [:intercept, :Exp],
    withinspec = :fisyr,
    yspec=:b1lRet3yrNetFF3Con)

  println(describe(??df))
  println(describe(??dfcross))

  #run the regressions
  ??dfs::Vector{DataFrame} = [??df, ??dfcross, ??df, ??dfcross, ??df, ??dfcross]
  computeFMLMresults!(??dfs, specs)

  return specs
end
