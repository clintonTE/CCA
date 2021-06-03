


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

  local βdf::DataFrame
  local βdfcross::DataFrame
  local focalfields::Vector{Symbol}
  #run the first pass regressions to get the betas
  #println("q expense complete: ",
  #  sum(completecases(df[[:lreturn, :Llcontributions, :lprogramexpense_avg]])))

  regressbynonprofit!(df,
    YField = :lreturn,
    XFields = [:Lpcontributions],
    βNames = [:b0temp1, :b1lRetLCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmsp5001yr,
    XFields = [:Lpcontributions],
    βNames = [:b0temp2, :b1lRetNetSP500LCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmff31yr,
    XFields = [:Lpcontributions],
    βNames = [:b0temp3, :b1lRetNetFF3LCon],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  deletecols!(df, [:b0temp1, :b0temp2, :b0temp3])

  βfields::Vector{Symbol} = [:b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

  (f->winsorize!(df[f], prop=winsorizeprop)).(βfields)
  focalfields = [:ein; βfields]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  βdf = unique(df[focalfields])
  for f ∈ rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    βdf = join(βdf, unique(subdf), on=:ein, kind=:left)
  end

  CSV.write("$WORKING_PATH\\betatest.csv", βdf)

  #make sure the joins went as planned
  @assert size(unique(βdf),1) == size(βdf,1)
  focalfields = [:ein; :fisyr; βfields]

  βdfcross = df[focalfields]

  for f ∈ rhs2ndpassfields
    subdf::DataFrame = view(df, completecases(df[:,[:ein,:fisyr,f]]), [:ein,:fisyr,f])
    βdfcross = join(βdfcross, unique(subdf), on=[:ein, :fisyr], kind=:left)
  end

  #make sure all of the betas within each entity are unique
  @assert size(unique(βdfcross),1) == size(βdfcross,1)

  #make the fisyrs into symbols
  rename!(βdfcross, [:fisyr=>:fisyrold])
  βdfcross[:fisyr] = (y::Int->Symbol("y", y)).(βdfcross[:fisyrold])
  deletecols!(βdfcross, :fisyrold)

  #Regress the betas on the focal variables
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRetLCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :fisyr,
    yspec=:b1lRetLCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRetNetSP500LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :fisyr,
    yspec=:b1lRetNetSP500LCon)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRetNetFF3LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :fisyr,
    yspec=:b1lRetNetFF3LCon)

  println(describe(βdf))
  println(describe(βdfcross))

  #run the regressions
  βdfs::Vector{DataFrame} = [βdf, βdfcross, βdf, βdfcross, βdf, βdfcross]
  computeFMLMresults!(βdfs, specs)

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
