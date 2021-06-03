



#runs the basic panel form regressions
function olscontributionregcomb(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Llreturn3yrwin"),
    xnames = [:intercept, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse("Lpcontributions"),
    xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnwin)

  push!(specs,
  xspec=Meta.parse("Llreturnnetbmsp5003yrwin"),
  xnames = [:intercept, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)

  push!(specs,
  xspec=Meta.parse("Llreturnnetbmff33yrwin"),
  xnames = [:intercept, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
  xspec=Meta.parse("Lpcontributions"),
  xnames = [:intercept, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)

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
    yspec=:lreturnwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnwin)=#

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnwin)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)=#

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)=#

  push!(specs,
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end



#runs the basic panel form regressions
function olssensitivitypanelregPEhist(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse(
      "pprogramexpense_avgXLlreturn3yrwin + Llreturn3yrwin"),
    xnames = [:intercept, :LagReturn3yrAvgExp, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse(
      "LpprogramexpenseXLlreturn3yrwin +Llreturn3yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yrExp, :LagReturn3yr, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse(
      "pprogramexpense_avgXLlreturnnetbmsp5003yrwin + Llreturnnetbmsp5003yrwin"),
    xnames = [:intercept, :LagReturn3yrAvgExp, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse(
      "LpprogramexpenseXLlreturnnetbmsp5003yrwin + Llreturnnetbmsp5003yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yrExp, :LagReturn3yr, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse(
      "pprogramexpense_avgXLlreturnnetbmff33yrwin + Llreturnnetbmff33yrwin"),
    xnames = [:intercept, :LagReturn3yrAvgExp, :LagReturn3yr],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  push!(specs,
    xspec=Meta.parse(
      "LpprogramexpenseXLlreturnnetbmff33yrwin + Llreturnnetbmff33yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yrExp, :LagReturn3yr, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:pcontributions)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end


#runs the basic panel form regressions
function olssensitivitypanelregPEcat(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnwin)=#

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnwin)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmsp5001yrwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)=#

  push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmsp5001yrwin)

  push!(specs,
    xspec=Meta.parse("pprogramexpense_avgXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmff31yrwin)

  #=push!(specs,
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)=#

  push!(specs,
    #xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions"),
    xspec=Meta.parse("LpprogramexpenseXLpcont + Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmff31yrwin)

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
  df = prefilterforolssensitivity(deepcopy(df), returnfield = :lreturnwin,
    essentialfields=[:Lpcontributions, :Lpprogramexpense])

  regressbynonprofit!(df,
    YField = :lreturnwin,
    XFields = [:Lpcontributions, :Lpprogramexpense],
    βNames = [:b0temp1, :b1lRetLCon, :b2temp1],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmsp5001yrwin,
    XFields = [:Lpcontributions, :Lpprogramexpense],
    βNames = [:b0temp2, :b1lRetNetSP500LCon, :b2temp2],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :lreturnnetbmff31yrwin,
    XFields = [:Lpcontributions, :Lpprogramexpense],
    βNames = [:b0temp3, :b1lRetNetFF3LCon, :b2temp3],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  deletecols!(df, [:b0temp1, :b0temp2, :b0temp3])

  βfields::Vector{Symbol} = [:b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

  (f->winsorize!(df[f], prop=winsorizeprop)).(βfields)
  focalfields = [:ein; βfields]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  βdf = unique(df[[focalfields; :categoryfiltered]], focalfields)
  for f ∈ rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    βdf = join(βdf, unique(subdf), on=:ein, kind=:left)
  end

  plotβtab!(βdf, tabfield = :b1lRetLCon)
  #CSV.write("$WORKING_PATH\\betatest.csv", βdf)

  #make sure the joins went as planned
  @assert size(unique(βdf),1) == size(βdf,1)
  focalfields = [:ein; :fisyr; βfields; :categoryfiltered]

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
    withinspec = :categoryfiltered,
    clusterspec = :categoryfiltered,
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
    clusterspec = :categoryfiltered,
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
    clusterspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    #withinspec = :fisyr,
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  println(describe(βdf))
  println(describe(βdfcross))

  #run the regressions
  βdfs::Vector{DataFrame} = [βdf, βdfcross, βdf, βdfcross, βdf, βdfcross]
  computeFMLMresults!(βdfs, specs)

  return specs
end

#runs the basic panel form regressions
function olssensitivity2passhist(df::DataFrame;
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
  df = prefilterforolssensitivity(deepcopy(df), returnfield = :L3Llreturn3yrwin,
    essentialfields=[:Lpcontributions, :Lpprogramexpense, :L3Llreturn3yrwin,
      :L3Llreturnnetbmff33yrwin, :L3Llreturnnetbmsp5003yrwin])

  regressbynonprofit!(df,
    YField = :Lpcontributions,
    XFields = [:L3Llreturn3yrwin, :Lpprogramexpense],
    βNames = [:b0temp1, :b1lRetLCon, :b2temp1],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :Lpcontributions,
    XFields = [:L3Llreturnnetbmsp5003yrwin, :Lpprogramexpense],
    βNames = [:b0temp2, :b1lRetNetSP500LCon, :b2temp2],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  regressbynonprofit!(df,
    YField = :Lpcontributions,
    XFields = [:L3Llreturnnetbmff33yrwin, :Lpprogramexpense],
    βNames = [:b0temp3, :b1lRetNetFF3LCon, :b2temp3],
    groupbyfield = :ein,
    predictwithresults = false,
    checkunique=false)

  deletecols!(df, [:b0temp1, :b0temp2, :b0temp3])

  βfields::Vector{Symbol} = [:b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

  (f->winsorize!(df[f], prop=winsorizeprop)).(βfields)
  focalfields = [:ein; βfields]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  βdf = unique(df[[focalfields; :categoryfiltered]], focalfields)
  for f ∈ rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    βdf = join(βdf, unique(subdf), on=:ein, kind=:left)
  end

  plotβtab!(βdf, tabfield = :b1lRetLCon)
  #CSV.write("$WORKING_PATH\\betatest.csv", βdf)

  #make sure the joins went as planned
  @assert size(unique(βdf),1) == size(βdf,1)
  focalfields = [:ein; :fisyr; βfields; :categoryfiltered]

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
    withinspec = :categoryfiltered,
    clusterspec = :categoryfiltered,
    yspec=:b1lRetLCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    withinspec = :categoryfiltered, #no clustered errors due to Newey-West SEs
    yspec=:b1lRetLCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    withinspec = :categoryfiltered,
    clusterspec = :categoryfiltered,
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
    clusterspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    #withinspec = :fisyr,
    withinspec = :categoryfiltered,
    yspec=:b1lRetNetFF3LCon)

  println(describe(βdf))
  println(describe(βdfcross))

  #run the regressions
  βdfs::Vector{DataFrame} = [βdf, βdfcross, βdf, βdfcross, βdf, βdfcross]
  computeFMLMresults!(βdfs, specs)

  return specs
end
