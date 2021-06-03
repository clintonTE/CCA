

#runs the basic panel form regressions
function olsretcontributionL13(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Llreturnwin"),
    xnames = [:intercept, :LagReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)
  push!(specs,
    xspec=Meta.parse("Llreturn3yrwin"),
    xnames = [:intercept, :LagReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)
  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5001yrwin"),
    xnames = [:intercept, :LagReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5003yrwin"),
    xnames = [:intercept, :LagReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff31yrwin"),
    xnames = [:intercept, :LagReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff33yrwin"),
    xnames = [:intercept, :LagReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end



#runs the basic panel form regressions
function olscontributionregretcon(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnwin)

  push!(specs,
    xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
    xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnwin)

  push!(specs,
  xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
  xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yrwin)

  push!(specs,
  xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
  xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmsp5001yrwin)

  push!(specs,
  xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
  xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yrwin)

  push!(specs,
  xspec=Meta.parse("Lpcontributions + Lpprogramexpense"),
  xnames = [:intercept, :LagContrib, :LagExp],
    clusterspec = :categoryfiltered,
    withinspec = :categoryfiltered,
    yspec=:lreturnnetbmff31yrwin)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end

#runs the basic panel form regressions
function olscontributionregconret(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Llreturnwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)
  push!(specs,
    xspec=Meta.parse("Llreturn3yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yr, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)
  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5001yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5003yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yr, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff31yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff33yrwin + Lpprogramexpense"),
    xnames = [:intercept, :LagReturn3yr, :LagExp],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end



function olsretcontributionLL13(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("Llreturnwin + LLlreturnwin"),
    xnames = [:intercept, :LagReturn, :Lag2xReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturn3yrwin + L3Llreturn3yrwin"),
    xnames = [:intercept, :LagReturn3yr, :Lag2xReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5001yrwin + LLlreturnnetbmsp5001yrwin"),
    xnames = [:intercept, :LagReturn, :Lag2xReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmsp5003yrwin + L3Llreturnnetbmsp5003yrwin"),
    xnames = [:intercept, :LagReturn3yr, :Lag2xReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff31yrwin + LLlreturnnetbmff31yrwin"),
    xnames = [:intercept, :LagReturn, :Lag2xReturn],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  push!(specs,
    xspec=Meta.parse("Llreturnnetbmff33yrwin + L3Llreturnnetbmff33yrwin"),
    xnames = [:intercept, :LagReturn3yr, :Lag2xReturn3yr],
    yspec=:pcontributions,
    clusterspec = :ein,
    withinspec = :ein)

  #run the regressions
  computeFMLMresults!(df, specs)

  return specs
end
