

#does any data processign that is needed
#generally called from ComputeNCCS since it is cacheable

function prepsensitivityregressions!(data::NCCSData;
    toosmalltolog::Float64 = TOO_SMALL_TO_LOG)

  #logs if possible or returns missing
  plog(f::MFloat64)::MFloat64 = ismissing(f) || toosmalltolog ≤ TOO_SMALL_TO_LOG ? missing : log(f)

  #cleans up a common-size ratio
  cleanratio(f::MFloat64) = (!ismissing(f)) && isfinite(f) && f ≤ 1.0  && f > 0.0 ? f : missing

  #log our focal variables
  data.df[:lcontributions] = (plog).(data.df[:contributions])
  data.df[:lprogramexpense] = (plog).(data.df[:programexpense])
  data.df[:lprogramexpense_avg] = (plog).(data.df[:programexpense_avg])

  #taking the average for control purposes only with valid contribution fields
  #data.df[:lregnetassets] = ((c::MFloat64, a::MFloat64)->
  #  ismissing(c) ? missing : a).(data.df[:contributions], data.df[:lnetassets])
  lagwithin!(data, :netassets)

  #scale everythign to lagged net assets
  data.df[:pcontributions] = data.df[:contributions] ./ data.df[:netassets]
  data.df[:pcontributions] .= (cleanratio).(data.df[:pcontributions])
  #data.df[:lpcontributions] .= (plog).(data.df[:lpcontributions])

  data.df[:pprogramexpense] = data.df[:programexpense] ./ data.df[:netassets]
  data.df[:pprogramexpense] .= (cleanratio).(data.df[:pprogramexpense])
  #data.df[:lpprogramexpense] .= (plog).(data.df[:pprogramexpense])

  #get the average expense
  aggdf::DataFrame = aggregatefielddf(data.df[[:ein, :pprogramexpense
    ]], rejoin=false)
  unique!(aggdf)
  data.df = join(data.df, aggdf, on=[:ein], kind=:left)

  #lag our focal variables
  lagwithin!(data, :lcontributions)
  lagwithin!(data, :lpcontributions)
  lagwithin!(data, :publicsupport)
  lagwithin!(data, :lprogramexpense)
  lagwithin!(data, :pprogramexpense)
  #lagwithin!(data, :lpprogramexpense)


  #println(describe(data.df[[:ein, :pprogramexpense, :pprogramexpense_avg, :Lpprogramexpense]]))

  #generate the interactions
  data.df[:LpublicsupportXLlcont] = data.df[:Lpublicsupport] .* data.df[:Llcontributions]
  data.df[:publicsupport_avgXLlcont] = data.df[:publicsupport_avg] .* data.df[:Llcontributions]

  data.df[:LlprogramexpenseXLlcont] = data.df[:Llprogramexpense] .* data.df[:Llcontributions]
  data.df[:lprogramexpense_avgXLlcont] = data.df[:lprogramexpense_avg] .* data.df[:Llcontributions]

  data.df[:LpprogramexpenseXLlcont] = data.df[:Lpprogramexpense] .* data.df[:Llcontributions]
  data.df[:pprogramexpense_avgXLlcont] = data.df[:pprogramexpense_avg] .* data.df[:Llcontributions]

  #generate the first differences
  #differencewithin!(data, :lcontributions)
  #differencewithin!(data, :publicsupport)
  #differencewithin!(data, :lprogramexpense)
  #differencewithin!(data, :lprogramexpense)


  #I think regressions work better with symbols instead of strings
  rename!(data.df, :ein=>:ein_old)
  data.df[:ein] = (Symbol).(data.df[:ein_old])
  deletecols!(data.df, :ein_old)

  return nothing
end


#runs the basic panel form regressions
function olssensitivitypanelregPE(df::AbstractDataFrame)
  specs::FMSpecs = FMSpecs(6)

  #set up the results
  push!(specs,
    xspec=Meta.parse("lprogramexpense_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("LlprogramexpenseXLlcont + Llcontributions + Llprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("lprogramexpense_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("LlprogramexpenseXLlcont + Llcontributions + Llprogramexpense"),
    xnames = [:intercept, :LagContribExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("lprogramexpense_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgExp, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  push!(specs,
    xspec=Meta.parse("LlprogramexpenseXLlcont + Llcontributions + Llprogramexpense"),
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
    xspec=Meta.parse("publicsupport_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    withinspec = :ein,
    clusterspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLlcont + Llcontributions + Lpublicsupport"),
    xnames = [:intercept, :LagContribSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturn)

  push!(specs,
    xspec=Meta.parse("publicsupport_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLlcont + Llcontributions + Lpublicsupport"),
    xnames = [:intercept, :LagContribSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmsp5001yr)

  push!(specs,
    xspec=Meta.parse("publicsupport_avgXLlcont + Llcontributions"),
    xnames = [:intercept, :LagContribAvgSup, :LagContrib],
    clusterspec = :ein,
    withinspec = :ein,
    yspec=:lreturnnetbmff31yr)

  push!(specs,
    xspec=Meta.parse("LpublicsupportXLlcont + Llcontributions + Lpublicsupport"),
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
  rhs2ndpassfields::Vector{Symbol} = [:Lpprogramexpense#=, :Lpublicsupport_avg=#])

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

  focalfields = [:ein, :b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

  #now fold in the 2nd pass regression terms through repeated inner joins and uniqueness filters
  βdf = unique(df[focalfields])
  for f ∈ rhs2ndpassfields_avg
    subdf::DataFrame = view(df, completecases(df[:,[:ein,f]]), [:ein, f])
    βdf = join(βdf, unique(subdf), on=:ein, kind=:left)
  end

  CSV.write("$WORKING_PATH\\betatest.csv", βdf)

  #make sure the joins went as planned
  @assert size(unique(βdf),1) == size(βdf,1)
  focalfields = [:ein, :fisyr, :b1lRetLCon, :b1lRetNetSP500LCon,
    :b1lRetNetFF3LCon]

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
    #withinspec = :fisyr,
    yspec=:b1lRetLCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRetNetSP500LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    #withinspec = :fisyr,
    yspec=:b1lRetNetSP500LCon)


  push!(specs,
    xspec=Meta.parse("pprogramexpense_avg"),
    xnames = [:intercept, :AvgExp],
    yspec=:b1lRetNetFF3LCon)

  push!(specs,
    xspec=Meta.parse("Lpprogramexpense"),
    xnames = [:intercept, :LagExp],
    #withinspec = :fisyr,
    yspec=:b1lRetNetFF3LCon)


  #run the regressions
  βdfs::Vector{DataFrame} = [βdf, βdfcross, βdf, βdfcross, βdf, βdfcross]
  computeFMLMresults!(βdfs, specs)

  return specs
end

#primarily generates the gregeression table
function olssensitivitypanel(df::AbstractDataFrame;
  decimals::Int=DECIMALS,
  suffix::String = "",
  displayedvars::Union{Vector{Symbol},Nothing} = nothing,
  contentRowNames::Union{Vector{String},Nothing}=nothing,
  displayedvarsscaling::Union{Vector{Float64},Nothing} = nothing,
  regressionfunction::Function = olssensitivitypanelreg,
  titleCaption::String = "TBD",
  outputpath::String = TABLE_PATH,
  errorfunction::Function = getClustered!)

  #first run the regressions
  local specs::FMSpecs = regressionfunction(df)

  #apply the defaults based on the sepcs
  isnothing(displayedvars) &&  (displayedvars = setdiff(unique([specs.xnames...;]), [:intercept]))
  isnothing(contentRowNames) &&  (contentRowNames = (string).(displayedvars))
  isnothing(displayedvarsscaling) && (displayedvarsscaling = fill(10000.0, length(displayedvars)))


  local descRowNames::Vector{String} = [
    "Nonprofit(NP) FE",
    "Year FE",
    "NP Clustered",
    "\$R^{2}\$ (within)",
    "N"]

  #allocate space for the descriptive rows
  local descContent::Vector{Vector{String}} =
    ((i::Int)->Vector{String}(undef, specs.N[])).(1:length(descRowNames))

    #this builds the descriptive rows. There is Probably a better way to do this,
  #but its fairly specific to the project.
  for i ∈ 1:specs.N[]
    r::Int = 0
    descContent[r+=1][i] = "$(specs.withinspecs[i]==:ein ||
        occursin("ein", string(specs.xspecs[i])) ? "X" : "")"
    descContent[r+=1][i] = "$(specs.withinspecs[i]==:fisyr ||
        occursin("fisyr", string(specs.xspecs[i])) ? "X" : "")"
    descContent[r+=1][i] = "$(specs.clusterspecs[i]==:ein ? "X" : "")"
    descContent[r+=1][i] =
      "$(num2Str(getR²(specs.results[i], adjusted=false),  decimals))"
    descContent[r+=1][i] = "$(specs.results[i].N)"
  end

  #NOTE delete this if the console gets too cluttered
  display(descContent)

  local colNames::Vector{Vector{String}} = [
    ["return", "SP500 Adj", "FF3 Adj"],
    specs.specnames]

  local widthColNames::Vector{Vector{Int}} = [[2,2,2], ones(Int, specs.N[])]

  local alignmentColNames::Vector{Vector{String}} =
    [fill("c",length(widthColNames[1])), fill("r", specs.N[])]

  #checks based on the dimensions of the table parameters
  @assert (minimum((sum).(widthColNames)) == specs.N[] &&
    maximum((sum).(widthColNames)) == specs.N[] &&
    length(contentRowNames) == length(displayedvars)
    )



  local tabletext::String = texTable(specs.results,
    errorfunction#=getNeweyWestFunc(5) (lm::CTLM)->getNeweyWestSlow(lm, 5)=#,
    displayedvars,
    titleCaption = titleCaption,
    colNames = colNames,
    widthColNames = widthColNames,
    contentRowNames = contentRowNames,
    descRowNames = descRowNames,
    descContent = descContent,
    decimalDigits = 1,
    scaling= displayedvarsscaling,
    stars=true,
    starStrings = OVERRIDE_STAR_STRINGS,
    #clearMem = USE_AGGRESSIVE_GC,
    caption = "to be written",
    nakedTable=true)


  writeNakedTable(tabletext, path=outputpath, outName="olspanel$suffix.tex")
  println("Regressions for olspanel$suffix.tex complete and recorded.")

  return tabletext
end

function prefilterforolssensitivity(df::DataFrame;
  minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA)::DataFrame

  local subdf::SubDataFrame
  local ssubdf::SubDataFrame


  #make sure we ahve enough unique points for each organization
  df[:useforregressions] = falses(size(df,1))
  checkfields::Vector{Symbol} = [:lreturn, :Llcontributions, :useforregressions]

  #loop over all unique eins
  for subdf ∈ groupby(df, :ein)

    ssubdf = view(subdf, completecases(subdf[checkfields]), checkfields)
    subN = size(unique(ssubdf[:Llcontributions]),1)


    #make sure we have enough points for the regression
    if subN ≥ minpointspluskforβ + 2 #add two for the intercept and Llcontributions
      ssubdf[:useforregressions] .= true
    end
  end

  return df[df[:useforregressions],:]
end

function makeolspaneltables(data::NCCSData)

  #this gets repeated a lot so we can save time doing it now
  regdf::DataFrame = prefilterforolssensitivity(data.df)

  olssensitivitypanel(regdf, suffix = "programexpense",
    #displayedvars = nothing,
    #contentRowNames = nothing,
    regressionfunction = olssensitivitypanelregPE,
    titleCaption = "Program Expense vs Contribution Return Sensitivity")

  #NOTE: Failed proxy
  olssensitivitypanel(regdf, suffix = "publicsupport",
    #displayedvars = nothing,
    #contentRowNames = nothing,
    regressionfunction = olssensitivitypanelregPS,
    titleCaption = "PublicSupport vs Contribution Return Sensitivity")

  olssensitivitypanel(regdf, suffix = "2pass",
    #displayedvars = nothing,
    #contentRowNames = nothing,
    regressionfunction = olssensitivity2pass,
    errorfunction = getWhiteΣ!,
    titleCaption = "Effect of public support and program expense on sensitivity")

end



function olssensitivitytables(data::NCCSData)
  makeolspaneltables(data)
end
