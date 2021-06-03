



#does any data processign that is needed
#generally called from ComputeNCCS since it is cacheable
function prepsensitivityregressions!(data::NCCSData;
    toosmalltolog::Float64 = TOO_SMALL_TO_LOG,
    period::Symbol = :fisyr, group=:ein,)

  N::Int = size(data.df, 1)
  #logs if possible or returns missing
  plog(f::MFloat64)::MFloat64 = ismissing(f) || f ≤ toosmalltolog ? missing : log(f)

  #cleans up a common-size ratio
  cleanratio(f::MFloat64) = ((!ismissing(f)) && isfinite(f) &&
    (f ≤ 1.0)  && (f > toosmalltolog) ? f : missing)

  #log our focal variables
  #data.df[:lcontributions] = (plog).(data.df[:contributions])
  #data.df[:lprogramexpense] = (plog).(data.df[:programexpense])
  #data.df[:lprogramexpense_avg] = (plog).(data.df[:programexpense_avg])

  #taking the average for control purposes only with valid contribution fields
  #data.df[:lregnetassets] = ((c::MFloat64, a::MFloat64)->
  #  ismissing(c) ? missing : a).(data.df[:contributions], data.df[:lnetassets])
  sort!(data.df, [group, period])
  lagwithin!(data, :netassets, sorted=true)
  #  lagwithin!(data, :totexp, sorted=true)

  #scale everythign to lagged net assets
  data.df[:pcontributions] = Vector{MFloat64}(undef, N)
  data.df[:pcontributions] .= data.df[:contributions] ./ data.df[:Lnetassets]
  data.df[:pcontributions] .= (cleanratio).(data.df[:pcontributions])
  #data.df[:lpcontributions] = (plog).(data.df[:pcontributions])

  data.df[:nonzeroprogramexpense] = (f->
    (!ismissing(f)) && (f>0.0) ? f : missing ).(data.df[:programexpense])
  #data.df[:nonzeromgmtexpense] = (f->
  #  (!ismissing(f)) && (f>0.0) ? f : missing ).(data.df[:mgmtexpense])

  data.df[:pprogramexpense] = Vector{MFloat64}(undef, N)
  data.df[:pprogramexpense] .= data.df[:nonzeroprogramexpense] ./ data.df[:Lnetassets]
    #data.df[:mgmtexpense] ./ data.df[:Lnetassets]
    #(data.df[:totexp] .- data.df[:nonzeroprogramexpense]) ./ data.df[:Lnetassets]
    #data.df[:programexpense] ./ data.df[:totexp] #./ data.df[:Lnetassets]
  data.df[:pprogramexpense] .= (cleanratio).(data.df[:pprogramexpense])
  #data.df[:lpprogramexpense] = (plog).(data.df[:pprogramexpense])

  #get the average expense
  fieldstoaggregate::Vector{Symbol} = [:pprogramexpense]
  aggdf::DataFrame = aggregatefielddf(data.df[[:ein;
    fieldstoaggregate]], rejoin=false)
  unique!(aggdf)

  #the setdiff below is to avoid duplicates
  fieldstonotmerge::Vector{Symbol} = (s->Symbol(s,"_avg")).(fieldstoaggregate)
  data.df = join(data.df[setdiff(names(data.df), fieldstonotmerge)],
    aggdf, on=[:ein], kind=:left)

  #resort after the join
  sort!(data.df, [group, period])
  #lag our focal variables
  #lagwithin!(data, :lcontributions)
  lagwithin!(data, :pcontributions, sorted=true)
  #lagwithin!(data, :lpcontributions, sorted=true)
  lagwithin!(data, :publicsupport, sorted=true)
  #lagwithin!(data, :lprogramexpense, sorted=true)
  lagwithin!(data, :pprogramexpense, sorted=true)
  #lagwithin!(data, :lpprogramexpense, sorted=true)
  lagwithin!(data, :lreturn3yr, sorted=true)
  lagwithin!(data, :lreturnnetbmsp5003yr, sorted=true)
  lagwithin!(data, :lreturnnetbmff33yr, sorted=true)

  createquantiles!(data.df, :pcontributions)
  createquantiles!(data.df, :pprogramexpense)
  createquantiles!(data.df, :lreturn3yr)
  createquantiles!(data.df, :lreturnnetbmsp5003yr)
  createquantiles!(data.df, :lreturnnetbmff33yr)

  differencewithin!(data, :pcontributions, sorted=true)
  createquantiles!(data.df, :Dpcontributions)

  #println(describe(data.df[[:ein, :lpprogramexpense, :lpprogramexpense_avg, :Llpprogramexpense]]))

  #generate the interactions
  data.df[:LpublicsupportXLpcont] = data.df[:Lpublicsupport] .* data.df[:Lpcontributions]
  data.df[:publicsupport_avgXLpcont] = data.df[:publicsupport_avg] .* data.df[:Lpcontributions]

  data.df[:LpprogramexpenseXLpcont] = data.df[:Lpprogramexpense] .* data.df[:Lpcontributions]
  data.df[:pprogramexpense_avgXLpcont] = data.df[:pprogramexpense_avg] .* data.df[:Lpcontributions]

  data.df[:publicsupportXpcont] = data.df[:publicsupport] .* data.df[:pcontributions]
  data.df[:publicsupport_avgXpcont] = data.df[:publicsupport_avg] .* data.df[:pcontributions]

  data.df[:pprogramexpenseXpcont] = data.df[:pprogramexpense] .* data.df[:pcontributions]
  data.df[:pprogramexpense_avgXpcont] = data.df[:pprogramexpense_avg] .* data.df[:pcontributions]

  #I think regressions work better with symbols instead of strings
  rename!(data.df, :ein=>:ein_old)
  data.df[:ein] = (Symbol).(data.df[:ein_old])
  deletecols!(data.df, :ein_old)

  return nothing
end

#makes a scatterplot of betas
function plotβ!(df::DataFrame;
  xfield::Symbol = :pprogramexpense_avg,
  yfield::Symbol = :b1lRet3yrNetSP500Con,#:b1lRet3yrCon,
  suffix::String = "",
  outputpath::String = GRAPH_PATH,
  shortnameindex::Dict = SHORT_NAME_INDEX,
  computequantiles::Bool = true
  )

  local plotnames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local title::String

  if computequantiles
    df = deepcopy(df)
    createquantiles!(df, yfield, groupfields = nothing)
    #    createquantiles!(df, xfield, groupfields = nothing)
    #    xfield = Symbol(:q, xfield)
    yfield = Symbol(:q, yfield)
    push!(plotnames, "betaplot-$(shortnameindex[yfield])V$(shortnameindex[xfield])$suffix")
    title = "Return-Contrib. Beta ($(shortnameindex[yfield]))
      vs. Prog. Expense Percentile ($(shortnameindex[xfield])$suffix)"

  else
    push!(plotnames, "betaplot-$(shortnameindex[yfield])V$(shortnameindex[xfield])$suffix")
    title = "Return-Contribution Beta ($(shortnameindex[yfield]))
      vs. Program Expense ($(shortnameindex[xfield])$suffix)"
  end

  push!(plots, plot(df[completecases(df[[xfield,yfield]]),:],
          x=xfield, y = yfield,
          Guide.title("Return-Contribution Beta ($(shortnameindex[yfield]))
            vs. Program Expense ($(shortnameindex[xfield])$suffix)"),
          Guide.xlabel("Governance ($(shortnameindex[xfield]))"),
          Guide.ylabel("density$(shortnameindex[yfield])"),
          Geom.hexbin(xbincount=25, ybincount=25),
          #Geom.point,
          #Coord.cartesian(xmin=0.0, xmax=1.0)
          ))


  for i ∈ 1:length(plotnames)
    draw(PDF("$outputpath\\$(plotnames[i]).pdf", 9inch, 7inch), plots[i])
  end

  return nothing
end


#primarily generates the regeression table
function olssensitivitypanel(df::AbstractDataFrame;
  decimals::Int=DECIMALS,
  suffix::String = "",
  displayedvars::Union{Vector{Symbol},Nothing} = nothing,
  contentRowNames::Union{Vector{String},Nothing}=nothing,
  displayedvarsscaling::Union{Vector{Float64},Nothing} = nothing,
  regressionfunction::Function = olssensitivitypanelreg,
  titleCaption::String = "TBD",
  outputpath::String = TABLE_PATH,
  errorfunction::T where T<:Union{Function, Vector{Function}} = getClustered!,
  scalingvalue::Float64 = 100.0,
  categorynames::Vector{String} = ["return", "SP500 Adj", "FF3 Adj"],
  categorywidths::Vector{Int} = [2,2,2])

  #first run the regressions
  local specs::FMSpecs = regressionfunction(df)

  #apply the defaults based on the sepcs
  isnothing(displayedvars) &&  (displayedvars = setdiff(unique([specs.xnames...;]), [:intercept]))
  isnothing(contentRowNames) &&  (contentRowNames = (string).(displayedvars))
  isnothing(displayedvarsscaling) && (displayedvarsscaling = fill(scalingvalue, length(displayedvars)))


  local descRowNames::Vector{String} = [
    "Nonprofit(NP) FE",
    "Category FE",
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
    descContent[r+=1][i] = (
        occursin("ein", string(specs.xspecs[i], specs.withinspecs[i])) ? "X" : "")
    descContent[r+=1][i] = (
        occursin("category", string(specs.xspecs[i], specs.withinspecs[i])) ? "X" : "")
    descContent[r+=1][i] = "$(specs.clusterspecs[i]==:ein ? "X" : "")"
    descContent[r+=1][i] =
      "$(num2Str(getR²(specs.results[i], adjusted=false),  decimals))"
    descContent[r+=1][i] = "$(specs.results[i].N)"
  end

  local colNames::Vector{Vector{String}} = [categorynames, specs.specnames]
  local widthColNames::Vector{Vector{Int}} = [categorywidths, ones(Int, specs.N[])]
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
    alignmentColNames = alignmentColNames,
    titleCaption = titleCaption,
    colNames = colNames,
    widthColNames = widthColNames,
    contentRowNames = contentRowNames,
    descRowNames = descRowNames,
    descContent = descContent,
    decimalDigits = 2,
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
  minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA,
  returnfield::Symbol = :lreturn,
  essentialfields::Vector{Symbol}=[#=:Llcontributions,=# :Lpcontributions])::DataFrame

  local subdf::SubDataFrame
  local ssubdf::SubDataFrame


  #make sure we ahve enough unique points for each organization
  df[:useforregressions] = trues(size(df,1))


  checkfields::Vector{Symbol} = [returnfield; essentialfields; :useforregressions]

  #loop over all unique eins
  for subdf ∈ groupby(df, :ein)

    ssubdf = view(subdf, completecases(subdf[checkfields]), checkfields)
    subN = minimum((f->size(unique(ssubdf[f]),1)).(essentialfields))


    #make sure we have enough points for the regression
    if subN < minpointspluskforβ + 2 #add two for the existing dof
      ssubdf[:useforregressions] .= false
    end
  end

  return df[df[:useforregressions],:]
end

function makeolspaneltables(data::NCCSData; rerunprep::Bool = true)

  #NOTE: This is also called in the compute function
  #if nothing is changed it does not need to be run
  rerunprep && prepsensitivityregressions!(data)

  #can disable the below

  local displayedvars::Vector{Symbol}

  #this gets repeated a lot so we can save time doing it now
  regdf::DataFrame = prefilterforolssensitivity(data.df)

  #Univariate Regressions
  #=displayedvars = [:intercept, :LagContrib, :LagExp]
  olssensitivitypanel(regdf, suffix = "univariate",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olscontributionreg,
    titleCaption = "Program Expense vs Contribution Return Sensitivity")

  displayedvars = [:intercept, :Contrib, :Exp]
  olssensitivitypanel(regdf, suffix = "univariatehist",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olscontributionreghist,
    titleCaption = "Program Expense vs Contribution Return Sensitivity")=#

  displayedvars = [:intercept, :LagContrib, :Contrib,]
  olssensitivitypanel(regdf, suffix = "univariatecomb",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olscontributionregcomb,
    titleCaption = "Program Expense vs Contribution Return Sensitivity")

  #expense regressions
  displayedvars = [:LagContribAvgExp, :LagContribExp, :LagContrib]
  olssensitivitypanel(regdf, suffix = "programexpense",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olssensitivitypanelregPE,
    titleCaption = "Program Expense vs Contribution Return Sensitivity")

  displayedvars = [:ContribAvgExp, :ContribExp, :Contrib]
  olssensitivitypanel(regdf, suffix = "programexpensehist",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olssensitivitypanelregPEhist,
    titleCaption = "Governance vs Contribution Lagged 3YReturn Sensitivity")
  #NOTE: Failed proxy: Public Support
  #=displayedvars = [:LagContribAvgSup, :LagContribSup, :LagContrib]
  olssensitivitypanel(regdf, suffix = "publicsupport",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olssensitivitypanelregPS,
    titleCaption = "PublicSupport vs Contribution Return Sensitivity")=#



  #pick the standard error functions
  nwerror::Function = getNeweyWestFunc(3)
  mwerror::Function = getWhiteΣ!
  errorfunctions::Vector{Function} = [mwerror, nwerror, mwerror, nwerror, mwerror, nwerror]
  olssensitivitypanel(regdf, suffix = "2pass",
    #displayedvars = nothing,
    #contentRowNames = nothing,
    regressionfunction = olssensitivity2pass,
    errorfunction = errorfunctions,
    titleCaption = "Effect of program expense on sensitivity to expected returns")
  #olssensitivitypanel(regdf,
  #  suffix = "2passhist",
  #  #displayedvars = nothing,
  #  #contentRowNames = nothing,
  #  regressionfunction = olssensitivity2passhist,
  #  errorfunction = errorfunctions,
  #  titleCaption = "Effect of program expense on sensitivity to historical returns")


  #tab::CrossTabTable = CrossTabTable(regdf, tabfield=:qDpcontributions)
  #writecrosstabtable(tab)

end



function olssensitivitytables(data::NCCSData)
  makeolspaneltables(data)
end
