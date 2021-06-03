


function prefilterforolscontributions(df::DataFrame;
  minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA,
  returnfield::Symbol = :lreturn,
  essentialfields::Vector{Symbol}=[:Llreturnwin, :pcontributions])::DataFrame

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

function makeolscontributiontables(data::NCCSData; rerunprep::Bool = true)

  #NOTE: This is also called in the compute function
  #if nothing is changed it does not need to be run

  #can disable the below

  local displayedvars::Vector{Symbol}

  #this gets repeated a lot so we can save time doing it now
  regdf::DataFrame = prefilterforolscontributions(data.df)

  #Univariate Regressions
  displayedvars = [:LagReturn, :LagReturn3yr]
  olspanel(regdf, suffix = "contribv2",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olsretcontributionL13,
    titleCaption = "Contributions vs Lagged Returns")


  displayedvars = [:LagContrib, :LagExp]
  olspanel(regdf, suffix = "retcon",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olscontributionregretcon,
    titleCaption = "Returns vs Lagged Contributions")

  displayedvars = [:LagReturn, :LagReturn3yr, :LagExp]
  olspanel(regdf, suffix = "conret",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olscontributionregconret,
    titleCaption = "Returns vs Lagged Contributions")

  #expense regressions
  #=displayedvars = [:LagReturn, :Lag2xReturn, :LagReturn3yr, :Lag2xReturn3yr]
  olspanel(regdf, suffix = "contrib2xlagv2",
    displayedvars = displayedvars,
    contentRowNames = (string).(displayedvars),
    regressionfunction = olsretcontributionLL13,
    titleCaption = "Contributions vs Twice-lagged Returns")=#



end



function olscontributiontables(data::NCCSData)
  makeolspaneltables(data)
end
