mcmcresultsname(::Val{:ols}) = "ans_$(PARAM[:crspdatalabel])" *
  "_$(PARAM[:mcmcresultssuffix])" *
  "_$(PARAM[:mcmcmodel])" *
  "_$(PARAM[:mcmcprior])" *
  "_$(Dates.format(now(),"yyyymmdd_HHMMSS"))"

#optresultsname() = optresultsname(Val(PARAM[:itermodel]))
#primary function for estimating the measure using optimizations
#note that this still uses the same volume parts construct as flux
function formmcmcregression(panel::DataFrame, ms::MeasureSpec,
  MCMCModel::Val = Val(PARAM[:mcmcmodel]))

  #issorted(panel, [:permno, :date]) || error("panel must be sorted")
  sort!(panel, :date)


  #in contrast to other methods, drop any column where we don't have
  #the lag weight (improves performance I think)

  #conditionforturing!(panel::SubDataFrame, ms::MeasureSpec)
  zs::ZSpec = condition!(panel, ms)

  #testing
  if PARAM[:testmeasure]
    testmcmc(panel, zs)
  end

  if PARAM[:refreshmeasure]
    results::DataFrame = mcmcmodel(panel, zs, MCMCModel)
    #destandardize!(results, zs, intercept=false)

    #this will be the name of the results file
    #resultname::String = mcmcresultsname(λ⁺)

    #results |> CSV.write("$(PARAM[:outputpath])\\$(resultname).csv")

  end

end
