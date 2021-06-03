
function computemeasure(;
    MeasureType::Val = Val(PARAM[:measuretype]),
    regressionmethod=REGRESSION_METHOD[PARAM[:regressionmethod]],
    testdate::Date = PARAM[:testdate],
    includecomp::Bool = PARAM[:compincludecomp]
  )

  #acquire the panel. be careful about short-circuiting refresh panel and changing progrma parameters
  panel = formpanel()

  #debugmsg("dbgfinal", panel)
  #throw("stop")

  (includecomp == (:gvkey ∈ propertynames(panel))) || error("includecomp = $includecomp
  yet :gvkey ∈ panel.propertynames = $(:gvkey ∈ propertynames(panel))!!! This implies that
  crsp or crsp and comp must be refreshed.")

  #this assigns the primary weights and sets parameters for estimating the measure
  #see Measures.jl for details
  (panel, ms::MeasureSpec) = measure!(panel, MeasureType)

  #used in calculations for the controls
  formmarketweights!(panel)

  #get the strategy returns for post estimation purposes
  returns = formreturns(panel, ms)
  funds = prepfunds(returns, ms)

  #form the controls seperate from the focal characteristics
  formcontrols!(panel, ms)

  #estimate the measure
  regressionmethod(panel, ms; returns, funds)

  @info "Measure $(PARAM[:refreshmeasure] ? "" : "NOT") computed.
    Panel size: $(size(panel)), uniquepermno: $(unique(panel.permno) |> length)"

  testregressioninputs(panel)

  return nothing
end


function testregressioninputs(panel::DataFrame)
  if PARAM[:testmeasureinputs]

    #time series test output
    testpermno::Int = PARAM[:testpermnomult]
    testview::SubDataFrame = view(panel, (panel.permno .% testpermno) .== 0, :)
    if size(testview,1)>20_000 #prevents scenarios where the file is larger than needed
      testview=view(testview, 1:20_000,:)
    end
    testview |> CSV.write("$(PARAM[:testpath])\\$(
      PARAM[:crspdaily])_$(PARAM[:crspfrequencysuffix])_measuredat.csv")

    #cross-sectional test output
    testdate::Date = PARAM[:testdate]
    testview = view(panel, panel.date .== testdate,:)
    testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdaily])_$(
        PARAM[:crspfrequencysuffix])_measuredat_$testdate.csv")
  end


  #special case where we replace the actuald ata with simulated data
  if PARAM[:measuretype] == :testsimulated
    @warn "Estimated simulated data only!!!"
    panel |> CSV.write("$(PARAM[:testpath])\\completesimulated_measuredat_$testdate.csv")
  end
end
