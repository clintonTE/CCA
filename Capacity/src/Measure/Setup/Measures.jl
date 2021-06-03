

#eliminates missing values
function cleanpanel(panel::DataFrame, ξs::Vector{Characteristic})
  essentialfields =  [ξ.Fξ for ξ ∈ ξs]


  pushmaybe(s::Symbol) = s ∈ essentialfields ? essentialfields : push!(essentialfields, s)
  pushmaybe(::Nothing) = essentialfields
  for ξ ∈ ξs
    if ξ.ws.specname === PARAM[:caraname]
      pushmaybe(PARAM[:iivol])
    elseif ξ.ws.specname === PARAM[:sortedname]
      pushmaybe(PARAM[:sorteddefaultconditioning])
    elseif ξ.ws.specname === PARAM[:doublesortedname]
      pushmaybe(PARAM[:doublesorteddefaultconditioning])
      pushmaybe(PARAM[:doublesorteddefaultsecondary])
    elseif ξ.ws.specname === PARAM[:zname]
      continue
    elseif ξ.ws.specname === PARAM[:fixedname]
      continue
    else
      @assert false
    end
  end

  for f ∈ essentialfields
    panel[!, f] .= finiteormissing.(panel[!, f])
  end

  panel = panel[completecases(panel, essentialfields), :]

  return panel
end

function basiccharacteristics!(panel::DataFrame;
  strategies=error("strategies is required"),
  strategyprefix::Symbol = error("strategyprefix is required"),
  weightspec::Symbol=error("weightspec is required"),
  preconditionertype::Symbol=error("preconditioner is required"))

  #form the weight spec
  local preconditioner
  if preconditionertype === :quantiledemean
    preconditioner = demeanξ ∘ quantiles
  elseif preconditionertype === :quantile
    preconditioner = quantiles
  elseif preconditionertype === :identity
    preconditioner = (x)->x
  else
    error("unrecognized preconditionertype $preconditionertype")
  end

  local ws::WeightSpec
  if weightspec === :z
    ws = zweightspec(; preconditioner)
  elseif weightspec === :cara
    ws = caraweightspec(; preconditioner)
  elseif weightspec === :sorted
    ws = sortedweightspec(; preconditioner)
  elseif weightspec === :sorted2
    ws = doublesortedweightspec(; preconditioner)
  elseif weightspec === :fixed
    ws = fixedweightspec(; preconditioner)
  else
    error("unrecognized weightspec $weightspec")
  end

  Fstrategies = Symbol.(strategyprefix, strategies)


  ξs = (s->Characteristic(ws, s)).(Fstrategies)

  return ξs
end
basiccharacteristics!(panel::DataFrame, ::Val{:momentum};
  strategies=PARAM[:momentumstrategiesused],
  strategyprefix::Symbol = PARAM[:momentumstrategyprefix],
  preconditionertype::Symbol = PARAM[:momentumpreconditionertype],
  weightspec::Symbol=PARAM[:momentumweightspec]) = basiccharacteristics!(panel;
    strategies, strategyprefix, preconditionertype, weightspec)

basiccharacteristics!(panel::DataFrame ,::Val{:placebo};
  strategies=PARAM[:placebostrategiesused],
  strategyprefix::Symbol = Symbol(),
  preconditionertype::Symbol = PARAM[:placebopreconditionertype],
  weightspec::Symbol=PARAM[:placeboweightspec]) = basiccharacteristics!(panel;
    strategies, strategyprefix, preconditionertype, weightspec)

basiccharacteristics!(panel::DataFrame, ::Val{:fundamental};
  strategies=PARAM[:fundamentalstrategiesused],
  strategyprefix::Symbol = PARAM[:fundamentalstrategyprefix],
  preconditionertype::Symbol = PARAM[:fundamentalpreconditionertype],
  weightspec::Symbol=PARAM[:fundamentalweightspec]) = basiccharacteristics!(panel;
    strategies, strategyprefix, preconditionertype, weightspec)


basiccharacteristics!(panel::DataFrame, ::Val{:fixed};
  strategies=PARAM[:fixedstrategiesused],
  strategyprefix::Symbol = PARAM[:fixedstrategyprefix],
  preconditionertype::Symbol = :identity,
  weightspec::Symbol=:fixed) = basiccharacteristics!(panel;
    strategies, strategyprefix, preconditionertype, weightspec)

#####################

#forms a measure from a combination of different characteristic categories
#this is the main function for forming the weights and measure spec,
#though we can make others as needed with the measuretype parameter
function measure!(panel::DataFrame, ::Val{:combinedbasic};
  strategycategories::Vector{Symbol} = PARAM[:combinedstrategycategories])

  ξs = reduce(vcat, (C->basiccharacteristics!(panel, Val(C))).(strategycategories))

  panel = cleanpanel(panel, ξs)

  #assign weights
  (ξ->assignweights!(panel, ξ)).(ξs)
  #error("Got through it! make sure all the testing code is reverted")

  #form the measure
  ms::MeasureSpec = MeasureSpec(ξs...)

  return (panel, ms)
end

#a simple version for grabbing the MeasureSpec
function MeasureSpec(::Val{:combinedbasic};
  strategycategories::Vector{Symbol} = PARAM[:combinedstrategycategories])

  ξs = reduce(vcat, (C->basiccharacteristics!(panel, Val(C))).(strategycategories))
  return MeasureSpec(ξs...)
end


##########################################
#Simulated data from excel
#replaces actual data with simulated data
function measure!(::Any, ::Val{:testsimulated};
  Fxis = [:xi1,:xi2],
  testpath::String = PARAM[:testpath],
  testdataname::String = PARAM[:testdataname],
  testdateformat=PARAM[:wrdsdateformat])

  #read in the test dataset
  panel = CSV.File("$testpath\\$testdataname.csv") |> DataFrame
  panel.date = (s->Dates.Date("$s", testdateformat)).(panel.date)
  sort!(panel, [:permno,:date])

  panel.Wdvol = panel.dvol #do not winsorize since we are trying to match excel
  #differencewithin2!(panel, [:Wdvol], :permno, date=:date)

  wsz::WeightSpec = zweightspec(preconditioner=demeanξ)
  #ξs = (Fxi->Characteristic(wsz, Symbol(:x_,Fxi))).(Fxis)
  ξs = (Fxi->Characteristic(wsz, Fxi)).(Fxis)
  (ξ->assignweights!(panel, ξ)).(ξs)
  ms::MeasureSpec = MeasureSpec(ξs...)


  return (panel, ms)
end
