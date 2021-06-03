
#iterate through the momentums trategies and put them into the dataframe
function momentumstrategies!(panel::DataFrame;
  strategies::Dict=PARAM[:momentumstrategies], Fconditional::NSymbol=nothing)
  for (Fmomentum, strategy) ∈ strategies
    #we will attatch this message to an error if the property names don't line up as expected
    if strategy.type ≡ :longshort
      longshortmomentumstrategy!(panel, strategy; Fmomentum, Fconditional)
    elseif strategy.type ≡ :longonly
      length(propertynames(strategy)) == 3 || error(errmsg)
      longonlymomentumstrategy!(panel, strategy;
        Fmomentum, Flret, Fconditional)
    else
      throw(errmsg)
    end
  end

end

#this function computes a long-short momentum value
function longshortmomentumstrategy!(panel::DataFrame, strategy;
  Fmomentum=throw("Fmomentum is required"),
  Fconditional::NSymbol,
  )
  @unpack shortperiod, longperiod, shortfraction, longfraction, Flret = strategy
  @assert longperiod > shortperiod

  datemin=minimum(skipmissing(panel.date))
  datemax=maximum(skipmissing(panel.date))
  numdates::Int=panel.date |> unique |> length


  Fmeasurelong::Symbol = Symbol(Flret,replace("$longperiod", " "=>""))
  Fmeasureshort::Symbol = Symbol(Flret,replace("$shortperiod", " "=>""))

  trailingret!(panel, [longperiod, shortperiod],
    [longfraction, shortfraction],
    Fmeasures=[Fmeasurelong, Fmeasureshort]; Flret, Fconditional)

  panel[!, Fmomentum] = panel[!, Fmeasurelong] .- panel[!, Fmeasureshort]

  if PARAM[:testmomentum]
    testpermno::Int = PARAM[:testpermnomult]
    testview::SubDataFrame = view(panel, (panel.permno .% testpermno) .== 0, :)
    if size(testview,1)>20_000 #prevents scenarios where the file is larger than needed
      testview=view(testview, 1:20_000,:)
    end
    testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(Fconditional)_$Fmomentum.csv")
  end


  select!(panel, Not([Fmeasurelong,Fmeasureshort]))


  return nothing
end

#this function computes a long only momentum value
function longonlymomentumstrategy!(panel::DataFrame, strategy;
  Fmomentum=throw("Fmomentum is required"),
  Fconditional::NSymbol,
  )

  #@info "$Fmomentum: minpoints=$minpoints"
  @unpack period, fraction
  Fmeasure::Symbol = Symbol(Flret,replace("$period", " "=>""))

  if strategy.strictindex
    @unpack Flret = strategy
    trailingret!(panel, [calendarinterval], [fraction], Fmeasures=[Fmomentum]; Flret, Fconditional)

  if PARAM[:testmomentum]
    testpermno::Int = PARAM[:testpermnomult]
    testview::SubDataFrame = view(panel, (panel.permno .% testpermno) .== 0, :)
    if size(testview,1)>20_000 #prevents scenarios where the file is larger than needed
      testview=view(testview, 1:20_000,:)
    end
    testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(Fconditional)_$Fmomentum.csv")
  end

  return nothing
end


#this function computes a price-based momentum value
function longonlymomentumstrategy!(panel::DataFrame, strategy;
  Fmomentum=throw("Fmomentum is required"),
  Fconditional::NSymbol,
  )

  #@info "$Fmomentum: minpoints=$minpoints"
  @unpack period, fraction, Findex = strategy
  Fmeasure::Symbol = Symbol(Flret,replace("prc$period", " "=>""))

  stricttrailingret!(panel, [period], [fraction], Fmeasures=[Fmomentum]; Findex, Fconditional)

  if PARAM[:testmomentum]
    testpermno::Int = PARAM[:testpermnomult]
    testview::SubDataFrame = view(panel, (panel.permno .% testpermno) .== 0, :)
    if size(testview,1)>20_000 #prevents scenarios where the file is larger than needed
      testview=view(testview, 1:20_000,:)
    end
    testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(Fconditional)_$(Fmomentum)_idx.csv")
  end

  return nothing
end

#############################trailingret########################
#used to compute the actual trailing returns
hpr(v::AbstractVector) = hpr(exp(sum(v))-1)
hpr(x::Real) = x
hpr(::Missing) = throw("missings found when computing trailing ret!!!")
function trailingret!(crsp::DataFrame, calendarintervals::Vector{<:DatePeriod},
  fractions::Vector{Float64};
  Flret::Symbol,
  Fconditional::NSymbol,
  Fmeasures::Vector{Symbol})

  for (calendarinterval, fraction, Fmeasure) ∈ zip(calendarintervals, fractions, Fmeasures)
    crsp[!, Fmeasure] = Vector{MFloat64}(undef, size(crsp,1))
    trailingmeasure!(hpr, crsp, calendarinterval, fraction; Ftarget=Flret, Fconditional, Fmeasure)
  end

  return nothing
end

idxret(v::AbstractVector{T})::T where T = idxret(idxret[end]/idxret[1] - one(T))
idxret(x::Real) = isfinite(x) ? x : idxret(missing)
idxret(::Missing) = throw("missings found when computing trailing index ret!!!")
function stricttrailingret!(crsp::DataFrame, calendarintervals::Vector{<:DatePeriod},
  fractions::Vector{Float64};
  Findex::Symbol,
  Fconditional::NSymbol,
  Fmeasures::Vector{Symbol})

  for (calendarinterval, fraction, Fmeasure) ∈ zip(calendarintervals, fractions, Fmeasures)
    crsp[!, Fmeasure] = Vector{MFloat64}(undef, size(crsp,1))
    trailingmeasure!(hpr, crsp, calendarinterval, fraction; Ftarget=Flret, Fconditional, Fmeasure)
  end

  return nothing
end


############post-processing of momentum strategies

function conditionmomentumstrategies!(panel::DataFrame; strategies::Dict=PARAM[:momentumstrategies])
  for (Fmomentum, strategy) ∈ strategies
    #winsorize the momentum
    FWmomentum::Symbol = Symbol(:W, Fmomentum)
    FWXmomentum::Symbol = Symbol(:WX, Fmomentum)
    panel[!,FWmomentum] = winsorizequantile(panel[!, Fmomentum], PARAM[:qwinsorize], twosided=true)
    winsorizewithin!(panel, Fsource=Fmomentum, prop=PARAM[:qwinsorize], Fgroup=:date,
      twosided=true, Fnew=FWXmomentum)
    lagwithin2!(panel, [Fmomentum, FWmomentum, FWXmomentum], :permno, date=:date) #this could be parallelized

    #these wont work in the current setup
    #FWiivolWmomentum = Symbol(:Wiivol, FWmomentum)
    #panel[!,Symbol(:Wiivol, FWmomentum)] =  panel[!, :Wiivol] .* panel[!, FWmomentum]
    #panel[!,Symbol(:LWiivolL, FWmomentum)] =  panel[!, :LWiivol] .* panel[!, Symbol(:L, FWmomentum)]
    #lagwithin2!(panel, [FWiivolWmomentum], :permno, date=:date) #this could be parallelized
  end



  return nothing
end
