#holds the parameters for forming the weights
#not meant to hold any meaningful data- just parameters
struct WeightSpec
  preconditioner::Function
  crosssectionweights::Function
  weightprefix::Symbol

  gross::Float64
  net::Float64
  long::Float64
  short::Float64

  specname::Symbol
end

#const WeightSpec = NamedTuple

#helper functions
@inline computelong(gross::Real, net::Real) = (gross+net)/2
@inline computeshort(gross::Real, net::Real) = (gross-net)/2
@inline maxx0(x::Real) = max(x, 0.0)
@inline minx0(x::Real) = min(x, 0.0)

function assignweights!(panel::DataFrame, ξ::AbstractCharacteristic)

  (ξ.X[:Fw] ∈ names(panel)) && error("$(ξ.X[:Fw]) already exists in dataframe")
  panel[!,ξ.X[:Fw]] = Vector{MFloat64}(undef, size(panel,1))
  @info "assigning weight for $(ξ.X[:Fw])"

  spanels::GroupedDataFrame = groupby(panel, :date)

  ####### NOTE can delete this after testing
  thresholds=PARAM[:sorteddefaultthresholds]
  Fw = ξ.X[:Fw]
  Fξ = ξ.Fξ
  ws = ξ.ws
  ######

  #=Threads.@threads =#for i ∈ 1:length(spanels)
    spanel = spanels[i]

    ####WARNING FOR TESTING ONLY!!!!
    #=ξv::Vector{Float64} = spanel[:,Fξ]#ws.preconditioner(spanel[:,Fξ])
    ranks::Vector{Float64} = ξv  |> quantiles

    Fwcol = zeros(length(ranks))
    Fwcol[ranks .≥ thresholds.upper] .=  1.0
    Fwcol[ranks .< thresholds.lower] .=  -1.0
    #Fwcol .= ((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower))
    if !(minimum(Fwcol) ≈ -1.0 * maximum(Fwcol))
      testpanel = spanel |> DataFrame
      testpanel.raw = ξv
      testpanel.unconditioned = (ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower)
      testpanel.ranks = ranks
      testpanel |> CSV.write("$(PARAM[:testpath])\\error_subpanel_$(panel[1,:date])_$(Fw).csv")

      error("Ranking Failed!!
        for $Fξ -> $Fw, $i, $(spanel.date[1])
        sort category: $((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower))
        ranks: \$ranks
        spanel[:,Fξ]: \$(spanel[:,Fξ])
        typeof:$(typeof((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower)))
        $(spanel[1, Fw] |> typeof)")
    end=#
    ###############################
    ξ.ws.crosssectionweights(spanels[i], ξ.ws, ξ.Fξ, ξ.X[:Fw])

    ####sanity checks
    if ξ.ws.specname === :fixedweight #this is a long only strategy
      @assert 1.0 ≈ sum(spanels[i][!, ξ.X[:Fw]]) ≈ sum(abs.(spanels[i][!, ξ.X[:Fw]]))
      continue
    end
    @assert 1.0 ≈ 1.0 + sum(spanels[i][!, ξ.X[:Fw]])
    @assert 2.0 ≈ sum(abs.(spanels[i][!, ξ.X[:Fw]]))
  end



  #we want the weight fields, and the lagged product of the weights with the returns
  #NOTE: the accounting identity requires contemporaneous return
  #the characteristic however can be defined as either lag, contemporaneous or soemthing else
  # its all definitional
  lagwithin2!(panel, [ξ.X[:Fw], :ret], :permno, date=:date)
  panel[!, ξ.X[:FRLw]] = panel[!,ξ.X[:FLw]] .* (1 .+ panel[!, :ret])

  lagwithin2!(panel, [ξ.X[:FLw]], :permno, date=:date)
  panel[!, ξ.X[:FRLLw]] = panel[!,ξ.X[:FLLw]] .* (1 .+ panel[!, :Lret])
  lagwithin2!(panel, [ξ.X[:FLLw]], :permno, date=:date)

  #compute the forward return weights
  transform!(groupby(panel, :date),
    [ξ.X[:FRLw], ξ.X[:Fw]] => ((RLw,w) -> sum(skipmissing(RLw)) .* w) => ξ.X[:FRtw])

#  @assert ismissing.(panel[!, ξ.X[:FRLw]]) .== (ismissing.(panel[!,ξ.X[:FLw]]) .& ismissing.(panel[!, :ret])) |> all

  #check this was done correctly
  for spanel ∈ groupby(view(panel, completecases(panel, ξ.X[:FLw]), :), :date)
    #@assert (!ismissing).(spanel[:, ξ.X[:FRtw]]) |> all

    @assert (spanel[:, ξ.X[:FRtw]] ≈ (
      sum(spanel[:,ξ.X[:FLw]] .* (1.0 .+ spanel.ret)) .* spanel[:, ξ.X[:Fw]]))
  end

    return nothing
end


###########################################CARA Specs
#returns a vector of weights according to the CARA methodology
function caraweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  Fiivol::Symbol = PARAM[:portfolioiivol], minweightmagnitude=PARAM[:minweightmagnitude])

  local longraw::Float64
  local shortraw::Float64

  ξ = ws.preconditioner(panel[:,Fξ])

  panel[:, Fw] .= ξ .* panel[!, Fiivol]
  panel[:, Fw] .-= mean(panel[:, Fw])

  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  #correct for small values
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    if wₙ ≥ 0 && wₙ/longraw*ws.long < minweightmagnitude
      panel[i, Fw] = 0.0
    elseif wₙ < 0 && abs(wₙ)/shortraw*ws.short < minweightmagnitude
      panel[i, Fw] = 0.0
    end
  end

  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  if (!isfinite(longraw)) || (longraw==0.0)
    println("ERROR!")
    println(panel[!,[:date, :permno, Fiivol, Fw]])
    error("Longraw ($longraw) is not non-zero finite.")
  elseif (!isfinite(shortraw)) || (shortraw==0.0)
    println(panel[!,[:date, :permno, Fiivol, Fw]])
    error("Shortraw ($shortraw) is not non-zero finite.")
  end

  #assigns the weight by first checking if the position is long or short
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    panel[i, Fw] = ifelse(wₙ ≥ 0, wₙ/longraw*ws.long, wₙ/shortraw*ws.short)
  end

  return nothing
end


#creates a specification for cara based weights
caraweightspec(;
  preconditioner = error("preconditioning function is required"),
  crosssectionweights = caraweights!,
  weightprefix::Symbol = PARAM[:caraweightprefix],
  gross=PARAM[:caragross],
  net=PARAM[:caranet],
  specificationname=PARAM[:caraname]) = WeightSpec(
    preconditioner,
    crosssectionweights,
    weightprefix,
    gross,
    net,
    computelong(gross,net),
    computeshort(gross,net),
    specificationname)


################################Z-weight spec (Equal conditional weights)

function zweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  minweightmagnitude=PARAM[:minweightmagnitude])

  ξ = ws.preconditioner(panel[:,Fξ])

  panel[:, Fw] .= ξ
  panel[:, Fw] .-= mean(panel[:, Fw])

  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  #correct for small values
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    if wₙ ≥ 0 && wₙ/longraw*ws.long < minweightmagnitude
      panel[i, Fw] = 0.0
    elseif wₙ < 0 && abs(wₙ)/shortraw*ws.short < minweightmagnitude
      panel[i, Fw] = 0.0
    end
  end

  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  if (!isfinite(longraw)) || (longraw==0.0)
    println("ERROR!")
    println(panel[!,[:date, :permno, Fw]])
    error("Longraw ($longraw) is not non-zero finite.")
  elseif (!isfinite(shortraw)) || (shortraw==0.0)
    println(panel[!,[:date, :permno, Fw]])
    error("Shortraw ($shortraw) is not non-zero finite.")
  end

  #assigns the weight by first checking if the position is long or short
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    panel[i, Fw] = ifelse(wₙ ≥ 0, wₙ/longraw*ws.long, wₙ/shortraw*ws.short)
  end

  #=panel[!, Fw] .= panel[!, Fξ] .- mean(panel[:, Fξ])
  panel[!, Fw] .= panel[!, Fξ] ./ sum((abs).(panel[!, Fw])) .* 2

  (sum((x->ismissing(x)||(!isfinite(x))).(panel[!,Fw])) == 0) || error(
    "INVALID values detected in weights for Fξ")=#

  @assert sum(panel[!, Fw]) + 1.0 ≈ 1.0
  @assert sum(abs.(panel[!, Fw])) ≈ 2.0


  return nothing
end

#creates a specification given equal conditional weights
zweightspec(;
  preconditioner = error("preconditioning function is required"),
  crosssectionweights = zweights!,
  weightprefix::Symbol = PARAM[:zweightprefix],
  gross=PARAM[:zgross],
  net=PARAM[:znet],
  specificationname=PARAM[:zname]) = WeightSpec(
    preconditioner,
    crosssectionweights,
    weightprefix,
    gross,
    net,
    computelong(gross,net),
    computeshort(gross,net),
    specificationname)


#####################sortspec
#Applies preconditioner to the characteristic -usually either demeaning, quantilizaiton, or both,
#The preconditioned characteristic is then multiplied by a conditioner (such as market cap)
#The conditioner weights the characteristic within a particular sort threshold
#Then the portfolios are formed within quantile thresholds, both long and short
#These are then finally normalized to add to a value (usually 1) per side
#note the conditioning vector can be overridden if customization is needed
function _sortedweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  thresholds::NamedTuple=error("thresholds is required"),
  conditioning::AbstractVector=error("conditioning is required"),
  minweightmagnitude=PARAM[:minweightmagnitude])

  ξ::Vector{Float64} = ws.preconditioner(panel[:,Fξ])
  ranks::Vector{Float64} = ξ  |> quantiles
  terciles = ((ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower))
  minimum(terciles) ≈ -1.0 * maximum(terciles) || error("
    (minimum(ξ.X[:Fw]) /≈ -1.0 * maximum(ξ.X[:Fw])): This probably
    means the sort failed to differentiate between the quantiles. Check winsorization
    and sorting thresholds to try again.")

  #conditioning is usually marketcap or a vector of ones for equal wieghting
  #see other sortwedweights! function
  panel[:, Fw] .= terciles .* conditioning



  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  #correct for small values
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    if wₙ ≥ 0 && wₙ/longraw*ws.long < minweightmagnitude
      panel[i, Fw] = 0.0
    elseif wₙ < 0 && abs(wₙ)/shortraw*ws.short < minweightmagnitude
      panel[i, Fw] = 0.0
    end
  end

  longraw = sum((maxx0).(panel[!, Fw]))
  shortraw = sum((minx0).(panel[!, Fw])) .* -1.0

  if (!isfinite(longraw)) || (longraw==0.0)
    #println("ERROR!")
    #println(panel[1:10,[:date, :permno, Fw]])
    error("Longraw ($longraw) is not non-zero finite for date/field $(panel[1,[:date, Fw]]).")
  elseif (!isfinite(shortraw)) || (shortraw==0.0)
    #println(panel[1:10,[:date, :permno, Fw]])
    #=testpanel = panel |> DataFrame
    testpanel.conditioning = conditioning
    testpanel.unconditioned = (ranks .≥ thresholds.upper) .- (ranks .< thresholds.lower)
    testpanel.ranks = ranks
    testpanel |> CSV.write("$(PARAM[:testpath])\\error_subpanel_$(panel[1,:date])_$(Fw).csv")=#

    error("Shortraw ($shortraw) is not non-zero finite for date/field $(panel[1,[:date, Fw]]).")
  end

  #assigns the weight by first checking if the position is long or short
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    panel[i, Fw] = ifelse(wₙ ≥ 0, wₙ/longraw*ws.long, wₙ/shortraw*ws.short)
  end

  #=panel[!, Fw] .= panel[!, Fξ] .- mean(panel[:, Fξ])
  panel[!, Fw] .= panel[!, Fξ] ./ sum((abs).(panel[!, Fw])) .* 2

  (sum((x->ismissing(x)||(!isfinite(x))).(panel[!,Fw])) == 0) || error(
    "INVALID values detected in weights for Fξ")=#

  @assert sum(panel[!, Fw]) + 1.0 ≈ 1.0
  @assert sum(abs.(panel[!, Fw])) ≈ 2.0


  return nothing
end


#override the kw arguments to get characteristic specific threshholds and conditioning
function sortedweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  thresholds=PARAM[:sorteddefaultthresholds],
  Fconditioning::Union{Nothing, Symbol} = PARAM[:sorteddefaultconditioning],)

  #extract the conditioning vector
  conditioning =  (Fconditioning === nothing ?
    ones(size(panel, 1)) |> Vector{MFloat64} : panel[!, Fconditioning])

  return _sortedweights!(panel, ws, Fξ, Fw; thresholds, conditioning)
end




#creates a specification given equal conditional weights
#WARNING- switch this back to zweights!
sortedweightspec(;
  preconditioner = error("preconditioning function is required"),
  crosssectionweights = sortedweights!,
  weightprefix::Symbol = PARAM[:sortedweightprefix],
  gross=PARAM[:sortedgross],
  net=PARAM[:sortednet],
  specificationname=PARAM[:sortedname]) = WeightSpec(
    preconditioner,
    crosssectionweights,
    weightprefix,
    gross,
    net,
    computelong(gross,net),
    computeshort(gross,net),
    specificationname)

#####################doublesortspec
#Applies preconditioner to the characteristic -usually either demeaning, quantilizaiton, or both,
#The preconditioned characteristic is then multiplied by a conditioner (such as market cap)
#The conditioner weights the characteristic within a particular sort threshold
#Then the portfolios are formed within quantile thresholds, both long and short
#These are then finally normalized to add to a value (usually 1) per side
#note the conditioning vector can be overridden if customization is needed
function _doublesortedweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  thresholds::NamedTuple=error("thresholds is required"),
  conditioning::AbstractVector=error("conditioning is required"),
  secondary::AbstractVector = error("secondary sort variable \"secondary\" is required"),
  minweightmagnitude=PARAM[:minweightmagnitude])

  ξ::Vector{Float64} = ws.preconditioner(panel[:,Fξ])
  primaryranks::Vector{Float64} = ξ  |> quantiles
  primaryterciles::Vector{Int} = ((primaryranks .≥ thresholds.primaryupper) .-
    (primaryranks .< thresholds.primarylower))
  minimum(primaryterciles) ≈ -1.0 * maximum(primaryterciles) || error("
    (minimum(ξ.X[:Fw]) /≈ -1.0 * maximum(ξ.X[:Fw])): This probably
    means the primary sort failed to differentiate between the quantiles. Check winsorization
    and sorting thresholds to try again.")

  secondaryranks::Vector{Float64} = ws.preconditioner(secondary) |> quantiles
  secondaryterciles::Vector{Int} = ((secondaryranks .≥ thresholds.secondaryupper) .-
    (secondaryranks .< thresholds.secondarylower))
  minimum(secondaryterciles) ≈ -1.0 * maximum(secondaryterciles) || error("
    (secondary) /≈ -1.0 * maximum(secondary)): This probably
    means the secondary sort failed to differentiate between the quantiles. Check winsorization
    and sorting thresholds to try again.")


  #conditioning is usually marketcap or a vector of ones for equal wieghting
  #see other sortwedweights! function
  w = (
    longupper = ((primaryterciles .== 1) .& (secondaryterciles .== 1)) .* conditioning,
    longlower = ((primaryterciles .== 1) .& (secondaryterciles .== -1)) .* conditioning,
    shortupper = -1 .* ((primaryterciles .== -1) .& (secondaryterciles .== 1)) .* conditioning,
    shortlower = -1 .* ((primaryterciles .== -1) .& (secondaryterciles .== -1)) .* conditioning)

  #below checks to make sure there are no collisions
  wcheck = reduce(+, w |> values)
  @assert (wcheck .+ 1.0 .≈ 1.0) .| (abs.(wcheck) .≈ conditioning) |> all
  @assert all(vcat(w.longupper, w.longlower, -1 .* w.shortupper, -1 .* w.shortlower) .≥ 0.0)
  longupperraw = sum(w.longupper)
  longlowerraw = sum(w.longlower)
  shortupperraw = sum(w.shortupper)  .* -1.0
  shortlowerraw = sum(w.shortlower)  .* -1.0


  #zero the small weights- note we divide the target leg leverage by two since we weight
  #the secondary upper and secondary lower portfolio components equally
  zerosmallweights(wₙ, raw, legleverage) = ifelse(
    abs(wₙ)/raw*legleverage<minweightmagnitude, 0.0, wₙ)
  w.longupper .= zerosmallweights.(w.longupper, longupperraw, ws.long/2)
  w.longlower .= zerosmallweights.(w.longlower, longlowerraw, ws.long/2)
  w.shortupper .= zerosmallweights.(w.shortupper, shortupperraw, ws.short/2)
  w.shortlower .= zerosmallweights.(w.shortlower, shortlowerraw, ws.short/2)

  #update weight totals
  longupperraw = sum(w.longupper)
  longlowerraw = sum(w.longlower)
  shortupperraw = sum(w.shortupper)  .* -1.0
  shortlowerraw = sum(w.shortlower)  .* -1.0

  @assert all(isfinite.(w |> values .|> sum))
  @assert !all((legtot -> abs(legtot) ≈ 1.0).(w |> values .|> sum))

  #assigns the weight by first checking if the position is long or short
  panel[:, Fw] .= 0.0
  panel[:, Fw] .+= w.longupper/longupperraw*ws.long/2
  panel[:, Fw] .+= w.longlower/longlowerraw*ws.long/2
  panel[:, Fw] .+= w.shortupper/shortupperraw*ws.short/2
  panel[:, Fw] .+= w.shortlower/shortlowerraw*ws.short/2

  #sanity check on total exposure
  @assert sum(panel[panel[!,Fw] .> 0,Fw]) ≈ 1.0
  @assert sum(panel[panel[!,Fw] .< 0,Fw]) ≈ -1.0



  return nothing
end


#override the kw arguments to get characteristic specific threshholds and conditioning
function doublesortedweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  thresholds=PARAM[:doublesorteddefaultthresholds],
  Fconditioning::Union{Nothing, Symbol} = PARAM[:doublesorteddefaultconditioning],
  Fsecondary = PARAM[:doublesorteddefaultsecondary],)

  #extract the conditioning vector
  conditioning =  (Fconditioning === nothing ?
    ones(size(panel, 1)) |> Vector{MFloat64} : panel[!, Fconditioning])

  secondary = panel[!, Fsecondary]

  return _doublesortedweights!(panel, ws, Fξ, Fw; thresholds, conditioning, secondary)
end




#creates a specification given equal conditional weights
#WARNING- switch this back to zweights!
doublesortedweightspec(;
  preconditioner = error("preconditioning function is required"),
  crosssectionweights = doublesortedweights!,
  weightprefix::Symbol = PARAM[:doublesortedweightprefix],
  gross=PARAM[:doublesortedgross],
  net=PARAM[:doublesortednet],
  specificationname=PARAM[:doublesortedname]) = WeightSpec(
    preconditioner,
    crosssectionweights,
    weightprefix,
    gross,
    net,
    computelong(gross,net),
    computeshort(gross,net),
    specificationname)

#############fixedweightspec

function fixedweights!(panel::AbstractDataFrame, ws::WeightSpec, Fξ::Symbol, Fw::Symbol;
  minweightmagnitude=PARAM[:minweightmagnitude])

  ξ = ws.preconditioner(panel[:,Fξ])

  panel[:, Fw] .= ξ
  longraw = sum((maxx0).(panel[!, Fw]))

  #correct for small values
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    if wₙ ≥ 0 && wₙ/longraw*ws.long < minweightmagnitude
      panel[i, Fw] = 0.0
    elseif wₙ < 0 && abs(wₙ)/shortraw*ws.short < minweightmagnitude
      error("negative weights not allowed for fixed weight spec")
      #panel[i, Fw] = 0.0
    end
  end

  longraw = sum((maxx0).(panel[!, Fw]))

  if (!isfinite(longraw)) || (longraw==0.0)
    println("ERROR!")
    println(panel[!,[:date, :permno, Fw]])
    error("Longraw ($longraw) is not non-zero finite.")
  end

  #assigns the weight by first checking if the position is long or short
  for (i,wₙ) ∈ enumerate(panel[!, Fw])
    panel[i, Fw] = wₙ/longraw*ws.long
  end


  @assert sum(panel[!, Fw]) ≈ sum(abs.(panel[!, Fw])) ≈ 1.0


  return nothing
end




#creates a specification given equal conditional weights
fixedweightspec(;
  preconditioner = error("preconditioning function is required"),
  crosssectionweights = fixedweights!,
  weightprefix::Symbol = PARAM[:fixedweightprefix],
  gross=1.0,
  net=1.0,
  specificationname=PARAM[:fixedname]) = WeightSpec(
    preconditioner,
    crosssectionweights,
    weightprefix,
    gross,
    net,
    computelong(gross,net),
    computeshort(gross,net),
    specificationname)
