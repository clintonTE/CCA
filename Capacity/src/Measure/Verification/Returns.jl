
function formreturns(panel::DataFrame, ms::MeasureSpec;
  returnsname::String = PARAM[:returnsname],
  savepath::String = PARAM[:workingpath],
  returntype::String = PARAM[:crspdatalabel],
  frequencysuffix::String = PARAM[:crspfrequencysuffix],
  Fvol = PARAM[:Fvol],
  volscale::Float64 = PARAM[:volscale]
  #testdate::Date = PARAM[:testdate],
  )::DataFrame

  #housekeeping
  local returnspath::String = "$savepath\\$(returntype)_$(returnsname).csv"
  local returns::DataFrame


  #The below block directs the panel creation
  if PARAM[:refreshreturns]
    #returns = DataFrame(date=sort!(unique(panel.date)))
    returns = vcat(
      DataFrame(date=Date[], Wvolume=MFloat64[], volume=MFloat64[]),
      combine(groupby(panel, :date),
        :date => last => :date,
        :Wdvol => ((v)->sum(skipmissing(v))) => :Wvolume,
        :dvol => ((v)->sum(skipmissing(v)) ./ volscale) => :volume,
        ))

    #pad the dates to get an extra data point of returns
    if frequencysuffix ≡ "m"
      push!(returns, (; date=minimum(returns.date) - Month(1), Wvolume=missing, volume=missing))
    elseif frequencysuffix ≡ "w"
      push!(returns, (; date=minimum(returns.date) - Day(7), Wvolume=missing, volume=missing))
    else
      throw("unrecognized date frequency")
    end

    sret = returns[completecases(returns, [:volume, :Wvolume]), :]
    println("initial: size: $(size(sret)) dates: $(minimum(sret.date)):$(maximum(sret.date))")
    sort!(returns, :date)

    returns = addffreturns!(returns)
    sret = returns[completecases(returns, [:volume, :Wvolume]), :]
    println("afterff: size: $(size(sret)) dates: $(minimum(sret.date)):$(maximum(sret.date))")
    sort!(returns, :date)

    returns = addstrategyreturns!(returns, panel, ms)
    sret = returns[completecases(returns, [:volume, :Wvolume]), :]
    println("afterstrategy: size: $(size(sret)) dates: $(minimum(sret.date)):$(maximum(sret.date))")
    sort!(returns, :date)

    returns |> CSV.write(returnspath)

    sum(completecases(returns)) == nrow(returns) || @warn "WARNING: Some dates in the returns file
      are missing benchmark returns. Make sure you know why before disabling this error- in
      particular, it could be a symptom of stale benchmark data. You can check the returns file
      to see what is missing. On the other hand, if its just a result of using a lag
      then its probably safe to ignore this warning."
  else
    returns = CSV.File(returnspath) |> DataFrame
  end


  return returns
end

getreturns() = CSV.File(
    "$(PARAM[:workingpath])\\$(PARAM[:crspdatalabel])_$(PARAM[:returnsname]).csv") |> DataFrame


#compute the strategy returns
function addstrategyreturns!(returns::DataFrame, panel::DataFrame, ms::MeasureSpec)

  @assert issorted(returns, :date)

  for ξ ∈ ms.ξs
    #acquire the weight field
    #NOTE- not sure whether to use Lw or w here- doesn't mkae a big difference
    Fw = ξ.X[:FLw]
    strategyreturns = combine(groupby(view(panel, completecases(panel, [Fw]),:), :date),
      :date => last => :date, #just record the date
      [Fw, :ret] => weightedreturn => ξ.Fξ,
      [Fw, :ret] => ((w,ret)->weightedreturn(w .|> (x)->max(x,0.0),ret)) => Symbol(ξ.Fξ,:_long),
      [Fw, :ret] => ((w,ret)->weightedreturn(w .|> (x)->min(x,0.0),ret)) => Symbol(ξ.Fξ,:_short),
      )

    #=the below code reverses the sign of the return to the same as the effect on leverage
    #note this is not the same as passive change in portfolio leverage- rather, its the
    #passive change in portfolio leverage caused by net gains and losses
    Basically, profits from the long side are leveraging, while profits from short side
    #are deleveraging
    Logic to create the measure is as follows: start with the leverage effect
      (long profits - short profits)
      Then if part of the leverage effect is due to comovement (e.g. one side
      created gains, the other losses) then subtract out the offsetting part of the positions
      as that is most likely the effect of the market. If both sides created the same
      sign P&L, then the leverage effect is entirely due to the P&L so don't subtract anything.
      (this is actually the logic of the check)
      =#

    Flong, Fshort, Fnet = Symbol(ξ.Fξ,:_long), Symbol(ξ.Fξ,:_short), ξ.Fξ
    transform!(strategyreturns,
      [Flong, Fshort] => ((long, short) ->  long .- short) => Symbol(ξ.Fξ,:_retlevgross),
      [Flong, Fshort, Fnet] => ByRow((long, short, net) ->
        ifelse(sign(long)==sign(short),long - short, ifelse(abs(long)>abs(short),net, -net))) =>
        Symbol(ξ.Fξ,:_retlev),
      [Flong, Fshort, Fnet] => ByRow((long, short, net) -> long - short -
        ifelse(sign(long)==sign(short),0.0, 2.0*sign(long - short)*min(abs(long),abs(short)))) =>
        Symbol(ξ.Fξ,:_retlevcheck))
    @assert strategyreturns[!, Symbol(ξ.Fξ,:_retlev)] ≈ strategyreturns[!, Symbol(ξ.Fξ,:_retlevcheck)] "
      zip(strategyreturns[!, Symbol(ξ.Fξ,:_retlev)],strategyreturns[!, Symbol(ξ.Fξ,:_retlevcheck)])
        $(zip(strategyreturns[!, Symbol(ξ.Fξ,:_retlev)],
          strategyreturns[!, Symbol(ξ.Fξ,:_retlevcheck)]) |> collect)"
    select!(strategyreturns, Not(Symbol(ξ.Fξ,:_retlevcheck)))


    #the below is another choice- using a leftjoin uses all available data, but
    #using an innerjoin keeps the time windows consistent
    #the below is a compromise- use a left join but bound how far it goes back
    startdate = returns.date[
      max(something(findfirst(isequal(minimum(strategyreturns.date)), returns.date),1)-1,1)]
    boundedreturns = view(returns,
      (returns.date .≥ startdate) .&
      (returns.date .≤ maximum(strategyreturns.date)), :)
    returns=leftjoin(boundedreturns, strategyreturns, on=:date)

    if !issorted(returns, :date)
      sort!(returns, :date)
    end
    #below is the old innerjoin code
    #=returns = innerjoin(returns, combine(groupby(strategyreturns, :date),
      :date => last => :date, #just record the date
      [Fw, :ret] => weightedreturn => ξ.Fξ), on=:date)=#
  end

  return returns
end

#note this natually will lead to higher/lower returns depending on the gross leverage employed
#e.g. if the long leg returns +10%, short leg returns -10%, and both sides have weights
#that add to 1.0/-1.0, then the resulting return will be 1.0*0.1 - 1.0*-0.1 = 20%
weightedreturn(weights, returns) = sum(weights .* returns)


#function addstrategy

#merges in the FF returns
function addffreturns!(returns::DataFrame)
  #this is less important when just combining the time series as opposed to creating a panel
  #but we might as well be consistent
  local datefrequency
  if PARAM[:crspfrequencysuffix] == "m"
    datefrequency = :month
  elseif PARAM[:crspfrequencysuffix] == "w"
    datefrequency = :week
  else
    @assert false
  end

  #do I break up mergeff or bootstrap a group into returns?
  returns = mergeff(returns; datefrequency, testlabel="FFreturns")

  return returns
end

#################################################
#the rest of these methods handle processing wide-form returns
#and merging the returns with a panel

#computes a cumulative return skipping missing values in place
function skippingcumret(rets::AbstractVector{T}) where T<:MFloat64
  base::Float64 = 1.0

  out::Vector{T} = rets .+ 1.0
  completeout = view(out,out .!== missing) #want to skip missing values while totalling
  completeout .= cumprod(completeout)

  return out
end

#creates a cumulative index from a return stream
function returns2index(dates::AbstractVector{Date},
  rets::AbstractVector{<:MFloat64})

  issorted(dates) || error("dates must be sorted to form index")

  @assert length(dates) == length(rets)
  index::Vector{Float64} = skippingcumret(rets)

  return :index, DataFrame(date=dates, index=index)
end


#creates a cumulative index from a return stream
function returns2index(df::AbstractDataFrame;
  Fdate::Symbol = :date,
  Frets::Vector{Symbol} = error("Frets is required"))

  issorted(df, Fdate) || error("df must be sorted by date")
  #create the cumulative index
  indexdf = DataFrame(date=df.date)
  Findexes = (Fret->Symbol(Fret,:_idx)).(Frets)
  for (Fret, Findex) ∈ zip(Frets, Findexes)
    indexdf[!, Findex] = skippingcumret(df[!, Fret])
  end
  @assert issetequal(propertynames(indexdf), [:date; Findexes])

  #WARNING the compelete cases is important here
  indexdf = indexdf[completecases(indexdf[!, Findexes]),:]

  return Findexes,indexdf
end

#merge in returns from a file of indexes
#1) merge indexdf with unique dates
#2) generate returns
#3) merge the returns with the provided df

#WARNING WARNING WARNING the join here does not necessitate a complete match between
#the df and the index- callers responsibility to check index dates are comprehensive
#WARNING WARNING WARNING order of rows is undefined
function index2returns!(df::DataFrame, indexdf::AbstractDataFrame;
  Findexes::Vector{Symbol} = throw("Findexes is required"),
  Frets::Vector{Symbol} = throw("Frets is required"),
  datefrequency::Symbol = throw("date frequency is required (:day, :month, :week)"),
  testlabel::String = throw("testlabel is required"),
  Fgroup::NSymbol = nothing, #only used for testing purposes if provided
  validateexactfrequency = false)

  if Fgroup === nothing
    allunique(df.date) || throw("Fgroup must be provided for testing purposes
      if dates in df are not unique")
  end

  @assert allunique(indexdf.date)
  @assert issorted(indexdf, :date)
  @assert all(indexdf.date .!== missing)
  @assert all(df.date .!== missing)
  if :indexdate ∈ propertynames(df)
    throw("Column indexdate should not be in df. Call select!(bm, Not(:indexdate))
      to delete the column prior to calling index2returns")
  end


  #join on the dates, treating months as a special case
  bm = DataFrame(date=sort!(unique(df.date)))
  if datefrequency === :month
    indexdf.dateym = (d->year(d) + month(d) / 100).(indexdf.date)
    indexdf = combine(last, groupby(indexdf, :dateym))
    bm.dateym = (d->year(d) + month(d) / 100).(bm.date)

    #avoid duplicate columns in the join
    rename!(indexdf, :date=>:indexdate)

    #throwing the error here from funds-mf:Merge key(s) in df1 are not unique
    bm = leftjoin(bm, indexdf, on=:dateym, validate=(true, true))
    sort!(bm, [:dateym])

    #the following creates a crosswalk for testing and validation
    indexdf = leftjoin(indexdf, DataFrame(date=bm.date, dateym = bm.dateym), on=:dateym)
    sort!(indexdf, :dateym)
  elseif datefrequency ∈ [:day, :week]
    bm = leftjoin(bm, indexdf, on=:date, validate=(true, true))
    sort!(bm, :date)
  elseif datefrequency === :year
    throw("year date frequency not supported (see month for an idea how to do this)")
  else
    throw("unrecognized date frequency $datefrequency")
  end

  #bm |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_bmlab1.csv")
  #@info "Wrote bm1. datefrequency=$datefrequency"
  #check for missing values
  #=for Findex ∈ Findexes
    indicesmissing = sum(bm[!, Findex] .≡ missing)
    if indicesmissing > 0
      @warn ("BM dates found with no corresponding index. Testlabel: $testlabel, " *
        "Field: $Findex, nummissing=$indicesmissing
        Check testfiles ending with $(testlabel)_ind2ret and $(testlabel)_indmerged to" *
        "get some details.")
    end
  end=#

  #another check on the above join check on the above join
  validatedates(bm.date |> Vector{Date}, frequency=datefrequency; validateexactfrequency)


  #validate dates by comparing the sets of the benchmark and the df,
  #ignoring mismatches beyond the endpoints
  indexdates = indexdf[completecases(indexdf, [Findexes; :date]), :date]
  dfnotindex = setdiff(bm.date, indexdates)
  dfnotindex = (!isempty(dfnotindex))  ? dfnotindex[(minimum(indexdates) .≤ dfnotindex) .&
    (dfnotindex .≤ maximum(indexdates))] : dfnotindex

  #=indexnotdf = setdiff(indexdates, bm.date)
  indexnotdf = (!isempty(indexnotdf)) ? indexnotdf[(minimum(bm.date) .≤ indexnotdf) .&
    (indexnotdf .≤ maximum(bm.date))] : indexnotdf=#

  (!isempty(dfnotindex)) && (@warn("df dates not in index (interior only): $dfnotindex" ))
  #(!isempty(indexnotdf)) && (@warn("index dates not in df (interior only): $indexnotdf") )

  for (Fret, Findex) ∈ zip(Frets, Findexes)
    bm[!, Fret] = [missing; bm[2:end, Findex] ./ bm[1:(end-1), Findex] .- 1]
  end

  #join the index with the panel data
  df = leftjoin(df, bm, on=:date, validate=(false, true))


  #write output for testing purposes
  if PARAM[:testindex2returns]
    bm |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(testlabel)_ind2ret.csv")
    local testview::SubDataFrame
    if Fgroup ≡ nothing
      testview = view(df, :, :)
    elseif Fgroup ≡ :permno
      testpermno::Int = PARAM[:testpermnomult]
      testview = view(df, (df[!, Fgroup] .% testpermno) .== 0, :)
    elseif Fgroup ≡ :fundid
      testids = values(PARAM[:fundbetatestfundid])
      testview = view(df, (g->g ∈ testids).(df[!, Fgroup]), :)
    end
    testview |> CSV.write(
      "$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(testlabel)_indmerged.csv")
  end

  if :dateym ∈ propertynames(df)
    select!(df, Not([:dateym]))
  end
  select!(df, Not(Findexes))
  return df
end

#merge in returns from the returns file
#core of code is in index2returns- just calls a generic function to generate the indexdf
#WARNING WARNING WARNING the joins here do not necessitate a complete match between
#the df and the index- callers responsibility to check index dates are comprehensive
#WARNING WARNING WARNING order of rows is undefined
function mergereturns(df::DataFrame, returns::AbstractDataFrame;
  Frets::Vector{Symbol} = throw("Frets is required"),
  datefrequency::Symbol = throw("date frequency is required (:day, :month, :week)"),
  testlabel::String = throw("testlabel is required"),
  Fgroup::NSymbol = nothing, #only used for testing purposes if provided
  validateexactfrequency = false)

  #form an index for each return stream
  Findexes, indexdf = returns2index(returns; Frets)

  #returns |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_returnstemp.csv")
  #indexdf |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_indextemp.csv")

  #WARNING WARNING WARNING order of rows is undefined
  return index2returns!(df, indexdf;
    Frets, Findexes, datefrequency, testlabel, Fgroup, validateexactfrequency)
end
