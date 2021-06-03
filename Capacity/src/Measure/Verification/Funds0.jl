
#load and process funds as desired
function prepfunds(returns::DataFrame,
    ms::Union{Nothing, MeasureSpec}=nothing,
    ;refreshfunds::Bool = PARAM[:refreshfunds],
    refreshcomom::Bool = PARAM[:refreshcomom],
    workingpath::String = PARAM[:workingpath],
    loadpath::NString = nothing,
    hfrname::String = PARAM[:hfrname],
    mfname::String = PARAM[:mfname],
    returntype::String = PARAM[:crspdatalabel],
    plname::String = PARAM[:plname],
    weeklycrsptype::String = string(PARAM[:crspweekly]),
    dtformat = PARAM[:wrdsdateformat])

  local hf::DataFrame
  local mf::DataFrame
  local comom::DataFrame

  if refreshfunds

    hf = prephf(returns::DataFrame, ms::MeasureSpec)
    hf |> CSV.write("$workingpath\\$(returntype)_$(hfrname).csv")

    mf = prepmf(returns::DataFrame, ms::MeasureSpec)
    mf |> CSV.write("$workingpath\\$(returntype)_$(mfname).csv")
  else
    hf = CSV.File("$workingpath\\$(returntype)_$(hfrname).csv") |> DataFrame
    mf = CSV.File("$workingpath\\$(returntype)_$(mfname).csv") |> DataFrame
  end
  if refreshcomom
    if PARAM[:crspfrequencysuffix] ≡ "m"
      comom = generatemonthlycomomentum()
      #comom = generatecomomentumbb()
    elseif PARAM[:crspfrequencysuffix] ≡ "w"
      comom = generateweeklycomomentum()
    else
      @assert false
    end
    comom |> CSV.write("$workingpath\\$(returntype)_$(plname).csv")
  else
    comom = CSV.File("$workingpath\\$(returntype)_$(plname).csv") |> DataFrame
  end

  return (;hf, mf, comom)
end

#=function betaweightsslice(df, extrabetaweights::Vector{Symbol}=PARAM[:fundextrabetaweights],
  Fbetacontrols=PARAM[:fundbetacontrols],
  betacorridorradius=PARAM[:fundbetacorridorradius])

  throw("THIS NEEDS TO BE COMPLETED")
  #this will hold the weights by fundid
  #the below version must work with missing data
  weights = DataFrame(fundid=unique(df.fundid))
  function betaw(t,Frhs)
    completeindices = sum()
    X = reduce(hcat, (F->t[F]).(Frhs)

    (svd(Xclean))\(t.netret))[1]
  end
  local Fweights = Vector{Symbol}()
  #now execute the regressions
  for Fbetaweight ∈ Fbetaweights
    Fweight = Symbol(:w_, Fbetaweight)
    Fabsweight = Symbol(:absw_, Fbetaweight)
    Fabscorweight = Symbol(:abscorw_, Fbetaweight)
    newweightfields = [Fweight,Fabsweight,Fabscorweight]

    #form a view with no missing values
    focalfields = [:date; :fundid; :some1s; Fbetaweight; :netret; Fbetacontrols]
    Frhs = [Fbetaweight; Fbetacontrols; :some1s;]
    cleandf = view(df, completecases(df, [Fbetaweight; Fbetacontrols; :netret]), focalfields)

    #run the regression for each fundid
    #we are primarily interested in the absolute vlaue of the weights
    #the corrodor version creates a buffer to guard against measurement error in the controls
    weight = combine(groupby(cleandf, :fundid), (:) |> AsTable =>
      ((t)->betaw(t,Frhs)) => Fweight)
    weight[!, Fabsweight] = abs.(weight[!, Fweight])
    weight[!, Fabscorweight] = (x->ifelse(x<betacorridorradius,0.0, x)).(weight[!, Fabsweight])
    weights = innerjoin(weights, weight[:, [:fundid; newweightfields]], on=[:fundid])
    append!(Fweights, newweightfields)
  end

end=#

function betaweights(df, returns, ms::MeasureSpec;
  extrabetaweights::Vector{Symbol}=PARAM[:fundextrabetaweights],
  Fbetacontrols=PARAM[:fundbetacontrols],
  betacorridorradius=PARAM[:fundbetacorridorradius],
  datefrequency::Symbol = throw("datefrequency is required"),
  testlabel=throw("testlabel is required"),
  Frfr::Symbol = PARAM[:fundrfr],
  )

  @assert issorted(df, [:fundid, :date])

  Fbetaweights::Vector{Symbol} = [ms.Fξs; extrabetaweights;]

  #just get the parts we want to merge from the returns dataframe
      #crashed here in the returns file
  #throwing the error here from funds-mf: Merge key(s) in df1 are not unique
  df = mergereturns(df, returns,
    Fgroup = :fundid,
    Frets = [Fbetaweights;Fbetacontrols;Frfr],
    validateexactfrequency=true;
    datefrequency, testlabel)

  #=WARNING- the below code uses all available data, but will use different fund-date
  #universes for different beta controls. I leave it disabled because its also very slow
  #also it shouldn't mkae a big difference
  for Fret ∈ [Fbetaweights;Fbetacontrols]
    if :indexdate ∈ propertynames(df)
      select!(df, Not([:indexdate]))
    end
    df = mergereturns(df, returns,
      Fgroup = :fundid,
      Frets = [Fret;],
      validateexactfrequency=true;
      datefrequency, testlabel)
  end=#

  #this is the intercept
  df.some1s = ones(nrow(df))
  df.netret = df.ret .- df[!, Frfr]

  #this will hold the weights by fundid
  weights = DataFrame(fundid=unique(df.fundid))

  #runs the regression
  #form the LHS matrix then
  #betaw(t,Fbetaweight) = (svd([t[Fbetaweight] t.some1s])\(t.ret))[1]
  betaw(t,Frhs) = (svd(reduce(hcat, (F->t[F]).(Frhs)))\(t.netret))[1]

  local Fweights = Vector{Symbol}()
  #now execute the regressions
  for Fbetaweight ∈ Fbetaweights
    Fweight = Symbol(:w_, Fbetaweight)
    Fabsweight = Symbol(:absw_, Fbetaweight)
    Fabscorweight = Symbol(:abscorw_, Fbetaweight)
    newweightfields = [Fweight,Fabsweight,Fabscorweight]

    #form a view with no missing values
    focalfields = [:date; :fundid; :some1s; Fbetaweight; :netret; Fbetacontrols]
    Frhs = [Fbetaweight; Fbetacontrols; :some1s;]
    cleandf = view(df, completecases(df, [Fbetaweight; Fbetacontrols; :netret]), focalfields)

    #run the regression for each fundid
    #we are primarily interested in the absolute vlaue of the weights
    #the corrodor version creates a buffer to guard against measurement error in the controls
    weight = combine(groupby(cleandf, :fundid), (:) |> AsTable =>
      ((t)->betaw(t,Frhs)) => Fweight)
    weight[!, Fabsweight] = abs.(weight[!, Fweight])
    weight[!, Fabscorweight] = (x->ifelse(x<betacorridorradius,0.0, x)).(weight[!, Fabsweight])
    weights = innerjoin(weights, weight[:, [:fundid; newweightfields]], on=[:fundid])
    append!(Fweights, newweightfields)
  end

  #merge in the new weight fields
  df = leftjoin(df, weights, on=:fundid, validate=(false, true))
  if !issorted(df, [:fundid, :date])
    sort!(df, [:fundid, :date])
  end
  select!(df, Not([:some1s]))

  for Fweight ∈ Fweights
    #=winsorizewithin!(df, Fsource=Fweight, prop=PARAM[:fundwinsorizebeta], Fgroup=:date,
      twosided=true, Fnew=Symbol(:temp, Fweight))=#
    df[!, Symbol(:temp, Fweight)] = winsorizequantile(df[!,Fweight],
      PARAM[:fundwinsorizebeta]; twosided=true)
  end




  if PARAM[:testfundbeta]
    testmonth = PARAM[:fundbetatestmonth]
    testview::SubDataFrame = view(df, (month.(df.date) .== testmonth.month) .&
      (year.(df.date) .== testmonth.year), :)
    testview[1:(min(20_000, nrow(testview))),:] |> CSV.write(
      "$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(Tuple(testmonth))_$testlabel.csv")

    testfundid = PARAM[:fundbetatestfundid][testlabel |> Symbol]
    testview = view(df, df.fundid .== testfundid, :)
    testview[1:(min(20_000, nrow(testview))),:] |> CSV.write(
      "$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_$(testfundid)_$testlabel.csv")
  end


  #@info propertynames(df)
  for Fweight ∈ Fweights
    select!(df, Not([Fweight]))

    rename!(df, Dict(Symbol(:temp, Fweight)=>Fweight))
  end

  return Fweights, df
end

hfrbetaweights(hfr, returns, ms::MeasureSpec) = betaweights(
    hfr, returns, ms, testlabel="hfr", datefrequency=:month)

function prephf(returns::AbstractDataFrame, ms::MeasureSpec;
  hfrnamealive::String = PARAM[:hfrnamealive],
  hfrnamedead::String = PARAM[:hfrnamedead],
  hfrpath::String = PARAM[:hfrpath],
  minfundmonths::Int = PARAM[:fundminfundmonths],
  dropfinalmonths::Int = PARAM[:hfrdropfinalmonths],
  lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
  incsvstream::Function = IN_CSV_STREAM,
  csvextension::String = CSV_EXTENSION)::DataFrame

  hfralive::DataFrame = incsvstream("$hfrpath\\$hfrnamealive.$csvextension") |> CSV.File |> DataFrame
  hfralive.alive = trues(size(hfralive,1))
  hfrdead::DataFrame = incsvstream("$hfrpath\\$hfrnamedead.$csvextension") |> CSV.File |> DataFrame
  hfrdead.alive = falses(size(hfrdead,1))

  hfr = vcat(hfralive, hfrdead)

  normalizename(s::String) = replace(s, "_"=>"") |> lowercase
  rename!(hfr, normalizename.(names(hfr)))

  #need to parse the dates
  for datefield ∈ [:date, :inception, :fundassetsasof, :firmassetsasof]
    hfr[!, datefield] = (d->parsehfr(Date, d)).(hfr[!, datefield])
  end

  sort!(hfr, [:fundid, :date,])
  hfr.firstdomdate = hfr.date .|> firstdayofmonth
  @assert sum(nonunique(hfr, [:fundid, :firstdomdate]))==0

  #classify the strategies- first read in a spreadsheet of the classifications
  stratcodes = CSV.File("$hfrpath\\classifications.csv") |>  DataFrame
  rename!(stratcodes, normalizename.(names(stratcodes)))
  stratcodes=stratcodes[:, Not([:mainstrategy, :substrategy])]

  #conform index classifications
  hfr.hfri = hfr.inhfri .=== "Yes"
  hfr.hfrx = hfr.inhfrx .=== "Yes"

  #then join the classifications with the funds
  hfr = innerjoin(hfr, stratcodes,on=[:strategycode, :substrategycode])
  if !issorted(hfr, [:fundid, :date])
    sort!(hfr, [:fundid, :date,])
  end
  @assert size(hfr,1) == size(hfrdead,1) + size(hfralive,1) #sanity check

  #keep only the relevant dates- is this the right approach? should I take the max range?
  #keep only dates within the return range
  #also drop the final month(s) due to incomplete data
  mindate=minimum(returns.date)|>firstdayofmonth
  maxdate = min(maximum(returns.date) |> lastdayofmonth,
    maximum(hfr.date) - Month(dropfinalmonths) |> lastdayofmonth)
  hfr = hfr[(mindate .< hfr.date) .& (hfr.date .< maxdate), :]
  hfr.ret = hfr.performance ./ 100
  select!(hfr, Not(:performance))

  #now compute AUMs by strategy bunches
  hfr = hfr[hfr.fundassetsdenomin .== "USD",:]
  hfr.assets ./= 1_000
  rename!(hfr, :assets=>:aum) #what about nav?

  #acquire the index of strategy returns
  hfr = hfr[completecases(hfr, [:aum, :ret]), :]
  transform!(groupby(hfr, :fundid), nrow => :fundmonths)
  hfr = hfr[hfr.fundmonths .≥ minfundmonths, :]

  println("Attempting to weight hfr fund aum by returns...")
  Fhfweights, hfr = hfrbetaweights(hfr, returns, ms)
  println("weighting successful.")

  #sums values in an array if the conditional is true and nothing is missing
  dotifs(v1,v2 = 1.0, v3 = 1.0) = sum((v1 .* v2 .* v3) .|> (x)->coalesce(x, 0.0))

  #now group the assets
  shfs = groupby(hfr, :date)
  hf = combine(shfs ,
    :aum => dotifs => :aum,
    [:aum, :equityonly] => dotifs => :A_aumequityonly,
    [:aum, :equity] => dotifs => :A_aumequity,
    [:aum, :quant] => dotifs => :A_aumquant,
    [:aum, :hfri] => dotifs => :A_aumhfri,
    [:aum, :hfrx] => dotifs => :A_aumhfrx,
  )
  sort!(hf, :date)

  #group the regression weights
  for Fhfweight ∈ Fhfweights
    hfweightdf = combine(shfs,
      [:aum, Fhfweight] => dotifs => Symbol(:A_,Fhfweight),
      [:aum, Fhfweight, :equity] => dotifs => Symbol(:A_aumequityX,Fhfweight),
      [:aum, Fhfweight, :equityonly] => dotifs => Symbol(:A_aumequityonlyX,Fhfweight),
    )
    hf = innerjoin(hf, hfweightdf, on=:date, validate=(true,true))
  end
  sort!(hf, :date)

  #create the log growwth columns
  logmaybe(::Missing) = missing
  logmaybe(x) = x ≤ 0 ? missing : log(x)
  for acol ∈ findcols(hf, :asset)
    @assert occursin(r"^A", string(acol))
    hf[!, Symbol(:lG_,acol)] = [missing; logmaybe.(hf[2:end, acol] ./ hf[1:(end-1), acol])]
    #hf[!, Symbol(:G_,acol)] = [missing; hf[2:end, acol] ./ hf[1:(end-1), acol] .- 1]
  end
  #average all cols via the LP filter
  hf = lpfilter(hf, focalcols=setdiff(propertynames(hf), [:date]); lpfilterperiods)

  #create normalized versions NOTE- I don't think this is the right approach- disabled
  #create normalized versions
  #=focalcols = setdiff(propertynames(hf), [:date, :aum])
  for f ∈ focalcols
    transform!(hf, [f, :aum] => ((x,y)->x ./ y) =>Symbol(:N, f))
  end

  for per ∈ lpfilterperiods
    (per === nothing) && continue

    #the suffix is how we will identify the lp columns during analysis
    hf = lpfilter(hf, lpfilterperiods=[per], focalcols=[focalcols; :aum])
    lpcols = focalcols .|> (f)->Symbol(f,lpsuffix(per))
    aumcol = Symbol(:aum, lpsuffix(per))
    lpcols = (s->Symbol(s, lpsuffix(per))).(focalcols)
    for lpcol ∈ lpcols
      transform!(hf, [lpcol, aumcol] => ((x,y)->x ./ y) =>Symbol(:N, lpcol),)
    end
  end=#

  return hf
end

parsehfr(::Type{Date}, s::String, hfrdateformat::DateFormat=dateformat"yyyymmdd"
  ) = Dates.Date(s, hfrdateformat)
parsehfr(::Type{<:Any}, ::Missing) = missing
parsehfr(::Type{Date}, i::Int, hfrdateformat::DateFormat=dateformat"yyyymmdd") = parsehfr(
  Date, "$i", hfrdateformat)


############################MF

mfbetaweights(mf, returns, ms::MeasureSpec;
  datefrequency=:month) = betaweights(
    mf, returns, ms;  testlabel="mf", datefrequency)

function prepmf(returns::AbstractDataFrame, ms::MeasureSpec;
  mfdatefrequency::Symbol=PARAM[:mfdatefrequency],
  mfname::String = PARAM[:mfname],
  mfnamedesc::String = PARAM[:mfnamedesc],
  mfpath::String = PARAM[:mfpath],
  minfundmonths::Int = PARAM[:fundminfundmonths],
  minfunds::Int = PARAM[:fundminfunds],
  lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
  incsvstream::Function = IN_CSV_STREAM,
  csvextension::String = CSV_EXTENSION,
  dtformat::DateFormat = PARAM[:wrdsdateformat])

  #scrubs the names and dates
  function initialclean(df::DataFrame; #=startdate::Date = PARAM[:mfstartdate]=#)
    rename!(df, (cleanname).(names(df)))
    df.date = (i::Int->Date("$i", dtformat)).(df.caldt)
    #df.year = (d->year(d)).(df.date)
    transform!(df, :crspfundno => ByRow((s)->Symbol(:mf, s)) => :fundid)
    select!(df, Not([:crspfundno, :caldt]))

    return df
  end

  #start by working with the description file
  mfdesc::DataFrame = CSV.File("$mfpath\\$mfnamedesc.csv") |> DataFrame
  mfdesc = initialclean(mfdesc)
  sort!(mfdesc, [:fundid, :date])
  #mfdesc = combine(groupby(mfdesc, [:fundid, :year]), propertynames(mfdesc) .=>
  #  (v)->v[end], renamecols=false)
  #dups = mfdesc[nonunique(mfdesc[:, [:fundid, :year]]),:]
  #dups.id = 1:nrow(dups)
  ##println(innerjoin(mfdesc, dups[:, [:fundid, :year, :id]], on=[:fundid, :year]))
  #@assert sum(nonunique(mfdesc[:, [:fundid, :year]])) == 0

  mfdesc.lipperclassname .= (s::MString->
    ismissing(s) ? missing : lowercase(s)).(mfdesc.lipperclassname)
  mfdesc.lipperobjname .= (s::MString->
    ismissing(s) ? missing : lowercase(s)).(mfdesc.lipperobjname)

  select!(mfdesc, Not([:date]))

  #note- the classifications are by hand- the below is ONLY for support of the by-hand classifications
  unique(mfdesc[:, [:lipperclassname, :lipperobjname]]) |> CSV.write("mfclassificationsraw.csv")

  #only interested in the most recent, non-missing row for each fund
  mfdesc = mfdesc[(mfdesc.lipperclassname .!== missing) .| (mfdesc.lipperobjname .!== missing), :]
  mfdesc = combine(last, groupby(mfdesc, :fundid), renamecols=false)
  @assert sum(nonunique(mfdesc[:, [:fundid]])) == 0

  #get the classifications from a seperate file
  stratcodes = CSV.File("$mfpath\\mfclassifications.csv") |>  DataFrame
  select!(stratcodes, [:lipperclassname, :lipperobjname, :equity,:equityonly, :domesticeo])
  oldnrows = nrow(mfdesc)
  mfdesc = innerjoin(mfdesc, stratcodes, on=[:lipperclassname, :lipperobjname], matchmissing=:equal)
  @assert nrow(mfdesc) == oldnrows
  @assert sum(nonunique(mfdesc[:, [:fundid]])) == 0

  #clean the data
  mf::DataFrame = incsvstream("$mfpath\\$mfname.$csvextension") |> CSV.File |> DataFrame
  mf = initialclean(mf)


  sort!(mf, [:fundid, :date])
  @assert sum(nonunique(mf[!, [:fundid, :date]])) == 0

  #parse relevant fields to a number- letter codes generally indicate missing
  for (Fold, Fnew) ∈ [:mret=>:ret, :mtna=>:aum, :mnav=>:nav]
    mf[!, Fnew] = (s->parsemf(Float64, s)).(mf[!, Fold])
  end

  #invalid values also generally represent missing
  mf.ret .= mf.ret .|> (x)-> (x ≡ missing ) || (x < -1.0) ? missing : x
  mf.aum .= mf.aum .|> (x)-> (x ≡ missing ) || (x < 0.0) ? missing : x
  mf.aum ./= 1000

  mf = dropmissing!(mf, :ret)

  #implmeent the reversal check of PST2015
  lagwithin2!(mf, [:aum], :fundid, date=:date)
  leadwithin2!(mf, [:aum], :fundid, date=:date)
  mf.Gaum = (mf.aum .- mf.Laum) ./ mf.Laum
  mf.rev = (mf.Naum .- mf.aum) ./ (mf.aum .- mf.Laum)
  transform!(groupby(mf, :fundid), [:Gaum, :rev, :Laum, :aum] => ByRow(
    (Gaum, rev, Laum, aum)->ifelse(
      (Gaum !== missing) && (rev !== missing) &&
      (abs(Gaum) > 0.5) && (-1.25 ≤ rev ≤ -0.75) && (Laum > 0.01), missing, aum))
    => :adjaum)

  @info "Dropped $(sum(mf.aum .!== mf.adjaum)) month-fund rows due to rapid  aum reversals"
  select!(mf, Not([:Gaum, :rev, :Laum, :aum, :Naum]))
  rename!(mf, :adjaum => :aum)



  mf = patchmfmonthlyreturngaps(mf)

  #no point in keeping a row if the assets field is missing
  mf = dropmissing!(mf, :aum)
  transform!(groupby(mf, :fundid), nrow => :fundmonths)
  mf = mf[mf.fundmonths .≥ minfundmonths, :]

  #trim down to the desired range with a buffer
  startdate = minimum(returns.date) - Year(1)
  enddate =  maximum(returns.date) + Year(1)
  mf = mf[(mf.date .≥ startdate) .&  (mf.date .≤ enddate), :]

  #if we want to patch in daily or weekly return data, this is where we do that
  select!(mf, [:fundid, :date, :aum, :ret])
  mf = patchmfreturns(mf, mfdatefrequency)

  #patch in the description info
  mf = leftjoin(mf, mfdesc[:, [:fundid, :equity, :equityonly, :domesticeo,
    :lipperclassname, :lipperobjname]], on=[:fundid])
  if !issorted(mf, [:fundid, :date])
    sort!(mf, [:fundid, :date])
  end

  transform!(groupby(mf, :date), nrow => :numfunds)
  mf = mf[mf.numfunds .≥  minfunds,:]

  #merge in the beta-based fund weights
  println("Attempting to weight mf fund aum by returns...")
  Fmfweights, mf = mfbetaweights(mf, returns, ms; datefrequency=mfdatefrequency)
  println("...weighting successful.")

  #sums values in an array if the conditional is true and nothing is missing
  zeromissings(x) =  coalesce(x, 0.0)
  dotifs(v1,v2 = 1.0, v3 = 1.0) = sum((v1 .* v2 .* v3) .|> zeromissings)
  desccount(aum,beta = 1.0, cond = 1.0) = sum((aum .* cond .* (beta .!==missing)) .|> zeromissings .|> (x)->x>0.0)
  descmin(beta,cond=true) = minimum((beta .* cond) |> skipmissing)
  descmax(beta,cond=true) = maximum((beta .* cond) |> skipmissing)
  descstd(beta,cond=true) = std((beta .* cond) |> skipmissing)
  descmedianabs(beta,cond=true) = median((beta .* cond) .|> abs |> skipmissing)
  descmeanabs(beta,cond=true) = mean((beta .* cond) .|> abs |> skipmissing)


  combine(groupby(view(mf, mf.date .== Date(2019,12,27),:), [:lipperclassname, :lipperobjname]),
    Fmfweights .=> ((v)->mean(skipmissing(abs.(v)))) .=> (s->Symbol(:meanabs_,s)).(Fmfweights),
    :aum => (v)->sum(skipmissing(v)) => :aum,
    zip([:aum for i ∈ 1:length(Fmfweights)], Fmfweights) .|> collect .=> dotifs .=> (
      s->Symbol(:A_,s)).(Fmfweights),
    :equity => last => :equity,
    :equityonly => last => :equityonly,
    :domesticeo => last => :domesticeo) |> CSV.write("$(PARAM[:testpath])\\mf_asset aggregates 20191227.csv")

  throw("stop")

  #now group the assets
  smfs = groupby(mf, :date)
  funds = combine(smfs,
    :aum => dotifs => :aum,
    [:aum, :equity] => dotifs => :aumequity,
    [:aum, :equityonly] => dotifs => :aumequityonly,
    [:aum, :domesticeo] => dotifs => :aumdomesticeo,
  )
  if !issorted(funds, :date)
    sort!(funds, :date)
  end

  #group the regression weights
  for Fmfweight ∈ Fmfweights
    oldnames = propertynames(funds)
    mfweightdf = combine(smfs,
      [:aum, Fmfweight] => dotifs => Symbol(:A_, Fmfweight),
      [:aum, Fmfweight, :equity] => dotifs => Symbol(:A_aumequityX,Fmfweight),
      [:aum, Fmfweight, :equityonly] => dotifs => Symbol(:A_aumequityonlyX,Fmfweight),
      [:aum, Fmfweight, :domesticeo] => dotifs => Symbol(:A_aumdomesticeoX,Fmfweight),


      [:aum, Fmfweight] => desccount => Symbol(:descn_, Fmfweight),
      [Fmfweight] => descmedianabs => Symbol(:descmedabs_, Fmfweight),
      [Fmfweight] => descmeanabs => Symbol(:descmeanabs_, Fmfweight),
      [Fmfweight] => descstd => Symbol(:descstd_, Fmfweight),
      [Fmfweight] => descmin => Symbol(:descmin_, Fmfweight),
      [Fmfweight] => descmax => Symbol(:descmax_, Fmfweight),


      [:aum, Fmfweight, :equity] => desccount => Symbol(:descn_aumequityX,Fmfweight),
      [Fmfweight, :equity] => descmedianabs => Symbol(:descmedabs_aumequityX,Fmfweight),
      [Fmfweight, :equity] => descmeanabs => Symbol(:descmeanabs_aumequityX,Fmfweight),
      [Fmfweight, :equity] => descstd => Symbol(:descstd_aumequityX,Fmfweight),
      [Fmfweight, :equity] => descmin => Symbol(:descmin_aumequityX,Fmfweight),
      [Fmfweight, :equity] => descmax => Symbol(:descmax_aumequityX,Fmfweight),


      [:aum, Fmfweight, :equityonly] => desccount => Symbol(:descn_aumequityonlyX,Fmfweight),
      [Fmfweight, :equityonly] => descmedianabs => Symbol(:descmedabs_aumequityonlyX,Fmfweight),
      [Fmfweight, :equityonly] => descmeanabs => Symbol(:descmeanabs_aumequityonlyX,Fmfweight),
      [Fmfweight, :equityonly] => descstd => Symbol(:descstd_aumequityonlyX,Fmfweight),
      [Fmfweight, :equityonly] => descmin => Symbol(:descmin_aumequityonlyX,Fmfweight),
      [Fmfweight, :equityonly] => descmax => Symbol(:descmax_aumequityonlyX,Fmfweight),

      [:aum, Fmfweight, :domesticeo] => desccount => Symbol(:descn_aumdomesticeoX,Fmfweight),
      [Fmfweight, :domesticeo] => descmedianabs => Symbol(:descmedabs_aumdomesticeoX,Fmfweight),
      [Fmfweight, :domesticeo] => descmeanabs => Symbol(:descmeanabs_aumdomesticeoX,Fmfweight),
      [Fmfweight, :domesticeo] => descstd => Symbol(:descstd_aumdomesticeoX,Fmfweight),
      [Fmfweight, :domesticeo] => descmin => Symbol(:descmin_aumdomesticeoX,Fmfweight),
      [Fmfweight, :domesticeo] => descmax => Symbol(:descmax_aumdomesticeoX,Fmfweight),

    )
    funds = innerjoin(funds, mfweightdf, on=:date, validate=(true,true))
  end
  sort!(funds, :date)

  logmaybe(::Missing) = missing
  logmaybe(x) = x ≤ 0 ? missing : log(x)
  lpfiltercols = findcols(funds, :asset)
  for acol ∈ findcols(funds, :asset)
    @assert occursin(r"^A", string(acol))

    push!(lpfiltercols, Symbol(:lG_,acol))
    funds[!, lpfiltercols[end]] = [missing; logmaybe.(funds[2:end, acol] ./ funds[1:(end-1), acol])]
    #funds[!, Symbol(:G_,acol)] = [missing; funds[2:end, acol] ./ funds[1:(end-1), acol] .- 1]
  end

  funds = lpfilter(funds, focalcols=lpfiltercols; lpfilterperiods)



  return funds
end

function patchmfmonthlyreturngaps(mf::AbstractDataFrame,
    maxpatchlength=PARAM[:mftimetostale])

  #sanity checks to make sure everything is as expected
  @assert issorted(mf, [:fundid, :date])
  @assert sum(completecases(mf, [:ret])) == nrow(mf)

  #temporary helper fields
  mf.ret1 = mf.ret .+ 1
  mf.maxvaliddate = (mf.date .+ maxpatchlength) .|> lastdayofmonth

  function patchwithret(dates, maxdate, aum, ret)

    sum(aum .≡ missing) == length(aum) && return aum
    T = searchsortedfirst(dates, maximum(dates[aum .!== missing]))

    maxvaliddate = Date(1900,1,1)
    for t ∈ 1:T

      #if the aum is not missing, then record the max valid date
      if aum[t] !== missing
        maxvaliddate = maxdate[t]
        continue
      end

      if dates[t] ≤ maxvaliddate
        aum[t] = (1+ret[t]) * aum[t-1]
      end
    end

    return aum
  end

  #correct the dates
  mf.aumcorrected = mf.aum |> deepcopy
  transform!(groupby(mf, :fundid),
    [:date, :maxvaliddate, :aumcorrected, :ret] => patchwithret => :aumcorrected)
  mf.corrected = (mf.aum .≡ missing ) .& (mf.aumcorrected .!== missing)

  #check
  uncorrectedcomplete = view(mf, completecases(mf, :aum), :)
  @assert (uncorrectedcomplete.aum .== uncorrectedcomplete.aumcorrected) |> all

  testfunds = mf[mf.corrected, :fundid] |> unique |> (v)->v[1:3]
  mf[
    (mf.fundid .== testfunds[1] ) .| (mf.fundid .== testfunds[2] ) .| (mf.fundid .== testfunds[3])
    , :] |> CSV.write("$(PARAM[:testpath])\\mf_monthlyaumpatch.csv")

  mf.aum = mf.aumcorrected |> deepcopy
  select!(mf, Not([:aumcorrected, :ret1, :maxvaliddate]))
  return mf
end






function patchmfreturns(mf::DataFrame, mfdatefrequency::Symbol;
    mfdailyname::String = PARAM[:mfdailyname],
    mfpath::String = PARAM[:mfpath],
    incsvstream::Function = IN_CSV_STREAM,
    csvextension::String = CSV_EXTENSION,
    dtformat = PARAM[:wrdsdateformat],
    interpolate::Bool = PARAM[:mfdailyinterpolate],
  )

  if mfdatefrequency ≡ :month
    dropmissing!(mf, [:ret])
    winsorizewithin!(mf, Fsource=:ret, prop=PARAM[:mfqwinsorize], Fgroup=:date,
      twosided=true, Fnew=:WXret)
    select!(mf, Not([:ret]))
    rename!(mf, :WXret=>:ret)
    return mf #nothing to patch in this case
  end

  #we may be able to drop this assertion if it is checked later
  @assert mfdatefrequency ∈ [:day, :week]
  oldmonthlyaum = mf[!, [:date, :fundid, :aum]]
  rename!(oldmonthlyaum, :aum=>:aumoldmonthly)

  #load and clean daily data
  mfd::DataFrame = incsvstream(
    "$mfpath\\$mfdailyname.$csvextension") |> CSV.File |> DataFrame
  rename!(mfd, (cleanname).(names(mfd)))
  mfd.date = (i::Int->Date("$i", dtformat)).(mfd.caldt)
  mfd.yearmonth = mfd.date .|> (d->year(d) + month(d)/100)
  #df.year = (d->year(d)).(df.date)
  transform!(mfd, :crspfundno => ByRow((s)->Symbol(:mf, s)) => :fundid)
  mfd.ret = mfd.dret .|> (r-> parsemf(Float64, r))
  dropmissing!(mfd, [:fundid, :date, :ret])

  #quick and easy performance enhancement based on the time period covered
  mfd = mfd[
    (mfd.date .≥ firstdayofmonth(minimum(mf.date))) .&
    (mfd.date .≤ lastdayofmonth(maximum(mf.date))),:]

  #prep mf
  @assert issorted(mf, [:fundid, :date])

  #we want the beginning of period AUM, not the end of period
  mfN = deepcopy(mf)
  mfN.Nyearmonth = mf.date .+ Month(1) .|> (d->year(d) + month(d)/100)
  mfN = select(mfN, [:fundid, :Nyearmonth, :aum])

  #we want to use the eom aum onthe last day of the month,
  #and then use returns to interpolate at other parts of the month
  mf.yearmonth = mf.date .|> (d->year(d) + month(d)/100)
  mf = select!(mf, [:fundid, :yearmonth, :aum])
  rename!(mf, :aum=>:aumeom)
  mfN = leftjoin(mfN, mf, on=[:fundid, :Nyearmonth => :yearmonth])


  #NOTE: the below restriction ensures that only the AUM is shifted as desired
  #WARNING: read above before adding columns to below

  select!(mfd, [:fundid, :date, :yearmonth, :ret])
  mf = innerjoin(mfN, mfd, on=[:fundid, :Nyearmonth => :yearmonth], validate=(true, false))

  #validate df and make sure everything is as expected
  @assert nonunique(mf, [:fundid, :date]) |> !all

  mf = sort(mf, [:fundid, :date])
  @assert issorted(mf, [:fundid, :Nyearmonth])
  dropmissing!(mf, [:ret, :aum])

  #winsorze for data quality
  winsorizewithin!(mf, Fsource=:ret, prop=PARAM[:mfqwinsorize], Fgroup=:Nyearmonth,
    twosided=true, Fnew=:WXret)
  select!(mf, Not([:ret]))
  rename!(mf, :WXret=>:ret)
  println(describe(mf))

  mf.aumbom = mf.aum |> deepcopy
  mfbyfund = groupby(mf, [:fundid, :Nyearmonth])


  #the interpolation routine backs out net deposits, then discounts
  #them back to the beginning of period. The discounted deposits are then added
  #pro-rata to the base beginning of period amount and discoutned back up to the
  #focal day. In equations,
  # Let CUMR(t) be the cumulative return for the month on day t. Then
  # totaldepbom = (AUM(T)-AUM(0)*(1+CUMR(T))/(1+CUMR(T))
  # AUM(t)=(AUM(0)+totaldepbom*t/T)*CUMR(t)
  function interpolatereturns(
      ret::AbstractVector,  aum::AbstractVector, aumeom::AbstractVector) where T

    cumret = cumprod(ret .+ 1.0) .- 1.0
    netgains = aumeom[end] - (1.0 + cumret[end]) * aum[1]
    netgainsbom = coalesce(netgains / (1.0+cumret[end]), 0.0)
    N::Int = length(ret)

    interpolated = [(t/N*netgainsbom+aum[1])*(1+cumr) for (t,cumr) ∈ enumerate(cumret)]

    #maybe disable the assertion if it slows things down
    @assert (interpolated[end] + 1.0 ≈ coalesce(aumeom[end], aum[1]*(1+cumret[end])) + 1.0) "
      aumeom[end]=$(aumeom[end])
      but interpolated[end]=$(interpolated[end])"
    return interpolated
  end

  @assert issorted(mf, [:fundid, :date])
  if interpolate

    #this will be correct for all but the last day of the month
    transform!(mfbyfund, [:ret, :aum, :aumeom] => interpolatereturns => :aum)
  else
    transform!(mfbyfund,
      [:ret, :aum] => ((ret,aum)->cumprod(ret .+ 1.0) .* aum) => :aum)
    #now patch in the eom aum at the end of the month
    transform!(mfbyfund,
      [:aum, :aumeom] => ((aum, aumeom)->aumeom[end] ≡ missing ?
        aum : [aum[1:(end-1)]; aumeom[end]]) => :aum)
  end






  if mfdatefrequency ≡ :week
    mf.weekid = assignweekid(mf.date)
    mfw = combine(groupby(mf, [:fundid, :weekid], sort=true),
      :ret => (r -> prod(1.0 .+ r) - 1.0) => :ret,
      :date => last => :date,
      :aum => last => :aum,
      )

    #we need to conform the weekly dates
    conformeddates = combine(groupby(mfw[!, [:weekid, :date]], :weekid),
      :date => maximum => :date)
    select!(mfw, Not([:date]))
    oldrows = nrow(mfw)
    mfw = innerjoin(mfw, conformeddates, on=:weekid, validate=(false,true))
    @assert oldrows == nrow(mfw)
    if !issorted(mfw, [:fundid, :weekid])
      sort!(mfw, [:fundid, :weekid])
    end

    #writes the output to a CSV for inspection
    if PARAM[:testfunds]
      testfundid = mfw[findmax(mfw.aum)[2], :fundid]
      testmfw = mfw[mfw.fundid .== testfundid,:]
      rename!(testmfw, Dict(:ret=>:retw, :aum=>:aumw, :date=>:datew))
      testmf = innerjoin(testmfw, mf, on=[:fundid, :weekid])

      testmf.ym = testmf.date .|> (d->year(d) + month(d)/100)
      oldmonthlyaum.ym = oldmonthlyaum.date .|> (d->year(d) + month(d)/100)
      rename!(oldmonthlyaum, :date=>:dateoldaum)
      testmf = leftjoin(testmf, oldmonthlyaum, on=[:fundid, :ym])
      sort!(testmf, [:fundid, :date])
      testmf |> CSV.write(
        "$(PARAM[:testpath])\\$(PARAM[:crsptypesuffix])dailyforweeklymf_$(mfdatefrequency).csv")

    end
    mf = mfw
  else
    @assert mfdatefrequency ≡ :day
  end

  #throw("got through it- check now!")

  return mf

end




parsemf(::Type{Date}, s::String, mfdateformat::DateFormat=dateformat"yyyymmdd"
  ) = Dates.Date(s, mfdateformat)
parsemf(::Type{<:Any}, ::Missing) = missing
parsemf(::Type{Date}, i::Int, mfdateformat::DateFormat=dateformat"yyyymmdd") = parsemf(
  Date, "$i", mfdateformat)

#parse to a number if possible
parsemf(T::Type{<:Real}, s::String) = something(tryparse(T, s), missing)
parsemf(::Type{T}, x::T) where T = x
