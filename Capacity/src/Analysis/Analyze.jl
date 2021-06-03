function analyze()

  measure, measureinfo = loadmeasure()
  @unpack prefix, analysispath, outputpath, resultid = measureinfo

  try
    PARAM[:refreshverify] && verifymeasure(; measure, measureinfo)
    PARAM[:refreshleadlag] && analyzeleadlag(; measure, measureinfo)

  catch err

    mv(outputpath, "$(analysispath)\\$prefix\\ERR_$resultid")
    #show(stdout, MIME"text/plain"(), stacktrace(catch_backtrace()))
    @error "Verificaiton failed" exception=(err, catch_backtrace())
    throw("Verification failed")
  end

end


#performs file io and initial setup for a new measure estimation
function loadmeasure(;
    burnin::DatePeriod = PARAM[:measureburnin],
    measurefilename::String = PARAM[:analysismeasurefilename],
    analysispath = PARAM[:analysispath],
    returntype::String = PARAM[:crspdatalabel],
    workingpath::String = PARAM[:workingpath],
    hfrname::String = PARAM[:hfrname],
    mfname::String = PARAM[:mfname],
    plname::String = PARAM[:plname],)

  datefrequency = Dict("w"=>:week, "m"=>:month)[PARAM[:crspfrequencysuffix]]

  prefix = (returntype === measurefilename[1:length(returntype)]
    ) ? measurefilename : "$(returntype)_$measurefilename"

  #create a directory to hold anlsyis results
  resultid = "rid$(Dates.format(now(),"yymmdd-HHMM"))"

  #load the measure file
  local measure::DataFrame
  inputpath = "$analysispath\\$prefix\\inputs"

  #each results set gets its own output path- do the below only once
  if !ispath("$analysispath\\$prefix")
    mkdir("$analysispath\\$prefix")
    mkdir("$analysispath\\$prefix\\inputs")
    mv("$analysispath\\$measurefilename.csv",  "$inputpath\\$measurefilename.csv")
    mv("$analysispath\\$measurefilename.txt",  "$inputpath\\$measurefilename.txt")
    cp("$workingpath\\$(returntype)_$(hfrname).csv", "$inputpath\\$(returntype)_$(hfrname).csv")
    cp("$workingpath\\$(returntype)_$(mfname).csv", "$inputpath\\$(returntype)_$(mfname).csv")
    #NOTE: momentum info is disabled for now- need to re-run data if going to try it
    if isfile("$workingpath\\$(returntype)_$(plname).csv")
      @warn "Momentum file detected! Note you will need to disable the error barrier in order
        to actually run the momentum tests"
      cp("$workingpath\\$(returntype)_$(plname).csv", "$inputpath\\$(returntype)_$(plname).csv")
    end
  else
    isfile("$analysispath\\$measurefilename.csv") && throw("Duplicate measure file
      $measurefilename found in analysis folder outside of input/results folder $prefix.
      This should not happen- measurefiles should correspond 1:1 with the input/results folder,
      unless the input/results folder does not exist (at which point it will be created
      immediately).")
  end


  #load the input data
  measure = CSV.File("$inputpath\\$measurefilename.csv") |> DataFrame

  #create a directory to hold anlsyis results
  resultid = "rid$(Dates.format(now(),"yymmdd-HHMM"))"
  outputpath = "$analysispath\\$prefix\\$resultid"


  @assert !ispath("$analysispath\\$prefix\\$resultid")
  mkdir(outputpath)
  cp(PARAMETER_FILE, "$outputpath\\$(resultid)_Parameters.xlsx")


  #now pre-process the columns
  #sums accross assets
  if datefrequency === :week && PARAM[:analysisaggregateweeklycols]
    assetmat = measure[!, findcols(measure, :asset)] |> Matrix{MFloat64}
    oldcols = [findcols(measure, :growth); findcols(measure, :asset);
      propertynames(measure[!, r"^V_"])]
    aggassets = sum(assetmat, dims=2) |> vec
    measure[!, Symbol(:A_, weeklyaggregatename)] = aggassets
    measure[!, Symbol(:G_, weeklyaggregatename)] = [
      missing; aggassets[2:end] ./ aggassets[1:(end-1)]; ]

    #also need to some accross volume
    volumemat = measure[!, propertynames(measure[!, r"^V_"])] |> Matrix{MFloat64}
    aggvolume = sum(volumemat, dims=2) |> vec
    measure[!, Symbol(:V_, weeklyaggregatename)] = aggvolume
    select!(measure, Not(oldcols))
  end

  #pre-process the measure data
  measure = measure[measure.date .≥ (minimum(measure.date) + burnin), :]
  derivemeasurecols!(measure; datefrequency)

  #setup neweywest SEs and moving block bootstrap
  #below was previously ceil(max(nrow(measure),nrow(funds))^0.25) |> Int
  neweylags::Int = nrow(measure)^0.25 |> ceil |> Int
  nwΣ = neweywestΣfunc(neweylags)

  bsaggfunc(v::AbstractVector) = (;se=std(skipmissing(v)), p=nonparametricp(v))
  #bsaggfunc(v) = (;se=std(skipmissing(v)), p=nonparametricp(v))

  bsmeasurecolsregex::Regex = Regex(
    PARAM[:bsmeasurecolsregex] === nothing ? "" : PARAM[:bsmeasurecolsregex])
  bsbenchmarkcolsregex::Regex = Regex(
    PARAM[:bsbenchmarkcolsregex] === nothing ? "" : PARAM[:bsbenchmarkcolsregex])
  bssamples = PARAM[:bssamples]

  #note- the below was previously #ceil(max(nrow(measure),nrow(funds))^(1/3)) |> Int
  bswidth = nrow(measure)^(1/3) |> ceil |> Int

  #we can reuse some of the bootstrap aggregation code for the null simulations
  nullar1aggfunc(v::AbstractVector,statname,stattestvals) = (;
    se=std(skipmissing(v)),
    p=nonparametricnullp(v, stattestvals[statname]))
  nullar1aggfunc(M::AbstractMatrix,statname,stattestvals) = (;
    se=std(M, dims=1) |> vec,
    p=nonparametricnullp(M, stattestvals[statname]))
  nullar1simulations = bssamples

  #lightly package the io data for use as needed
  measureinfo = (; resultid, inputpath, outputpath, returntype,
    analysispath, datefrequency, prefix,
    neweylags, nwΣ,
    bsaggfunc, bsmeasurecolsregex, bsbenchmarkcolsregex, bssamples, bswidth,
    nullar1aggfunc, nullar1simulations)

  return (measure, measureinfo)
end


#helpful for grabbing the column name sans prefix
stripcolprefix(s) = replace(s, r"^l?[A-Z0-9]+_"=>"") |> Symbol
stripcolprefix(s::Symbol) = stripcolprefix(string(s))

#compute the strategy pressure measure- AUM×turnover/aggvol
function derivemeasurecols!(measure;
    Fassets::Vector{Symbol}=findcols(measure, :asset, includelp=false),
    datefrequency::Symbol,
    lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
    Ftotalvolume::Symbol = Symbol(:desc_, PARAM[:analysistotalvolume]),
    lladditionalreturncolsidx = PARAM[:lladditionalreturncolsidx],
    lladditionalnormalizers = PARAM[:lladditionalnormalizers],
    llmarketcol = PARAM[:llmarketcol])

    #other helper methods
    simpledif(v::AbstractVector) = [missing; v[2:end] .- v[1:(end-1)]]
    _logdif(v) = [missing; log.(v[2:end]) .- log.(v[1:(end-1)])]
    logdif(v) = _logdif(abs.(v)) .* sign.(v)

    #generate log growth columns
    for gcol ∈ findcols(measure, :growth)
      measure[!, gcol] = [missing; measure[2:end, gcol]]
      measure[!, Symbol(:l,gcol)] = log.(abs.(measure[!, gcol])) .* sign.(measure[!, gcol])
      measure[!, Symbol(:g,string(gcol)[2:end])] = measure[!, gcol] .- 1

      rcol = Symbol(:R_, stripcolprefix(gcol))
      acol = Symbol(:A_, stripcolprefix(gcol))
      agcol = Symbol(:ag, string(gcol)[2:end])
      ngcol = Symbol(:ng, string(gcol)[2:end])
      cggcol = Symbol(:cgg, string(gcol)[2:end])
      measure[!, agcol] = measure[!, gcol] .- measure[!, rcol] .- 1.0
      measure[!, ngcol] = measure[!, gcol] .- measure[!, llmarketcol] .- 1.0
      measure[!, cggcol] = measure[!, rcol] |> deepcopy

      #active asset growth
      measure[!, Symbol(:aA, string(gcol)[2:end])] = [missing; measure[1:(end-1), acol] .* measure[2:end, agcol]]
      measure[!, Symbol(:nA, string(gcol)[2:end])] = [missing; measure[1:(end-1), acol] .* measure[2:end, ngcol]]
    end


    #Fassetcols = assetcols(measure)
    #make growth cols
    measure = lpfilter(measure,focalcols = [
      findcols(measure, :loggrowth); findcols(measure, :asset); :mctotal])
    colstodif = findcols(measure, :loggrowth)
    transform!(measure,  colstodif .=>
      (v -> [missing; v[2:end] .- v[1:(end-1)]])
      .=> (s->Symbol(:G, s)).(colstodif))

    Fvols = propertynames(measure[!, r"^V_"])
    Fpurchases = propertynames(measure[!, r"^Vs_"])
    Fpressures = (s->Symbol(:P_, stripcolprefix(s))).(Fassets)
    Fparticipations = (s->Symbol(:p_, stripcolprefix(s))).(Fassets)
    Flgpressures = (s->Symbol(:lG, s)).(Fpressures)
    totalvolume = measure[!, Ftotalvolume]
    totalcapitalization = measure[!, :mctotal]

    Fpurchaseparticipations = (s->Symbol(:pp_, stripcolprefix(s))).(Fassets)
    Fpurchasefractions = (s->Symbol(:VsA_, stripcolprefix(s))).(Fassets)
    Fpurchasecaps = (s->Symbol(:VsMC_, stripcolprefix(s))).(Fassets)

    #compute the presure characteristic
    transform!(measure, Fvols .=>
      ((vol)-> vol ./ totalvolume) .=> Fparticipations)
    transform!(measure, Fpurchases .=>
      ((pur)-> pur ./ totalvolume) .=> Fpurchaseparticipations)
    transform!(measure, zip(Fassets, Fpurchases) .|> collect .=>
      ((asset,pur)-> pur ./ asset) .=> Fpurchasefractions)
    transform!(measure,  Fpurchases .=>
      ((pur)-> pur ./ totalcapitalization) .=> Fpurchasecaps)
    transform!(measure, Fvols .=>
      ((vol)-> log.(vol)) .=> ((s->Symbol(:l,s)).(Fvols)))
    transform!(measure, zip(Fassets, Fvols) .|> collect .=>
      ((asset, vol)->asset .* vol ./ totalvolume) .=> Fpressures)
    transform!(measure, Fpressures .=> logdif .=> Flgpressures)
    lpfilter(measure, focalcols=Flgpressures)
    Flgpressures = [Flgpressures; lpcols.(Flgpressures)...;]

    #compute a version just of the volume-pressure
    Fvpressures = (s->Symbol(:VP_, stripcolprefix(s))).(Fassets)
    Flgvpressures = (s->Symbol(:lGVP_, stripcolprefix(s))).(Fassets)
    transform!(measure, Fvols .=> (v -> logdif(v ./ totalvolume)) .=> Flgvpressures)
    lpfilter(measure, focalcols=Flgvpressures)
    Flgvpressures = [Flgvpressures; lpcols.(Flgvpressures)...]

    #compute a measure of dollar return
    Frets = (s->Symbol(:R_,s)).(stripcolprefix.(Fassets))
    Fdolrets = (s->Symbol(:X,s)).(Frets)
    transform!(measure, zip(Fassets, Frets) .|> collect .=> (
      (a,r)-> a .* r) .=> Fdolrets)


    #somewhat different calculation methedology here for the non-logged version
    #in order to avoid various operations on ratios
    lpfilter(measure, focalcols=Fvols;) #[Fvols; Ftotalvolume;])
    Fvpressures = [Fvpressures; lpcols.(Fvpressures)...;]
    Fgvpressures = (s->Symbol(:G, s)).(Fvpressures)
    Fvolsforvpressures = [Fvols; lpcols.(Fvols)...;]
    Ftotalvolumes = [Ftotalvolume for i ∈ 1:length(Fvolsforvpressures)]
    transform!(measure, zip(Fvolsforvpressures, Ftotalvolumes) .|> collect .=>
      ((v,tv) -> simpledif(v) ./ tv) .=> Fgvpressures)


    #generate the appropriate column labels
    Fmcassetcols = [Fassets; lpcols.(Fassets)...;] #all asset cols
    Fmccols = findcols(measure, :mc)
    Fmccols = [Fmccols; (Fmccol->[Fmccol for i in 1:length(lpcols(Fmccol))]).(Fmccols)...;]
    Famccols = Fassets .|>  stripcolprefix .|> (s)->Symbol(:M_, s)
    Famccols = [Famccols; lpcols.(Famccols)...;]
    Fgamccols = Famccols .|> s-> Symbol(:G, s)

    #these will be used for return predictability
        #WARNING: quick and dirty test
    #measure.mctotal = [missing; measure.mctotal[1:(end-1)]]
    lpfilter(measure, focalcols=[:mctotal]) #need to smooth the mctotal
    Fmctotcols = [:mctotal for i ∈ 1:length(Fassets)]
    Fmctotcols = [Fmctotcols; (Fmctotcol->
      [Fmctotcol for i in 1:length(lpcols(Fmctotcol))]).(Fmctotcols)...;]
    Famctotcols = Fassets .|>  stripcolprefix .|> (s)->Symbol(:MT_, s)
    Famctotcols = [Famctotcols; lpcols.(Famctotcols)...;]
    Flamctotcols = Famctotcols .|> s-> Symbol(:l, s)
    Fgamctotcols = Famctotcols .|> s-> Symbol(:G, s)


    Faamcactiveassetcols = Fassets .|> (s)->Symbol(:a, s)
    Faamcactiveassetcols = [Faamcactiveassetcols; lpcols.(Faamcactiveassetcols)...;]
    Fagamctotcols = Famctotcols .|> s-> Symbol(:aG, s)

    Fnamcactiveassetcols = Fassets .|> (s)->Symbol(:n, s)
    Fnamcactiveassetcols = [Fnamcactiveassetcols; lpcols.(Fnamcactiveassetcols)...;]
    Fngamctotcols = Famctotcols .|> s-> Symbol(:nG, s)

    #NOTE we take the ratio of the average rather than the average of the ratio,
    #while also following the Welch 2020 RoC approach

    transform!(measure,
      zip(Fmcassetcols, Fmccols) .|> collect .=>
        ((asset, mc)->(asset ./ mc)) .=> Famccols,
      zip(Fmcassetcols, Fmccols) .|> collect .=>
        ((asset, mc)->(simpledif(asset) ./ mc)) .=> Fgamccols,

      zip(Fmcassetcols, Fmctotcols) .|> collect .=>
        ((asset, mctot)->(asset ./ mctot)) .=> Famctotcols,
      zip(Fmcassetcols, Fmctotcols) .|> collect .=>
        ((asset, mctot)->log.(abs.(asset ./ mctot)) .* sign.(asset)) .=> Flamctotcols,
      zip(Fmcassetcols, Fmctotcols) .|> collect .=>
        ((asset, mctot)->(simpledif(asset) ./ mctot)) .=> Fgamctotcols,

      zip(Faamcactiveassetcols, Fmctotcols) .|> collect .=>
        ((activeasset, mctot)->(activeasset ./ mctot)) .=> Fagamctotcols,
      zip(Fnamcactiveassetcols, Fmctotcols) .|> collect .=>
        ((netasset, mctot)->(netasset ./ mctot)) .=> Fngamctotcols,
        )


    #the below algorithm normalizes for an arbitrary denominator
    for (prefix, Fdenom) ∈ zip(keys(lladditionalnormalizers), values(lladditionalnormalizers))
      @assert Fdenom ∈ propertynames(measure)
      lpfilter(measure, focalcols=[Fdenom]) #need to smooth the denom

      #form the denominator names
      Fdenoms = [Fdenom for i ∈ 1:length(Fassets)]
      Fdenoms = [Fdenoms; (Fdenomcol->
        [Fdenomcol for i in 1:length(lpcols(Fdenomcol))]).(Fdenoms)...;]

      #form the names for the new fields
      Fadenoms = Fassets .|>  stripcolprefix .|> (s)->Symbol(prefix, s)
      Fadenoms = [Fadenoms; lpcols.(Fadenoms)...;]
      Fladenoms = Fadenoms .|> s-> Symbol(:l, s)
      Fgadenoms = Fadenoms .|> s-> Symbol(:G, s)
      Fagadenoms = Fadenoms .|> s-> Symbol(:aG, s)
      Faladenoms = Fadenoms .|> s-> Symbol(:al, s)

      #normalize as necesary
      transform!(measure,
        zip(Fmcassetcols, Fdenoms) .|> collect .=>
          ((asset, denom)->(asset ./ denom)) .=> Fadenoms,
        zip(Fmcassetcols, Fdenoms) .|> collect .=>
          ((asset, denom)->log.(abs.(asset ./ denom)) .* sign.(asset)) .=> Fladenoms,
        zip(Fmcassetcols, Fdenoms) .|> collect .=>
          ((asset, denom)->(simpledif(asset) ./ denom)) .=> Fgadenoms,
        zip(Faamcactiveassetcols, Fdenoms) .|> collect .=>
          ((activeasset, denom)->(activeasset ./ denom)) .=> Fagadenoms,)

    end


    #lgvpressures doesn't include lpcols, but Fgvpressures does
    allcols = [Flgpressures; Fgvpressures; Flgvpressures]
    doubledifcols = allcols .|> s->Symbol(:G, s)
    transform!(measure, allcols .=> simpledif .=> doubledifcols)


    measure |> CSV.write(
      "$(PARAM[:testpath])\\$(PARAM[:analysismeasurefilename])_derive.csv")

    return measure
    #Fpressures =
end


#helper functions for acquiring data from the measure file
function findcols(df, coltype::Symbol; includelp=true, args...)::Vector{Symbol}

  simplepatterns = Dict(
    :loggrowth => r"^lG_",
    :growth => r"^G_",
    :growthM1 => r"^g_",
    :activegrowth => r"^ag_",
    :netgrowth => r"^ng_",
    :capitalgainsgrowth => r"^cgg_",
    :activeassetalloc => r"^aA_",
    :z => r"^Z_",
    :volume => r"^V_",
    :lvolume => r"^lV_",
    :purchases => r"^Vs_",
    :purchasefraction => r"^VsA_",
    :purchasecapitalization => r"^VsMC_",
    :return => r"^R_",
    :ff => r"^F_",
    :logpressure => r"^lP_",
    :growthpressure => r"^GP_",
    :loggrowthpressure => r"^lGP_",
    :participation => r"^p_",
    :purchaseparticipation => r"^pp_",
    :sumloggrowthpressure => r"^lGSP_",
    :growthvolumepressure => r"^GVP_",
    :loggrowthvolumepressure => r"^lGVP_",
    :assetmarketcap => r"^M_",
    :assetmarketcapgrowth => r"^GM_",

    :assetmf => r"^AN_",
    :lassetmf => r"^lAN_",
    :assetmfgrowth => r"^GAN_",
    :activelassetmf => r"^alAN_",
    :activeassetmfgrowth => r"^aGAN_",

    :assetmfequity => r"^EN_",
    :lassetmfequity => r"^lEN_",
    :assetmfequitygrowth => r"^GEN_",
    :activeassetmfequitygrowth => r"^aGEN_",

    :assetmfequityonly => r"^EON_",
    :lassetmfequityonly => r"^lEON_",
    :assetmfequityonlygrowth => r"^GEON_",
    :activeassetmfequityonlygrowth => r"^aGEON_",

    :assetmfdomesticeo => r"^DEON_",
    :lassetmfdomesticeo => r"^lDEON_",
    :assetmfdomesticeogrowth => r"^GDEON_",
    :activeassetmfdomesticeogrowth => r"^aGDEON_",

    :assetmarketcaptotal => r"^MT_",
    :lassetmarketcaptotal => r"^lMT_",
    :assetmarketcaptotalgrowth => r"^GMT_",
    :activeassetmarketcaptotalgrowth => r"^aGMT_",
    :netassetmarketcaptotalgrowth => r"^nGMT_",
    )

  local foundcols
  if coltype ∈ keys(simplepatterns)
    foundcols = propertynames(df[!, simplepatterns[coltype]])
  else
    foundcols = findcols(df, Val{coltype}(); args...)
  end

  if !includelp #sometimes we may want to exclude the filtered cols
    foundcols = setdiff(foundcols, propertynames(df[!, r"_LP[0-9]+d$"]))
  end
  return foundcols
end


findcols(measure::AbstractDataFrame, ::Val{:asset}; includelp=true) = setdiff(
  propertynames(measure[!, r"^A_"]), includelp ? [] : propertynames(measure[!, r"_LP[0-9]+d$"]))
findcols(measure::AbstractDataFrame, ::Val{:mc};
  mctype=:gross) = propertynames(measure[!, Regex("^mc$(mctype)_")])

#generic fallback
findcols(measure::AbstractDataFrame, reg::Regex; includelp=true) = setdiff(
  propertynames(measure[!, reg]), includelp ? [] : propertynames(measure[!, r"_LP[0-9]+d$"]))
findcols(measure::AbstractDataFrame, reg::String; args...) = findcols(
  measure, Regex(reg), args...)

lpcols(F::Symbol;
  lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods]) = (
    (p->Symbol(F, lpsuffix(p))).(setdiff(lpfilterperiods, [nothing])))

differencecol(s::Symbol) = Symbol(:G, s)
