function verifymeasure(;
    measure::DataFrame,
    measureinfo::NamedTuple,
    lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
    benchmarkfunds::Bool = PARAM[:analysisbenchmarkfunds],
    benchmarkcomomentum::Bool=PARAM[:analysisbenchmarkcomomentum],
    comomcols::Union{Vector{Symbol}, Nothing} = PARAM[:analysiscomomcols],
    hfrname::String = PARAM[:hfrname],
    mfname::String = PARAM[:mfname],
    plname::String = PARAM[:plname],
    weeklyaggregatename = PARAM[:analysisweeklyaggregatename]
  )

  @unpack inputpath, outputpath, resultid, returntype, datefrequency = measureinfo


  hf = CSV.File("$inputpath\\$(returntype)_$(hfrname).csv") |> DataFrame
  mf = CSV.File("$inputpath\\$(returntype)_$(mfname).csv") |> DataFrame


    #compare the measure to particular benchmarks
    if datefrequency === :month
      validatedates(measure.date, frequency=:month, validateexactfrequency=true)

      if benchmarkfunds
        hfresults = comparetofundassets(measure, measureinfo, hf;
            fundslabel=:hf,
            datefrequency=:month,
            measurecols = findcols(measure, :loggrowth),
            benchmarkcols = findcols(hf, :loggrowth),
            usefocalcolheuristics=true,
            matchonlpfilter=true)

        mfresults = comparetofundassets(measure, measureinfo, mf;
            fundslabel=:mf,
            datefrequency=:month,
            measurecols = findcols(measure, :loggrowth),
            benchmarkcols = findcols(mf, :loggrowth),
            usefocalcolheuristics=true,
            matchonlpfilter=true)
        fundresults = [hfresults.results; mfresults.results;]
        sort!(fundresults, :cor)
        hfresults.data |> CSV.write("$outputpath\\$(resultid)_$(hfrname)_series.csv")
        mfresults.data |> CSV.write("$outputpath\\$(resultid)_$(mfname)_series.csv")
        fundresults |> CSV.write("$outputpath\\$(resultid)_fund_corresults.csv")

    elseif datefrequency === :week
      validatedates(measure.date, frequency=:week, validateexactfrequency=true)
      benchmarkfunds && throw("Only monthly data is supported when comparetofunds=true")
    end

    else

      @info "Note: fund benchmark analysis not run as analysisbenchmarkfunds=false"

    end

    if benchmarkcomomentum
      throw("Momentum not set up!!! Need to re-run everything before disabling this error.
        To set up, re-run/refresh momentum soup-2-nuts, then disable this error
        (can also disable the load warning in the analyze file)")

      comom = CSV.File("$inputpath\\$(returntype)_$(plname).csv") |> DataFrame

      benchmarkcols = isnothing(comomcols) ? setdiff(propertynames(comom), [:date]) : comomcols
      measurecols = [findcols(measure, :loggrowth); findcols(measure, :loggrowthpressure);
        findcols(measure, :assetmarketcap); findcols(measure, :loggrowthvolumepressure);
        findcols(measure, :growthvolumepressure)]

      println("benchmarkcols: $benchmarkcols")
      measurecols = [measurecols; measurecols .|> differencecol]
      comomresults = comparetofundassets(measure, measureinfo, comom,
        fundslabel = :pl,
        matchonlpfilter=false,
        usefocalcolheuristics=false;
        datefrequency, benchmarkcols, measurecols)
      sort!(comomresults.results, :cor)
      comomresults.data |> CSV.write("$outputpath\\$(resultid)_$(plname)_series.csv")
      comomresults.results |> CSV.write("$outputpath\\$(resultid)_comom_corresults.csv")
    else
      @info "Note: comomentum not run as analysisbenchmarkcomomentum=false"
    end




end


#useful function to translate monthyl dates between dataframes
function conformmonthlydates(old::AbstractVector{Date}, valid::AbstractVector{Date};
  datejoin = leftjoin)

  validatedates(old, frequency=:month, validateexactfrequency=true)
  validatedates(valid, frequency=:month, validateexactfrequency=true)

  #form a lookup table out of the available dates
  olddf = DataFrame(idx=1:length(old) |> collect, old=old)
  olddf.year = year.(old)
  olddf.month = month.(old)

  validdf = DataFrame(valid=valid)
  validdf.year = year.(valid)
  validdf.month = month.(valid)
  @assert sum(nonunique(validdf, [:year, :month])) == 0


  joint = datejoin(olddf, validdf, on=[:year, :month], validate=(false, true))
  @assert sum(nonunique(joint, [:year, :month])) == 0
  sort!(joint, :idx)

  return joint.valid
end

#useful function to translate weekly dates between dataframes
function conformweeklydates(old::AbstractVector{Date}, valid::AbstractVector{Date};
  datejoin = leftjoin)

  validatedates(old, frequency=:week, validateexactfrequency=true)
  validatedates(valid, frequency=:week, validateexactfrequency=true)

  olddf = DataFrame(idx=1:length(old) |> collect, old=old, lastdow = lastdayofweek.(old))
  validdf = DataFrame(valid=valid, lastdow = lastdayofweek.(valid))
  @assert sum(nonunique(validdf, [:lastdow])) == 0

  joint = datejoin(olddf, validdf, on=:lastdow, validate=(false, true))
  @assert sum(nonunique(joint, :lastdow)) == 0
  sort!(joint, :idx)

  return joint.valid
end



#compute correlations and statistics
function comparetofundassets(measure, measureinfo, funds;
    fundslabel::Symbol = throw("results label is required"),
    datefrequency::Symbol = throw("datefrequency is required"),
    measurecols = setdiff(propertynames(measure), [:date]),
    matchonlpfilter=true,
    benchmarkcols = setdiff(propertynames(funds), [:date]), #benchmarkcols
    extrabetaweights = PARAM[:fundextrabetaweights],
    extrabetaweightsidx = PARAM[:fundextrabetaweightsidx],
    minσ = eps(Float64)^0.5,
    usefocalcolheuristics = true
    )

  function bsstatfunc(sim)
    b,m = sim[1], sim[2]

    ρ = cor(b,m)
    return Dict(
      :cor  => ρ,
      :kend =>corkendall(b,m),
      :spear =>corspearman(b,m),
      :betabm => ρ * std(b)/std(m),
      :betamb => ρ * std(m)/std(b))
  end

  nullar1statfunc(b,m) = bsstatfunc((b,m))

  @unpack (neweylags, nwΣ,
    bsaggfunc, bsmeasurecolsregex, bsbenchmarkcolsregex, bssamples, bswidth,
      nullar1aggfunc, nullar1simulations) = measureinfo

  if datefrequency ≡ :month
    funds.measuredate = conformmonthlydates(funds.date, measure.date)
  elseif datefrequency ≡ :week
    funds.measuredate = conformweeklydates(funds.date, measure.date)
  else
    throw("unrecognized date frequency")
  end

  funds = funds[completecases(funds, :measuredate),:]
  select!(funds, Not(:date))
  joint = innerjoin(measure, funds,
    on=:date=>:measuredate, validate=(true,true))
  sort!(joint, :date)

  #compute the statistics pairwise in a semi-long format
  results = DataFrame(
    label=Symbol[], #label in case we concatenate this df with others
    measure=Symbol[], #measure field that we are benchmarking
    benchmark=Symbol[], #benchmark field

    cor=MFloat64[],
    cort=MFloat64[],
    corhac=MFloat64[],
    cormbbz = MFloat64[],
    cormbbp = MFloat64[],
    corar1h0z = MFloat64[],
    corar1h0p = MFloat64[],

    spear=MFloat64[],
    speart=MFloat64[],
    spearmbbz=MFloat64[],
    spearmbbp=MFloat64[],
    spearar1h0z=MFloat64[],
    spearar1h0p=MFloat64[],
    #spearz=MFloat64[],
    kend=MFloat64[],
    kendz=MFloat64[],
    kendmbbz=MFloat64[],
    kendmbbp=MFloat64[],
    kendar1h0z=MFloat64[],
    kendar1h0p=MFloat64[],

    betabm=MFloat64[],
    betabmmw=MFloat64[],
    betabmhac=MFloat64[],
    betabmmbbp=MFloat64[],
    betabmmbbz=MFloat64[],
    betabmar1h0z=MFloat64[],
    betabmar1h0p=MFloat64[],

    betamb=MFloat64[],
    betambmw=MFloat64[],
    betambhac=MFloat64[],
    betambmbbp=MFloat64[],
    betambmbbz=MFloat64[],
    betambar1h0z=MFloat64[],
    betambar1h0p=MFloat64[],
    N=MInt[])


  resultsbythread = [deepcopy(results) for t ∈ 1:Threads.nthreads()]

  #this subfunction determines which measure cols will be tested against a bcol
  function formfocalcols(bcol)
    matchedcols = matchonlpfilter ? matchcolsonlpsuffix(bcol, measurecols) : measurecols
    (!usefocalcolheuristics) && return measurecols

    focalmeasurecols = Symbol[]

    #=if (strippedbcol |> Symbol) ∉ extrabetaweights
      #along with maybe the s
      strippedbcol = stripcolprefix(bcol)=#

    for mcol ∈ matchedcols
      local isfocal::Bool
      strippedmcol = stripcolprefix(mcol) |> string

      if occursin(strippedmcol, bcol |> string)
        isfocal=true

      #see if we explicitly want to include it as a focal col
      else
        isfocal = false
        for k ∈ extrabetaweights
          #in this case, the specified extra column is embeded in the name
          # of the bcol
          if occursin(k|> string, bcol |> string) && (
            any((m->
              occursin(m|>string, strippedmcol)).(extrabetaweightsidx[k])))

            isfocal=true
          end
        end
      end

      if isfocal
        push!(focalmeasurecols, mcol)
      end
    end

    return focalmeasurecols
  end
  for bcol ∈ benchmarkcols
    tresults = resultsbythread[Threads.threadid()]

    focalmeasurecols = formfocalcols(bcol)
    #focalmeasurecols = matchonlpfilter ? matchcolsonlpsuffix(bcol, measurecols) : measurecols
    #we want to match the smoothing/lp filter level between the measure and the benchmark

    for mcol ∈ focalmeasurecols



      jointnomissing = view(joint, completecases(joint, [mcol,bcol]), [:date, mcol, bcol])
      m::Vector{Float64} = jointnomissing[:,mcol]
      b::Vector{Float64} = jointnomissing[:,bcol]
      N=nrow(jointnomissing)

      #try
        #stats all from wikipedia, so probably want better sources if used
        pear = cor(m,b) #pearson correlation
        peart = pear*((N-2)/(1-pear^2))^(0.5)
        spear = StatsBase.corspearman(m, b) #spearman correlation
        #spearz = atanh(spear)*((N-3)/1.06)^0.5
        speart = spear*((N-2)/(1-spear^2))^(0.5)
        kend = StatsBase.corkendall(m, b) #kendall's tau
        kendz = kend/((2*(2*N+5))/(9*N*(N-1)))^0.5 #not sure about this- investigate if using

        #=Xb = hcat(b, ones(N))
        lm = FMLM(Xb, m)
        beta=lm.β[1]
        @assert lm.β[1] ≈ pear * std(m)/std(b)=#
        Xm = hcat(m, ones(N))
        lmbm = FMLM(Xm, b)
        betabm=lmbm.β[1]
        @assert (lmbm.β[1] ≈ pear * std(b)/std(m)) || (lmbm.β[1] ≡ pear * std(b)/std(m)) "
          lmbm.β[1]= $(lmbm.β[1]) while pear * std(b)/std(m) = $(pear * std(b)/std(m))!!"
        betabmmw = betabm/modifiedwhiteΣ!(lmbm)[1,1]^0.5
        betabmhac = betabm/nwΣ(lmbm)[1,1]^0.5

        Xb = hcat(b, ones(N))
        lmmb = FMLM(Xb, m)
        betamb=lmmb.β[1]
        @assert ((lmmb.β[1] ≈ pear * std(m)/std(b)) || (lmmb.β[1] ≡ pear * std(m)/std(b)))
        betambmw = betamb/modifiedwhiteΣ!(lmmb)[1,1]^0.5
        betambhac = betamb/nwΣ(lmmb)[1,1]^0.5

        #compute the pearson correlation as a regression
        #this is kinda hacky- might want to replace with a block bootstrap
        #=Xbcor = hcat(b .* std(m) ./ std(b), ones(N))
        lmcor = FMLM(Xbcor, m)=#
        Xmcor = hcat(m .* std(b) ./ std(m), ones(N))
        lmcor = FMLM(Xmcor, b)
        @assert ((lmcor.β[1] ≈ pear) || (lmcor.β[1] ≡ pear))
        pearhacz = pear/nwΣ(lmcor)[1,1]^0.5

        #perform moving block bootstrap if desired
        if (occursin(bsmeasurecolsregex, string(mcol)) &&
            occursin(bsbenchmarkcolsregex, string(bcol)) &&
            (std(b) > minσ) && (std(m) > minσ))
          #@info "bcol:$bcol, mcol:$mcol"
          bs = mbbootstrap([b,m], bsstatfunc, bsaggfunc, B=bssamples, w=bswidth)
          cormbbz = pear/bs.cor.se
          spearmbbz = spear/bs.spear.se
          kendmbbz =  kend/bs.kend.se
          betabmmbbz =  betabm/bs.betabm.se
          betambmbbz =  betamb/bs.betamb.se

          cormbbp = bs.cor.p
          spearmbbp = bs.spear.p
          kendmbbp =  bs.kend.p
          betabmmbbp =  bs.betabm.p
          betambmbbp =  bs.betamb.p

          ar1h0teststats = Dict(:cor=>pear, :kend=>kend, :spear=>spear,
            :betabm => betabm, :betamb => betamb)
          ar1h0 = comparear1undernull(b,m, nullar1statfunc,
            (v,s)->nullar1aggfunc(v,s,ar1h0teststats),
            Nsimulations=nullar1simulations)
          #these are the alternate hypotheses to test



          corar1h0z =  pear/ar1h0.cor.se
          spearar1h0z = spear/ar1h0.spear.se
          kendar1h0z =  kend/ar1h0.kend.se
          betabmar1h0z =  betabm/ar1h0.betabm.se
          betambar1h0z =  betamb/ar1h0.betamb.se

          corar1h0p = ar1h0.cor.p
          spearar1h0p = ar1h0.spear.p
          kendar1h0p =  ar1h0.kend.p
          betabmar1h0p =  ar1h0.betabm.p
          betambar1h0p =  ar1h0.betamb.p
        else
          cormbbp=missing; spearmbbp=missing; kendmbbp=missing; betabmmbbp=missing; betambmbbp=missing;
          cormbbz=missing; spearmbbz=missing; kendmbbz=missing; betabmmbbz=missing; betambmbbz=missing;

          corar1h0p=missing; spearar1h0p=missing; kendar1h0p=missing; betabmar1h0p=missing; betambar1h0p=missing;
          corar1h0z=missing; spearar1h0z=missing; kendar1h0z=missing; betabmar1h0z=missing; betambar1h0z=missing;
        end

        push!(tresults,(;
          label=fundslabel,
          measure=mcol,
          benchmark=bcol,
          cor=pear,
          cort=peart,
          corhac=pearhacz,
          cormbbz, cormbbp, corar1h0z, corar1h0p,
          spear, speart, spearmbbz, spearmbbp, spearar1h0z, spearar1h0p,
          kend, kendz, kendmbbp, kendmbbz, kendar1h0p, kendar1h0z,
          betabm, betabmmw, betabmhac, betabmmbbz,  betabmmbbp, betabmar1h0z,  betabmar1h0p,
          betamb, betambmw, betambhac, betambmbbz,  betambmbbp, betambar1h0z,  betambar1h0p,
          N))

    end
  end

  joint |> CSV.write(
    "$(PARAM[:testpath])\\$(PARAM[:analysismeasurefilename])_$(fundslabel)_jointcor.csv")

  #collect the results in a single df
  append!(results, reduce(append!, resultsbythread))
  sort!(results, [:benchmark, :measure])

  return (;results, data=joint)
end



#normalizedassetcols(measure) = propertynames(measure[!, r"^NA_"])
