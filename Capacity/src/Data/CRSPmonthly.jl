
#takes a crsp daily file and makes it monthly
function makemonthly(selectionop::Function = (noops(x)=x),
    monthlyop::Function = (noopm(x)=x);
    crspnameout=PARAM[:crspmonthly],
    mainpanel::Bool=true,
    binprefix = "m_",
    computeivol::Bool = PARAM[:computeivol])

  local crsp::DataFrame


  #WARNING- crspd is NOT valid, below is only to create crspm
  crsp = prepcrsp(crsptype=PARAM[:crspdaily],
    csvextension = CSV_EXTENSION,
    incsvstream = IN_CSV_STREAM,
    binprefix = binprefix,
    refreshcrsp=PARAM[:refreshmonthly],
    crspcolumns=[PARAM[:monthlycrspcols]; keys(PARAM[:momentumstrategies]) |> collect],
    crsprequiredcolumns=PARAM[:dailycrsprequired]) do crspd

    @assert (PARAM[:crspfrequencysuffix] == "m") || (!mainpanel)
    sort!(crspd, [:permno, :date]) #cache the sort

    #computerangemeasures!(crspd, sorted=true) #not used, mainly for kyle's lamda
    crspd.dvol = computedollarvol(crspd)
    rename!(crspd, :vol=>:svol)
    #crspd.shares = computeshares(crspd)

    #NOTE: testing note- change to trailingsigmatest! to test vol routine
    #println(" checksum: $(sum(skipmissing(crspd.sigma)))")

    crspd = mergeff(crspd, Fgroup=:permno, datefrequency=:day)
    #crspd = dailyop(crspd)

    if !issorted(crspd, [:permno, :date])
      sort!(crspd, [:permno, :date])
    end

    crspm::DataFrame = makemonthly(selectionop, crspd; computeivol)
    validatedates(crspm.date |> unique |> sort!; frequency=:month, validateexactfrequency=true)
    sort!(crspm, [:permno, :date])
    validatesummaryintegrity(crspm, Fid=:monthid)

    if PARAM[:testmonthlycsv]

      #time series test output
      testpermno::Int = PARAM[:testpermnomult]
      #testcrspm::DataFrame = crspm[crspm.permno .% testpermno .== 0, :]
      testcrspm::DataFrame = crspm[crspm.permno .== 14541, :]
      #if size(testcrspm,1)>1_000 #prevents scenarios where the file is larger than needed
      #  testcrspm=testcrspm[1:1_000]
      #end

      #=permnos::Vector{Int} = unique(crspm.permno)[1:400:min(7000, end)]
      crspm.totest = falses(size(crspm,1))
      for scrspm ∈ groupby(crspm, :permno)
        (first(scrspm.permno) ∈ permnos) && (scrspm.totest .= true)
      end
      testcrspm = crspm[crspm.totest,:]
      select!(crspm, Not(:totest))=#
      rename!(testcrspm, :date=>:datem, :dvol=>:dvolm, :ivol=>:ivolm, :ret=>:retm, :price=>:pricem,
        :tridx=>:tridxm, :shares=>:sharesm, :sigma=>:sigmam, :mc=>:mcm, :plret=>:plretm,
        :exchcd=>:exchcdm, :siccd=>:siccdm,)
      rename!(testcrspm, (k->k=>Symbol(k, :m)).(keys(PARAM[:momentumstrategies])))
      #make sure we merge using the same df as used to create the monthly df
      testcrspout = innerjoin(testcrspm, crspd, on=[:permno, :monthid])
      sort!(testcrspout, [:permno, :monthid])
      @assert size(testcrspout,1) ≤ 30_000

      testcrspout |> CSV.write(
        "$(PARAM[:testpath])\\$(binprefix)dailyformonthly_$(PARAM[:crspdaily]).csv")
    end

    crspm = monthlyop(crspm)

    return crspm
  end #end do block

  #println(describe(crsp))
  return crsp
end


function debugmsg(hdr, df)
  #=meanret = mean(skipmissing(df.ret)) |> (x) -> round(x,sigdigits=12)
  meanprice = mean(skipmissing(df.price)) |> (x) -> round(x,sigdigits=12)
  meanvol = mean(skipmissing(df.dvol)) |> (x) -> round(x,sigdigits=12)
  meanpermno = mean(df.permno) |> (x) -> round(x,sigdigits=12)=#
  print("$hdr hdr: size(df): $(size(df))")
  for s ∈ propertynames(df)
    (!(eltype(df[!,s]) <: Union{Real, Missing})) && continue
    all(df[!,s] .=== missing) && continue
    occursin("placebo", string(s)) && continue
    print(" $s $(mean(skipmissing(df[!, s]))),")
  end
  print("\n")
end


#aggregates a daily dataframe into a monthly dataframe
function makemonthly(selectionop::Function, crspd::AbstractDataFrame;
  computeivol::Bool)

  #Assign each date to a monthid
  crspd.monthid = assignmonthid(crspd.date)
  crspd.retP1 = 1.0 .+ crspd.ret #convert to gross returns for geometric compounding
  crspd.absret = (abs).(crspd.ret)

  issorted(crspd, [:permno, :date]) || error("crspd must be sorted in makemonthly")
  computemonthend!(crspd)
  trailingsigma!(crspd, Fendpoint=:monthend)
  select!(crspd, Not(:monthend))

  #println(describe(crspd))

  computemonthend!(crspd)
  if computeivol
    regressatperiodend!(crspd,
      minpoints=PARAM[:sigmaminpoints],
      maxpoints=PARAM[:sigmamaxpoints], #this is mainly for technical reasons
      calendarinterval = PARAM[:sigmacalendarinterval],
      FXs = PARAM[:ivolfields],
      Fagg = :ivol,
      aggfunc = std,
      Fendpoint = :monthend,
      #=aggfunc = (lm::FMLM)->std(lm.ε)=#)

  else
    @info "ivol not computed"
    crspd.ivol = missings(Float64, nrow(crspd))
  end


  #sanity checks
  #debugmsg("dbg1", crspd)
  #drop rows with too much missing data
  @assert issorted(crspd, [:permno, :date])
  oldnrow = nrow(crspd)
  crspd::DataFrame = crspd[completecases(crspd, PARAM[:monthlycrsprequired]), :]
  @assert (size(crspd,1) ≥ oldnrow * 0.5) "monthlycrsprequired dropped more than 50% of rows!!!
    (oldnrow=$oldnrow, size(crspd)=$(size(crspd)))"
  #debugmsg("dbg2", crspd)
  valid = selectionop(crspd, Ftimegroup=:monthid, datefrequency=:month)
  #debugmsg("dbg2.5", valid)
  #valid = combine(groupby(dvalid[!, [:permno, :monthid, ]], [:permno, :monthid]), nrow => :ndays)
  #valid = ensurecompleteness(valid, Ftimegroup=:monthid, datefrequency=:month)
  crspd = innerjoin(crspd, valid, on=[:permno, :monthid], validate=(false,true))
  #debugmsg("dbg3", crspd)
  if !issorted(crspd, [:permno, :date])
    sort!(crspd, [:permno, :date])
  end

  #need to adjust for industry for certain versions of momentum
  crspd.plret = adjustindustry(crspd)
  #crspd.plret = Vector{Float64}(crspd.plret)
  crspd.lplret = (log).(crspd.plret .+ 1)



  print("Time to compute monthly momentum: ")
  @time momentumstrategies!(crspd, Fconditional=:monthend)
#  throw("Stop-inspect momentum")
  select!(crspd, Not(:monthend))

  if !issorted(crspd, [:permno, :date])
    sort!(crspd, [:permno, :date])
  end


  #scrspd::SubDataFrame = view(crspd, completecases(crspd, PARAM[:monthlycrsprequired]), :)
  #@assert size(crspd,1) .* 0.5 .≤ size(scrspd,1) "
  #  monthlycrsprequired dropped more than 50% of rows!!!"
  #debugmsg("dbg4", crspd)
  #combine by month
  Fmomentums = keys(PARAM[:momentumstrategies]) |> collect
  s2crspds::GroupedDataFrame = groupby(crspd, [:permno, :monthid])
  crspm::DataFrame = combine(s2crspds,
    #:permno => length=> :ndays,
    :ndays => last => :ndays,
    :date => last=> :date,
    :sigma=>last=> :sigma,
    :shares=>last=> :shares,

    :dvol => sum=> :dvol,
    :price => last=> :price,
    :tridx => last=> :tridx,

    :mc => last => :mc,
    :retP1 => prod=> :ret, #WARNING- to be adjusted later

    :ivol => last=> :ivol,
    :exchcd => last => :exchcd,
    :siccd => last => :siccd,

    :lplret => sum => :plret, #WARNING- to be adjusted later
    Fmomentums .=> last .=> Fmomentums,
    #  :ret => std=> :std,
    #yzstd = [:hi, :lo, :open, :price, :Ldistadjclose]=>yzstd
    )

  crspm.ret .-= 1.0
  crspm.plret .= exp.(crspm.plret) .- 1.0
  #debugmsg("dbg5", crspm)
  #sanity check
  testcrspms = groupby(crspm, [:permno, :monthid])
  Threads.@threads for k ∈ keys(s2crspds)
    testscrspm = testcrspms[(permno=k.permno,monthid=k.monthid)]
    @assert sum(log.(1 .+ s2crspds[k].ret)) |> exp ≈ testscrspm.ret[1] + 1.0
  end
  #conform the dates
  scrspms::GroupedDataFrame = groupby(crspm, :monthid)
  Threads.@threads for i ∈ 1:length(scrspms)
    scrspm::SubDataFrame = scrspms[i]
    @assert (maximum(scrspm.date) - minimum(scrspm.date)).value < 31
    scrspm.date .= maximum(scrspm.date)
  end

  #debugmsg("dbg6", crspm)
  return crspm
end


#recommend calling this function whenever we want only the ends of the month, and deleting afterwar
function computemonthend!(crspd::DataFrame; createnew::Bool = true)
  #compute change in month
  (:monthend ∉ names(crspd)) || (error(":monthend already in crspd.
    You probably should have deleted it the last time you used it."))
  crspd.monthend = falses(size(crspd,1))
  @assert issorted(crspd, [:permno, :date])

  scrspds = groupby(crspd, [:permno, :monthid])
  Threads.@threads for i::Int ∈ 1:length(scrspds)
    scrspd::SubDataFrame = scrspds[i]
    scrspd[end, :monthend] = true
  end

  return nothing
end

#assigns a monthid for each date
function assignmonthid(dates::AbstractVector{Date}, mindate::Date=Date(1900,1,1))::Vector{Int}
  day(mindate) == 1 || (
    @warn "day of $mindate ≠ 1- mindate must be first day of month")

  #first compute all possible first days of the month
  maxdate = maximum(dates)
  allmonths = mindate:Month(1):firstdayofmonth(maxdate) |> collect

  #create an index to look up the monthid
  monthidindex = Dict(allmonths[mid] => mid for mid ∈ 1:length(allmonths))

  monthid::Vector{Int} = (d::Date->monthidindex[d |> firstdayofmonth]).(dates)

  return monthid
end
