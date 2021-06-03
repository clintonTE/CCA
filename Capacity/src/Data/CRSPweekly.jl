
#takes a crsp daily file and makes it weekly
function makeweekly(selectionop::Function = (noopd(x)=x),
    weeklyop::Function = (noopw(x)=x);
    crspnameout=PARAM[:crspweekly],
    binprefix="w_",
    mainpanel::Bool = true,
    computeivol::Bool = PARAM[:computeivol])

  local crsp::DataFrame


  #WARNING- crspd is NOT valid, below is only to create crspw
  crsp = prepcrsp(crsptype=PARAM[:crspdaily],
    csvextension = CSV_EXTENSION,
    incsvstream = IN_CSV_STREAM,
    binprefix = binprefix,
    refreshcrsp=PARAM[:refreshweekly],
    crspcolumns=[PARAM[:weeklycrspcols]; keys(PARAM[:momentumstrategies]) |> collect],
    crsprequiredcolumns=PARAM[:dailycrsprequired]) do crspd

    @assert (PARAM[:crspfrequencysuffix] == "w") || (!mainpanel)
    sort!(crspd, [:permno, :date]) #cache the sort

    #computerangemeasures!(crspd, sorted=true) #not used, mainly for kyle's lamda
    crspd.dvol = computedollarvol(crspd)
    rename!(crspd, :vol=>:svol)
    #crspd.shares = computeshares(crspd)

    #NOTE: testing note- change to trailingsigmatest! to test vol routine
    #println(" checksum: $(sum(skipmissing(crspd.sigma)))")

    crspd = mergeff(crspd, Fgroup=:permno, datefrequency=:day)
    if !issorted(crspd, [:permno, :date])
      sort!(crspd, [:permno, :date])
    end
    #crspd = dailyop(crspd)

    crspw::DataFrame = makeweekly(selectionop, crspd; computeivol)
    validatedates(crspw.date |> unique |> sort!; frequency=:week, validateexactfrequency=true)

    sort!(crspw, [:permno, :date])
    validatesummaryintegrity(crspw, Fid=:weekid)

    if PARAM[:testweeklycsv]
      testpermno::Int = PARAM[:testpermnomult]
      #testcrspw::DataFrame = crspw[crspw.permno .% testpermno .== 0, :]
      testcrspw::DataFrame = crspw[crspw.permno .== 14541, :] |> deepcopy
      #if size(testcrspw,1)>1_000 #prevents scenarios where the file is larger than needed
      #  testcrspw=testcrspw[1:1_000]
      #end

      rename!(testcrspw, :date=>:datew, :dvol=>:dvolw, :ivol=>:ivolw,
        :ret=>:retw, :price=>:pricew, :tridx=>:tridxw,
        :shares=>:sharesw, :mc=>:mcw, :sigma=>:sigmaw,
        :plret=>:plretw, :exchcd=>:exchcdw, :siccd=>:siccdw,)
      rename!(testcrspw, (k->k=>Symbol(k, :w)).(keys(PARAM[:momentumstrategies])))


      #make sure we merge using the same df as used to create the weekyl df
      testcrspout = innerjoin(
        testcrspw, crspd,#view(crspd,completecases(crspd, PARAM[:weeklycrsprequired]),:),
        on=[:permno, :weekid])
      sort!(testcrspout, [:permno, :weekid])
      testcrspout |> CSV.write(
        "$(PARAM[:testpath])\\$(binprefix)dailyforweekly_$(PARAM[:crspdaily]).csv")

    end

    crspw = weeklyop(crspw)

    return crspw
  end #end do block

  #println(describe(crsp))
  return crsp
end

#aggregates a daily dataframe into a weekly dataframe
function makeweekly(selectionop::Function, crspd::DataFrame;
  computeivol::Bool)

  #Assign each date to a weekid
  crspd.weekid = assignweekid(crspd.date)
  crspd.retP1 = 1.0 .+ crspd.ret #convert to gross returns for geometric compounding
  crspd.absret = (abs).(crspd.ret)

  issorted(crspd, [:permno, :date]) || error("crspd must be sorted in makeweekly")
  computeweekend!(crspd)
  trailingsigma!(crspd, Fendpoint=:weekend)
  select!(crspd, Not(:weekend))

  #compute change in week
  computeweekend!(crspd)
  if computeivol
    regressatperiodend!(crspd,
      minpoints=PARAM[:sigmaminpoints],
      maxpoints=PARAM[:sigmamaxpoints], #this is mainly for technical reasons
      calendarinterval = PARAM[:sigmacalendarinterval],
      FXs = PARAM[:ivolfields],
      Fagg = :ivol,
      aggfunc = std,
      Fendpoint = :weekend,
      #=aggfunc = (lm::FMLM)->std(lm.ε)=#)
  else
    @info "ivol not computed"
    crspd.ivol = missings(Float64, nrow(crspd))
  end



  #drop rows with too much missing data and apply selection criteria (e.g. min price, etc)
  @assert issorted(crspd, [:permno, :date])
  oldnrow = nrow(crspd)
  crspd::DataFrame = crspd[completecases(crspd, PARAM[:weeklycrsprequired]), :]
  @assert (size(crspd,1) ≥ oldnrow * 0.5) "weeklycrsprequired dropped more than 50% of rows!!!"

  valid = selectionop(crspd, Ftimegroup=:weekid, datefrequency=:week)
  #valid = combine(groupby(dvalid[!, [:permno, :weekid, ]], [:permno, :weekid]), nrow => :ndays)
  #valid = ensurecompleteness(valid, Ftimegroup=:weekid, datefrequency=:week)
  crspd = innerjoin(crspd, valid, on=[:permno, :weekid], validate=(false,true))

  crspd.plret = adjustindustry(crspd)
  #crspd.plret = Vector{Float64}(crspd.plret)
  crspd.lplret = (log).(crspd.plret .+ 1.0)
  if !issorted(crspd, [:permno, :date])
    sort!(crspd, [:permno, :date])
  end

  momentumstrategies!(crspd, Fconditional=:weekend)
  select!(crspd, Not(:weekend))

  if !issorted(crspd, [:permno, :date])
    sort!(crspd, [:permno, :date])
  end

  #combine by week
  Fmomentums = keys(PARAM[:momentumstrategies])  |> collect
  s2crspds::GroupedDataFrame = groupby(crspd, [:permno, :weekid])
  crspw::DataFrame = combine(s2crspds,
    #:permno => length=> :ndays,
    :ndays => last => :ndays,
    :date => last=> :date,
    :sigma=>last=> :sigma,
    :shares=>last=> :shares,
    #:svol => sum=> :svol,
    :dvol => sum=> :dvol,
    :price => last=> :price,
    :tridx => last => :tridx,
    #:price => mean=> :priceavg,
    :mc => last => :mc,
    :retP1 => prod=> :ret, #WARNING- to be adjusted later
    #:absret=>sum=> :sumabsret,
    :ivol => last=> :ivol,
    :exchcd => last => :exchcd,
    :siccd => last => :siccd,
    :lplret => sum => :plret, #WARNING- to be adjusted later
    Fmomentums .=> last .=> Fmomentums,
    #  :ret => std=> :std,
    #yzstd = [:hi, :lo, :open, :price, :Ldistadjclose]=>yzstd
    )

  crspw.ret .-= 1.0
  crspw.plret .= exp.(crspw.plret) .- 1.0

  testcrspws = groupby(crspw, [:permno, :weekid])
  Threads.@threads for k ∈ keys(s2crspds)
    testscrspw = testcrspws[(permno=k.permno,weekid=k.weekid)]
    @assert sum(log.(1 .+ s2crspds[k].ret)) |> exp ≈ testscrspw.ret[1] + 1.0
  end

  #conform the dates
  scrspws::GroupedDataFrame = groupby(crspw, :weekid)
  Threads.@threads for i ∈ 1:length(scrspws)
    scrspw::SubDataFrame = scrspws[i]
    @assert (maximum(scrspw.date) - minimum(scrspw.date)).value < 7
    scrspw.date .= maximum(scrspw.date)
  end


  return crspw
end


#recommend calling this function whenever we want only the ends of the week, and deleting afterwar
function computeweekend!(crspd::DataFrame; createnew::Bool = true)
  #compute change in week
  (:weekend ∉ names(crspd)) || (error(":weekend already in crspd.
    You probably should have deleted it the last time you used it."))
  crspd.weekend = falses(size(crspd,1))
  @assert issorted(crspd, [:permno, :date])

  scrspds = groupby(crspd, [:permno, :weekid])
  Threads.@threads for i::Int ∈ 1:length(scrspds)
    scrspd::SubDataFrame = scrspds[i]
    scrspd[end, :weekend] = true
  end

  return nothing
end

#assigns a weekid for each date
function assignweekid(dates::AbstractVector{Date}, mindate::Date=Date(1900,1,1))::Vector{Int}
  dayofweek(mindate) == 1 || (
    throw("day of $mindate ≠ 1, so day of week will not align with dayofweek function"))

  maxdate = maximum(dates)
  weekid::Vector{Int} = (d::Date->(d-mindate).value ÷ 7).(dates)

  return weekid
end

avgabs(v::Vector{MFloat64}) = mean((abs).(v))

function validatesummaryintegrity(crsp::DataFrame;
    Fid::Symbol=error("summary id Fid is required"))
  local scrsps::GroupedDataFrame

  #check weekid=>date
  scrsps = groupby(crsp, Fid)
  Threads.@threads for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]
    (maximum(scrsp.date) ≠ minimum(scrsp.date)) && (
      error("multiple dates found for 1 Fid $Fid"))
  end

  scrsps = groupby(crsp, :date)
  Threads.@threads for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]
    (maximum(scrsp[!, Fid]) ≠ minimum(scrsp[!, Fid])) && (
      error("multiple Fid $Fids found for 1 date"))
  end

  @info "summary data validated using summary id $Fid"
  return nothing
end


function regressatperiodend!(crspd::AbstractDataFrame;
  minpoints::Int = error("minpoints is required"),
  maxpoints::Int = error("maxpoints is required"),
  calendarinterval::DatePeriod = error("calendarinterval is required"),
  Fret::Symbol = :ret,
  Fendpoint::Symbol = error("Fendpoint is required"),
  FXs::Vector{Symbol} =error("FXs are required"),
  aggfunc::Function = error("aggfunc is required"),
  Fagg::Symbol = error("Fagg is required"))

  #local declarations
  local scrspds::GroupedDataFrame

  #this will hold the results
  crspd[!, Fagg] =  Vector{MFloat64}(undef, size(crspd,1))

  #form the regression formula and variable names
  Xexpr::FMExpr = Meta.parse(join((string).(FXs), " + "))
  Xnames::Vector{Symbol} = [:intercept; FXs]
  Yname::Symbol = Fret

  #avoid missing values
  completecrsp::SubDataFrame = view(crspd, completecases(crspd, [Fret; FXs]), :)

  #pre-allocate space before the threaded regressions
  #Since we have variable sized arrays and matrixes, we pre-allocate for each size

  K::Int = length(Xnames)
  Xs::Vector{Matrix{Float64}} = (t::Int->
    Matrix{Float64}(undef, maxpoints, K)).(1:Threads.nthreads())
  ys::Vector{Vector{Float64}} = (t::Int->
    Vector{Float64}(undef, maxpoints)).(1:Threads.nthreads())
  εs::Vector{Vector{Float64}} = (t::Int->
    Vector{Float64}(undef, maxpoints)).(1:Threads.nthreads())

  #first column is the intercept
  for t ∈ 1:Threads.nthreads()
    Xs[t][:, 1] .= ones(maxpoints)
  end


  #helper function to get the residuals
  @inline function resid!(X::AbstractMatrix, y::AbstractVector, m::Vector{Float64}, ε::Vector{Float64})
    ε .= X*m .- y
    return nothing
  end
  #this function gets the beginning of the interval
  @inline getj₀(d::Date, dates::AbstractVector{Date}) = searchsortedfirst(
    dates, d + Day(1) - calendarinterval)

  scrspds = groupby(completecrsp, :permno)
  regressionctr::Vector{Int} = zeros(Int, Threads.nthreads())
  print("Running ivol regressions... ")
  #=@time Threads.@threads =#for i ∈ 1:length(scrspds)
    local scrspd::SubDataFrame = scrspds[i]
    local t::Int = Threads.threadid()


    nrows::Int = size(scrspd,1)
    (nrows < minpoints) && continue

    #loop over possible bounds
    for j::Int ∈ minpoints:nrows
      regressionctr[t] += 1
      #only running regressions at the end of the week
      scrspd[j, Fendpoint] || continue

      #check for stale windows
      j₀::Int = getj₀(scrspd[j, :date], scrspd.date)
      N::Int = j-j₀ +1
      (N < minpoints)  && continue #elliminates windows without enough points

      #get appropriately sized matrices and arrays
      sX::SubArray = view(Xs[t], 1:N, :)
      resize!(ys[t], N)
      resize!(εs[t], N)

      @assert length(ys[t]) == length(εs[t]) == size(sX,1)

      #populate Y vector and X matrix
      ys[t] .= scrspd[j₀:j, Fret]
      for (c, FX) ∈ enumerate(FXs)
        sX[:, c+1] .= scrspd[j₀:j, FX]
      end

      resid!(sX, ys[t], GLM.coef(lm(sX, ys[t])), εs[t])
      scrspd[j, Fagg] = aggfunc(εs[t])
      #scrspd[j, Fagg] = aggfunc(FMLM(Xs[t], ys[t], Xnames=Xnames, Yname=Yname))
    end
  end
  println("$(sum(regressionctr)) regressions complete.")

  return nothing
end
