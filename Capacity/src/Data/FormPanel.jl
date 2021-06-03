
function formpanel(;
  panelname::String = PARAM[:panelname],
  savepath::String = PARAM[:workingpath],
  testdate::Date = PARAM[:testdate],
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  binextension::String = BIN_EXTENSION,
  includecomp::Bool = PARAM[:compincludecomp]
  )::DataFrame

  #housekeeping
  local panel::DataFrame
  local crsp::DataFrame
  local panelpath::String = "$savepath\\$panelname.$binextension"

  if PARAM[:crspfrequencysuffix] == "w"
    @info "loading weekly data"
  elseif PARAM[:crspfrequencysuffix] == "m"
    @warn "Using monthly data"
  else
    @assert false
  end

  #The below block directs the panel creation
  if PARAM[:refreshpanel]

    #quick way to filter the dates
    if PARAM[:crspfrequencysuffix] == "w"
      crsp = makeweekly(selectsummarydata)
    elseif PARAM[:crspfrequencysuffix] == "m"
      crsp = makemonthly(selectsummarydata)
    end

    if includecomp #allows us to short circuit the merge if not using comp data
      comp = prepcomp()
      panel = mergecrspcomp(;crsp, comp)
    else
      panel = crsp
    end
    @assert includecomp == (:gvkey ∈ propertynames(panel))
    @info "$(includecomp ? "comp rows included! unmerged rows dropped" :
      "comp excluded, so rows that would have been unmerged were kept.")"

    panel = processpanel(panel)

    if PARAM[:testpanel]

      #time series test output
      testpermno::Int = PARAM[:testpermnomult]
      testview::SubDataFrame = view(panel, (panel.permno .% testpermno) .== 0, :)
      if size(testview,1)>20_000 #prevents scenarios where the file is larger than needed
        testview=view(testview, 1:20_000,:)
      end
      testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_panel.csv")

      #cross-sectional test output
      testview = view(panel, panel.date .== testdate,:)
      testview |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_panel_$testdate.csv")

    end

    @info "Panel formed of size $(size(panel))."

    outbinstream(panelpath, panel)
  else

    panel = inbinstream(panelpath)
  end

  return panel
end

function boundpanelbydate!(panel::DataFrame;
  lowerdatebound::NDate=PARAM[:lowerdatebound],
  upperdatebound::NDate=PARAM[:upperdatebound])
  if !(lowerdatebound === nothing)
    @info "filtering panel to date≥$(lowerdatebound) (current lowest: $(minimum(panel.date)))"
    filter!(:date=>d->d≥lowerdatebound, panel)
  end
  if !(upperdatebound === nothing)
    @info "filtering panel to date≤$(upperdatebound) (current highest: $(maximum(panel.date)))"
    filter!(:date=>d->d≤upperdatebound, panel)
  end

end

function placebostrategies!(panel::DataFrame)
  N=size(panel,1)
  panel.placebo = rand(N) |> Vector{MFloat64}
  panel.placebofixed = Vector{MFloat64}(undef, N)
  panel.placebofixeddiscrete = Vector{MFloat64}(undef, N)

  discretevector::Vector{Float64} = collect(0.5:0.5:5.0)
  for spanel ∈ groupby(panel, :permno)
    spanel.placebofixed .= rand()
    spanel.placebofixeddiscrete .= rand(discretevector)
  end

  return nothing
end


#iterate through the momentums trategies and put them into the dataframe
function fundamentalstrategies!(panel::DataFrame;
    strategies=PARAM[:fundamentalstrategies],
    bookequityfloor::Float64 = PARAM[:compbookequityfloor])


    #book to market
    panel.bm = finiteormissing.(panel.bkequity .* 1_000_000 ./ panel.mc)


    ######testing code for comparison to Lewellen
    assetfloor(x) = max(x,bookequityfloor)#max(x,bookequityfloor)
    assetfloor(::Missing) = missing
    panel.fyear = year.(panel.fdate)
    transform!(groupby(panel, [:permno, :fyear]), [:Lmc, :bkequity] |> AsTable =>
      ((t)->log(assetfloor(t.bkequity[1]) * 1_000_000 / t.Lmc[1])) =>
      :lbLmquarterly)
    @assert length(panel.lbLmquarterly) == size(panel,1)
    panel.lbLm = log.(panel.bkequity .* 1_000_000 ./ panel.Lmc) .|> finiteormissing
    #######


    winsorizewithin!(panel, Fsource=:ret, prop=PARAM[:qwinsorize], Fgroup=:date,
      twosided=true, Fnew=:WXret)

    for Fstrategy ∈ strategies
      FWstrategy::Symbol = Symbol(:W, Fstrategy)
      FWXstrategy::Symbol = Symbol(:WX, Fstrategy)
      panel[!,FWstrategy] = winsorizequantile(panel[!, Fstrategy], PARAM[:qwinsorize], twosided=true)


      winsorizewithin!(panel, Fsource=Fstrategy, prop=PARAM[:qwinsorize], Fgroup=:date,
        twosided=true, Fnew=FWXstrategy)
      lagwithin2!(panel, [Fstrategy, FWstrategy, FWXstrategy], :permno, date=:date) #this could be parallelized

      #NOTE- fairly sure we don't need the below
      #FWiivolWstrategy = Symbol(:Wiivol, FWstrategy)
      #panel[!,Symbol(:Wiivol, FWstrategy)] =  panel[!, :Wiivol] .* panel[!, FWstrategy]
      #panel[!,Symbol(:LWiivol, FWstrategy)] =  panel[!, :LWiivol] .* panel[!, Symbol(:L, FWstrategy)]
    end

end

function processpanel(panel; momentumstrategies=PARAM[:momentumstrategies],
    includecomp::Bool = PARAM[:compincludecomp])
  #winsorize and compute the nomalized volume
  panel.Wdvol = (winsorizequantile(
    panel.dvol, PARAM[:qwinsorize], twosided=true) ./ PARAM[:volscale])

  panel.dvolscaled = panel.dvol ./ PARAM[:volscale]
  differencewithin2!(panel, [:Wdvol], :permno, date=:date)
  #computenvol!(panel, :Wdvol)

  #this forms a market cap control and the basis of the size characteristic
  panel.lmc = log.(panel.mc) .|> finiteormissing
  panel.mcscaled = panel.mc ./ PARAM[:volscale]
  lagwithin2!(panel, [:mc, :lmc, :mcscaled], :permno, date=:date)
  lagwithin2!(panel, [:Lmc, :Llmc, :Lmcscaled], :permno, date=:date)
  #transform!(groupby(panel, :permno),
  #  )
  transform!(groupby(panel, :permno),
    :lmc => ((lmc)->mean(skipmissing(lmc))) => :fixedlmc,
    :mc => ((mc)->mean(skipmissing(mc))) => :fixedmc)

  panel.vecof1s = ones(size(panel,1)) #used for fixed (constant)

  #winsorize idiosyncratic volatiltiy
  if (panel.ivol .!== missing ) |> all
    panel.Wivol = winsorizequantile(panel.ivol, PARAM[:qwinsorize], twosided=true)
    panel.Wiivol = 1.0 ./ panel.Wivol
    lagwithin2!(panel, [:Wivol, :Wiivol], :permno, date=:date)
  else
    select!(panel, Not(:ivol))
  end

  #form the momentum characteristic
  #panel=panel[panel.ndays .≥ 10, :] #WARNING- comment this
  if PARAM[:momentumstrategyfrequency] === :aggregate
    select!(panel, Not([:plret; keys(PARAM[:momentumstrategies])|>collect]))
    panel.lret = (log).(panel.ret .+ 1)
    panel.plret = adjustindustry(panel, Fmc=:mc)
    panel.lplret = (log).(panel.plret .+ 1)
    momentumstrategies!(panel)
  elseif PARAM[:momentumstrategyfrequency] === :day
    @assert setdiff([:plret; keys(PARAM[:momentumstrategies])|>collect], propertynames(panel)) |> isempty
  elseif PARAM[:momentumstrategyfrequency] === :test
    rename!(panel, Dict(s=>Symbol(:d_,s) for s ∈ [:plret; keys(PARAM[:momentumstrategies])|>collect]))
    panel.lret = (log).(panel.ret .+ 1)
    panel.plret = adjustindustry(panel, Fmc=:mc)
    panel.lplret = (log).(panel.plret .+ 1)
    momentumstrategies!(panel)
    throw("aggregate momentum strategies computed")
  else
    throw("unrecognized momentum strategy frequency (:day, :aggregate) are acceptable")
  end

  #use this for within-winsorization
  panel.yearmonth = Symbol.(:Y, year.(panel.date), :M, month.(panel.date))
  #momentumstrategies!(panel)
  conditionmomentumstrategies!(panel)


  if includecomp
    fundamentalstrategies!(panel)
  end

  #form characteristics normalized by market value


  #form market value characteristic
  panel.mcthorn = finiteormissing.(panel.mc) ./ 10.0^9
  lagwithin2!(panel, [:mcthorn], :permno, date=:date)
  lagwithin2!(panel, [:Lmcthorn], :permno, date=:date)

  boundpanelbydate!(panel) #restrict the dates
  placebostrategies!(panel) #create placebo characteristics


  return panel
end


#=function ensurecompleteness(crsp::DataFrame; Ftimegroup::Symbol,
  datefrequency::Symbol)

  #drop all time groups (month or week) with insufficient data
  #use a loose absolute criteria and a stricter criteria based on a fraction of the median
  mindaysmedianfrac::Float64 = PARAM[:mindaysmedianfrac][datefrequency]
  absolutemindays::Int = PARAM[:absolutemindays][datefrequency]
  crsp.tokeep = trues(size(crsp,1))
  scrsps::GroupedDataFrame = groupby(crsp, Ftimegroup)
  Threads.@threads for i ∈ 1:length(scrsps)
    scrsp::SubDataFrame = scrsps[i]
    #now employ some heuristics to make sure we don't have too much missing data
    medianndays::Int = median(scrsp.ndays)
    #minndays::Int = medianndays ≥ 2 ? medianndays - 1 : medianndays
    minndays::Int = max(ceil(medianndays*mindaysmedianfrac) |> Int, absolutemindays)
    scrsp.tokeep .= scrsp.ndays .≥ minndays
  end

  crsp = crsp[crsp.tokeep, :]
  select!(crsp, Not(:tokeep))

  return crsp
end=#


function selectsummarydata(crsp::DataFrame;
  #NOTE: Probably ok for now, but consider doing these in quantiles by cross-section instead
  minprice::Float64 = PARAM[:minprice],
  minmarketcap::Float64 = PARAM[:minmarketcap],
  minnysequantile::Float64 = PARAM[:minnysequantile],
  Ftimegroup::Symbol,
  datefrequency::Symbol)

  #@assert all((!ismissing).(crsp.mc))
  #@assert all((!ismissing).(crsp.exchcd))
  crsp.mc = Vector{Float64}(crsp.mc)
  crsp.exchcd = Vector{Int}(crsp.exchcd)

  function assignnysequantile(t)
    nysemc = t.mc[t.exchcd .== 1]
    nysecdf = ecdf(nysemc)
    return nysecdf.(t.mc)
  end

  @assert issorted(crsp, [:permno, :date])
  @assert issorted(crsp, [:permno, Ftimegroup])
  filtercols = [:permno, Ftimegroup, :mc, :exchcd, :price]
  @assert setdiff(filtercols, propertynames(crsp)) |> isempty
  crspeop::DataFrame = combine(groupby(crsp[!,filtercols], [:permno, Ftimegroup]),
    filtercols .=> last,
    nrow => :ndays, renamecols=false)
  #debugmsg("dbg2.01", crspeop)
  @assert issorted(crspeop, [:permno, Ftimegroup])
  transform!(groupby(crspeop, Ftimegroup), [:mc, :exchcd] |> AsTable =>
    assignnysequantile => :nysequantile)
  #debugmsg("dbg2.02", crspeop)
  nysebreakpoint(t) = minimum(t.mc[t.nysequantile .≥ minnysequantile])
  transform!(groupby(crspeop, Ftimegroup), [:mc, :nysequantile] |> AsTable =>
    nysebreakpoint => :breakpoints)
  #debugmsg("dbg2.03", crspeop)
  #compute quantiles based on nyse thresholds
  averagebreakpoint::Float64 = mean(combine(groupby(crspeop, Ftimegroup), :breakpoints =>
    last, renamecols=false).breakpoints)

  @info "Average NYSE breakpoint: $averagebreakpoint"

  #question- should I filter out firms
  filter!(AsTable([:price, :mc, :nysequantile]) =>
    (r ->
    (r.price ≥ minprice) &&
    (r.mc ≥ minmarketcap) &&
    (r.nysequantile ≥ minnysequantile)), crspeop)

  @assert all(crspeop.mc .≥ crspeop.breakpoints)
  crspeop = crspeop[!, [:permno, Ftimegroup, :ndays, :nysequantile]]
  #debugmsg("dbg2.04", crspeop)
  #drop all time groups (month or week) with insufficient data
  #use a loose absolute criteria and a stricter criteria based on a fraction of the median
  #mindaysmedianfrac::Float64 = PARAM[:mindaysmedianfrac][datefrequency]
  absolutemindays::Int = PARAM[:absolutemindays][datefrequency]


  medianndays(ndays) = max(floor(median(ndays)) - 1 |> Int,1)
  transform!(groupby(crspeop, Ftimegroup), :ndays => medianndays => :medianndays)
  crspeop = crspeop[(crspeop.ndays .≥ crspeop.medianndays) .& (crspeop.medianndays .≥ absolutemindays), :]
  select!(crspeop, Not([:medianndays]))

  #debugmsg("dbg2.05", crspeop)


  return crspeop
end
