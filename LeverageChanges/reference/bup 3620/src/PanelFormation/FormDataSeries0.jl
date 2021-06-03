
function formdataseries(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  datapath::String=DATA_PATH,
  levname::String = LEV_NAME,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  binextension::String = BIN_EXTENSION,
  outpath::String = OUT_PATH,
  refreshdataseries::Bool = true,
  trialrunonly::Bool = false,
  retaineddataseriescolumns = RETAINED_COLS_DATA_SERIES)::DataFrame

  local lev::DataFrame
  local fyear::Int

  local dailyunlesstrial::Symbol = trialrunonly ? :monthly : :daily
  refreshdataseries = (refreshdataseries || refreshcomp)


  levname = "$levname-$dailyunlesstrial"

  if refreshdataseries
    lev = prepcomp(refreshcomp=refreshcomp)

    #supplements the market value and constructs the trailing daily vol
    lev = maketvolmcapseries!(
      lev, refreshtvolmcapcrsp=refreshtvolmcapcrsp, crsptype=dailyunlesstrial)

    lev = marketretseries!(lev, refreshtretcrsp = refreshtretcrsp, crsptype=:monthly)

    issorted(lev, [:permno, :date]) || error("lev should be sorted here")
    CSV.write("output\\$(levname)_ts.csv", lev[1:20_000, :])

    winsorizekeydata!(lev)
    constructportfolios!(lev)
    doubleincreasegroups!(lev)

    issorted(lev, [:permno, :date]) || sort!(lev, [:permno, :date])
    if TEST_OUTPUT
      fyear = 2015
      CSV.write("output\\$(levname)_$fyear.csv", lev[(y->y==fyear).(lev.fyear), :])

      fyear = 2014
      CSV.write("output\\$(levname)_$fyear.csv", lev[(y->y==fyear).(lev.fyear), :])
    end

    select!(lev, RETAINED_COLS_DATA_SERIES)

    outbinstream("$datapath\\$levname.$binextension", lev)

  else

    lev = inbinstream("$datapath\\$levname.$binextension")

  end

  return lev
end


function winsorizekeydata!(lev::DataFrame;
    fields2winsorize::AbstractArray{Symbol}=[CQUART_FVALS; CQUART_FVALS_D; CQUART_FVALS_BRS_D],
    winsorizeprop=WINSORIZE_PROP, groupfield=:fyear)

  #dups will lead to race errors
  @assert allunique(fields2winsorize)

  @mpar for s ∈ fields2winsorize
    local nbefore::Int = length(unique(lev[!, s]))
    winsorize!(lev[!, s], prop = winsorizeprop)
    (length(unique(lev[!, s])) == nbefore) && (@warn "winsorization didn't do anything for field $s.
      Something is probably wrong.")
  end

  return nothing
end


#takes the last merged record for each merged comp record
function condensewithlast!(panel::DataFrame;
    condense::Bool = true, writepartout::Bool = false)::Nothing
  issorted(panel, [:permno, :date]) || error("panel must be sorted")
  panel.annflag = falses(size(panel,1))

  #now collapse the return data
  panelgrp::GroupedDataFrame = groupby(panel, :permno)
  @mpar for i ∈ 1:length(panelgrp)
    spanel::SubDataFrame = panelgrp[i]

    for sspanel ∈ groupby(spanel, :compid)
      sspanel.annflag[end] = true
    end
  end

  if writepartout
    panel[1:10_000,:] |> CSV.write("output\\condensewlast.csv")
  end

  #now drop all other rows
  #can either free the memory or do this in place
  if condense
    filter!(r->r.annflag, panel)
    select!(panel, Not([:annflag]))
  end

  return nothing
end

function marketretseries!(lev::DataFrame;
  refreshtretcrsp::Bool = true,
  crsptype::Symbol = :daily,
  months2stale::Int = MONTHS_2_STALE_SHORT)

  crspcolumns::Vector{Symbol} = [:date, :permno, :ret, :lret,
    :ret12m, :lret12m#=, :net12m, :lnet12m=#]

  periods2stale = Month(months2stale)

  validpermnos::Vector{Int} = unique(lev.lpermno)

  crsp::DataFrame = prepcrsp(
    crspcolumns = crspcolumns,
    crsptype = :monthly,
    refreshcrsp=refreshtretcrsp,
    validpermnos = validpermnos,
    prefix = "tret") do df::DataFrame
    computenetreturns!(df)
    computetrailingreturns!(df, 12, Flret = :lret, Fnewret = :ret12m,
      months2stale = MONTHS_2_STALE_SHORT + 12,
      minpointsperwindow=REQUIRED_MONTHLY_TRAILING)
    #computetrailingreturns!(df, 12, Flret = :lnet, Fnewret = :net12m,
    #  months2stale = MONTHS_2_STALE_SHORT + 12)
  end

  lev = mergecrsp(crsp, lev, refreshmerge=true, cachedata = false, crsptype=crsptype, prefix="tret")
  obliterate!(crsp)

  #now collapse the return data
  #***NOTE***: We are collecting the latest entry for each period
  #if periods are at irregular intervals, this is technically different than what is stated in the paper
  #but it fully avoids contaminated date announcements
  #the alternative might be to drop all such overlapping entries, or make the intervals irregular
  sort!(lev, [:permno, :date])

  if TEST_OUTPUT
    CSV.write("output\\marketretseries.csv", lev[(end-10_000):end, :])
    @info "marketretseries.csv written."
  end

  #compute the returns over the fiscal year
  lev.fret = Vector{MFloat64}(missing, size(lev,1))
  slevs::GroupedDataFrame = groupby(lev, [:permno])
  @mpar for i ∈ 1:length(slevs)
    slev::SubDataFrame = slevs[i]
    for sslev ∈ groupby(slev, :compid)
      fdate::Date = sslev.fdate[end]
      Lfdate::MDate = sslev.Lfdate[end]
      ismissing(Lfdate) && continue
      s3lev::SubDataFrame = view(slev, (r::DataFrameRow->
        (r.date ≤ fdate) && (r.date > Lfdate)).(eachrow(slev)), :)
      sslev.fret .= exp(sum(skipmissing(s3lev.lret))) - 1.0
    end
  end

  condensewithlast!(lev, writepartout = true)

  #now sort and lag
  sort!(lev, [:permno, :fyear])
  lagwithin!(lev, [:ret12m], :permno, :fyear, sorted=true)
  lagwithin!(lev, [:Lret12m], :permno, :fyear, sorted=true)

  #println(describe(lev))
  ##now make the table 5 leverage versions net of stock returns

  #this is the stock return constribution of leverage
  lev.RLmkequity = lev.Lmkequity .* (1 .+ lev.fret)
  lev.RDmkflev = lev.Lnetdebt ./ (lev.Lnetdebt .+ lev.RLmkequity) .- lev.Lmkflev
  lev.RDmkliab = lev.Lnetliab ./ (lev.Llt .+ lev.RLmkequity) .- lev.Lmkliab


  #the residual is the managerial component
  lev.MDmkflev = lev.Dmkflev .- lev.RDmkflev
  lev.MDmkliab = lev.Dmkliab .- lev.RDmkliab



  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\levuniv_$d.csv", lev[(cd->cd==d).(lev.date), :])
  end


  return lev
end

function maketvolmcapseries!(lev::DataFrame;
  refreshtvolmcapcrsp::Bool = true,
  crsptype::Symbol = :daily,
  months2stale::Int = MONTHS_2_STALE_SHORT,
  requiredannualdatacols::Vector{Symbol} = REQUIRED_ANNUAL_DATA_COLS,
  lagmonthrange::UnitRange = LAG_MONTH_RANGE)

  crspcolumns::Vector{Symbol} = [:date, :permno, :mktcap, :vol252dnet, :vol252dret]
  periods2stale = Month(months2stale)

  validpermnos::Vector{Int} = unique(lev.lpermno)

  print("time to load vol: ")
  @time crsp::DataFrame = prepcrsp(
    crspcolumns = crspcolumns,
    crsptype = crsptype,
    refreshcrsp=refreshtvolmcapcrsp,
    validpermnos = validpermnos,
    prefix = "tvolmcap") do df::DataFrame
    computecrspmktcap!(df)
    computenetreturns!(df)
    #winsorize!(df.net, prop=0.999)
    computetrailingvol!(df, crsptype==:daily ? 252 : 12, [:net, :ret],
      Fvols = [:vol252dnet, :vol252dret],
      months2stale = MONTHS_2_STALE_SHORT + 12,
      minpointsperwindow = crsptype==:daily ? REQUIRED_DAILY_TRAILING : REQUIRED_MONTHLY_TRAILING)

    #winsorize per Ivo
    df.vol252dnet .= (f::MFloat64->winsorizelevel(f, 0., 0.1)).(df.vol252dnet)
    df.vol252dret .= (f::MFloat64->winsorizelevel(f, 0., 0.1)).(df.vol252dret)
    #computetrailingvol!(df, crsptype==:daily ? 252 : 12,
    #  Fret = :ret, Fvol = :vol252dret, months2stale = MONTHS_2_STALE_SHORT + 12)

  end

  lev = mergecrsp(crsp, lev,
    refreshmerge=true,
    cachedata = false,
    crsptype=crsptype,
    prefix="tvolmcap")

  obliterate!(crsp)

  #supplement the market value
  sort!(lev, [:permno, :date])
  lev.last = falses(size(lev,1))

  #gets the row in a dataframe given a date or the msot recent past row
  @inline function rowfromdatesorted(df::AbstractDataFrame, dt::Date)::Union{Nothing, DataFrameRow}
    local r::Union{DataFrameRow,Nothing}

    i::Int = searchsortedfirst(df.date, dt)

    #check if the index is valid for our purposes
    r = ((i>size(df,1)) #this means the date or a larger date wasn't found
      ? nothing : df[i,:])

    return r
  end

  slevs::GroupedDataFrame = groupby(lev, :permno)
  mkequityadded::Int = sum((ismissing).(lev.mkequity))
  @mpar for i ∈ 1:length(slevs)
    slev::SubDataFrame = slevs[i]

    for sslev ∈ groupby(slev, :compid)
      sslev.last[end] = true

      #now supplement the market equity
      local Nmissing::Int = sum((ismissing).(sslev.mkequity))
      if Nmissing > 1
        (Nmissing < size(sslev,1)) && error("Inconsistent presence of mkequity")

        dt::Date = eom(sslev.fyrmonth[1])

        #find the first date on or before the end of the fiscal year if available
        local r::Union{Nothing, DataFrameRow} = rowfromdatesorted(slev, dt)

        #note- need to adjust for different scaling between crsp and compustat
        #NOTE: In PreprocessCRSP we scaled mktcap to dollars, now we scale it down to millions
        (!isnothing(r)) && ( dt - periods2stale ≤ r.date) && (sslev.mkequity .= r.mktcap/1_000_000)
      end
    end
  end

  mkequityadded -= sum((ismissing).(lev.mkequity))
  @info "Supplemented $mkequityadded mkequity records using $crsptype data"

  #now that we have patched the market value, drop all other rows
  #can either free the memory or do this in place
  lev = filter!(r->r.last, lev)

  #improve mkequity quality, and make sure we have the field
  lev.mkequity .= (floororx).(lev.mkequity)
  lev = lev[(m::MFloat64->(!ismissing(m)) && (m>MIN_MKEQUITY)).(lev.mkequity),:]

  #make the cash field
  lev.mkat = lev.at .- lev.bkequity .+ lev.mkequity #old: lev.lt .+ lev.mkequity
  lev.mkat .= (floororx).(lev.mkat)
  lev.mknegcash = 1.0 .- (lev.cash ./ lev.mkat)
  lev.mknegcash .= (boundunit).(lev.mknegcash)

  #compute liab
  lev.mkliab = lev.netliab ./ lev.mkat #netliab is compa.lt .- compa.cash from preprocesscomp
  lev.mkliab .= (boundunit).(lev.mkliab)

  #compute #flev field
  lev.mkflev = lev.netdebt ./ (lev.netdebt .+ lev.mkequity)
  lev.mkflev .= (boundunit).(lev.mkflev)

  #make finite (note we winsorize later on)
  @mpar for s ∈ [:mknegcash, :mkflev, :mkliab]
    lev[!,s] .= (f::MFloat64->coalesce(isfinite(f),false) ? f : missing).(lev[!,s])
  end

  sort!(lev, [:permno, :fyear]) #not completely sure we need this sort

  lev.bm = lev.bkequity ./ lev.mkequity
  finiteormissing!(lev, :bm)
  lev.bm .= ((x::MFloat64) -> winsorizelevel(x, 0.0, 5.0)).(lev.bm)

  #make lagged change in assets as a control variable
  #println(describe(lev))
  #println("Busy work: $(mean(rand(5000,5000)))")
  lev.lat = (log).(lev.at)
  finiteormissing!(lev, :lat)

  #maybe I should lag the book stuff earlier? maybe not
  differencewithin!(lev, [:mknegcash, :bknegcash, :mkflev,
    :bkflev, :mkliab, :bkliab, :lat, :vol252dnet, :mkequity], :permno, :fyear, sorted=true)

  lev.Dat = (exp).(lev.Dlat)
  lev.lmkequity = (log).(lev.mkequity*1_000_000) #rescale

  finiteormissing!(lev, :lmkequity)
  lagwithin!(lev, [:vol252dnet, :Dvol252dnet,
    :mkat, :mkequity, :lmkequity, :bm, :Dlat, :nop,
    :mknegcash, :bknegcash, :mkflev, :bkflev, :mkliab, :bkliab,
    :Dmknegcash, :Dbknegcash, :Dmkflev, :Dbkflev, :Dmkliab, :Dbkliab,
    :netdebt, :lt, :at, :netliab],
    :permno, :fyear, sorted=true)
  lagwithin!(lev, :Lvol252dnet, :permno, :fyear, sorted=true)


  #WARNING: FILTERS
  #now we only allow records in Ivo's range
  #WARNING this drops a lot of records. Implemented 2152020
  filter!(lev) do r
    (!ismissing(r[:fyrmonth])) &&
    (!ismissing(r[:Lfyrmonth])) &&
    (r[:fyrmonth] - r[:Lfyrmonth] ∈ lagmonthrange)
  end

  dropmissing!(lev, requiredannualdatacols)

  ######END IVO FILTERS
  rename!(lev, [:date=>:dailyenddate, :compid=>:compidtvolmcap])
  select!(lev, Not([:permno, :last])) #dup of lpermno

  if TEST_OUTPUT
    TEST_OUTPUT && CSV.write(
      "output\\volmcap_$(crsptype).csv",
        lev[(end-20000):end, :])
  end

  return lev
end

#creates the control variables with appropriate lagging
function makepanelmonthly(panel::DataFrame;
  crsptype::Symbol = :daily)

  local crsptvol::DataFrame
  local crsptret::DataFrame
  local crsp::DataFrame

  #load the previous crsp fields
  crsptvol = prepcrsp(
    crspcolumns = nothing,
    crsptype = crsptype,
    refreshcrsp=false,
    validpermnos = nothing,
    prefix = "tvolmcap") do df::DataFrame
      @assert false #do not refresh the crsp file here
    end
  select!(crsptvol, [:permno, :date, :vol252dnet, ])

  crsptret = prepcrsp(
    crspcolumns = nothing,
    crsptype = :monthly,
    refreshcrsp=false,
    validpermnos = nothing,
    prefix = "tret") do df::DataFrame
      @assert false #do not refresh the crsp file here
    end
  select!(crsptret, [:permno, :date, :ret12m, :ret, :lret])

  #use the dialy data to supplement the existing monthly data (new dataset is monthly)
  crsp = join(crsptvol, crsptret, on=[:permno, :date])
  obliterate!(crsptvol)
  obliterate!(crsptret)
  crsp[1:20_000,:] |> CSV.write("$OUT_PATH\\crspret_$(crsptype)_ts.csv")

  #this will make the lookups easier
  crsp.dyrmonth = (d->YearMonth(Year(d).value, Month(d).value)).(crsp.date)

  select!(panel, Not([:vol252dnet, :ret12m])) #drop these since we will work with the monthly versions
  panel = mergecrsp(crsp, panel, refreshmerge=true, cachedata = false,
    crsptype=crsptype, prefix="panel")

  #preallocate
  N::Int = size(panel,1)
  panel.L12ret12m = Vector{MFloat64}(undef, N)
  panel.L24ret12m = Vector{MFloat64}(undef, N)
  panel.ayret12m = Vector{MFloat64}(undef, N)

  panel.L12vol252dnet = Vector{MFloat64}(undef, N)
  panel.L24vol252dnet= Vector{MFloat64}(undef, N)
  panel.ayvol252dnet= Vector{MFloat64}(undef, N)

  panel.ayoffset = panel.dyrmonth .- panel.ayrmonth

  sort!(panel, [:permno, :dyrmonth])
  @info "beginning panel creation"
  spanels::GroupedDataFrame = groupby(panel, [:permno])
  #=@mpar =#for i ∈ 1:length(spanels)
    local spanel::SubDataFrame = spanels[i]
    local Nspanel::Int = size(spanel,1)
    local ptr::Int
    local target::YearMonth

    #this should be true always
    (allunique(spanel.dyrmonth)) || error("Dates for dyrmonth are not unique w/in permno")


    for (j, r) ∈ enumerate(eachrow(spanel))

      if j > 12 #compute teh 1yr lagged controls
        target = r.dyrmonth - 12
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L12ret12m = spanel[ptr, :ret12m]
          r.L12vol252dnet = spanel[ptr, :vol252dnet]
        end
      end

      if j > 24 #compute teh 2yr lagged controls
        target = r.dyrmonth - 24
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L24ret12m = spanel[ptr, :ret12m]
          r.L24vol252dnet = spanel[ptr, :vol252dnet]
        end
      end

      if j > r.ayoffset
        target = r.dyrmonth - r.ayoffset
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.ayret12m = spanel[ptr, :ret12m]
          r.ayvol252dnet = spanel[ptr, :vol252dnet]
        end
      end
    end
  end

  #println(unique(panel.fyoffset))
  #error("stop here")

  #println(describe(panel))
  @info "panel creation complete"

  return panel
end
