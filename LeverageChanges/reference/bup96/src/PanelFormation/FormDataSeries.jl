
function formdataseries(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  datapath::String=DATA_PATH,
  levname::String = LEV_NAME,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  outpath::String = OUT_PATH,
  refreshdataseries::Bool = true,
  trialrunonly::Bool = false,
  archive::Bool = true,
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

    lev = maketretseries!(lev, refreshtretcrsp = refreshtretcrsp, crsptype=:monthly)

    winsorizekeydata!(lev)
    constructportfolios!(lev)

    if TEST_OUTPUT
      fyear = 2015
      CSV.write("output\\$(levname)_$fyear.csv", lev[(y->y==fyear).(lev.fyear), :])
    end

    if TEST_OUTPUT
      fyear = 2014
      CSV.write("output\\$(levname)_$fyear.csv", lev[(y->y==fyear).(lev.fyear), :])
    end

    select!(lev, RETAINED_COLS_DATA_SERIES)

    outjlsstream("$datapath\\$levname.jls.lz4") do s
      serialize(s, lev)
    end

    archive && archivefile("$levname.jls.lz4") #never hurts to have a backup
  else

    try
      injlsstream("$datapath\\$levname.jls.lz4") do s
        lev = deserialize(s) #NOTE: DEBUG ONLY
      end
    catch err
      @warn "Could not load $levname. Error: $err"
      @info "Attempting to restore back-up and load from archive..."
      unarchivefile("$levname.jls.lz4")
      injlsstream("$datapath\\$levname.jls.lz4") do s
        lev = deserialize(s) #NOTE: DEBUG ONLY
      end
      @info "Successfully restored from archive."
    end
  end

  return lev
end

function winsorizekeydata!(lev::DataFrame;
    fields2winsorize::AbstractArray{Symbol}=[CQUART_FVALS; ],
    winsorizeprop=WINSORIZE_PROP, groupfield=:fyear)

  #dups will lead to race errors
  @assert length(fields2winsorize) == length(unique(fields2winsorize))

  @mpar for s ∈ fields2winsorize
    local nbefore::Int = length(unique(lev[!, s]))
    winsorize!(lev[!, s], prop = winsorizeprop)
    (length(unique(lev[!, s])) == nbefore) && (@warn "winsorization didn't do anything for field $s.
      Something is probably wrong.")
  end

  return nothing
end


# a couple of helper functions
@inline function rowbeforedatesorted(df::AbstractDataFrame, dt::Date)::Union{Nothing, DataFrameRow}
  local r::Union{DataFrameRow,Nothing}

  i::Int = searchsortedfirst(df.date, dt)

  #check if the index is valid for our purposes
  r = ((i>size(df,1)) || #this means the date or a larger date wasn't found
    (i==1)
    ? nothing : df[i-1,:])

  return r
end

#takes the last merged record for each merged comp record
function condensewithlast!(panel::DataFrame; sortedpermnodate=false,
    condense::Bool = true)::Nothing
  (!sortedpermnodate) && sort!(panel, [:permno, :date])
  panel.annflag = falses(size(panel,1))

  #now collapse the return data
  panelgrp::GroupedDataFrame = groupby(panel, :permno)
  @mpar for i ∈ 1:length(panelgrp)
    spanel::SubDataFrame = panelgrp[i]

    for sspanel ∈ groupby(spanel, :compid)
      sspanel.annflag[end] = true
    end
  end

  #now drop all other rows
  #can either free the memory or do this in place
  if condense
    filter!(r->r.annflag, panel)
    select!(panel, Not([:annflag]))
  end

  return nothing
end

function maketretseries!(lev::DataFrame;
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
      months2stale = MONTHS_2_STALE_SHORT + 12)
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

  #=if TEST_OUTPUT
    CSV.write("output\\maketretseries.csv", lev[(end-10_000):end, :])
    @info "maketretseries.csv written."
  end=#

  condensewithlast!(lev, sortedpermnodate=true)

  #now sort and lag
  sort!(lev, [:permno, :fyear])
  lagwithin!(lev, [:ret12m], :permno, :fyear, sorted=true)
  lagwithin!(lev, [:Lret12m], :permno, :fyear, sorted=true)

  #=if TEST_OUTPUT
    CSV.write("output\\maketretseries-last.csv", lev[(end-10_000):end, :])
    @info "maketretseries-last.csv written."
  end=#

  return lev
end

function maketvolmcapseries!(lev::DataFrame;
  refreshtvolmcapcrsp::Bool = true,
  crsptype::Symbol = :daily,
  months2stale::Int = MONTHS_2_STALE_SHORT)

  crspcolumns::Vector{Symbol} = [:date, :permno, :mktcap, :vol252d]
  periods2stale = Month(months2stale)

  validpermnos::Vector{Int} = unique(lev.lpermno)

  crsp::DataFrame = prepcrsp(
    crspcolumns = crspcolumns,
    crsptype = crsptype,
    refreshcrsp=refreshtvolmcapcrsp,
    validpermnos = validpermnos,
    prefix = "tvolmcap") do df::DataFrame
    computecrspmktcap!(df)
    computenetreturns!(df)
    computetrailingvol!(df, crsptype==:daily ? 252 : 12,
      Fret = :net, Fvol = :vol252d, months2stale = MONTHS_2_STALE_SHORT + 12)
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

  #now supplement the market data
  levgrp::GroupedDataFrame = groupby(lev, :permno)
  mkequityadded::Int = sum((ismissing).(lev.mkequity))
  @mpar for i ∈ 1:length(levgrp)
    slev::SubDataFrame = levgrp[i]

    for sslev ∈ groupby(slev, :compid)
      sslev.last[end] = true

      local Nmissing::Int = sum((ismissing).(sslev.mkequity))
      if Nmissing > 1
        (Nmissing < size(sslev,1)) && error("Inconsistent presence of mkequity")

        dt::Date = eom(sslev.fyrmonth[1])

        #find the first date before the end of the fiscal year if available
        #need to get this from slev
        local r::Union{Nothing, DataFrameRow} = rowbeforedatesorted(slev, dt)
        (!isnothing(r)) && ( dt - periods2stale ≤ r.date) && (sslev.mkequity .= r.mktcap)
      end
    end
  end

  mkequityadded -= sum((ismissing).(lev.mkequity))
  @info "Supplemented $mkequityadded mkequity records using $crsptype data"



  #now that we have patched the market value, drop all other rows
  #can either free the memory or do this in place
  lev = filter!(r->r.last, lev)

  #make the cash field
  lev.mkat = lev.lt .+ lev.mkequity
  lev.mknegcash = 1.0 .- (lev.cash ./ lev.mkat)

  #compute #flev field
  lev.mkflev = lev.netdebt ./ (lev.netdebt .+ lev.mkequity)

  #compute liab
  lev.mkliab = lev.lt ./ lev.mkat

  #make finite (note we winsorize later on)
  @mpar for s ∈ [:mknegcash, :mkflev, :mkliab]
    lev[!,s] = (f::MFloat64->coalesce(isfinite(f),false) ? f : missing).(lev[!,s])
  end

  sort!(lev, [:permno, :fyear]) #not completely sure we need this sort

  #make lagged change in assets as a control variable
  lev.lat = (log).(lev.at)
  differencewithin!(lev, :lat, :gvkey, :fyear, sorted=true)
  lev.Dat = (exp).(lev.Dlat)

  lev.bm = lev.bkequity ./ lev.mkequity

  lagwithin!(lev, [:vol252d, :mkat, :mkequity, :bm, :Dat, :op],
    :permno, :fyear, sorted=true)
  lagwithin!(lev, :Lvol252d, :permno, :fyear, sorted=true)
  #lagwithin!(lev, [:LLvol252d, :LLvol504d], :fyear, sorted=true)



  rename!(lev, [:date=>:dailyenddate, :compid=>:compidtvolmcap])
  select!(lev, Not([:permno, :last])) #dup of lpermno

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
  select!(crsptvol, [:permno, :date, :vol252d])

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

  #this will make the lookups easier
  crsp.dyrmonth = (d->YearMonth(Year(d).value, Month(d).value)).(crsp.date)

  select!(panel, Not([:vol252d, :ret12m])) #drop these since we will work with the monthly versions
  panel = mergecrsp(crsp, panel, refreshmerge=true, cachedata = false,
    crsptype=crsptype, prefix="panel")

  #preallocate
  N::Int = size(panel,1)
  panel.L12ret12m = Vector{MFloat64}(undef, N)
  panel.L24ret12m = Vector{MFloat64}(undef, N)
  panel.fyret12m = Vector{MFloat64}(undef, N)

  panel.L12vol252d = Vector{MFloat64}(undef, N)
  panel.L24vol252d= Vector{MFloat64}(undef, N)
  panel.fyvol252d= Vector{MFloat64}(undef, N)

  panel.fyoffset = panel.dyrmonth .- panel.fyrmonth

  sort!(panel, [:permno, :dyrmonth])
  @info "beginning panel creation"
  spanels::GroupedDataFrame = groupby(panel, [:permno])
  #=@mpar =#for i ∈ 1:length(spanels)
    local spanel::SubDataFrame = spanels[i]
    local Nspanel::Int = size(spanel,1)
    local ptr::Int
    local target::YearMonth

    #this should be true always
    (Nspanel == length(Set(spanel.dyrmonth))) || error("Dates for dyrmonth are not unique w/in permno")


    for (j, r) ∈ enumerate(eachrow(spanel))

      if j > 12 #compute teh 1yr lagged controls
        target = r.dyrmonth - 12
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L12ret12m = spanel[ptr, :ret12m]
          r.L12vol252d = spanel[ptr, :vol252d]
        end
      end

      if j > 24 #compute teh 2yr lagged controls
        target = r.dyrmonth - 24
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L24ret12m = spanel[ptr, :ret12m]
          r.L24vol252d = spanel[ptr, :vol252d]
        end
      end

      if j > r.fyoffset #compute teh 2yr lagged controls
        target = r.dyrmonth - r.fyoffset
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.fyret12m = spanel[ptr, :ret12m]
          r.fyvol252d = spanel[ptr, :vol252d]
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
