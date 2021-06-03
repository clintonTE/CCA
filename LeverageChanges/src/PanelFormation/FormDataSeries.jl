
function formdataseries(;
  refreshcomp::Bool = true,
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

    #lev = marketretseries!(lev, refreshtretcrsp = refreshtretcrsp, crsptype=:monthly)

    issorted(lev, [:lpermno, :fdate]) || error("lev should be sorted here")
    CSV.write("output\\$(levname)_ts.csv", lev[1:20_000, :])

    winsorizekeydata!(lev)
    constructportfolios!(lev)
    doubleincreasegroups!(lev)

    issorted(lev, [:lpermno, :fdate]) || sort!(lev, [:lpermno, :fdate])
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
    winsorizeprop=WINSORIZE_PROP)

  #dups will lead to race errors
  @assert allunique(fields2winsorize)

  @mpar for s ∈ fields2winsorize
    local nbefore::Int = length(unique(lev[!, s]))
    lev[!, s] .= winsorizequantile(lev[!, s], winsorizeprop, twosided=true)
    (length(unique(lev[!, s])) == nbefore) && (@warn "winsorization didn't do anything for field $s.
      Something is probably wrong.")
  end

  return nothing
end

function maketvolmcapseries!(lev::DataFrame;
  refreshtvolmcapcrsp::Bool = true,
  crsptype::Symbol = :daily,
  months2stale::Int = MONTHS_2_STALE_SHORT,
  requiredannualdatacols::Vector{Symbol} = REQUIRED_ANNUAL_DATA_COLS,
  lagmonthrange::UnitRange = LAG_MONTH_RANGE)

  crspcolumns::Vector{Symbol} = [:date, :permno, :mktcap,
    :vol1yearnet, :vol1yearret, :ret1year, :lret, :ret]
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

    #it will be useful to flag end of months (or end of periods if cut short)
    issorted(df, [:permno, :date]) || error("crsp must be sorted")
    df.eom = Vector{Bool}(undef, size(df,1))
    sdfs::GroupedDataFrame = groupby(df, :permno)
    @mpar for i ∈ 1:length(sdfs)
      sdf::SubDataFrame = sdfs[i]
      months::Vector{Int} = (month).(sdf.date)
      sdf.eom .=  [months[1:(end-1)] .≠ months[2:end]; true]
    end
    #compute the trailing rets
    computetrailingreturns!(df, Year(1),  Flret = :lret, Fnewret = :ret1year,
      minpointsperwindow = crsptype==:daily ? REQUIRED_DAILY_TRAILING : REQUIRED_MONTHLY_TRAILING,
      eomonly = true)

    #compute the vols
    computetrailingvol!(df, Year(1), [:net, :ret], Fvols = [:vol1yearnet, :vol1yearret],
      minpointsperwindow = crsptype==:daily ? REQUIRED_DAILY_TRAILING : REQUIRED_MONTHLY_TRAILING,
      eomonly = true)



    #winsorize per Ivo
    df.vol1yearnet .= (f::MFloat64->winsorizelevel(f, 0., 0.1)).(df.vol1yearnet)
    df.vol1yearret .= (f::MFloat64->winsorizelevel(f, 0., 0.1)).(df.vol1yearret)

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


  #compute the returns and standard deviation over the fiscal year
  lev.fret = Vector{MFloat64}(missing, size(lev,1))
  lev.fvol = Vector{MFloat64}(missing, size(lev,1))
  periods2stalelong = Month(MONTHS_2_STALE_LONG)
  slevs = groupby(lev, [:permno]) #refresh the grouped df, not sure if this is necessary
  #=@mpar =#for i ∈ 1:length(slevs)
    slev::SubDataFrame = slevs[i]
    sslevs::GroupedDataFrame = groupby(slev, :compid)
    for j ∈ 1:length(sslevs)
      sslev = sslevs[j]
      fdate::Date = sslev.fdate[end]
      Lfdate::MDate = sslev.Lfdate[end]
      if ismissing(Lfdate)
        #as a check, if the lag is missing look for fdates within the acceptable window for lags
        if (j≠1) && (sslev[1,:fdate] - periods2stalelong ≤ sslevs[j-1][1, :fdate])
          println(view(slev, slev.last, [:permno, :date, :fdate, :Lfdate, :fyear]))
          error("unexpected missing lagged fdate")
        end
        Lfdate = sslev[1,:date] - Day(1)
      end

      #make sure the fdate interval is acceptable
      ((Lfdate) ≥ firstdayofmonth(fdate) - Month(18)) || continue
      ((Lfdate) ≤ lastdayofmonth(fdate) - Month(6)) || continue

      #ismissing(Lfdate) && continue
      s3lev::SubDataFrame = view(slev, (r::DataFrameRow->
        (r.date ≤ fdate) && (r.date > Lfdate)).(eachrow(slev)), :)
      sslev.fret .= exp(sum(skipmissing(s3lev.lret))) - 1.0
      sslev.fvol .= std(skipmissing(s3lev.ret))
    end
  end

  lev[1:20_000,:] |> CSV.write("output\\fretvollevcheck.csv")


  #now that we have patched the market value, drop all other rows
  #can either free the memory or do this in place
  @info "rows before dropping all but last in retvol: $(size(lev,1))"
  filter!(r->r.last, lev)
  @info "rows after dropping all but last in retvol: $(size(lev,1))"

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

  #sort!(lev, [:permno, :fdate]) #not completely sure we need this sort

  lev.bm = lev.bkequity ./ lev.mkequity
  finiteormissing!(lev, :bm)
  lev.bm .= ((x::MFloat64) -> winsorizelevel(x, 0.0, 5.0)).(lev.bm)

  #make lagged change in assets as a control variable
  #println(describe(lev))
  #println("Busy work: $(mean(rand(5000,5000)))")
  lev.lat = (log).(lev.at)
  finiteormissing!(lev, :lat)

  #maybe I should lag the book stuff earlier? maybe not
  issorted(lev, [:permno, :fdate]) || error("lev should be sorted on permno, fdate")
  differencewithin2!(lev, [:mknegcash, :bknegcash, :mkflev,
    :bkflev, :mkliab, :bkliab, :lat, :vol1yearnet, :vol1yearret,:mkequity], :permno, date=:fdate)

  lev.Dat = (exp).(lev.Dlat)
  lev.lmkequity = (log).(lev.mkequity*1_000_000) #rescale

  finiteormissing!(lev, :lmkequity)
  lagwithin2!(lev, [:vol1yearnet, :Dvol1yearnet,:vol1yearret, :Dvol1yearret,
    :mkat, :mkequity, :lmkequity, :bm, :Dlat, :nop,
    :mknegcash, :bknegcash, :mkflev, :bkflev, :mkliab, :bkliab,
    :Dmknegcash, :Dbknegcash, :Dmkflev, :Dbkflev, :Dmkliab, :Dbkliab,
    :netdebt, :lt, :at, :netliab],
    :permno, date=:fdate)
  lagwithin2!(lev, [:Lvol1yearnet, :Lvol1yearret], :permno, date=:fdate)


  ############################ return-related procedures
  #this is remaining content from the old return buildup procedure
  issorted(lev, [:permno, :fdate]) || error("lev should be sorted on permno, fdate")
  lagwithin2!(lev, [:ret1year], :permno, date=:fdate)
  lagwithin2!(lev, [:Lret1year], :permno, date=:fdate)

  #println(describe(lev))
  ##now make the table 5 leverage versions net of stock returns

  #this is the stock return constribution of leverage
  lev.RLmkequity = lev.Lmkequity .* (1 .+ lev.fret)
  lev.RDmkflev = lev.Lnetdebt ./ (lev.Lnetdebt .+ lev.RLmkequity) .- lev.Lmkflev
  lev.RDmkliab = lev.Lnetliab ./ (lev.Llt .+ lev.RLmkequity) .- lev.Lmkliab


  #the residual is the managerial component
  lev.MDmkflev = lev.Dmkflev .- lev.RDmkflev
  lev.MDmkliab = lev.Dmkliab .- lev.RDmkliab

  issorted(lev, [:permno, :fdate]) || error("lev should be sorted on permno, fdate")
  #WARNING: FILTERS
  #now we only allow records in Ivo's range
  #WARNING this drops a lot of records. Implemented 2152020
  #unlike Ivo I drop both next and lagged (C-L and L-LL)
  filter!(lev) do r
    (!ismissing(r[:fyrmonth])) &&
    (!ismissing(r[:Lfyrmonth])) &&
    (!ismissing(r[:LLfyrmonth])) &&
    (r[:fyrmonth] - r[:Lfyrmonth] ∈ lagmonthrange) &&
    (r[:Lfyrmonth] - r[:LLfyrmonth] ∈ lagmonthrange)
  end

  dropmissing!(lev, requiredannualdatacols)

  ######END IVO FILTERS
  rename!(lev, [:date=>:dailyenddate, :compid=>:compidtvolmcap])
  select!(lev, Not([:permno, :last])) #dup of lpermno

  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\levuniv_$d.csv", lev[(cd->cd==d).(lev.dailyenddate), :])
  end

  if TEST_OUTPUT
    TEST_OUTPUT && CSV.write(
      "output\\volretmcap_$(crsptype).csv",
        lev[(end-20000):end, :])
  end

  return lev
end

#creates the control variables with appropriate lagging
function makepanelmonthly(panel::DataFrame;
  crsptype::Symbol = :daily,
  refreshtretcrsp::Bool = true)

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
  select!(crsptvol, [:permno, :date, :vol1yearnet, :vol1yearret,:ret1year])

  crsptret = prepcrsp(
    crspcolumns = nothing,
    crsptype = :monthly,
    refreshcrsp=false,
    validpermnos = nothing,
    prefix = "tret") do df::DataFrame
    computenetreturns!(df)
    end
  select!(crsptret, [:permno, :date, :ret, :lret])

  #use the dialy data to supplement the existing monthly data (new dataset is monthly)
  crsp = innerjoin(crsptvol, crsptret, on=[:permno, :date])
  obliterate!(crsptvol)
  obliterate!(crsptret)
  crsp[1:20_000,:] |> CSV.write("$OUT_PATH\\crspret_$(crsptype)_ts.csv")

  #this will make the lookups easier
  crsp.dyrmonth = (d->YearMonth(Year(d).value, Month(d).value)).(crsp.date)

  select!(panel, Not([:vol1yearnet, :ret1year])) #drop these since we will work with the monthly versions
  panel = mergecrsp(crsp, panel, refreshmerge=true, cachedata = false,
    crsptype=crsptype, prefix="panel")

  #preallocate
  N::Int = size(panel,1)
  panel.L12ret1year = Vector{MFloat64}(undef, N)
  panel.L24ret1year = Vector{MFloat64}(undef, N)
  #panel.ayret12m = Vector{MFloat64}(undef, N)

  panel.L12vol1yearret = Vector{MFloat64}(undef, N)
  panel.L24vol1yearret= Vector{MFloat64}(undef, N)
  #panel.ayvol1yearnet= Vector{MFloat64}(undef, N)

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
    if (!allunique(spanel.dyrmonth))
      println(spanel[!, [:permno, :gvkey, :adate, :begindate, :enddate, :dyrmonth]])
      error("Dates for dyrmonth are not unique w/in permno")
    end


    for (j, r) ∈ enumerate(eachrow(spanel))

      if j > 12 #compute teh 1yr lagged controls
        target = r.dyrmonth - 12
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L12ret1year = spanel[ptr, :ret1year]
          r.L12vol1yearret = spanel[ptr, :vol1yearret]
        end
      end

      if j > 24 #compute teh 2yr lagged controls
        target = r.dyrmonth - 24
        ptr = searchsortedfirst(spanel.dyrmonth, target)
        if spanel.dyrmonth[ptr] == target
          r.L24ret1year = spanel[ptr, :ret1year]
          r.L24vol1yearret = spanel[ptr, :vol1yearret]
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
