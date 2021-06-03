

function constructportfolios(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  refreshdataseries::Bool = true,
  refreshportfolios::Bool = true,
  portname::String = PORT_NAME,
  panelname::String = PANEL_NAME,
  datapath::String=DATA_PATH,
  inbinstream::Function = IN_BIN_STREAM,
  outbinstream::Function = OUT_BIN_STREAM,
  binextension::String = BIN_EXTENSION,
  outpath::String = OUT_PATH,
  trialrunonly::Bool = false)::DataFrame

  local panel::DataFrame
  local port::DataFrame

  #load the leverage data set
  local dailyunlesstrial::Symbol = trialrunonly ? :monthly : :daily
  local cqports::CQuartileConstructions = leverageconstructions()
  panelname = "$panelname-$dailyunlesstrial"
  portname = "$portname-$dailyunlesstrial"

  if refreshportfolios
    local panel = formdataseries(
      refreshcomp = refreshcomp,
      refreshtvolmcapcrsp = refreshtvolmcapcrsp,
      refreshtretcrsp = refreshtretcrsp,
      refreshdataseries=refreshdataseries,
      trialrunonly = trialrunonly)

    panel = makepanelmonthly(panel, crsptype=dailyunlesstrial)
    port = formfactors(panel, crsptype=dailyunlesstrial)
    cqports = merge(cqports, filteredleverageconstructions())
    #println(describe(port))
    port = mergeoutside(port)
    #println(describe(port))

    #merge and verify integrity
    println("Size of panel: $(size(panel))")
    panel = join(panel, port, on=:date, kind=:inner)
    (s->finiteormissing!(panel, s)).([:ret; :lret; T1_CONTROL10])
    #winsorizekeydata!(panel, fields2winsorize = [:ret; :lret; :ret12m; :LLret12m; :vol252d; :LLvol252d; T1_CONTROL10])

    checkdfkey(panel, [:date, :lpermno])
    println("Size of panel: $(size(panel))")

    checksum::Float64 = computechecksum(panel) + computechecksum(port)
    @info "Panel constructed with checksum $(checksum)"

    condensewithlast!(panel, condense=false)
    if TEST_OUTPUT
      d = Date(2015,12,31)
      CSV.write("output\\$(panelname)_$d.csv", panel[(cd->cd==d).(panel.date), :])
      CSV.write("output\\$portname.csv", port[:,
        filter(s::Symbol->string(s)[1:2] ≠ "S_", names(port))])
    end

    #WARNING These rows save space, but are not essential
    select!(port, Not(Fns(cqports)))
    retainedpanelcols = unique([:gvkey; :permno; :date; :fyear; :linkeffdate; :linkenddate;
      :adate; :begindate; :enddate; names(port); Fgroups(cqports); Fvals(cqports);
      :ret; :lret; T1_CONTROL10; :ret12m; :LLret12m; :vol252dnet; :Lvol252dnet; :LLvol252dnet;
      :annflag; (s->Symbol(:LD, s)).(CQUART_FVALS); CQUART_FVALS_BRS_D;])
    select!(panel, retainedpanelcols)

    sort!(panel, [:permno, :date])
    CSV.write("output\\$(panelname)_ts_$d.csv", panel[1:20_000, :])    
    outbinstream("$datapath\\$panelname.$binextension", panel)
    #CSV.write("$datapath\\$portname.csv")
  #error("A-OK!")
  else
    panel = inbinstream("$datapath\\$panelname.$binextension")
  end

  return panel
end

#=@inline function finiteormissing!(df::AbstractDataFrame, s::Symbol)::Nothing
  df[!,s] .= ((f::Union{Missing, Real})->((!ismissing(f)) && isfinite(f)) ? f : missing).(df[!,s])

  return nothing
end=#

#use this to create the default portoflio constructions
function leverageconstructions()
  local cqports::CQuartileConstructions

  local cqports = CQuartileConstructions(CQUART_FVALS, CQUART_FCONTROLS,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES)
  local cqportsalt = CQuartileConstructions(CQUART_FVALS, CQUART_FCONTROLS_ALT,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES_ALT)
  cqports = merge(cqports, cqportsalt)

  local cqportsd = CQuartileConstructions(CQUART_FVALS_D, CQUART_FCONTROLS_D,
    Fweightons=CQUART_FWEIGHTONS_D, namevec=CQUART_FNAMES_D)
  cqports = merge(cqports, cqportsd)
  local cqportsaltd = CQuartileConstructions(CQUART_FVALS_D, CQUART_FCONTROLS_ALT_D,
    Fweightons=CQUART_FWEIGHTONS_D, namevec=CQUART_FNAMES_ALT_D)
  cqports = merge(cqports, cqportsaltd)

  local cqportsrmd = CQuartileConstructions(CQUART_FVALS_RMD, CQUART_FCONTROLS_RMD,
    Fweightons=CQUART_FWEIGHTONS_RMD, namevec=CQUART_FNAMES_RMD)
  cqports = merge(cqports, cqportsrmd)
  local cqportsaltrmd = CQuartileConstructions(CQUART_FVALS_RMD, CQUART_FCONTROLS_ALT_RMD,
    Fweightons=CQUART_FWEIGHTONS_RMD, namevec=CQUART_FNAMES_ALT_RMD)
  cqports = merge(cqports, cqportsaltrmd)


  local cqportsbrs = CQuartileConstructions(CQUART_FVALS_BRS_D, CQUART_FCONTROLS_BRS_D,
    Fweightons=CQUART_FWEIGHTONS_BRS_D, namevec=CQUART_FNAMES_BRS_D)
  cqports = merge(cqports, cqportsbrs)

  return cqports
end

function filteredleverageconstructions()
  local cqfiltered::CQuartileConstructions

  cqfiltered = CQuartileConstructions(CQUART_FVALS_D, CQUART_FCONTROLS_D,
    Fweightons=CQUART_FWEIGHTONS_D, namevec=CQUART_FNAMES_FD)

  return cqfiltered
end

#computes a checksum of a dataframe. Used to check consistency.
function computechecksum(df::DataFrame; checksum::Float64 = 0.0)::Float64
  for c ∈ eachcol(df)
    if eltype(c) <: Union{Missing, Real}
      checksum += sum((x->ismissing(x) ? 0.0 : isfinite(x) ? x : 0.0).(c))
    end
  end
  checksum = log(checksum)

  return checksum
end

#this set of functions computes the returns from the panel
#accomplished via dispath
const NandRet = Tuple{Int, MFloat64}

function NandRet(panel::AbstractDataFrame, Fweighton::Symbol; Fret::Symbol=FRET)::NandRet
  spanel::SubDataFrame = view(
    panel, completecases(panel[!, [Fret, Fweighton]]), [Fret, Fweighton])

  N::Int = size(spanel,1)
  (N==0) && return (0,missing)
  return NandRet(N, spanel[!,Fret], spanel[!,Fweighton])
end

function NandRet(panel::AbstractDataFrame, ::Nothing; Fret::Symbol=FRET)::NandRet

  rets::T where T<:SubArray = view(panel[!, Fret], (!ismissing).(panel[!, Fret]))
  N::Int = length(rets)
  (N==0) && (return (0,missing))
  return NandRet(N, rets)
end

NandRet(N::Int, rets::AbstractVector)::NandRet = (N, sum(rets)/N)
NandRet(N::Int, rets::AbstractVector, weights::AbstractVector
  )::NandRet = (N, sum(rets .* weights)/sum(weights))

#gets a flat return file
#assumes the existance of a return file
function formfactors(panel::DataFrame; Fret::Symbol=FRET, crsptype::Symbol)::DataFrame

  local port::DataFrame
  local cqports::CQuartileConstructions = merge(leverageconstructions(), filteredleverageconstructions())

  local fgrps::Vector{Symbol} = Fgroups(cqports)
  local debugsuffix="$crsptype"

  #get the monthly returns

  #this will hold all of the portfolios
  #get the returns
  port = DataFrame(date = unique(panel.date))
  sort!(port, :date)

  Nport::Int = size(port, 1)

  #form the necessary columns
  for pc ∈ cqports
    port[!, pc.name] = missings(MFloat64, Nport)
    port[!, pc.Fn] = missings(MInt, Nport)
  end
  #index the return dates
  portrows::Dict = Dict(r.date=>r for r ∈ eachrow(port))

  #good point to check our answers
  if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\panel$(debugsuffix)_$d-mid.csv", panel[(cd->cd==d).(panel.date), :])
  end

  #form the portfolio returns
  #print("PROFILING: portfolio returns loop - ")
  spanels::GroupedDataFrame = groupby(panel, :date)
  @mpar for i ∈ 1:length(spanels)
    spanel::SubDataFrame = spanels[i]
    r::DataFrameRow = portrows[spanel.date[1]]

    #iterate through our portfolio constructions
    for pc ∈ cqports
      r[pc.name] = 0.0
      r[pc.Fn] = 0
      for grp::Int ∈ 1:pc.Ngroups #add the contribution of each group to the factor return
        sspanel = view(spanel, (g::MInt->(!ismissing(g)) && (g == grp)).(spanel[!,pc.Fgroup]),:)

        local N::Int
        local val::MFloat64
        (N, val) = NandRet(sspanel, pc.Fweighton)
        r[pc.Fn] += N
        r[pc.name] += val * pc.groupweights[grp]
      end
    end
  end


  return port
end

function constructportfolios!(lev::DataFrame)
  #first construct the portfolios
  local cqports::CQuartileConstructions = leverageconstructions()

  N::Int = size(lev,1)

  for pc ∈ cqports
    lev[!, pc.Fgroup] = Vector{MInt}(undef, N)
    lev[!, Symbol(:c, pc.Fgroup)] = Vector{MInt}(undef, N) #this is a temporary field holding the 4 within
  end

  #assign the gorups and weights
  lastcontrol::NSymbol = nothing
  one2four::Vector{Int} = collect(1:4)

  #use this to avoid allocation within the innerloop
  sortperms::Vector{Vector{Int}} = [Vector{Int}(undef, 4) for i ∈ 1:Threads.nthreads()]
  for pc ∈ cqports

    csym::Symbol = Symbol(:c, pc.Fgroup)
    #first do the control sort
    if lastcontrol ≠ pc.Fcontrol
      sort!(lev, pc.Fcontrol)
      lastcontrol = pc.Fcontrol
    end

    #complete cases only to make the logic easier
    local essentialfields::Vector{Symbol} = filter!(
      s->!isnothing(s), [pc.Fval, pc.Fcontrol, pc.Fweighton])
    local slev::SubDataFrame = view(lev,completecases(lev[:,essentialfields]), :)
    sslevs = groupby(slev, :fyear)
    for i ∈ 1:length(sslevs) #for each fiscal year
      local sslev::SubDataFrame = sslevs[i]
      local N::Int = size(sslev,1)

      sslev[1:N, csym] .= (i->i÷4).(0:(N-1)) #create the four-group blocks
      s3levs = groupby(sslev, csym)
      @mpar for j ∈ 1:length(s3levs) #for each block
        local s3lev::SubDataFrame = s3levs[j]
        (size(s3lev,1) ≠ 4) && continue

        tid::Int = Threads.threadid()
        sortperm!(sortperms[tid], s3lev[!, pc.Fval])
        #println("vals: $(s3lev[!, pc.Fval])")
        #println("order: $(sortperms[tid])\n")
        s3lev[sortperms[tid], pc.Fgroup] .= one2four
      end
      #N>10 && error("look at the results")
    end
  end



  if TEST_OUTPUT
    y::Int = 2015
    CSV.write("output\\levtemptest_$y.csv", lev[(cd->cd==y).(lev.fyear), :])
  end
  #println("*****THIS IS IT*****:\n", describe(lev))
end




#replaces certain group entries with missing
function filterportfoliogroup!(lev::DataFrame, Fgrp::Symbol, filterfunction::Function,
  Ffiltered::Symbol = Symbol(:F, Fgrp))

  lev[!,Ffiltered] = Vector{MInt}(undef, size(lev,1))

  for (i,r) ∈ enumerate(eachrow(lev))
    r[Ffiltered] = filterfunction(r, Fgrp)
  end

  return nothing
end

const DVOL_FILTER_FIELD = :LDvol252dnet
#only allows for grp 1 if leverage decreases and group 4 if leverage increases
@inline function filterbyvol(r::DataFrameRow, Fgrp::Symbol, FDvol::Symbol = DVOL_FILTER_FIELD)::MInt
  grp::MInt = r[Fgrp]

  Dvol::MFloat64 = r[FDvol]
  if ismissing(Dvol) || ismissing(grp)
    return missing
  end

  #only care about groups 1 and 4
  if (grp == 2 ||
    grp == 3 ||
    (grp == 1 && Dvol < 0) ||
    (grp == 4 && Dvol > 0))

    return grp
  end

  return missing
end

function doubleincreasegroups!(lev::DataFrame;
  grpfields::Vector{Symbol} = CQUART_FGROUPS_D,
  filteredgrpfields::Vector{Symbol} = CQUART_FGROUPS_FD,
  filterfunction::Function = filterbyvol)::Nothing

  ((Fgrp::Symbol, Ffiltered::Symbol)->
    filterportfoliogroup!(lev, Fgrp, filterfunction, Ffiltered)).(grpfields, filteredgrpfields)

  return nothing
end
