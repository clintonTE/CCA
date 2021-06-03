

function constructportfolios(;
  refreshcomp::Bool = true,
  refreshtretcrsp::Bool = true,
  refreshtvolmcapcrsp::Bool = true,
  refreshdataseries::Bool = true,
  refreshportfolios::Bool = true,
  portname::String = PORT_NAME,
  panelname::String = PANEL_NAME,
  datapath::String=DATA_PATH,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
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
    port = mergeoutside(port)
    println(describe(port))

    #merge and verify integrity
    println("Size of panel: $(size(panel))")
    #panel = join(panel, port, on=:date, kind=:inner)
    checkdfkey(panel, [:date, :lpermno])
    println("Size of panel: $(size(panel))")

    #WARNING These rows save space, but are not essential
    retainedpanelcols = unique([:gvkey; :lpermno; :date; :linkeffdate; :linkenddate;
      :adate; :begindate; :enddate; names(port); Fvals(cqport);
      :ret; :lret; T1_CONTROL10;])
    #select!(panel, retainedpanelcols)

    checksum::Float64 = computechecksum(panel) + computechecksum(port)
    @info "Panel constructed with checksum $(checksum)"

    if TEST_OUTPUT
      d = Date(2015,12,31)
      CSV.write("output\\$(panelname)_$d.csv", panel[(cd->cd==d).(panel.date), :])
      CSV.write("output\\$portname.csv", port[:,
        filter(s::Symbol->string(s)[1:2] ≠ "S_", names(port))])
    end

    outjlsstream("$datapath\\$panelname.jls.lz4") do s
      serialize(s, panel)
    end
    CSV.write("$datapath\\$portname.csv")
  #error("A-OK!")
  else
    injlsstream("$datapath\\$panelname.jls.lz4") do s
      panel = deserialize(s) #NOTE: DEBUG ONLY
    end
  end

  return panel
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

function NandRet(panel::AbstractDataFrame, weighton::Symbol; Fret::Symbol=FRET)::NandRet
  spanel::SubDataFrame = view(
    panel, completecases(spanel[!, [Fret, Fweighton]]), [Fret, Fweighton])

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
  local cqports::CQuartileConstructions = leverageconstructions()
  local fgrps::Vector{Symbol} = Fgroups(cqports)
  local debugsuffix="$crsptype"

  #get the monthly returns
  #print("PROFILING: make panel monthly - ")

  #this will hold all of the portfolios
  #get the returns
  #print("PROFILING: unstack returns - ")
  #port = unstack(view(panel, :, [:permno, :date, Fret]), :date, :permno,
  #  Fret, renamecols=(s->Symbol("S_", s)))
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

  #flatten the quartiles too
  #print("PROFILING: quartile flattening and joining - ")
  #=for pc ∈ cqports
    local portQ::DataFrame = unstack(view(panel, :, [:permno, :date, pc.Fgroup]), :date, :permno,
      pc.Fgroup, renamecols=(s->Symbol("Q_", pc.Fgroup, "_", s)))
    port = join(port, portQ, on=:date)
  end=#

  return port
end

#assigns the group for a portfolio
function assigngroups!(vals::AbstractVector{T},
  groups::AbstractVector{MInt},
  breakpoints::NTuple) where T

  N::Int = length(vals)
  quantiles::Vector{MFloat64} = Vector{MFloat64}(undef, N)

  #compute the quantiles
  qdist = StatsBase.ecdf((val::T->val).(vals))
  @inbounds @simd for i ∈ 1:N
    f::Float64 = vals[i]
    quantiles[i] = qdist(f)
  end

  #compute the groups
  groups .= (q::Float64-> findfirst(b::Float64->
      q≤b, breakpoints)).(quantiles)


  return nothing
end


#assigns the group for a portfolio
function assigngroups!(vals::AbstractVector{T},
  control::AbstractVector{V},
  groups::AbstractVector{MInt},
  breakpoints::NTuple) where {T,V}

  N::Int = length(vals)
  cgroups::Vector{MInt} = deepcopy(groups)

  assigngroups!(control, cgroups, breakpoints)

  for i ∈ 1:length(breakpoints)
    local svals::SubArray{T} = view(vals, cgroups .== i)
    local sgroups::SubArray{MInt} = view(groups, cgroups .== i)
    assigngroups!(svals, sgroups, breakpoints)
  end

  return nothing
end

#use this to create the default portoflio constructions
function leverageconstructions()
  local cqports::CQuartileConstructions
  local cqportsalt::CQuartileConstructions

  local cqports = CQuartileConstructions(CQUART_FVALS, CQUART_FCONTROLS,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES)
  local cqportsalt = CQuartileConstructions(CQUART_FVALS, CQUART_FCONTROLS_ALT,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES_ALT)

  cqports = merge(cqports, cqportsalt)

  return cqports
end


function constructportfolios!(lev::DataFrame)
  #first construct the portfolios
  local cqports::CQuartileConstructions = leverageconstructions()

  N::Int = size(lev,1)

  for pc ∈ cqports
    lev[!, pc.Fgroup] = Vector{MInt}(undef, N)
    #lev[!, pc.Fweight] = Vector{MFloat64}(undef, N)
  end

  #assign the gorups and weights
  slevs = groupby(lev, :fyear)
  @mpar for i ∈ 1:length(slevs)
    local slev::SubDataFrame = slevs[i]
    for pc ∈ cqports #for each portfolio construction

      #don't want to worry about missings here
      local essentialfields::Vector{Symbol} = filter!(
        s->!isnothing(s), [pc.Fval, pc.Fcontrol, pc.Fweighton])
      local sslev::SubDataFrame = view(slev,completecases(slev[:,essentialfields]), :)
      if size(sslev,1)>length(pc.breakpoints) #check that we have enough complete cases

        #assign the groups
        assigngroups!(sslev[!, pc.Fval], sslev[!, pc.Fcontrol], sslev[!, pc.Fgroup], pc.breakpoints)

        #=now assign the weights
        for s3lev ∈ groupby(sslev, pc.Fgroup)
          assignweightswithin!(s3lev, pc.Fweight)
        end=#
      end
    end
  end

end


#=function assignportfolios()
  returncontrolled::CQuartileConstructions = CQuartileConstructions(
    CQUART_FVALS, CQUART_FCONTROLS,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES)
  assetcontrolled::CQuartileConstructions = CQuartileConstructions(
    CQUART_FVALS, CQUART_FCONTROLS_ALT,
    Fweightons=CQUART_FWEIGHTONS, namevec=CQUART_FNAMES_ALT)
end=#
