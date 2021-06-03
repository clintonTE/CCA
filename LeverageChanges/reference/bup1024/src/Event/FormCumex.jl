
#makes a dataframe with one row per event
#useful for tables 9 and 10 and figure 3
function formcumex(dist::DataFrame;
  refreshcumex::Bool= true,
  cumexpath::String = EVENT_PATH,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  cumexname::String = CUMEX_NAME)

  local Fret::Symbol = DIST_RET
  local cumex::DataFrame

  cumexname = "$cumexname-$(REPLICATION_TYPE[])-$(DIST_TYPE[])"

  #get the appropriate table
  if refreshcumex
    cumex = formcumex(dist,Fret) #generates the table
    outjlsstream("$cumexpath\\$cumexname.jls.lz4") do s
      serialize(s, cumex)
    end
  cumex[cumex.exyear.==2015,:] |> CSV.write(
    "output\\cumex-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-2015.csv")
  println(describe(cumex))
  else #acquire the cached version
    injlsstream("$cumexpath\\$cumexname.jls.lz4") do s
      cumex = deserialize(s)
    end
  end


  return cumex
end

#make the actual table
function formcumex(dist::DataFrame, Fret::Symbol)::DataFrame
  local cumex::DataFrame
  local comp::DataFrame
  local cum::DataFrame
  local ex::DataFrame
  local cumexcols::Vector{Symbol}

  dist.abs = (abs).(dist[!,Fret])
  sdist::SubDataFrame = view(dist, dist.match .* dist.oneex .* dist.dclrexok,:)
  #sdist::SubDataFrame = view(dist, dist.match .* dist.dclrexok,:)

  #build the dataframe of yields

  #use this to aggregate over the days where appropriate
  @inline function aggfuncret(v::AbstractVector{<:MFloat64})
    (sum((ismissing).(v)) == length(v)) && return missing

    return prod(skipmissing(1.0 .+ v)) - 1.0
  end

  @inline function aggfuncabs(v::AbstractVector{<:MFloat64})
    (sum((ismissing).(v)) == length(v)) && return missing

    return sum(skipmissing(v))
  end

  #get the cumex returns over the respective windows
  #cumexcols = [:eid, :abs, Fret]
  #cum = aggregateoverdays(aggfunc, sdist, -7:-3, selectcols=cumexcols,
  #  colnames=[:eid, :cumabs, :cumret])
  #ex = aggregateoverdays(aggfunc, sdist, 3:7, selectcols=cumexcols,
  #  colnames=[:eid, :exabs, :exret])
  #cumex = join(cum, ex, on=:eid)

  cumexcols = [:eid, :abs]
  cum = aggregateoverdays(aggfuncabs, sdist, -7:-3, selectcols=cumexcols,
    colnames=[:eid, :cumabs])
  ex = aggregateoverdays(aggfuncabs, sdist, 3:7, selectcols=cumexcols,
    colnames=[:eid, :exabs])
  cumexabs = join(cum, ex, on=:eid)

  cumexcols = [:eid, Fret]
  cum = aggregateoverdays(aggfuncret, sdist, -7:-3, selectcols=cumexcols,
    colnames=[:eid, :cumret])
  ex = aggregateoverdays(aggfuncret, sdist, 3:7, selectcols=cumexcols,
    colnames=[:eid, :exret])
  cumexret = join(cum, ex, on=:eid)
  cumex = join(cumexabs, cumexret, on=:eid)

  #merge in the event-specific fields
  aggnames::Vector{Symbol} = [:permno, :eid, :yield, :exdate, :exyear,
    :primary, :distcd, :days2lastdiv, :ordinary]
  ssdist::SubDataFrame = view(sdist, (!ismissing).(sdist.yield), aggnames)
  agg::DataFrame = aggregate(groupby(ssdist, :eid), v->v[1])

  #fix the names
  @inline nofunction(s::Symbol) = Symbol(replace(string(s), "_function"=>""))
  rename!(agg, (s::Symbol-> s=>nofunction(s)).(names(agg)))
  #println(describe(DataFrame(agg)))
  cumex = join(cumex, agg, on=:eid)


  ssdist = view(sdist, ((r::DataFrameRow)->
    (r.day==0) && (!ismissing(r[Fret]))).(eachrow(sdist)), [:eid, Fret])
  cumex = join(cumex, ssdist, on=:eid)
  rename!(cumex, Fret=>:exdayret)


  #make the remaining time series that we need
  cumex.Dret = cumex.exret .- cumex.cumret
  cumex.Dabs = cumex.exabs .- cumex.cumabs
  cumex.yieldd1d = cumex.yield ./ (1.0 .- cumex.yield)
  finiteormissing!(cumex, :yieldd1d)

  @assert sum(cumex.yieldd1d .> .25 + 10^-8) == 0 #should be true due to prior winsorization

  #use this to create the unit-standardized fields
  function unitize(v::AbstractVector{T}) where T<:Union{Real,Missing}
    local low::Float64
    local high::Float64

    if Missing<:T
      low = minimum(skipmissing(v))
      high = maximum(skipmissing(v))
    else
      low = minimum(v)
      high = maximum(v)
    end
    return (v .- low) ./ (high - low)
  end


  cumex.UDret = unitize(cumex.Dret)
  cumex.UDabs = unitize(cumex.Dabs)
  cumex.Uyield = unitize(cumex.yield)
  cumex.Uyieldd1d = unitize(cumex.yieldd1d)
  cumex.Uexdayret = unitize(cumex.exdayret)

  requiredfields::Vector{Symbol} = [:yield, :Dret, :Dabs]
  cumex = cumex[completecases(view(cumex,:,requiredfields)), :]

  ##now we need the compustat data
  comp = prepcomp(refreshcomp=false)
  cumex = mergecrspcomp!(cumex, comp, joinkind=:left, Fcrspdate = :exdate, allowdatecollisions=true)
  excludedcompfields::Vector{Symbol} = setdiff(names(comp), [:gvkey, :lpermno, :bkliab, :at])
  select!(cumex, Not(excludedcompfields))
  #println(describe(cumex))

  #***now do the various table 10 classifications
  #start with housekeeping
  N::Int = size(cumex,1)
  cumex.Gyield = Vector{MSymbol}(undef, N)
  cumex.Gtaxdiv = Vector{MSymbol}(undef, N)
  cumex.Gfreqdiv = Vector{MSymbol}(undef, N)
  cumex.Gspecialdiv = Vector{MSymbol}(undef, N)
  cumex.Glastdiv = Vector{MSymbol}(undef, N)
  cumex.Gcompustat = Vector{MSymbol}(undef, N)
  cumex.Gbkliab = Vector{MSymbol}(undef, N)
  cumex.Gassets = Vector{MSymbol}(undef, N)

  #The following code handles dividend classifications

  #start with calssifciations based on distributions
  distcds::Vector{Vector{Int}} = (cd::Int->digits(cd)[4:-1:1]).(cumex.distcd) #(Reverse the order)
  #println(distcds)
  cumex.distcd1 = Vector{MInt}(missing, N)
  cumex.distcd2 = Vector{MInt}(missing, N)
  cumex.distcd3 = Vector{MInt}(missing, N)
  cumex.distcd4 = Vector{MInt}(missing, N)
  @mpar for i::Int ∈ 1:N
    (length(distcds[i])≠4) && continue
    cumex[i, :distcd1] = distcds[i][1]
    cumex[i, :distcd2] = distcds[i][2]
    cumex[i, :distcd3] = distcds[i][3]
    cumex[i, :distcd4] = distcds[i][4]
  end

  ###Do the classifications by row
  @mpar for r::DataFrameRow ∈ eachrow(cumex)
    r.Gyield = Gyield(r)
    r.Gtaxdiv = Gtaxdiv(r)
    r.Gfreqdiv = Gfreqdiv(r)
    r.Gspecialdiv = Gspecialdiv(r)
    r.Glastdiv = Glastdiv(r)
    r.Gcompustat = Gcompustat(r)
    r.Gbkliab = Gbkliab(r)
    r.Gassets = Gassets(r)
  end

  #drop some extraneous fields
  select!(cumex, Not([:distcd1, :distcd2, :distcd3, :distcd4]))

  return cumex
end


#checks asset thresholds
@inline function Gassets(r::DataFrameRow)::MSymbol

  (ismissing(r.at)) && (return missing)
  (r.at < 100.0) && (return :assets100)
  #=(100.0 ≤ r.at) && =#(r.at ≤ 3000.0) && (return :assets1003000)
  ( 3000.0 < r.at) && (return :assets3000)

  @assert false #should be MECE
end

#checks bkliab thresholds: bkliab13, bkliab1323, bkliab23
@inline function Gbkliab(r::DataFrameRow)::MSymbol

  (ismissing(r.bkliab)) && (return missing)
  (r.bkliab < 0.33333333) && (return :bkliab13)
  #=(0.33333333 ≤ r.bkliab) && =#(r.bkliab ≤ 0.66666667) && (return :bkliab1323)
  ( 0.66666667 < r.bkliab) && (return :bkliab23)

  @assert false #should be MECE
end

#checks gvkey for as evidence of the compustat merge: noncompustat, :compustat
@inline function Gcompustat(r::DataFrameRow)::MSymbol

  (ismissing(r.gvkey)) && (return :noncompustat)
  (!ismissing(r.gvkey)) && (return :compustat)

  @assert false #should be MECE
end

#classifies the lastdividend date as :lastdiv40, :lastdiv4080, :lastdiv80
@inline function Glastdiv(r::DataFrameRow)::MSymbol

  (ismissing(r.days2lastdiv)) && (return missing)
  (r.days2lastdiv < 40) && (return :lastdiv40)
  #=(40 ≤ r.days2lastdiv) && =#(r.days2lastdiv ≤ 80) && (return :lastdiv4080)
  (80 < r.days2lastdiv) && (return :lastdiv80)

  @assert false #should be MECE
end

#Gets the dividend frequency: :nonspecialdiv, :specialdiv
const SPECIAL_DIV = Dict{Int, MSymbol}(
  0=>missing,
  1=>missing,
  2=>:nonspecialdiv,
  3=>:nonspecialdiv,
  4=>:nonspecialdiv,
  5=>:nonspecialdiv, #for classification purposes
  6=>:nonspecialdiv, #assume yearend is same as annual
  7=>:specialdiv,
  8=>:specialdiv,
  9=>:specialdiv)
@inline function Gspecialdiv(r::DataFrameRow)::MSymbol

  (ismissing(r.distcd3) || ismissing(r.distcd1) || (r.distcd1 ≠ 1)) && (return missing)
  return SPECIAL_DIV[r.distcd3]


  @assert false #should be MECE
end

#Gets the dividend frequency: :monthlydiv, :quarterlydiv, :semiandannualdiv
const FREQ_DIV = Dict{Int, MSymbol}(
  0=>missing,
  1=>missing,
  2=>:monthlydiv,
  3=>:quarterlydiv,
  4=>:semiandannualdiv,
  5=>:semiandannualdiv, #for classification purposes
  6=>:semiandannualdiv, #assume yearend is same as annual
  7=>missing,
  8=>missing,
  9=>missing)
@inline function Gfreqdiv(r::DataFrameRow)::MSymbol

  (ismissing(r.distcd3) || ismissing(r.distcd1) || (r.distcd1 ≠ 1)) && (return missing)
  return FREQ_DIV[r.distcd3]


  @assert false #should be MECE
end

#Gets the dividend tax status: :taxdiv, :nontaxdiv
const TAX_STATUS = Dict{Int, MSymbol}(
  0=>missing,
  1=>missing,
  2=>:taxdiv,
  3=>:nontaxdiv,
  4=>:nontaxdiv, #per Ivo pg 27, return of capital dividends are considered non-taxable
  5=>:taxdiv,
  6=>:taxdiv,
  7=>:taxdiv,
  8=>:taxdiv,
  9=>:nontaxdiv)
@inline function Gtaxdiv(r::DataFrameRow)::MSymbol

  (ismissing(r.distcd4) || ismissing(r.distcd1) || (r.distcd1 ≠ 1)) && (return missing)
  return TAX_STATUS[r.distcd4]


  @assert false #should be MECE
end


#Yield grouping: :lowyield, :medyield, :highyield
@inline function Gyield(r::DataFrameRow)::MSymbol

  (ismissing(r.yield) || (!r.primary)) && (return missing)
  (r.yield < 0.0075) && (return :lowyield)
  #=(0.0075≤r.yield) && =#(r.yield≤0.015) && (return :medyield)
  (0.015 < r.yield) && (return :highyield)

  @assert false #should be MECE
end
