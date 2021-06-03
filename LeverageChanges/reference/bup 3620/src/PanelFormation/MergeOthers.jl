
#loads the value weighted crsp index
#also does some light preprocessing
function loadidxcrsp(idxcrspname::String, retainedidxcrspcols::Vector{Symbol};
  idxcrsppath::String = IDXCRSP_PATH,
  idxcrspdateformat::DateFormat = IDXCRSP_DATE_FORMAT,
  idxdatecol = IDX_DATE_COL,
  #=yearrange::UnitRange{Int} = YEAR_RANGE[]=#)

  idxcrsp::DataFrame = CSV.read("$idxcrsppath\\$idxcrspname.csv", threaded=CSV_THREADED) |> DataFrame
  idxcrsp.date = (i::Int->Date("$i", idxcrspdateformat)).(idxcrsp[!,idxdatecol])
  #filter!(r::DataFrameRow->year(r.date) ∈ yearrange, idxcrsp) #restrict to desired range

  #TEST_OUTPUT && (idxcrsp |> CSV.write("$OUT_PATH\\$idxcrspname.csv"))

  select!(idxcrsp, retainedidxcrspcols)
  rename(idxcrsp, Dict(s => Symbol(:P_,s) for s ∈ retainedidxcrspcols))


  return idxcrsp
end

#load the ff data
function loadff(ffname::String;
  ffpath::String = FF_PATH,
  ffdateformat::DateFormat = FF_DATE_FORMAT,
  fieldprefix::String = lowercase("$ffname"),
  rescale::Bool = true,
  yearrange::UnitRange{Int} = YEAR_RANGE[])

  ff::DataFrame = CSV.read("$ffpath\\$ffname.csv", threaded=CSV_THREADED) |> DataFrame

  #rename the columns
  ffnames::Vector{Symbol} = names(ff)
  ffnames[1] = :date
  for (i,s) ∈ enumerate(ffnames)
    local str::String = string(s)
    if i==1 #date column is a special case
      ffnames[1] = :date
      continue
    end

    str = replace(str, "-"=>"")
    str = replace(str, " "=>"")
    str = replace(str, "_"=>"")
    str = lowercase(str)
    ffnames[i] = Symbol("P_$(fieldprefix)_", str)
  end
  rename!(ff, ffnames)
  ff.date = (i::Int->Date("$i", ffdateformat)).(ff.date)

  filter!(r::DataFrameRow->year(r.date) ∈ yearrange, ff) #restrict to desired years

  #now make the colulmns for the index
  idxnames::Vector{Symbol} = setdiff(ffnames, [:date])
  idxnames .= (s->Symbol(s,"_idx")).(idxnames)
  for n ∈ idxnames
    ff[!,n] = Vector{MFloat64}(undef, size(ff,1))
  end

  @mpar for i ∈ 1:length(idxnames) #iterate the focal columns
    idx::Symbol = idxnames[i]
    ret::Symbol = Symbol(replace(string(idx), "_idx"=>""))

    base::Float64 = 1.0
    for r ∈ eachrow(ff) #now iterate down the rows
      if (!ismissing(r[ret]))
        if (r[ret] == -99.99) || (r[ret] < -100.) #missing codes for ff
          r[ret] = missing
        else
          rescale && (r[ret] /= 100)
          base *= 1.0 + r[ret]
          r[idx] = base
        end
      end
    end
  end

  TEST_OUTPUT && (ff |> CSV.write("$OUT_PATH\\$ffname.csv"))
  select!(ff, [:date; idxnames;])

  rename!(ff, Dict(n=>Symbol(replace(string(n), "_idx"=>"")) for n ∈ idxnames))

  return ff
end

#interpolates the index to every day in an interval using the appropriate compounding
function interpolateindex(idxcol::AbstractVector{T},
  datecol::AbstractVector{Date}; tol::Float64 = maximum(idxcol) *10^-13)::Dict{MDate,T} where T<:Any

  local idf::DataFrame
  local ccret::Vector{MFloat64}

  #check a couple basic types of errors
  (sum((ismissing).(idxcol))>0) && (error("interpolateindex does not support missing values.
    Consider supplying a view w/out missing vals."))
  @assert length(idxcol) == length(datecol)

  #accumulate dates
  begindate::Date = minimum(datecol)
  enddate::Date = maximum(datecol)
  alldates::Vector{Date} = collect(begindate:Day(1):enddate)
  actuals::Dict{Date,Union{Missing,T}} = Dict(datecol[i]=>idxcol[i] for i ∈ 1:length(datecol))

  #set the values to be interpolated between
  #first allocate space to work
  allN = length(alldates)
  allvals::Vector{Union{Missing,T}} = Vector{Union{Missing,T}}(undef, allN)
  nextvalids = Vector{Union{Missing,T}}(undef, allN)
  prevvalids = Vector{Union{Missing,T}}(undef, allN)
  nextvaliddates = Vector{MDate}(undef, allN)
  prevvaliddates = Vector{MDate}(undef, allN)


  #identify the most recent previous values
  for i ∈ 1:allN
    (i==1) && continue

    dprev::Date = alldates[i-1]
    if haskey(actuals, dprev)
      prevvalids[i] = actuals[dprev]
      prevvaliddates[i] = dprev
    else
      prevvalids[i] = prevvalids[i-1]
      prevvaliddates[i] = prevvaliddates[i-1]
    end
  end


  #identify the nextvalues by counting backwards
  for i ∈ allN:-1:1
    dnext::Date = alldates[i]
    if haskey(actuals, dnext)
      nextvalids[i] = actuals[dnext]
      nextvaliddates[i] = dnext
    else
      nextvalids[i] = nextvalids[i+1]
      nextvaliddates[i] = nextvaliddates[i+1]
    end
  end

  intervalperiods::Vector{MInt} = ((prev,next)->
    (!ismissing(prev)) && (!ismissing(next)) ?
    (next - prev).value : missing).(prevvaliddates, nextvaliddates)
  try
    ccret = (log).(nextvalids ./ prevvalids) ./ intervalperiods
  catch err
    println("nextvalids: $nextvalids")
    println("prevvalids: $prevvalids")
    error("Error: $err")
  end

  #finally do the interpolation
  for i ∈ 1:allN
    d::Date = alldates[i]
    if (!haskey(actuals, d)) #interpolation case
      allvals[i] = exp(ccret[i])*allvals[i-1]
    else #no need to interpolate if we have an actual value....
      allvals[i] = actuals[d]
      if i > 1 #but still check the answer if possible
        Δ::Float64 = abs(exp(ccret[i])*allvals[i-1] - allvals[i])
        (Δ > tol) && @warn("Tolerance $tol exceeded with $Δ on row $i.
          Value(t): $(allvals[i]), date: $d, return: $(ccret[i])")
      end
    end
  end

  interpolated::Dict{Date,T} = Dict(alldates[i] => allvals[i] for i ∈ 1:allN)
  return interpolated
end

#derives the returns from an outside cumulative return index
function returnsfromindex!(df::DataFrame, Fidxcol::Symbol; Fretcol::Symbol=Fidxcol,
  sortdf::Bool = true,
  Fdatecol::Symbol = :date, uselogs::Bool = false)::Nothing

  #pre-allocate and sort if needed
  N = size(df,1)
  sortdf && (sort!(df, Fdatecol))
  Fretcol_temp_ = Symbol(Fretcol,:_temp_)
  df[!,Fretcol_temp_] = Vector{MFloat64}(undef, N)

  local base::MFloat64 = missing
  for (i,r) ∈ enumerate(eachrow(df))
    if ismissing(base)
      base = df[i, Fidxcol]
      continue #otherwise we are just adding a 0 return row
    end

    base = ismissing(df[i-1, Fidxcol]) ? base : df[i-1, Fidxcol]

    #compute the return
    if !ismissing(r[Fidxcol])
      r[Fretcol_temp_] = r[Fidxcol]/base - 1.
      (uselogs) && (r[Fretcol_temp_] = log(1. + r[Fretcol_temp_]))
    end
  end

  #most likely scenario for the below is the idx and ret col have the same name
  if Fretcol ∈ names(df)
    select!(df, Not(Fretcol))
  end

  rename!(df, Fretcol_temp_ => Fretcol)

  return nothing
end

#helper function for multiple values
function returnsfromindex!(df::DataFrame, Fidxcols::Vector{Symbol};
  Fretcols::Vector{Symbol} = deepcopy(Fidxcols),
  sortdf::Bool = true, Fdatecol::Symbol = :date, uselogs::Bool = false)::Nothing

  sortdf && (sort!(df, Fdatecol)) #no need to sort multiple times

  for i ∈ 1:length(Fidxcols)
    returnsfromindex!(df, Fidxcols[i], Fretcol=Fretcols[i],
      sortdf = false,  Fdatecol = Fdatecol, uselogs=uselogs)
  end

  return nothing

end

function mergeoutsidedf!(::Type{T}, port::DataFrame, source::AbstractDataFrame,
  Fidxcol::Symbol; Fdatecol::Symbol = :date, interpolate::Bool = true)::Nothing where T<:Any

  local indexvals::Dict{Date,T}

  #make the index column in the output dataframe
  port[!,Fidxcol] = Vector{Union{Missing, T}}(undef, size(port,1))

  #avoid missing data in the input
  ssource::SubDataFrame = view(source, (!ismissing).(source[!, Fidxcol]),:)

  if interpolate
    indexvals = interpolateindex(ssource[!,Fidxcol], ssource[!,Fdatecol])
    maxdate = maximum(ssource[!,Fdatecol])
    mindate = minimum(ssource[!,Fdatecol])
    for r ∈ eachrow(port) #now update port
      if (r[Fdatecol] ≥ mindate) && (r[Fdatecol] ≤ maxdate)
        r[Fidxcol] = indexvals[r[Fdatecol]]
      end
    end
  else
    indexvals = Dict(r[Fdatecol]=>r[Fidxcol] for r ∈ eachrow(ssource))
    for r ∈ eachrow(port) #now update port but check that each entry exists
      if haskey(indexvals, r[Fdatecol])
        r[Fidxcol] = indexvals[r[Fdatecol]]
      end
    end
  end

  return nothing
end

#helper function that works with an aaray of columns
function mergeoutsidedf!(port::DataFrame, source::AbstractDataFrame,
  indexsyms = setdiff!(names(source), [:date]); interpolate::Bool = true)

  for s ∈ indexsyms
    T::Type = eltype(source[!,s])
    mergeoutsidedf!(T, port, source, s, interpolate=interpolate)
  end

  return nothing
end

#pulls the equal weighted return index values and interpolates for all dates
function mergevwcrspindex!(crsp::DataFrame; Findex::Symbol = VWCRSP_INDEX)
    vwcrsp = loadidxcrsp(VWCRSP_NAME, RETAINED_VWCRSP_COLS)

  idx::Dict{Date, Float64} = interpolateindex(vwcrsp[!, Findex], vwcrsp[!, :date])

  #testdf = DataFrame(date = collect(minimum(vwcrsp.date):Day(1):maximum(vwcrsp.date)))
  #testdf = join(testdf, vwcrsp, on=:date, kind=:left)
  #testdf.testidx = (d->ewidx[d]).(testdf.date) #DELETE THIS ROW
  #CSV.write("output\\vwcrspewindd.csv", testdf)

  crsp[!, Findex] = Vector{MFloat64}(undef, size(crsp,1))

  @assert length(crsp.date) == length(crsp[!, Findex])
  @mpar for i ∈ 1:length(crsp.date)
    crsp[i, Findex] = idx[crsp.date[i]]
  end
  #crsp.ewindd .= (d::Date->ewidx[d]).(crsp.date)
  crsp[!, Symbol(:l, Findex)]= (log).(crsp[!, Findex])

  #error("check now")
  return nothing

end

function mergeoutside(port::DataFrame)
  local vwcrsp::DataFrame
  local rfcrsp::DataFrame
  local ff5::DataFrame
  local ff3m::DataFrame

  #form an intermediate portfolio consisting of the unique panel dates
  outside::DataFrame = DataFrame(date=unique(port.date))

  #load the index portfolios
  vwcrsp = loadidxcrsp(VWCRSP_NAME, RETAINED_VWCRSP_COLS)
  rfcrsp = loadidxcrsp(RFCRSP_NAME, RETAINED_RFCRSP_COLS)
  ff3m = loadff(FF3M_NAME, rescale=false) #use the CRSP version, so no rescale
  ff5 = loadff(FF5_NAME)

  indicesadded::Vector{Symbol} = setdiff([names(vwcrsp); names(rfcrsp);
    names(ff3m); names(ff5)], [:date])

  mergectr::Int = 0
  for df::DataFrame ∈ [vwcrsp, rfcrsp, ff3m, ff5]
    mergectr+=1
    mergeoutsidedf!(outside, df)
  end

  #if TEST_OUTPUT
  #  CSV.write("output\\mergeoutsidetest.csv", port[:,
  #    filter(s::Symbol->string(s)[1:2] ≠ "S_", names(port))])
  #end

  #compute the index returns
  returnsfromindex!(outside, indicesadded, uselogs=false)

  #some housekeeping with the market index
  rename!(outside, Dict(:vwindd=>:mktgross, :t30ind=>:rfr))
  outside.mkt = outside.mktgross .- outside.rfr

  port = join(port, outside, on=:date)

  return port
end
