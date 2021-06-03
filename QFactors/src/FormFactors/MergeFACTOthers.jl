
#loads the value weighted crsp index
#also does some light preprocessing
function loadidxcrsp(idxcrspname::String, retainedidxcrspcols::Vector{Symbol};
  idxcrsppath::String = IDXCRSP_PATH,
  idxcrspdateformat::DateFormat = IDXCRSP_DATE_FORMAT,
  idxdatecol = IDX_DATE_COL)

  idxcrsp::DataFrame = CSV.read("$idxcrsppath\\$idxcrspname.csv") |> DataFrame
  idxcrsp.date = (i::Int->Date("$i", idxcrspdateformat)).(idxcrsp[!,idxdatecol])

  TEST_OUTPUT && (idxcrsp |> CSV.write("$OUT_PATH\\$idxcrspname.csv"))
  select!(idxcrsp, retainedidxcrspcols)

  return idxcrsp
end

#load the ff data
function loadffdailyports(ffportname::String = FF_PORT_NAME,
  ffpath::String = FF_PATH,
  ffportrename::Dict = FF_PORT_RENAME,
  ffportvals::Vector{Symbol} = FF_PORT_VALS,
  ffdateformat::DateFormat = FF_DATE_FORMAT)

  ff::DataFrame = CSV.read("$ffpath\\$ffportname.csv") |> DataFrame
  ff.date = (i::Int->Date("$i", ffdateformat)).(ff.crspdate)

  names!(ff, (s::Symbol->Symbol(lowercase(replace(string(s)," "=>"")))).(names(ff)))
  rename!(ff, ffportrename)
  TEST_OUTPUT && (ff |> CSV.write("$OUT_PATH\\$ffportname.csv"))
  #now convert the returns to indices
  for (i,r) ∈ enumerate(eachrow(ff))
    if i == 1 #start the index
      for s ∈ ffportvals
        r[s] = 1.0
      end
      continue
    end

    #=Threads.@threads =#
    for s ∈ ffportvals
      v = r[s]

      #handle missing values- assume no change in the index if missing
      ((v == -99.99) || (v < -100.)) && (v = 0.0)
      r[s] = (1. + v/100.) * ff[i-1,s]
    end
  end

  #TEST_OUTPUT && (ff |> CSV.write("$OUT_PATH\\$ffportname.csv"))
  select!(ff, [:date; ffportvals;])

  #rename index fields to avoid confusion with actual return fields
  rename!(ff, Dict(s=>Symbol(s,:idx) for s ∈ ffportvals))

  return ff
end

#interpolates the index to every day in an interval using the appropriate compounding
function interpolateindex(df::DataFrame, idxcol::Symbol,
  datecol::Symbol = :date, tol::Float64 = 10^-10)

  local idf::DataFrame

  df = df[completecases(df[:, [idxcol, datecol]]), :]

  begindate::Date = minimum(df[!,datecol])
  enddate::Date = maximum(df[!,datecol])
  alldates::Vector{Date} = collect(begindate:Day(1):enddate)

  #now create the interpolated dataframe
  Nidf::Int = length(alldates)
  idf = DataFrame(rid = collect(1:Nidf))
  idf[!,datecol] = alldates
  idf = join(idf, df, on=datecol, kind=:left)
  sort!(idf, datecol)

  #set the values to be interpolated between
  T::Type = eltype(idf[!,idxcol])
  idf.nextvalid = Vector{T}(undef, Nidf)
  idf.prevvalid = Vector{T}(undef, Nidf)
  idf.nextvaliddate = Vector{MDate}(undef, Nidf)
  idf.prevvaliddate = Vector{MDate}(undef, Nidf)

  #previous values
  for (i,r) ∈ enumerate(eachrow(idf))
    (i==1) && continue

    if  ismissing(idf[i-1, idxcol])
      r.prevvalid = idf[i-1, :prevvalid]
      r.prevvaliddate = idf[i-1, :prevvaliddate]
    else
      r.prevvalid = idf[i-1, idxcol]
      r.prevvaliddate = idf[i-1, datecol]
    end
  end


  #nextvalues
  for i ∈ Nidf:-1:1
    if ismissing(idf[i,idxcol])
      idf[i, :nextvalid] = idf[i+1, :nextvalid]
      idf[i, :nextvaliddate] = idf[i+1, :nextvaliddate]
    else
      idf[i, :nextvalid] = idf[i, idxcol]
      idf[i, :nextvaliddate] = idf[i, datecol]
    end
  end

  idf.period = ((prev,next)->(!ismissing(prev)) && (!ismissing(next)) ?
    (next - prev).value : missing).(idf.prevvaliddate, idf.nextvaliddate)
  idf.ccret = (log).(idf.nextvalid ./ idf.prevvalid) ./ idf.period

  #finally do the interpolation
  for (i,r) ∈ enumerate(eachrow(idf))
    if ismissing(r[idxcol])
      r[idxcol] = exp(r.ccret)*idf[i-1, idxcol]
    elseif i > 1 #do an error check
      Δ::Float64 = abs(exp(r.ccret)*idf[i-1, idxcol] - r[idxcol])
      (Δ > tol) && @warn("Tolerance $tol exceeded with $Δ on row $i. Row: $r")
    end
  end

  idf |> CSV.write("$OUT_PATH\\idf-$idxcol.csv")

  #drop extraneous values
  select!(idf, Not([:nextvalid, :nextvaliddate, :prevvalid,
    :prevvaliddate, :ccret, :period]))

  return idf
end

function indexreturns!(df::DataFrame, idxcol::Symbol, retcol::Symbol;
  logret::Bool = true, sortdf::Bool = true,
  datecol::Symbol = :date, uselogs::Bool = USE_LOGS)

  #pre-allocate and sort if needed
  N = size(df,1)
  sortdf && (sort!(df, datecol))
  df[!,retcol] = Vector{MFloat64}(undef, N)

  for (i,r) ∈ enumerate(eachrow(df))
    (i==1) && continue #special case since we are working with lags

    #compute the return
    r[retcol] = r[idxcol]/df[i-1, idxcol] - 1.
    (uselogs) && (r[retcol] = log(1. + r[retcol]))
  end

  return df
end



function mergeoutside!(fact::DataFrame; ffportvals::Vector{Symbol} = FF_PORT_VALS)
  local vwcrsp::DataFrame
  local rfcrsp::DataFrame
  local ff::DataFrame

  #first, load the index portfolios
  vwcrsp = loadidxcrsp(VWCRSP_NAME, RETAINED_VWCRSP_COLS)
  rfcrsp = loadidxcrsp(RFCRSP_NAME, RETAINED_RFCRSP_COLS)
  ff = loadffdailyports()

  #interpolate monthly indices where appropriate
  rfcrsp = interpolateindex(rfcrsp, :t30ind)

  #join the indices
  N = size(fact,1)
  fact = join(fact, vwcrsp, on=:date, kind=:inner)
  fact = join(fact, rfcrsp, on=:date, kind=:inner)
  fact = join(fact, ff, on=:date, kind=:inner)

  (N ≠ size(fact,1)) && @warn("Index dates do not match factor dates.
    Factor df size reduced from $N to $(size(fact,1))")



  #compute the index returns
  indexreturns!(fact, :vwindd, :mktgross)
  indexreturns!(fact, :t30ind, :rfr)
  fact.mkt = fact.mktgross .- fact.rfr
  for s ∈
    (s->indexreturns!(fact, Symbol(s,:idx), s)).(ffportvals)
  end

  return fact
end
