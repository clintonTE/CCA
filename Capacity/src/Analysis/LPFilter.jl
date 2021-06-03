
lpsuffix(per::Day)= "_LP$(per.value)d"
finiteor0(x) = isfinite(x) ? x : 0.0
finiteor0(::Missing) = 0.0
function lpfilter(df;
    focalcols::Vector{Symbol},
    lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
    buf=PARAM[:analysislpfilterbuffer],
    op=:mean
  )

  @assert issorted(df, :date)
  #these are the columns we will filter
  #originalgrowthcols::Vector{Symbol} = growthcols(df)
  @assert allunique(lpfilterperiods)

  #used for the start heuristic
  averageinterval = mean((df.date[2:end] .- df.date[1:(end-1)]) .|> (d)->d.value)
  for per ∈ lpfilterperiods
    (per === nothing) && continue #indicates the unfiltered columns

    #create the new cols
    lpcols = (s->Symbol(s,lpsuffix(per))).(focalcols)
    for lpcol ∈ lpcols
      df[!, lpcol] = missings(Float64, nrow(df))
    end

    startrow = (per.value / averageinterval) |> round |> Int
    for (i,r) ∈ enumerate(eachrow(df))
      #we don't want to start the series too early, but allow for some tolerance via buf
      #note buf needs to calibrated based on 1) the size of each period and 2) noting
      #that the dates are at the end of each period
      (i < startrow) && continue
      window = view(df,
        (df.date .≥ (r.date - per + buf)) .& (df.date .≤ r.date), focalcols)

      if op === :mean
        for (Fold, Fnew) ∈ zip(focalcols, lpcols)
          r[Fnew] = mean(skipmissing(window[!, Fold])) |> finiteormissing
        end

        @assert (finiteor0.(Vector(r[lpcols])) .≈
          finiteor0.(vec(sum(x->coalesce(x,0.0), Matrix(window),dims=1)) ./
          vec(sum(x->!ismissing(x), Matrix(window), dims=1)))) |> all
      elseif op === :sum
        for (Fold, Fnew) ∈ zip(focalcols, lpcols)
          r[Fnew] = sum(skipmissing(window[!, Fold])) |> finiteormissing
        end

        @assert (finiteor0.(Vector(r[lpcols])) .≈
          finiteor0.(vec(sum(x->coalesce(x,0.0), Matrix(window),dims=1)))) |> all
      else
        throw("unrecognized op $op")
      end
    end

    #adjust for growth rows having a missing first row
    for (Fold, Fnew) ∈ zip(focalcols, lpcols)
      if df[1,Fold] === missing #heuristic- if this is true, assume a growth row
        df[startrow, Fnew] = missing
      end
    end

  end


  df |> CSV.write("$(PARAM[:testpath])\\$(PARAM[:analysismeasurefilename])_testlp.csv")
  return df
end

#helper function to match smoothing levels bewteen measures and benchmarks
function matchcolsonlpsuffix(col, allcols)::Vector{Symbol}
  local strfocalcols
  strallcols = string.(allcols)
  if occursin(r"_LP[0-9]+d$", string(col))
    suffix = match(r"LP[0-9]+d$", string(col)).match
    strfocalcols = allcols[(s->occursin(Regex("$(suffix)\$"),s)).(strallcols)]
  else #no smoothing case
    strfocalcols = setdiff(allcols, allcols[(s->occursin(r"LP[0-9]+d$",s)).(strallcols)])
  end

  return Symbol.(strfocalcols)
end
