


#load the ff data
function loadff(ffname::String;
  ffpath::String = PARAM[:ffpath],
  ffdateformat::DateFormat = PARAM[:wrdsdateformat],
  fieldprefix::String = "F_$(lowercase(ffname))",
  ffscale::Bool = PARAM[:ffscale][ffname])

  ff::DataFrame = CSV.File("$ffpath\\$ffname.csv") |> DataFrame

  #rename the columns
  ffnames::Vector{String} = names(ff)
  for (i,str) ∈ enumerate(ffnames)
    if i==1 #date column is a special case
      ffnames[1] = "date"
      continue
    end

    str = replace(str, "-"=>"")
    str = replace(str, " "=>"")
    str = replace(str, "_"=>"")
    str = lowercase(str)
    ffnames[i] = length(fieldprefix) ≥ 1 ? string("$(fieldprefix)_", str) : str
  end
  rename!(ff, ffnames)
  ff.date = (i::Int->Date("$i", ffdateformat)).(ff.date)

  #check for missing data and scale
  Frets::Vector{String} = setdiff(ffnames, ["date"])
  Threads.@threads for Fret ∈ Frets #iterate the focal columns
    for r ∈ eachrow(ff) #now iterate down the rows
      if (r[Fret] == -99.99) || (r[Fret] < -100.) #missing codes for ff
        r[Fret] = missing
      elseif ffscale
        r[Fret] /= 100.
      end
    end
  end

  return ff
end



#loads the ff returns and converts the returns to an index
function loadffindexes(ffname::String)
  rawff::DataFrame = loadff(ffname)
  sort!(rawff, :date)

  #want the return names which include the file prefix
  retnames::Vector{Symbol} = (s->Symbol("F_$(ffname)_", s)).(PARAM[:ffretfieldidx][ffname])
  Findexes, indexdf::DataFrame = returns2index(rawff,Frets=retnames)

  #writes out the entire df and index
  if PARAM[:testff]
    testdf = leftjoin(rawff, indexdf, on=:date)
    sort!(testdf, :date)
    testdf |> CSV.write("$(PARAM[:testpath])\\$(ffname)_idxs.csv")
  end

  return Findexes, indexdf
end

#loads ff data and merges it into the existing dataframe
function mergeff(df::DataFrame;
  datefrequency::Symbol = throw("date frequency is required (:day, :month, :week)"),
  testlabel::String = "FF",
  Fgroup::NSymbol = nothing, #only used for testing purposes if provided
  validateexactfrequency = datefrequency != :day)

  local ffnames::Vector{String}=PARAM[:ffnames]
  local ffretfieldidx::Dict = PARAM[:ffretfieldidx]

  #load the idnex values, perform a basic validation and merge them into the main df
  for (i,ffname) ∈ enumerate(ffnames)
    Findexes, indexdf::DataFrame = loadffindexes(ffname)

    #make sure we have the right fields
    Frets = (s->Symbol("F_$(ffname)_", s)).(ffretfieldidx[ffname])
    @assert issetequal(Findexes, (Fret->Symbol(Fret,:_idx)).(Frets))
    @assert setdiff(Findexes, propertynames(indexdf)) |> isempty

    df = index2returns!(df, indexdf;
      Frets, Findexes, datefrequency, testlabel, Fgroup, validateexactfrequency)

    if datefrequency === :month
      select!(df, Not([:indexdate]))
    end
  end

  return df
end
