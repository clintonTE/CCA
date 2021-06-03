
function fixfiscal!(data::NCCSData; fiscaldatestr::String = FISCAL_DATE_FORMAT)
  #convert fiscal to dates
  rename!(data.df, :datefiscal=>:datefiscalstr)
  data.df[:datefiscal] =
    (s->ismissing(s) ? missing : Date(s,fiscaldatestr)).(data.df[:datefiscalstr])

  #correct bad fiscal dates
  data.df[:datefiscal] = ((fisyr,datefiscal)->
    (!ismissing(datefiscal)) && (Year(datefiscal) != fisyr) ?
      Date(fisyr, Dates.value(Month(datefiscal)), 1) : datefiscal
      ).(data.df[:fisyr], data.df[:datefiscal])

end

#this merges in a field or set of fields into data
#assumes the data being merged is logged returns
function mergelreturnfields!(data::NCCSData, indf::DataFrame, fieldstomerge::Vector{Symbol};
  fieldnames::Dict = Dict(f=>f for f ∈ fieldstomerge), #allows for renaming fields
  fieldtypes::Dict = Dict(f=>Union{Missing, eltype(indf[f])} for f ∈ fieldstomerge),
  fixfiscal::Bool = true, sortportfolio::Bool = true,
  lagfields::Bool = false) #allows for different field types

  #some time intensive activities
  fixfiscal && fixfiscal!(data) #only need to do this once
  sortportfolio && sort!(data.df, [:ein, :datefiscal])
  indf = sort(indf, :date)#sort!(indf, :date)

  #build an index for performance
  local indfindex::Dict{Union{Date,Missing},MInt} =
    Dict{Union{Date,Missing},MInt}(indf[i,:date]=>i for i ∈ 1:(size(indf,1)))
  indfindex[missing] = missing #provide a missing value if missing
  local idxbegin::Union{Int,Missing}
  local idxend::Union{Int,Missing}

  N::Int = size(data.df, 1)

  #pre-allocate the space
  for f ∈ fieldstomerge
    data.df[fieldnames[f]] = Vector{fieldtypes[f]}(undef, N)
    lagfields && (data.df[Symbol(:L, fieldnames[f])] = Vector{fieldtypes[f]}(undef, N))
  end

  fixfiscal && (data.df[:monthscovered] = Vector{MInt}(undef, N))

  #do merge at the NP level
  for subdf in groupby(data.df,:ein)
    subN::Int = size(subdf,1)

    for i::Int ∈ 1:subN
      idxend =  haskey(indfindex, subdf[i,:datefiscal]) ? indfindex[subdf[i,:datefiscal]] : missing
      if !ismissing(idxend)
        #assume eom dates (fiscal date given as yyyy-mm)
        #uses 12 month lag if no prior date is available

        if (((i==1) || (ismissing(subdf[i-1,:datefiscal]))) ||
          (subdf[i-1,:fisyr] != subdf[i,:fisyr] - 1))
            idxbegin = (idxend ≥ 12) ? idxend - 11 : missing
        elseif haskey(indfindex, subdf[i-1,:datefiscal])
            idxbegin = indfindex[subdf[i-1,:datefiscal]] + 1
        else
           idxbegin = missing
        end

        fixfiscal && (subdf[i, :monthscovered] = idxend - idxbegin + 1) #capture this, but only once

        #merge the data if it exists
        if !ismissing(idxbegin)
          for f ∈ fieldstomerge #merge the fields and the lag field if desired
            subdf[i,fieldnames[f]] = sum(indf[idxbegin:idxend, f])
            lagfields && (i≥2) && (subdf[i,Symbol(:L, fieldnames[f])] = subdf[i-1,fieldnames[f]])
          end
        end

      end
    end
  end

  return nothing
end

function mergeSP500!(data::NCCSData;
    refreshwrds::Bool = REFRESH_WRDS,
    minReturnCutoff = MIN_RETURN_CUTOFF, maxReturnCutoff = MAX_RETURN_CUTOFF)

  #get the S&P500 data and index it
  local spdf = getWRDSSP500(refreshwrds = REFRESH_WRDS)


  #preallocate and make a lookup table for efficiency
  data.df[:sp500ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:Lsp500ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:t30ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:lsp500t30] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:excess] = Vector{MFloat64}(missing, size(data.df,1))

  mergelreturnfields!(data, spdf, [:lvwretd, :lt30ret];
    fieldnames = Dict(:lvwretd=>:lsp500ret, :lt30ret=>:lt30ret),
    fixfiscal = true, sortportfolio=true, lagfields=true)

  data.df[:sp500ret] .= (exp).(data.df[:lsp500ret]) .- 1
  data.df[:Lsp500ret] .= (exp).(data.df[:Llsp500ret]) .- 1

  data.df[:t30ret] .= (exp).(data.df[:lt30ret]) .- 1
  data.df[:excess] .= data.df[:return] .- data.df[:t30ret]
  data.df[:lsp500t30] .= data.df[:lsp500ret] .- data.df[:lt30ret]


  #NOTE: Move these to another method???
  data.df[:pcontributions] = (data.df[:contributions]) ./ data.df[:lagnetassets]
  data.df[:pexpenses] = (data.df[:totexp]) ./ data.df[:lagnetassets]

  #trim rediculous values
  for v ∈ [data.df[:excess], data.df[:pcontributions], data.df[:pexpenses]]
    v .= (f->((ismissing(f)) || (f>maxReturnCutoff) || (f<minReturnCutoff)) ? missing : f).(v)
  end

  return data
end

#loads the Fama French Factors from a file
function loadFF(ffnamesuffix::String; refreshoutsidedata::Bool = true,
  ffpath::String = FF_PATH, ffname::String = FF_NAME,
  ffdateformat::String = FF_DATE_FORMAT)

  local ffdf::DataFrame
  local fffile::String = "$ffpath\\$(ffname)$(ffnamesuffix)"

  #load from a jls file if possible and desired
  if refreshoutsidedata || (!isfile("$fffile.csv"))
    ffdf= CSV.read("$fffile.csv")
    rename!(ffdf, :date=>:olddate)
    ffdf[:date] = (d->Date("$d",ffdateformat)).(ffdf[:olddate])
    deletecols!(ffdf, :olddate)


    serialize("$fffile.jls", ffdf)

  else
    ffdf = deserialize("$fffile.jls")
  end

  return ffdf
end

function mergeFF!(data::NCCSData, FFmodelsuffix::String; refreshoutsidedata::Bool = true)::Nothing
  N::Int = size(data.df,1)

  ffdf::DataFrame = loadFF(FFmodelsuffix, refreshoutsidedata=refreshoutsidedata)

  #merge in the raw data (fields from 3 and 5 factor models)
  fieldstomerge::Vector{Symbol} = FF_SUFFIX_FIELDS[FFmodelsuffix]
  lfieldstomerge::Vector{Symbol} = (s->Symbol(:l,s,FFmodelsuffix)).(fieldstomerge)

  ffN::Int = size(ffdf, 1)

  for i ∈ 1:length(fieldstomerge)
    ffdf[lfieldstomerge[i]] = Vector{MFloat64}(undef, ffN)
    ffdf[lfieldstomerge[i]] .= (f->ismissing(f) ? missing : log(1+f/100.0)).(ffdf[fieldstomerge[i]])
  end

  mergelreturnfields!(data, ffdf, lfieldstomerge, fixfiscal=false, lagfields=true)

  #readjust the mkt benchmark
  lmktt30sym::Symbol = Symbol("lmkt", FFmodelsuffix,"t30")
  data.df[lmktt30sym] = Vector{MFloat64}(undef, N)
  data.df[lmktt30sym] .= data.df[Symbol("lMkt_RF", FFmodelsuffix)] .+
    data.df[Symbol("lRF", FFmodelsuffix)] .- data.df[:lt30ret] #change rf benchmark

  return nothing
end


function mergeoutsidedatasources!(data::NCCSData; refreshoutsidedata::Bool = false,
  FFmodelsuffixes::Vector{String} = ["3", "5"])

  mergeSP500!(data::NCCSData)
  (FFmodelsuffix->mergeFF!(data, FFmodelsuffix)).(FFmodelsuffixes)
end
