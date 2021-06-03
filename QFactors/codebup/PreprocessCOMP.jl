
prepcomp(;refreshcomp=true) = loadcompdata(refreshcomp=refreshcomp)

function loadcompdata(;companame=COMP_A_NAME, compqname = COMP_Q_NAME, compname = COMP_NAME,
   datapath=DATA_PATH, refreshcomp=true,
   incsvstream::Function = IN_CSV_STREAM,
   injlsstream::Function = IN_JLS_STREAM,
   outjlsstream::Function = OUT_JLS_STREAM)::DataFrame

  local compa::DataFrame
  local compq::DataFrame
  local comp::DataFrame

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshcomp || (!isfile("$datapath\\$compname.jls.lz4"))
    compa = incsvstream("$datapath\\$companame.csv.lz4") |> CSV.read |> DataFrame
    compq = incsvstream("$datapath\\$compqname.csv.lz4") |> CSV.read |> DataFrame

    comp = preprocesscomp(compa, compq)

    outjlsstream("$datapath\\$compname.jls.lz4") do s
      serialize(s, comp) #NOTE: DEBUG ONLY
    end

    #NOTE: FINAL VERSION
    #preprocesscomp(compa, compq)
    #serialize("$datapath\\$compname.jls", comp)
  else
    injlsstream("$datapath\\$compname.jls.lz4") do s
      comp = deserialize(s) #NOTE: DEBUG ONLY
    end

    # comp = deserialize("$datapath\\$compname.jls") #NOTE: FINAL VERSION
  end



  return comp
end

#helper function to clean up names in the quarterly file
function properqname(s::Symbol)
  str::String = string(s)

  if str[end] == 'q' || (length(str) ≥ 3 && str[(end-2):end] == "qtr")
    return s
  end

  s = Symbol(s,:q)
  return s
end

#helper function to convert a fiscal year-quarter to a date
function fyearqtr2date(yq::YearQuarter, fyr::Int)

  dt::Date = Date(yq.y, fyr, 1)
  offsetmonths::Int = -(4 - yq.q)*3

  dt += Month(offsetmonths)

  return dt
end


#this processes both the q and the a comp frames, with as much commonly used code as possible
function preprocesscomp(compa::DataFrame, compq::DataFrame;
  compdateformat::DateFormat = COMP_DATE_FORMAT,
  usedcompafields::Vector{Symbol} = USED_COMPA_VALUE_FIELDS,
  usedcompqfields::Vector{Symbol} = USED_COMPQ_VALUE_FIELDS)

  #do operations that are common to both the annual and the quarterly data
  compa = preprocesscomp!(compa)
  compq = preprocesscomp!(compq)

  #use a consistent naming convention
  names!(compq, (properqname).(names(compq)))

  #Compute the quarters
  #no missing data for fiscal years
  compa = compa[(!ismissing).(compa.fyear), :]

  compq = compq[((cq, fq)->(!ismissing(cq)) && (!ismissing(fq))).(
    compq.datacqtr, compq.datafqtr),:]
  compq.fyearqtr =  (s->parsecomp(YearQuarter, s)).(compq.datafqtr)
  compq.cyearqtr =  (s->parsecomp(YearQuarter, s)).(compq.datacqtr)

  #convert the end of each fiscal period into the beginning of the appropriate calendar month
  compa.fdate = ((y::Int,m::Int)->Date(y,m,1)).(compa.fyear, compa.fyr)
  compq.fdateq = ((yq::YearQuarter,fyr::Int)-> fyearqtr2date(yq,fyr)).(compq.fyearqtr, compq.fyrq)

  #don't let the fiscal date be too far from the data date
  compa = compa[(r->(firstdayofmonth(r.ddate) - Month(6) ≤ r.fdate)).(eachrow(compa)),:]
  compq = compq[(r->(firstdayofmonth(r.ddateq) - Month(6) ≤ r.fdateq)).(eachrow(compq)),:]

  #delete rows with missing data
  compa = compa[(!ismissing).(compa.at), :]

  #now process cases where firms have duplicate date-firm key pairs
  #println("compa size before dedup:", size(compa))
  #println(describe(compa))
  compa = dedupcomp!(compa, [:gvkey, :fyear], usedcompafields, :indfmt, :ddate)
  compq = dedupcomp!(compq, [:gvkeyq, :fyearqtr], usedcompqfields, :indfmtq, :ddateq)

  #sort the data
  sort!(compa, [:gvkey, :fyear])
  sort!(compq, [:gvkeyq, :fyearqtr])

  ####now make the I/A amount, both in a quarterly and an annual form
  compa.ia = Vector{MFloat64}(undef, size(compa,1))
  compq.iaq = Vector{MFloat64}(undef, size(compq,1))

  lagwithin!(compa, :at, :gvkey, :fyear, sorted=true)
  compa.ia .= (compa.at .- compa.Lat) ./ compa.Lat
  compa.ia .= (f::MFloat64 -> (!ismissing(f)) && isfinite(f) ? f : missing).(compa.ia)

  lagwithin!(compq, :atq, :gvkeyq, :fyearqtr, sorted=true)
  compq.iaq .= (compq.atq .- compq.Latq) ./ compq.Latq
  compq.iaq .= (f::MFloat64 -> (!ismissing(f)) && isfinite(f) ? f : missing).(compq.iaq)

  ###derive ROE
  #first get shareholders equity
  compq.seq = (r::DataFrameRow-> coalesce(r.seqq, r.ceqq+r.pstkq, r.atq - r.ltq)).(eachrow(compq))
  #now add in deferred taxes and inv credits to get book equity, and subtract out preferred stock
  compq.beq = ((r::DataFrameRow)->
    r.seq + coalesce(r.txditcq,0.0) - coalesce(r.pstkrq, r.pstkq,0.0)).(eachrow(compq))
  lagwithin!(compq, :beq, :gvkeyq, :fyearqtr, sorted=true)
  compq.roeq = compq.ibq ./ compq.Lbeq
  compq.roeq .= (f::MFloat64->(!ismissing(f)) && isfinite(f) ? f : missing).(compq.roeq)

  comp::DataFrame = join(compa, compq, on=[:gvkey=>:gvkeyq, :fyear=>:fyearq], kind=:inner)

  #println(describe(comp))

  #CSV.write("output\\comp.csv", comp[comp.dyear .> 2010, :])


  return comp
end

#operations that are common to both the annual and the quarterly data
function preprocesscomp!(comp::DataFrame)
  #make all column names lower case
  names!(comp, (s::Symbol->Symbol(lowercase(string(s)))).(names(comp)))

  #format dates
  comp.ddate = (s->parsecomp(Date, s)).(comp.datadate)
  comp.dyear = (d::MDate-> ismissing(d) ? missing : Dates.year(d)).(comp.ddate)

  comp = comp[(!ismissing).(comp.dyear), :] #require the calendar date


  return comp
end

function dedupcomp!(comp::DataFrame,
  groupbyfields::Vector{Symbol}, usedcompfields::Vector{Symbol},
  formatfield::Symbol, ddatefield::Symbol)

  #keepflags
  comp.keep = trues(size(comp,1))

  #use this for deduping
  comp.valuesmissings = Vector{MFloat64}(undef, size(comp,1))

  ncols::Int = size(comp,2)
  for subdf ∈ groupby(comp, groupbyfields)
    if size(subdf,1) > 1

      #keep the most recent version
      maxdate::Date = maximum(subdf[!,ddatefield])
      subdf.keep .= (subdf[!,ddatefield] .== maxdate)

      #special case where we have two records of two different formats (indl and fs)
      #in this case, keep the format with the least missing values
      if sum(subdf.keep) == 2
        ssubdf::SubDataFrame = view(subdf, subdf.keep, :)
        if (length(unique(ssubdf.indfmt)) == 2)

          ssubdf.valuesmissings[1] = sum((c->ismissing(ssubdf[1,c])).(usedcompfields))
          ssubdf.valuesmissings[2] = sum((c->ismissing(ssubdf[2,c])).(usedcompfields))

          ssubdf.keep .= (ssubdf.valuesmissings .== minimum(ssubdf.valuesmissings))

          #tie breaker is to keep the INDL record (industrial since we are looking at non-financials)
          if sum(ssubdf.keep) > 1
            ssubdf.keep .= (ssubdf[!, formatfield] .== :INDL)
          end
        end
      end

      #if we still have a dup, delete all the records
      if sum(subdf.keep) > 1
        @warn "dedup failed group fields:
          $groupbyfields with values $((f->subdf[1,f]).(groupbyfields)) Records removed."
        subdf.keep .= false
      end
    end
  end

  comp = comp[comp.keep, :] #bandaid until filter is fixed

  select!(comp, Not([:keep, :valuesmissings]))

  return comp
end

parsecomp(::Type{Date}, s::String;
    compdateformat::DateFormat = COMP_DATE_FORMAT) = Dates.Date(s, compdateformat)

#helper methods and special cases
parsecomp(::Type{<:Any}, v::Missing) = missing
parsecomp(::Type{Date}, i::Int) = parsecomp(Date, "$i")

#parse a value where the type is known: Generic case
function parsecomp(::Type{T}, s::String;
    ##NOTE: uncomment when we have the missings parsedmissings = COMP_PARSED_MISSINGS
    ) where T

  #if it doesn't parse to the right type, set to missing
  v::Union{T,Missing} = something(tryparse(T, s), missing)

  #check against the list of missing codes
  #NOTE: uncomment when we have the missings
  #(!ismissing(v)) && (v ∈ parsedmissings) && (v=missing)

  return v
end

function parsecomp(::Type{YearQuarter}, s::String)
  local yq::MYearQuarter

  if length(s) ≠ 6
    yq = missing
  else
    y::MInt = parsecomp(Int, "$(s[1:4])")
    q::MInt = parsecomp(Int, "$(s[6])")
    (q==1) || (q==2) || (q==3) || (q==4) || (q=missing)

    yq = (!ismissing(y)) && (!ismissing(q)) ? YearQuarter(y,q) : missing
  end


  return yq
end

#=creates a type of quarter, which is sort of a hacked integer
function parsecomp(::Type{QQuarter}, s::String)
  local q::Union{Int,Missing}

  if length(s) ≠ 6
    q = missing
  else
    i::MInt = parsecomp(Int, "$(s[6])")
    q = ((!ismissing(i)) && (i ∈ 1:4)) ? QQuarter(i) : missing
  end

  ismissing(q) && @warn "date missing: $s"

  return q
end

#creates a type of QYear, which is sort of a hacked integer
function parsecomp(::Type{QYear}, s::String)
  local q::Union{Int,Missing}

  if length(s) ≠ 6
    q = missing
  else
    i::MInt = parsecomp(Int, "$(s[1:4])")
    q = (!ismissing(i)) ? Int(i) : missing
  end

  ismissing(q) && @warn "date missing: $s"

  return q
end=#
