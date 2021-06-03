

PARAM[:capacity] = Dict(
  :futuresfile=>"BBFutures",
  :futuresfileext=>"csv",
  :futurescolumns => [:fullticker, :description, :exchange,
    :volume10d,	:openinterest10d,	:subcategory],
  :periodicity => ["DAILY", "WEEKLY", "MONTHLY", "QUARTERLY", "YEARLY"],
  :pullfields => ["FUT_AGGTE_VOL", "PX_VOLUME", "PX_LAST", "PX_SETTLE",
    "CONTRACT_VALUE", "OPEN_INT"])

function futuresdata(;
  inputpath::String = PARAM[:inputpath],
  futuresfile::String = PARAM[:capacity][:futuresfile],
  futuresfileext::String = PARAM[:capacity][:futuresfileext],
  )

  futures = CSV.read("$inputpath\\$futuresfile.$futuresfileext") |> DataFrame!
  select!(futures, PARAM[:capacity][:futurescolumns])

  parsemfloat(f)=something(tryparse(Float64, f), missing)
  (:test ∈ propertynames(futures)) && (futures.test = parsemfloat.(futures.test))

  #display(describe(futures))
  return futures
end

#a function for pulling lots of tickers and automating common parameter choices
function multipullold(S;
  livepull::Bool="live pull is required",
  tickers::Vector{String}=error("tickers kw are required"),
  fields::Vector{String}=error("fields kw are required"),
  startdate::Date=Date(1900,1,1),
  enddate::Date=today(),
  periodicity::Vector{String}=error("periodicity kw are required"),
  CshAdjAbnormal::Bool =	true,
  CshAdjNormal::Bool =	true,
  currency::String=	"USD",
  DateFormat::String = "D",
  pullname::String=error("pullname is required"),
  inbinstream=IN_BIN_STREAM,
  outbinstream=OUT_BIN_STREAM,
  binext = BIN_EXTENSION)

  local results
  pullpath::String = "$(PARAM[:inputpath])\\$pullname.$binext"

  #WARNING- what to do about dates?

  #set up a parameter dict

  #this will hold each of the result dictionaries
  #allunique(tickers) || error("Tickers are not unique!")

  #batch the tickers for performance
  Ntickers::Int = length(tickers)
  Nperiods::Int = length(periodicity)
  tickerparts=Vector.(Iterators.partition(tickers, Threads.nthreads()) |> collect)
  Nparts::Int = length(tickerparts)

  sessions = (_->startbbsession()).(1:Threads.nthreads())

  if livepull
    try
      resultparts::Vector{Vector{Dict}} = [Vector{Dict}(undef, Nparts) for _ ∈ 1:Nperiods]
      @sync for (i, p) ∈ enumerate(periodicity)
        options::Dict{String, Any} = Dict(
          "CshAdjAbnormal"=>CshAdjAbnormal,
          "CshAdjNormal"=>CshAdjNormal,
          "currency"=>currency,
          #"DateFormat"=>DateFormat,
          )
        if p=="DAILY"
          #options["CDR"] = "7D"
          options["overrideOption"] = "OVERRIDE_OPTION_CLOSE"
          options["nonTradingDayFillOption"] = "ACTIVE_DAYS_ONLY"
          options["nonTradingDayFillMethod"] = "NIL_VALUE"
        end

        #query BB for each partition
        Threads.@threads for j ∈ 1:Nparts
          @async resultparts[i][j] = BLPData.bdh(sessions[Threads.threadid()], tickerparts[j], fields, startdate, enddate,
            periodicity=p, options=options)
        end
      end

    finally
      endbbsession.(sessions)
    end
    #combine the results
    @assert sum([(resultperiod->length.(resultperiod)).(resultparts)...;]) == Ntickers*Nperiods
    results = (resultperiod->reduce(merge, resultperiod)).(resultparts)
    @assert length(results) == Nperiods

    #write out the results
    outbinstream(pullpath, results)
  else
    results = inbinstream(pullpath)
  end

  #println(results)
  #println(results)
  results
end

function multipull(S;
  livepull::Bool="live pull is required",
  tickers::Vector{String}=error("tickers kw are required"),
  fields::Vector{String}=error("fields kw are required"),
  startdate::Date=Date(1900,1,1),
  enddate::Date=today(),
  periodicity::Vector{String}=error("periodicity kw are required"),
  CshAdjAbnormal::Bool =	true,
  CshAdjNormal::Bool =	true,
  currency::String=	"USD",
  DateFormat::String = "D",
  pullname::String=error("pullname is required"),
  inbinstream=IN_BIN_STREAM,
  outbinstream=OUT_BIN_STREAM,
  binext = BIN_EXTENSION)

  local results
  pullpath::String = "$(PARAM[:inputpath])\\$pullname.$binext"

  #WARNING- what to do about dates?

  #set up a parameter dict

  #this will hold each of the result dictionaries
  #allunique(tickers) || error("Tickers are not unique!")

  #batch the tickers for performance
  Ntickers::Int = length(tickers)
  Nperiods::Int = length(periodicity)

  resultparts::Vector{Vector{DataFrame}} = [Vector{DataFrame}(undef, Ntickers) for _ ∈ 1:Nperiods]
  if livepull
    @sync for (i, p) ∈ enumerate(periodicity)
      options::Dict{String, Any} = Dict(
        "CshAdjAbnormal"=>CshAdjAbnormal,
        "CshAdjNormal"=>CshAdjNormal,
        "currency"=>currency,
        #"DateFormat"=>DateFormat,
        )
      if p=="DAILY"
        #options["CDR"] = "7D"
        options["overrideOption"] = "OVERRIDE_OPTION_CLOSE"
        options["nonTradingDayFillOption"] = "ACTIVE_DAYS_ONLY"
        options["nonTradingDayFillMethod"] = "NIL_VALUE"
      end

      #query BB for each partition
      for j ∈ 1:Ntickers
        @async resultparts[i][j] = BLPData.bdh(S, tickers[j], fields, startdate, enddate,
          periodicity=p, options=options) |> DataFrame!
        #=for n ∈ names(resultparts[i][j])
          oldtype = eltype(resultparts[i][j][!,n])
          if !(Missing <: oldtype)
            resultparts[i][j][!,n] = Vector{Union{oldtype,Missing}}=#
      end
    end

    #label each df
    badtickers::Vector{String} = Vector{String}()
    for i ∈ 1:Nperiods
      for j ∈ 1:Ntickers
        if ncol(resultparts[i][j]) > 1#length(setdiff(fields, names(resultparts[i][j]))) == 0
          resultparts[i][j][!,:ticker] .= tickers[j]
        else
          push!(badtickers, "$(tickers[j]):$(periodicity[i])")
          resultparts[i][j] = DataFrame()
        end
      end
    end

    @info "Individual dataframes formed. $(length(badtickers)) bad tickers:\n$(badtickers)"

    #combine the results
    @assert sum([(resultperiod->length(resultperiod)).(resultparts)...;]) == Ntickers*Nperiods
    vcatunion(df1,df2) = vcat(df1,df2, cols=:union)
    results = (resultperiod->reduce(vcatunion, resultperiod)).(resultparts)
    @assert length(results) == Nperiods

    #write out the results
    outbinstream(pullpath, results)
  else
    results = inbinstream(pullpath)
  end

  #println(results)
  #println(results)
  results
end


function writefuturespull(df::AbstractDataFrame, outname::String;
  compress::Bool = size(df,1) > 150_000, #simple heuristic for compressing the CSV
  outputpath::String = PARAM[:outputpath],
  futuresfile::String = PARAM[:capacity][:futuresfile],
  outcsvstream = compress ? OUT_CSV_STREAM : OUT_SMALL_CSV_STREAM,
  csvext = compress ? CSV_EXTENSION : SMALL_CSV_EXTENSION
  )

  p::String = "$outputpath\\$outname.$csvext"
  outcsvstream(p, df)
end

#main entry point
function bbscript(S, ::Val{:futurespull};
  livepull=error("Live pull is required"),
  periodicity=PARAM[:capacity][:periodicity],
  fields=PARAM[:capacity][:pullfields])

  futures::DataFrame = futuresdata()
  tickers::Vector{String}=futures.fullticker

  #=tickers = ["AC1 Comdty", "C 1 Comdty", "DCS1 Comdty", "CZT1 Comdty", "WZ1 Comdty",
  "YW1 Comdty", "CRD1 Comdty", "IBW1 Comdty", "YC1 Comdty", "EP1 Comdty",
  "VB1 Comdty", "VBE1 Comdty", "MIE1 Comdty", "JC1 Comdty", "CBS1 Comdty",
  "JMS1 Comdty", "CXY1 Comdty", "VV1 Comdty", "C71 Comdty", "CT1 Comdty",
  #="QCY1 Comdty", "CTT1 Comdty", "CCL1 Comdty", "SHA1 Comdty", "P71 Comdty",
  "PFC1 Comdty", "PFR1 Comdty", "PAL1 Comdty", "CB1 Comdty", "DCE1 Comdty",
  "SB1 Comdty", "APW1 Comdty", "KO1 Comdty", "CC1 Comdty", "KC1 Comdty",
  "CFD1 Comdty", "DF1 Comdty", "REE1 Comdty", "QC1 Comdty", "QW1 Comdty",
  "SBT1 Comdty", "JCI1 Comdty", "UDR1 Comdty", "C21 Comdty", "KCT1 Comdty",
  "S91 Comdty", "CCT1 Comdty", "POC1 Comdty", "DA1 Comdty", "PMA1 Comdty",
  "JO1 Comdty", "MMR1 Comdty", "CHE1 Comdty", "LE1 Comdty", "V61 Comdty",
  "SF1 Comdty", "AX1 Comdty", "KV1 Comdty", "FEP1 Comdty", "MKP1 Comdty",
  "QSR1 Comdty", "OMP1 Comdty", "CPI1 Comdty", "POI1 Comdty", "CSO1 Comdty",
  "FSP1 Comdty", "BSC1 Comdty", "CSE1 Comdty", "QCR1 Comdty", "P81 Comdty",
  "MDS1 Comdty", "BUT1 Comdty", "JOT1 Comdty", "DRW1 Comdty", "QCE1 Comdty",
  "WWC1 Comdty", "WQS1 Comdty", "SAA1 Comdty", "LC1 Comdty", "LH1 Comdty",
  "FC1 Comdty", "JHT1 Comdty", "LS1 Comdty", "LET1 Comdty", "TFE1 Comdty",
  "OVO1 Comdty", "ZRR1 Comdty", "ZRO1 Comdty", "ANA1 Comdty", "M11 Comdty",
  "RS1 Comdty", "C01 Comdty", "DPR1 Comdty", "Q61 Comdty", "IJ1 Comdty",
  "RR1 Comdty", "O 1 Comdty", "SU1 Comdty", "QRI1 Comdty", "G01 Comdty",
  "FU1 Comdty", "MYDJ1 Comdty", "FY1 Comdty", "AKOJ1 Comdty", "IRL1 Comdty",=#
  #"ZVB1 Comdty", "IY1 Comdty", "IRI1 Comdty", "ZRC1 Comdty", "AE1 Comdty",
  #"SH1 Comdty", "AK1 Comdty", "S 1 Comdty", "SM1 Comdty", "BO1 Comdty",
  #"Q81 Comdty", "BP1 Comdty", "M71 Comdty", "SAR1 Comdty", "SAY1 Comdty",
  #"SAW1 Comdty", "SY1 Comdty", "YK1 Comdty", "FNS1 Comdty", "MSI1 Comdty",
  "VT1 Comdty", "SXY1 Comdty", "LIS1 Comdty", ]=#

  #=results = @btime multipull($S,
    livepull=$livepull,
    tickers=$tickers,
    fields=$fields,
    periodicity=$periodicity,
    pullname="ES1test")=#
  results = multipull(S,
    livepull=livepull,
    tickers=tickers,
    fields=fields,
    periodicity=periodicity,
    pullname="futures")

  results .= (result->leftjoin(result, futures, on=:ticker=>:fullticker)).(results)

  for (i,p) ∈ enumerate(periodicity)
    writefuturespull(results[i], "$(PARAM[:capacity][:futuresfile])_$p")
  end

end
