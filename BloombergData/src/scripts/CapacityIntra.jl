#`=@BDH("GCNY10Yr Index", "TRADE, BID, ASK", "1/1/2000 00:00:00", "7/7/2020 05:30:00",
#"IntrRw", "True", "Size", "S", "Type", "S","TimeZone=New York",)

PARAM[:geoff] = Dict{Symbol, Any}(
  :tickers => ["GCNY30YR Index", "GCNY10YR Index", "GCNY5YR Index", "GCNY2YR Index", "GCNY1YR Index"],
  :pullfields => ["TRADE", "BID", "ASK"],
  :file => "geoffresults")

function geoffpull(S;
  livepull::Bool="live pull is required",
  tickers::Vector{String}=error("tickers kw are required"),
  fields::Vector{String}=error("fields kw are required"),
  startdatetime = DateTime(2000,1,1,0,0,0),
  enddatetime = DateTime(2020,7,8,0,0,0), #WARNING- time zone not adjusted
  pullname::String=error("pullname is required"),
  outbinstream=OUT_BIN_STREAM,
  binext = BIN_EXTENSION)

  local results
  pullpath::String = "$(PARAM[:inputpath])\\$pullname.$binext"

  #set up a parameter dict

  #this will hold each of the result dictionaries
  #allunique(tickers) || error("Tickers are not unique!")

  #batch the tickers for performance
  Ntickers::Int = length(tickers)

  resultparts::Vector{DataFrame} = Vector{DataFrame}(undef, Ntickers)
  if livepull
    options::Dict{String, Any} = Dict(
      #"Size"=>"S",
      #"Type"=>"S",
      #"TimeZone" => "New York",
      #"DateFormat"=>DateFormat,
      )

    #query BB for each partition
    @sync for j ∈ 1:Ntickers
      @async resultparts[j] = BLPData.bdh_intraday_ticks(S, tickers[j], fields,
        startdatetime, enddatetime,
        options=options) |> DataFrame!
    end

    #label each df
    badtickers::Vector{String} = Vector{String}()
    for j ∈ 1:Ntickers
      if ncol(resultparts[j]) > 1#length(setdiff(fields, names(resultparts[i][j]))) == 0
        resultparts[j][!,:ticker] .= tickers[j]
      else
        push!(badtickers, "$(tickers[j])")
        resultparts[j] = DataFrame()
      end
    end

    @info "Individual dataframes formed. $(length(badtickers)) bad tickers:\n$(badtickers)"

    #combine the results
    @assert length(resultparts) == Ntickers
    vcatunion(df1,df2) = vcat(df1,df2, cols=:union)
    results = reduce(vcatunion, resultparts)

    #write out the results
    outbinstream(pullpath, results)
  else
    results = inbinstream(pullpath)
  end

  #println(results)
  #println(results)
  results
end


function writepull(df::AbstractDataFrame, outname::String;
  compress::Bool = false,#size(df,1) > 150_000, #simple heuristic for compressing the CSV
  outputpath::String = PARAM[:outputpath],
  futuresfile::String = PARAM[:geoff][:file],
  outcsvstream = compress ? OUT_CSV_STREAM : OUT_SMALL_CSV_STREAM,
  csvext = compress ? CSV_EXTENSION : SMALL_CSV_EXTENSION
  )

  p::String = "$outputpath\\$outname.$csvext"
  outcsvstream(p, df)
end

#main entry point
function bbscript(S, ::Val{:geoffpull};
  livepull=error("Live pull is required"))

  fields=PARAM[:geoff][:pullfields]
  tickers=PARAM[:geoff][:tickers]


  results = geoffpull(S,
    livepull=livepull,
    tickers=tickers,
    fields=fields,
    pullname="geoff")

  writepull(results, "$(PARAM[:geoff][:file])")

end
