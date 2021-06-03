

function futuresdata(;
  inputpath::String = PARAM[:inputpath],
  futuresfile::String = PARAM[:capacity][:futuresfile],
  futuresfileext::String = PARAM[:capacity][:futuresfileext],
  )

  futures = CSV.read("$inputpath\\$futuresfile.$futuresfileext") |> DataFrame
  select!(futures, PARAM[:capacity][:futurescolumns])

  parsemfloat(f)=something(tryparse(Float64, f), missing)
  futures.test = parsemfloat.(futures.test)

  #display(describe(futures))
  return futures
end

#a function for pulling lots of tickers and automating common parameter choices
function multipull(S;livepull::Bool,
  tickers::Vector{String}=error("tickers kw are required"),
  fields::Vector{String}=error("fields kw are required"),
  startdate::Date=Date(1900,1,1),
  enddate::Date=today(),
  periodicity::String=error("periodicity kw are required"),
  adjustmentAbnormal::Bool = true,
  adjustmentNormal::Bool=	true,
  currency::String=	"USD",
  DateFormat::String = "D",) #is this right?

  #WARNING- what to do about dates?

  #set up a parameter dict
  multipullparams::Dict{Symbol, Any} = Dict(:periodicity=>periodicity,
    :adjustmentAbnormal=>adjustmentAbnormal,
    :adjustmentNormal=>adjustmentNormal,
    :currency=>currency,
    #:DateFormat=>DateFormat,
    )
  if Period=="DAILY"
    multipullparams[:CDR] = "7D" #are tehse right? where is teh correct list of parameters?
    multipullparams[:Quote] = "Close"
  end

  if livepull
    results = BLPData.bdh(S, tickers, fields, startdate, enddate;
      multipullparams...) |> DataFrame
    display(results)
  else
    #serialization- compressed?
  end

  results
end

#main entry point
function bbscript(S, ::Val{:futurespull}; livepull=error("Live pull is required"))
  futures::DataFrame = futuresdata()
  tickers::Vector{String}=futures.fullticker

  pull = multipull(S,
    livepull=livepull,
    tickers=["ES1 Index"],
    fields=["px_last"],
    periodicity="Yearly")

end
