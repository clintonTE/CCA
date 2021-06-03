
const TAIL_MIN_ASSETS = 0.
const TAIL_NUM_YEARS = 5

function performancebyein(data::NCCSData; tailminassets = TAIL_MIN_ASSETS,
  numyears = TAIL_NUM_YEARS, lastyear = 2015,
  outputpath::String = TABLE_PATH,
  identifierfields::Vector{Symbol} = [:ein, :name, :fisyr, :datefiscal],
  return1y::Symbol = :lreturn,
  wealthfield::Symbol = :adjnetassets,
  sumfields::Vector{Symbol} = [:monthscovered, :lt30ret, :lsp500ret,
    :bmff3, :bmff5, :bmsp500],
  annualizefields::Bool = true)::DataFrame

  #make some assertions on what must be present in the field lists
  @assert :monthscovered ∈ sumfields
  @assert :fisyr ∈ identifierfields

  local firstyear::Int = lastyear - numyears
  local fields2annualize::Vector{Symbol}

  #partition off the part we want
  validfirmsdf = view(data.df, (data.df[wealthfield] .> tailminassets) .&
    (data.df[:fisyr] .≤ lastyear) .& (data.df[:fisyr] .> firstyear),
    [identifierfields; return1y; wealthfield; sumfields;]
    )

  #form the consolidated df

  taildf::DataFrame = DataFrame(ein=validfirmsdf[:ein])
  N::Int = size(taildf, 1)

  #pre-allocate the columns
  for f ∈ setdiff(names(validfirmsdf), [:ein])
    local newtype::Type
    local origtype::Type = eltype(data.df[f])

    if Float64 <: origtype
      newtype = MFloat64
    elseif Int <: origtype
      newtype = MInt
    elseif Date <: origtype
      newtype = MDate
    elseif String <: origtype
      newtype = MString
    else
      newtype = MSymbol
    end

    println("field $f: oldtype= $(origtype) newtype = $newtype")

      taildf[f] = Vector{Union{Missing, newtype}}(missing, N)
  end
  taildf[:N] = Vector{MInt}(undef, N)

  #create the same set of fields but annualized
  #NOTE: Annualized fields are not logged
  if annualizefields
    for f ∈ [return1y; sumfields]
      f_ann::Symbol = Symbol(f, "_ann")
      taildf[f_ann] = similar(taildf[f])
    end

    fields2annualize = setdiff([return1y; sumfields;], [:monthscovered])
  end

  taildfindex::Dict = Dict(taildf[i,:ein] => taildf[i,:] for i::Int  ∈ 1:N)


  #consolidate the eprformacne data
  for subdf ∈ groupby(validfirmsdf, :ein)
    tailrow = taildfindex[subdf[1,:ein]]

    #identifier fields receive the value from the first row
    #NOTE: No verification that identifier fields are the same for a given org
    for f ∈ identifierfields
      tailrow[f] = subdf[1,f]
    end

    tailrow[return1y] = sum(subdf[return1y])

    #subtotal the sum fields
    for f ∈ sumfields
      tailrow[f] = sum(subdf[f])
    end

    #if we want annualized copies
    if annualizefields
      yearscoveredinv::MFloat64 = 1.0 / (tailrow[:monthscovered] / 12.0 )

      for f ∈ fields2annualize
        tailrow[Symbol(f,"_ann")] = exp(tailrow[f] * yearscoveredinv) - 1.0
      end
    end

    #get average assets
    tailrow[wealthfield] = mean(subdf[wealthfield]) ./ 1_000_000

    tailrow[:N] = size(subdf,1)
  end

  taildf = taildf[completecases(taildf[[identifierfields; return1y; wealthfield]]),:]
  taildf = taildf[taildf[:N] .== numyears,:]

  return taildf
end


#writes out the data to the tails
function summarizesimpletails(data::NCCSData;
  numyears = TAIL_NUM_YEARS, lastyear = 2015,
  return1y::Symbol = :lreturn,
  numberpertail::Int = 200,
  outputpath::String = TABLE_PATH,
  datasuffix::String = DATA_SUFFIX)::Nothing

  firstyear::Int = lastyear - numyears
  outpath = "$outputpath\\tails_$(firstyear)_$(lastyear)_$return1y"

  #get a frame with cumulative returns by year
  taildf::DataFrame = performancebyein(data, numyears=numyears, lastyear=lastyear,
    return1y=return1y)

  #organize the return dataframe
  sort!(taildf, [return1y])
  N::Int = size(taildf, 1)

  taildfbottom::SubDataFrame = view(taildf, 1:numberpertail,:)
  taildftop::SubDataFrame = view(taildf, (N-numberpertail + 1):N,:)

  #finally write out all the tail data
  taildftop |> CSV.write("$(outpath)_top_$datasuffix.csv")
  taildfbottom |> CSV.write("$(outpath)_bottom_$datasuffix.csv")

  return nothing
end

#extract the tails
function extractpersistencetails(df::DataFrame, period::Int, endyear::Int, lbenchmark::NSymbol=nothing;
    outputpath::String = WORKING_PATH, threshold::Float64 = TAIL_THRESHOLD,
    minentriesperfile::Int = TAIL_MIN_ENTRIES_PER_FILE,
    consecutive::Bool = false, datasuffix::String=DATA_SUFFIX)::Nothing

  local subdf::SubDataFrame
  local topdf::SubDataFrame
  local bottomdf::SubDataFrame
  local topdfout::DataFrame
  local bottomdfout::DataFrame

  #get the field symbols given the period and benchmark
  local (percentilesym, Lpercentilesym) = buildpercentilesym(period, lbenchmark)
  local returnsym = Symbol("lreturn$(period)yr")

  #make the windows non-overlapping
  subdf = view(df, df[:fisyr] .≤ endyear, :)
  subdf = view(subdf, ((endyear .- subdf[:fisyr]) .% period) .== 0, :)

  #no missing data
  subdf = view(subdf, (!ismissing).(subdf[percentilesym]), :)

  #split out views for the perisstent tails
  topdf = view(subdf, (subdf[percentilesym] .≥ (1.0 - threshold)), :)
  bottomdf = view(subdf, (subdf[percentilesym] .≤ threshold), :)

  #use a different function for consecutive tails

  topdfout = aggregatetails(topdf, period, returnsym, percentilesym,
    minentriesperfile=minentriesperfile, consecutive=consecutive)
  bottomdfout = aggregatetails(bottomdf, period, returnsym, percentilesym,
    minentriesperfile=minentriesperfile, consecutive=consecutive)

  #sort to enable easy ranking
  topdfout = sort(topdfout,
    (order(:numappearances, rev=true),
    order(:avgreturnannual, rev=true)))

  bottomdfout = sort(bottomdfout,
    (order(:numappearances, rev=true),
    order(:avgreturnannual, rev=false)))

  #build the file outname
  outname::String = "perstails_per$(period)_minrec$(minentriesperfile)_thresh$threshold"
  consecutive && (outname = "$(outname)_consec")

  topdfout |> CSV.write("$outputpath\\$(outname)_top_$datasuffix.csv")
  bottomdfout |> CSV.write("$outputpath\\$(outname)_bottom_$datasuffix.csv")

  return nothing
end

#constructs the output dataframe
function aggregatetails(df::AbstractDataFrame, period::Int,
  returnsym::Symbol, percentilesym::Symbol;
  minentriesperfile::Int = TAIL_MIN_ENTRIES_PER_FILE,
  consecutive::Bool = false)

  local outdf::DataFrame
  local aggregatefunction!::Function
  local validrows::Vector{Bool}
  local appearancescutoff::Int

  #we want info at the non-profit level, so create a way to that that info
  outdf = DataFrame(ein = unique(df[:ein]))
  outN::Int = size(outdf,1)

  #this contains the information we want to look at
  outdf[:name] = Vector{MString}(undef, outN)
  outdf[:categorynamefiltered] = Vector{MSymbol}(undef, outN)
  outdf[:numappearances] = Vector{MInt}(undef, outN)
  outdf[:avgreturn] = Vector{MFloat64}(undef, outN)
  outdf[:avgreturnannual] = Vector{MFloat64}(undef, outN)
  outdf[:avgquantile] = Vector{MFloat64}(undef, outN)
  outdf[:appearances] = Vector{MString}(undef, outN)


  outindex::Dict = Dict(outdf[i,:ein] => outdf[i,:] for i ∈ 1:size(outdf,1))

  #NOTE: the kernel function called to aggreagte for each ein is set here
  if consecutive
    aggregatefunction! = ((subdf::AbstractDataFrame, r::DataFrameRow)->
      aggregateconsecutivetails!(subdf,r,period, returnsym, percentilesym)::Nothing)
  else
    aggregatefunction! = ((subdf::AbstractDataFrame, r::DataFrameRow)->
      aggregatetotaltails!(subdf,r,period, returnsym, percentilesym)::Nothing)
  end


  for ssubdf::SubDataFrame ∈ groupby(df, :ein)
    outrow::DataFrameRow = outindex[ssubdf[1,:ein]]
    aggregatefunction!(ssubdf, outrow)
  end


  validrows = falses(size(outdf, 1))
  appearancescutoff = maximum(outdf[:numappearances])

  while (sum(validrows) < minentriesperfile) && (appearancescutoff ≥ 1)
    appearancescutoff -= 1
    validrows .= ((i::MInt) -> (!ismissing(i)) && (i > appearancescutoff)).(outdf[:numappearances])
  end

   #outdf = DataFrame(view(outdf, validrows, :))

  return outdf[validrows, :]
end


#does the processing of the dataframe in the event we want a consecutive sequence
function aggregateconsecutivetails!(df::AbstractDataFrame,
    outrow::DataFrameRow, period::Int, returnsym::Symbol,
    percentilesym::Symbol)::Nothing

  N::Int = size(df, 1)

  local currentctr::Int = 0
  local maxconsecutive = 0

  local runningquantiletotal::Float64 = 0.0
  local runningreturntotal::Float64 = 0.0
  local runningstring::String = ""

  for i ∈ 1:N
    #first check if its a consecutive appearance
    currentctr += 1
    year::Int = df[i,:fisyr]

    #currentctr==1 is a special case
    if (i == 1) || (df[i-1,:fisyr] == year - period)
      runningreturntotal += df[i, returnsym]
      runningquantiletotal += df[i, percentilesym]
      runningstring = "$(runningstring)$(year), "

      if (maxconsecutive < currentctr)
        #record data on the max consecutive tail run
        maxconsecutive = currentctr
        outrow[:avgreturn] = runningreturntotal
        outrow[:avgquantile] = runningquantiletotal
        outrow[:appearances] = runningstring
      end
    else #otherwise reset the counter and the running totals
      runningreturntotal = df[i, returnsym]
      runningquantiletotal = df[i, percentilesym]
      runningstring = "$year, "
      currentctr = 1
    end
    #println("year: $year currentctr: $currentctr maxconsecutive: $maxconsecutive")
  end

  #=if maxconsecutive* 6 ≠ length(outrow[:appearances])
    println("output: $outrow")
    println(df[[:ein, :fisyr, returnsym, percentilesym]])
    println("period: $period, maxconsecutive: $maxconsecutive")
  end=#
  #@assert maxconsecutive* 6 == length(outrow[:appearances])

  outrow[:avgreturn] *= 1/maxconsecutive
  outrow[:avgreturnannual] = outrow[:avgreturn] / period
  outrow[:avgquantile] *= 1/maxconsecutive
  outrow[:numappearances] = maxconsecutive
  outrow[:name] = df[1,:name]
  outrow[:categorynamefiltered] = df[1,:categorynamefiltered]

  return nothing
end

#in this simpler case, we care only about tail appearances not consecutive appearances
function aggregatetotaltails!(df::AbstractDataFrame,
  outrow::DataFrameRow, period::Int, returnsym::Symbol, percentilesym::Symbol)::Nothing

  outrow[:numappearances] = size(df,1)
  outrow[:avgreturn] = mean(df[returnsym])
  outrow[:avgreturnannual] = outrow[:avgreturn] / period
  outrow[:avgquantile] = mean(df[percentilesym])
  outrow[:appearances] = join((string).(df[:fisyr]), ", ")
  outrow[:name] = df[1,:name]
  outrow[:categorynamefiltered] = df[1,:categorynamefiltered]

  return nothing
end

function extractpersistencetails(data::NCCSData;
  periods::Vector{Int} = PERSISTENCE_PERIODS,
  persistenceyears::Vector{Int} = PERSISTENCE_YEARS
  )::Nothing

  for p ∈ periods, y ∈ persistenceyears
    extractpersistencetails(data.df, p, y, consecutive = false)
    extractpersistencetails(data.df, p, y, consecutive = true)
  end

  return nothing
end

function summarizetails(data::NCCSData)
  #NOTE: Uncomment the below when we want simple tails
  #summarizesimpletails(data)
  extractpersistencetails(data)
end
