
const TAIL_MIN_ASSETS = 0.
const TAIL_NUM_YEARS = 5

function performancebyein(data::NCCSData; tailminassets = TAIL_MIN_ASSETS,
  numyears = TAIL_NUM_YEARS, lastyear = 2015,
  outputpath::String = TABLE_PATH,
  identifierfields::Vector{Symbol} = [:ein, :name, :fisyr, :datefiscal],
  return1y::Symbol = :lreturnnetbmsp5001yr,
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
function summarizetails(data::NCCSData;
  numyears = TAIL_NUM_YEARS, lastyear = 2015,
  return1y::Symbol = :lreturnnetbmsp5001yr,
  numberpertail::Int = 200,
  outputpath::String = TABLE_PATH)::Nothing

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
  taildftop |> CSV.write("$(outpath)_top.csv")
  taildfbottom |> CSV.write("$(outpath)_bottom.csv")

  return nothing
end
