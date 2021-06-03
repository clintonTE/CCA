

#does an initial filter to get down the file size
function filterToEndowments!(data::NCCSData)

  assetCutoff::Float64 = ASSET_CUTOFF
  assetLowerBound::Float64 = ASSET_LOWER_BOUND
  largeDifferenceCutoff::Float64 = LARGE_DIFFERENCE_CUTOFF
  #try to selection the desired endowments
  #remove missing data
  data.df = data.df[completecases(data.df[:,
    [:category, :netassets, :totrev, :totexp, :invrealizeddirect, :invrealizedindirect]]), :]
  data.df = data.df[data.df[:category] .== :B, :]
  data.df = data.df[data.df[:charitytype] .== :pc, :]
  N::Int = size(data.df,1)

  println("N1: $N")
  #filter out organizations which are too small
  data.df[:tokeep] = trues(size(data.df,1))
  sort!(data.df, [:ein, :fisyr])
  by(data.df, :ein) do subdf::SubDataFrame

    subN::Int = size(subdf,1)
    maxAssets::Float64 = maximum(subdf[:netassets])
    if maxAssets < assetCutoff || subN ≠ size(unique(subdf[:,[:ein, :fisyr]]),1)
      subdf[:tokeep] = falses(subN)
    end
  end

  #delete the invalid rows
  data.df = data.df[data.df[:tokeep] .== true, :]
  data.df = data.df[data.df[:netassets] .≥ assetLowerBound, :]
  N = size(data.df,1)
  println("N2: $N")

  data.df[:returnamt] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:return] = Vector{Union{Float64,Missing}}(missings(N))

  data.df[:returncheckamt] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:returncheck] = Vector{Union{Float64,Missing}}(missings(N))

  data.df[:largeDifference] = Vector{Union{Bool,Missing}}(missings(N))


  by(data.df, :ein) do subdf::SubDataFrame
    subN = size(subdf,1)

    @inbounds @simd for i ∈ 2:subN
      if subdf[i-1, :fisyr] == subdf[i, :fisyr] - 1

        #NOTE: returns = unrealized + inv income s.t. unrealized = Δassets - netincome
        returnCandidateAmt::Float64 = (subdf[i, :netassets] .- subdf[i-1, :netassets]
          ) .+ subdf[i, :invrealizeddirect] .- subdf[i, :netincome]
        returnCandidate::Float64 = returnCandidateAmt / subdf[i-1, :netassets]

        returnCheckCandidateAmt::Float64 = (subdf[i, :netassets] .- subdf[i-1, :netassets]
          ) .+ subdf[i, :invrealizedindirect] .- subdf[i, :netincome]
        returnCheckCandidate::Float64 = returnCheckCandidateAmt / subdf[i-1, :netassets]

        if (returnCandidate < 2. && returnCandidate > -0.9 &&
            abs(returnCandidate - returnCheckCandidate) ≤ largeDifferenceCutoff)# && abs(subdf[i, :totinv]) ≥ 1. #make sure we have valid numbers here
          subdf[i, :returnamt] = returnCandidateAmt
          subdf[i, :return] = returnCandidate

          subdf[i, :returncheckamt] = returnCheckCandidateAmt
          subdf[i, :returncheck] = returnCheckCandidate
        end
      end
    end
  end

  sort!(data.df, [:ein, :fisyr])

  return data
end


function makeWideEndowmentReturns(data::NCCSData)::DataFrame
  #form a wideDF
  widedf::DataFrame = DataFrame(ein = unique(data.df[:ein]))
  allYears::Dict = Dict(i=>Symbol(i) for i::Int ∈ collect(1989:2016))

  wideN::Int = size(widedf,1)
  N::Int = size(data.df,1)

  #make an index for the new dataframe
  einIndex::Dict = Dict(widedf[i, :ein]=>i for i::Int ∈ 1:wideN)

  #fill the rows for the wideDf
  for y::Int ∈ keys(allYears)
    widedf[allYears[y]] = Vector{Union{Missing, Float64}}(missings(wideN))
  end

  #this will hold the latest names
  widedf[:name] = Vector{Union{Missing, String}}(missings(wideN))

  sort!(data.df, [:ein, :fisyr])

  by(data.df, :ein) do subdf::SubDataFrame
    wideRow::Int = einIndex[subdf[1,:ein]]
    subN::Int = size(subdf,1)

    for i::Int ∈ 2:subN
      if !ismissing(subdf[i,:return])
        widedf[wideRow, allYears[subdf[i,:fisyr]]] = subdf[i,:return]
      end
    end

    #record the most recent name
    nameCandidates::Vector{String} = subdf[(!ismissing).(subdf[:,:name]), :name]
    if length(nameCandidates) > 0
      widedf[wideRow, :name] = nameCandidates[end]
    end
  end


  println(summary(widedf))

  return widedf
end

function computeEndowmentReturns(; refreshFilter::Bool = true,
    endowmentReturnsName::String = ENDOWMENT_RETURNS_NAME,
    writeEndowmentReturnsCSV::Bool = true,
    dataPath::String = DATA_PATH)::Nothing

  if refreshFilter
    data::NCCSData = deserializeNCCSData()
    data = filterToEndowments!(data)
    serializeNCCSData(data, dataName = endowmentReturnsName)
  else
    data = deserializeNCCSData(dataName = endowmentReturnsName)
  end

  widedf::DataFrame = makeWideEndowmentReturns(data)
  println("widedf size: $(size(widedf))")

  if writeEndowmentReturnsCSV
    outPath::String = "$dataPath\\$(endowmentReturnsName).csv"
    uCSV.write(outPath, widedf)

    outPath = "$dataPath\\$(endowmentReturnsName)_long.csv"
    uCSV.write(outPath, data.df)
  end



  return nothing

end
