
const DECIMALS = 2
const DECIMALS_RETURN = 2
const ASSETS_SUMMARY_PERIOD = 1994:2015
const RETURNS_SUMMARY_PERIOD = 2000:2015
const BY_CHARITY_FIELD = :charitytype
const BY_CATEGORY_FIELD = :categoryname



function assetsByYearDF(data::NCCSData; summaryPeriod::UnitRange = ASSETS_SUMMARY_PERIOD,
  focalField::Symbol = :adjnetassets)::DataFrame

  local outdf::DataFrame = DataFrame(fisyr = collect(summaryPeriod)) #this will hold the summary stats
  local subdf::SubDataFrame =  view(data.df,  (y::Int-> y ∈ summaryPeriod).(data.df[:fisyr]),
    [:ein, :fisyr, :charitytype, focalField])
  local nyears::Int = length(summaryPeriod) # of years

  local N::Vector{MInt} = Vector{MInt}(undef, nyears)
  local total::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local mean::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local skew::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local kurt::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local fracFoundation::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local wealth80::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local wealth95::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)

  #construct the info for each year

  #make an index so we can use the high performance by function
  local outdfIndex::Dict = Dict(outdf[i, :fisyr]=>i for i::Int ∈ 1:nyears)

  #this format is much faster since it does not imply recombination
  subdfByfisyr::GroupedDataFrame = groupby(subdf, :fisyr)
  for ssubdf::SubDataFrame in subdfByfisyr
  #by(subdf, :fisyr) do ssubdf
    r::Int = outdfIndex[ssubdf[1, :fisyr]]
    focal::Vector = ssubdf[:, focalField]
    nfirms::Int = length(focal)
    #println("Vector Length nfirms: $nfirms")
    focal = processmissing(focal) #handle missing values

    #println("Vector Length: $(length(focal))")
    N[r] = length(focal)
    total[r] = sum(focal)
    mean[r] = StatsBase.mean(focal)
    stddev[r] = StatsBase.std(focal)
    skew[r] = StatsBase.skewness(focal)
    kurt[r] = StatsBase.kurtosis(focal)
    p5[r] = StatsBase.quantile(focal, .05)
    p25[r] = StatsBase.quantile(focal, .25)
    p50[r] = StatsBase.quantile(focal, 0.5)
    p75[r] = StatsBase.quantile(focal, 0.75)
    p95[r] = StatsBase.quantile(focal, 0.95)
    fracFoundation[r] = sum(ssubdf[:charitytype] .== :pf) / N[r]

    #now get a couple distributional values
    cumv::Vector{Float64} = similar(focal) #will hold the running total
    cumv[1] = focal[1]
    for j::Int ∈ 2:length(focal)
      cumv[j] = cumv[j-1] + focal[j]
    end
    cumv .= cumv ./ cumv[end] #normalize

    wealth80[r] = sum(cumv .≤ 0.8) / N[r]
    wealth95[r] = sum(cumv .≤ 0.95) / N[r]
  end

  #now put the data into the dataframe
  outdf[:N] = N
  outdf[:total] = total
  outdf[:mean] = mean
  outdf[:stddev] = stddev
  outdf[:skew] = skew
  outdf[:kurt] = kurt
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95
  outdf[:fracFoundation] = fracFoundation
  outdf[:wealth80] = wealth80
  outdf[:wealth95] = wealth95

  return outdf
end

function assetsByYear(outdf::DataFrame; summaryPeriod::UnitRange = ASSETS_SUMMARY_PERIOD,
  focalField::Symbol = :adjnetassets, decimals::Int = DECIMALS)::String

  columnSymbols::Vector{Symbol} = setdiff(names(outdf), [:fisyr])

  numColumns = length(columnSymbols) #number of columns excluding the year
  columnIndex::Dict = Dict(columnSymbols[i]=>i for i ∈ 1:numColumns)

  #include all columns in the names
  columnNames::Vector{Vector{String}} = [(string).(columnSymbols)]
  descRowNames::Vector{String} = (string).(outdf[:fisyr])

  numRows = size(outdf,1)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(undef, numColumns)).(1:length(descRowNames))

  for r::Int ∈ 1:numRows
    haskey(columnIndex, :N) && (
      descContent[r][columnIndex[:N]] = "$(outdf[r, :N])")
    haskey(columnIndex, :total) && (
      descContent[r][columnIndex[:total]] = num2Str(outdf[r, :total], 0, scaleFactor=10^-6))
    haskey(columnIndex, :mean) && (
      descContent[r][columnIndex[:mean]] = num2Str(outdf[r, :mean], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :stddev) && (
      descContent[r][columnIndex[:stddev]] = num2Str(outdf[r, :stddev], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :skew) && (
      descContent[r][columnIndex[:skew]] = num2Str(outdf[r, :skew], decimals, scaleFactor=1.0))
    haskey(columnIndex, :kurt) && (
      descContent[r][columnIndex[:kurt]] = num2Str(outdf[r, :kurt], decimals, scaleFactor=1.0))
    haskey(columnIndex, :p5) && (
      descContent[r][columnIndex[:p5]] = num2Str(outdf[r, :p5], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p25) && (
      descContent[r][columnIndex[:p25]] = num2Str(outdf[r, :p25], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p50) && (
      descContent[r][columnIndex[:p50]] = num2Str(outdf[r, :p50], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p75) && (
      descContent[r][columnIndex[:p75]] = num2Str(outdf[r, :p75], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p95) && (
      descContent[r][columnIndex[:p95]] = num2Str(outdf[r, :p95], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :fracFoundation) && (
      descContent[r][columnIndex[:fracFoundation]] = num2Str(outdf[r, :fracFoundation], 2))
    haskey(columnIndex, :wealth80) && (
      descContent[r][columnIndex[:wealth80]] = num2Str(outdf[r, :wealth80], 2))
    haskey(columnIndex, :wealth95) && (
      descContent[r][columnIndex[:wealth95]] = num2Str(outdf[r, :wealth95], 2))
  end

  assetsByYearTable::String = texTable(titleCaption="Descriptive Data by Year",
      caption="""See Tex File""", #caption
      colNames= columnNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true
    )

    return assetsByYearTable

end

function assetsSummaryDF(data::NCCSData; summaryPeriod::UnitRange = ASSETS_SUMMARY_PERIOD,
  focalField::Symbol = :adjnetassets, decimals::Int = DECIMALS,
  byCharityField::Symbol = BY_CHARITY_FIELD,
  byCategoryField::Symbol = BY_CATEGORY_FIELD)::DataFrame


  local subdf::SubDataFrame =  view(data.df,  (y::Int-> y ∈ summaryPeriod).(data.df[:fisyr]),
    [:ein; :fisyr; :fullsample; byCharityField; byCategoryField; focalField])
  local byCharityValues::Vector{Symbol} = unique(processmissing(subdf[byCharityField]))
  local byCategoryValues::Vector{Symbol} = unique(processmissing(subdf[byCategoryField]))
  local rowFields::Vector{Symbol} = [byCharityValues; byCategoryValues; :fullsample]
  local nrows::Int = length(rowFields) # of rows we will have

  local N::Vector{MInt} = Vector{MInt}(undef, nrows)
  local minyear::Vector{MInt} = Vector{MInt}(undef, nrows)
  local maxyear::Vector{MInt} = Vector{MInt}(undef, nrows)
  local total::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local mean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local skew::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local kurt::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)

  #performance optimized way to split up the dataframe
  local outdfIndex::Dict = Dict(rowFields[i]=>i for i::Int ∈ 1:nrows)
  local rowNames::Vector{Symbol} =
    (s::Symbol->haskey(charitytype, s) ? charitytype[s] : s).(rowFields) #clarify some of the names
  local outdf::DataFrame = DataFrame(rownames = rowNames) #this will hold the summary stats

  #for each field we are considering
  #println("keys: $(collect(keys(outdfIndex)))")
  for f::Symbol ∈ [byCharityField; byCategoryField; :fullsample]
    #println("field $f")
    subdfByf::GroupedDataFrame = groupby(subdf, f)
    for ssubdf::SubDataFrame ∈ subdfByf
    #by(subdf, f) do ssubdf::SubDataFrame
      if !ismissing(ssubdf[1, f])
        r::Int = outdfIndex[ssubdf[1, f]]

        focal::Vector = ssubdf[:, focalField]
        nfirms::Int = length(focal)
        focal = processmissing(focal) #handle missing values

        N[r] = length(focal)
        minyear[r] = minimum(processmissing(ssubdf[:, :fisyr]))
        maxyear[r] = maximum(processmissing(ssubdf[:, :fisyr]))
        total[r] = StatsBase.sum(focal)
        mean[r] = StatsBase.mean(focal)
        stddev[r] = StatsBase.std(focal)
        skew[r] = StatsBase.skewness(focal)
        kurt[r] = StatsBase.kurtosis(focal)
        p5[r] = StatsBase.quantile(focal, .05)
        p25[r] = StatsBase.quantile(focal, .25)
        p50[r] = StatsBase.quantile(focal, 0.5)
        p75[r] = StatsBase.quantile(focal, 0.75)
        p95[r] = StatsBase.quantile(focal, 0.95)
      end
    end
  end

  #now put the data into the dataframe
  outdf[:N] = N
  outdf[:minyear] = minyear
  outdf[:maxyear] = maxyear
  outdf[:total] = total
  outdf[:mean] = mean
  outdf[:stddev] = stddev
  outdf[:skew] = skew
  outdf[:kurt] = kurt
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95

  return outdf
end

function assetsSummary(outdf::DataFrame; summaryPeriod::UnitRange = ASSETS_SUMMARY_PERIOD,
  focalField::Symbol = :adjnetassets, decimals::Int = DECIMALS,
  byCharityField::Symbol = BY_CHARITY_FIELD,
  byCategoryField::Symbol = BY_CATEGORY_FIELD)::String

  columnSymbols::Vector{Symbol} = setdiff(names(outdf), [:rownames])

  numColumns = length(columnSymbols) #number of columns excluding the year
  columnIndex::Dict = Dict(columnSymbols[i]=>i for i ∈ 1:numColumns)

  #include all columns in the names
  columnNames::Vector{Vector{String}} = [(string).(columnSymbols)]
  descRowNames::Vector{String} = (string).(outdf[:rownames])

  #NOTE: a couple of hacks to improve readability
  descRowNames[length(keys(charitytype))+1] = "\\midrule $(descRowNames[length(keys(charitytype))+1])"
  descRowNames[end] = "\\midrule $(descRowNames[end])"

  numRows = size(outdf,1)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = (
    (i::Int)->Vector{String}(undef, numColumns)).(1:length(descRowNames))

  for r::Int ∈ 1:numRows
    haskey(columnIndex, :N) && (
      descContent[r][columnIndex[:N]] = "$(outdf[r, :N])")
    haskey(columnIndex, :minyear) && (
      descContent[r][columnIndex[:minyear]] = "$(outdf[r, :minyear])")
    haskey(columnIndex, :maxyear) && (
      descContent[r][columnIndex[:maxyear]] = "$(outdf[r, :maxyear])")
    haskey(columnIndex, :total) && (
      descContent[r][columnIndex[:total]] = num2Str(outdf[r, :total], 0, scaleFactor=10^-6))
    haskey(columnIndex, :mean) && (
      descContent[r][columnIndex[:mean]] = num2Str(outdf[r, :mean], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :stddev) && (
      descContent[r][columnIndex[:stddev]] = num2Str(outdf[r, :stddev], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :skew) && (
      descContent[r][columnIndex[:skew]] = num2Str(outdf[r, :skew], decimals, scaleFactor=1.0))
    haskey(columnIndex, :kurt) && (
      descContent[r][columnIndex[:kurt]] = num2Str(outdf[r, :kurt], decimals, scaleFactor=1.0))
    haskey(columnIndex, :p5) && (
      descContent[r][columnIndex[:p5]] = num2Str(outdf[r, :p5], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p25) && (
      descContent[r][columnIndex[:p25]] = num2Str(outdf[r, :p25], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p50) && (
      descContent[r][columnIndex[:p50]] = num2Str(outdf[r, :p50], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p75) && (
      descContent[r][columnIndex[:p75]] = num2Str(outdf[r, :p75], decimals, scaleFactor=10^-6))
    haskey(columnIndex, :p95) && (
      descContent[r][columnIndex[:p95]] = num2Str(outdf[r, :p95], decimals, scaleFactor=10^-6))
  end

  assetsByYearTable::String = texTable(titleCaption="Descriptive Data by Charity Type and Category",
      caption="""See Tex File""", #caption
      colNames= columnNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true
    )

    return assetsByYearTable
end

function returnsByYearDF(data::NCCSData; summaryPeriod::UnitRange = RETURNS_SUMMARY_PERIOD,
  focalField::Symbol = :return, assetsField::Symbol = :adjnetassets)::DataFrame

  local outdf::DataFrame = DataFrame(fisyr = collect(summaryPeriod)) #this will hold the summary stats
  local subdf::SubDataFrame =  view(data.df,  (y::Int-> y ∈ summaryPeriod).(data.df[:fisyr]),
    [:ein, :fisyr, :charitytype, focalField, assetsField])
  local nyears::Int = length(summaryPeriod) # of years

  local N::Vector{MInt} = Vector{MInt}(undef, nyears)
  local equalmean::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local valuemean::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local skew::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local kurt::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nyears)

  #construct the info for each year
  local outdfIndex::Dict = Dict(outdf[i, :fisyr]=>i for i::Int ∈ 1:nyears)

  subdfByfisyr::GroupedDataFrame = groupby(subdf, :fisyr)
  for ssubdf::SubDataFrame ∈ subdfByfisyr
    r::Int = outdfIndex[ssubdf[1, :fisyr]]
    focal::Vector = ssubdf[:, focalField]
    nfirms::Int = length(focal)
    #println("Vector Length nfirms: $nfirms")
    focal = processmissing(focal) #handle missing values

    #println("Vector Length: $(length(focal))")
    N[r] = length(focal)
    equalmean[r] = StatsBase.mean(focal)
    valuemean[r] = #compute the value-weighted average
      sum(processmissing(ssubdf[:, focalField] .* ssubdf[:, assetsField])) /
      sum(processmissing(ssubdf[(!ismissing).(ssubdf[:, focalField]), assetsField]))
    stddev[r] = StatsBase.std(focal)
    skew[r] = StatsBase.skewness(focal)
    kurt[r] = StatsBase.kurtosis(focal)
    p5[r] = StatsBase.quantile(focal, .05)
    p25[r] = StatsBase.quantile(focal, .25)
    p50[r] = StatsBase.quantile(focal, 0.5)
    p75[r] = StatsBase.quantile(focal, 0.75)
    p95[r] = StatsBase.quantile(focal, 0.95)
  end

  #now put the data into the dataframe
  outdf[:N] = N
  outdf[:equalmean] = equalmean
  outdf[:valuemean] = valuemean
  outdf[:stddev] = stddev
  outdf[:skew] = skew
  outdf[:kurt] = kurt
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95

  return outdf
end

function returnsByYear(outdf::DataFrame; summaryPeriod::UnitRange = RETURNS_SUMMARY_PERIOD,
  focalField::Symbol = :return, decimals::Int = DECIMALS_RETURN)::String

  columnSymbols::Vector{Symbol} = setdiff(names(outdf), [:fisyr])

  numColumns = length(columnSymbols) #number of columns excluding the year
  columnIndex::Dict = Dict(columnSymbols[i]=>i for i ∈ 1:numColumns)

  #include all columns in the names
  columnNames::Vector{Vector{String}} = [(string).(columnSymbols)]
  descRowNames::Vector{String} = (string).(outdf[:fisyr])

  numRows = size(outdf,1)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(undef, numColumns)).(1:length(descRowNames))

  for r::Int ∈ 1:numRows
    haskey(columnIndex, :N) && (
      descContent[r][columnIndex[:N]] = "$(outdf[r, :N])")
    haskey(columnIndex, :equalmean) && (
      descContent[r][columnIndex[:equalmean]] = num2Str(outdf[r, :equalmean], decimals, scaleFactor=100.0))
    haskey(columnIndex, :valuemean) && (
      descContent[r][columnIndex[:valuemean]] = num2Str(outdf[r, :valuemean], decimals, scaleFactor=100.0))
    haskey(columnIndex, :stddev) && (
      descContent[r][columnIndex[:stddev]] = num2Str(outdf[r, :stddev], decimals, scaleFactor=100.0))
    haskey(columnIndex, :skew) && (
      descContent[r][columnIndex[:skew]] = num2Str(outdf[r, :skew], decimals, scaleFactor=1.0))
    haskey(columnIndex, :kurt) && (
      descContent[r][columnIndex[:kurt]] = num2Str(outdf[r, :kurt], decimals, scaleFactor=1.0))
    haskey(columnIndex, :p5) && (
      descContent[r][columnIndex[:p5]] = num2Str(outdf[r, :p5], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p25) && (
      descContent[r][columnIndex[:p25]] = num2Str(outdf[r, :p25], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p50) && (
      descContent[r][columnIndex[:p50]] = num2Str(outdf[r, :p50], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p75) && (
      descContent[r][columnIndex[:p75]] = num2Str(outdf[r, :p75], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p95) && (
      descContent[r][columnIndex[:p95]] = num2Str(outdf[r, :p95], decimals, scaleFactor=100.0))
  end

  returnByYearTable::String = texTable(titleCaption="Descriptive Return Data by Year",
      caption="""See Tex File""", #caption
      colNames= columnNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true
    )

    return returnByYearTable
end

function returnsSummaryDF(data::NCCSData; summaryPeriod::UnitRange = RETURNS_SUMMARY_PERIOD,
  focalField::Symbol = :return, decimals::Int = DECIMALS_RETURN,
  byCharityField::Symbol = BY_CHARITY_FIELD,
  byCategoryField::Symbol = BY_CATEGORY_FIELD, assetsField::Symbol = :adjnetassets)::DataFrame


  local subdf::SubDataFrame =  view(data.df,  (y::Int-> y ∈ summaryPeriod).(data.df[:fisyr]),
    [:ein; :fisyr; :fullsample; byCharityField; byCategoryField; focalField; assetsField])
  local byCharityValues::Vector{Symbol} = unique(processmissing(subdf[byCharityField]))
  local byCategoryValues::Vector{Symbol} = unique(processmissing(subdf[byCategoryField]))
  local rowFields::Vector{Symbol} = [byCharityValues; byCategoryValues; :fullsample]
  local nrows::Int = length(rowFields) # of rows we will have

  local N::Vector{MInt} = Vector{MInt}(undef, nrows)
  local minyear::Vector{MInt} = Vector{MInt}(undef, nrows)
  local maxyear::Vector{MInt} = Vector{MInt}(undef, nrows)
  local equalmean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local valuemean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local skew::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local kurt::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)

  #performance optimized way to split up the dataframe
  local outdfIndex::Dict = Dict(rowFields[i]=>i for i::Int ∈ 1:nrows)
  local rowNames::Vector{Symbol} =
    (s::Symbol->haskey(charitytype, s) ? charitytype[s] : s).(rowFields) #clarify some of the names
  local outdf::DataFrame = DataFrame(rownames = rowNames) #this will hold the summary stats

  #for each field we are considering
  for f::Symbol ∈ [byCharityField; byCategoryField; :fullsample]
    subdfByf::GroupedDataFrame = groupby(subdf, f)
    for ssubdf::SubDataFrame ∈ subdfByf
    #by(subdf, f) do ssubdf::SubDataFrame
      if !ismissing(ssubdf[1, f])
        r::Int = outdfIndex[ssubdf[1, f]]

        focal::Vector = ssubdf[:, focalField]
        nfirms::Int = length(focal)
        focal = processmissing(focal) #handle missing values

        N[r] = length(focal)
        minyear[r] = minimum(processmissing(ssubdf[:fisyr]))
        maxyear[r] = maximum(processmissing(ssubdf[:fisyr]))
        equalmean[r] = StatsBase.mean(focal)
        valuemean[r] = #compute the value-weighted average
          sum(processmissing(ssubdf[:, focalField] .* ssubdf[:, assetsField])) /
          sum(processmissing(ssubdf[(!ismissing).(ssubdf[:, focalField]), assetsField]))
        stddev[r] = StatsBase.std(focal)
        skew[r] = StatsBase.skewness(focal)
        kurt[r] = StatsBase.kurtosis(focal)
        p5[r] = StatsBase.quantile(focal, .05)
        p25[r] = StatsBase.quantile(focal, .25)
        p50[r] = StatsBase.quantile(focal, 0.5)
        p75[r] = StatsBase.quantile(focal, 0.75)
        p95[r] = StatsBase.quantile(focal, 0.95)
      end
    end
  end

  #now put the data into the dataframe
  outdf[:N] = N
  outdf[:minyear] = minyear
  outdf[:maxyear] = maxyear
  outdf[:equalmean] = equalmean
  outdf[:valuemean] = valuemean
  outdf[:stddev] = stddev
  outdf[:skew] = skew
  outdf[:kurt] = kurt
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95

  return outdf
end

function returnsSummary(outdf::DataFrame; summaryPeriod::UnitRange = RETURNS_SUMMARY_PERIOD,
  focalField::Symbol = :return, decimals::Int = DECIMALS_RETURN,
  byCharityField::Symbol = BY_CHARITY_FIELD,
  byCategoryField::Symbol = BY_CATEGORY_FIELD)::String

  columnSymbols::Vector{Symbol} = setdiff(names(outdf), [:rownames])

  numColumns = length(columnSymbols) #number of columns excluding the year
  columnIndex::Dict = Dict(columnSymbols[i]=>i for i ∈ 1:numColumns)

  #include all columns in the names
  columnNames::Vector{Vector{String}} = [(string).(columnSymbols)]
  descRowNames::Vector{String} = (string).(outdf[:rownames])

  #NOTE: a couple of hacks to improve readability
  descRowNames[length(keys(charitytype))+1] = "\\midrule $(descRowNames[length(keys(charitytype))+1])"
  descRowNames[end] = "\\midrule $(descRowNames[end])"

  numRows = size(outdf,1)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = ((i::Int)->Vector{String}(undef, numColumns)).(1:length(descRowNames))

  for r::Int ∈ 1:numRows
    haskey(columnIndex, :N) && (
      descContent[r][columnIndex[:N]] = "$(outdf[r, :N])")
    haskey(columnIndex, :minyear) && (
      descContent[r][columnIndex[:minyear]] = "$(outdf[r, :minyear])")
    haskey(columnIndex, :maxyear) && (
      descContent[r][columnIndex[:maxyear]] = "$(outdf[r, :maxyear])")
    haskey(columnIndex, :equalmean) && (
      descContent[r][columnIndex[:equalmean]] = num2Str(outdf[r, :equalmean], decimals, scaleFactor=100.0))
    haskey(columnIndex, :valuemean) && (
      descContent[r][columnIndex[:valuemean]] = num2Str(outdf[r, :valuemean], decimals, scaleFactor=100.0))
    haskey(columnIndex, :stddev) && (
      descContent[r][columnIndex[:stddev]] = num2Str(outdf[r, :stddev], decimals, scaleFactor=100.0))
    haskey(columnIndex, :skew) && (
      descContent[r][columnIndex[:skew]] = num2Str(outdf[r, :skew], decimals, scaleFactor=1.0))
    haskey(columnIndex, :kurt) && (
      descContent[r][columnIndex[:kurt]] = num2Str(outdf[r, :kurt], decimals, scaleFactor=1.0))
    haskey(columnIndex, :p5) && (
      descContent[r][columnIndex[:p5]] = num2Str(outdf[r, :p5], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p25) && (
      descContent[r][columnIndex[:p25]] = num2Str(outdf[r, :p25], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p50) && (
      descContent[r][columnIndex[:p50]] = num2Str(outdf[r, :p50], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p75) && (
      descContent[r][columnIndex[:p75]] = num2Str(outdf[r, :p75], decimals, scaleFactor=100.0))
    haskey(columnIndex, :p95) && (
      descContent[r][columnIndex[:p95]] = num2Str(outdf[r, :p95], decimals, scaleFactor=100.0))
  end

  assetsByYearTable::String = texTable(titleCaption="Descriptive Return Data by Charity Type and Category",
      caption="""See Tex File""", #caption
      colNames= columnNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true
    )

    return assetsByYearTable
end


function getAssetsTables(data::NCCSData; summaryPeriod::UnitRange = ASSETS_SUMMARY_PERIOD,
  focalField::Symbol = :adjnetassets, sorted::Bool=false, outputPath::String = TABLE_PATH,
  headerName::String = HEADER_NAME, footerName::String = FOOTER_NAME)::Nothing

  #sort if we need to
  (!sorted) && sort!(data.df, [:fisyr, :adjnetassets], (:fisyr, order(:adjnetassets, rev=true)))

  #we get a full table and a shorter table for the decks
  #start with getting the wealth info by year
  outassetsyeardf::DataFrame = assetsByYearDF(data, focalField=focalField, summaryPeriod=summaryPeriod)
  assetsByYearTable::String = assetsByYear(outassetsyeardf, focalField=focalField, summaryPeriod=summaryPeriod)
  assetsByYearTableShort::String = assetsByYear(
    outassetsyeardf[[:fisyr, :N, :total, :mean, :stddev, :p5, :p50, :p95, :wealth80, :wealth95]],
    focalField=focalField, summaryPeriod=summaryPeriod)

  #now get the wealth info by category
  outassetssummarydf::DataFrame = assetsSummaryDF(data, focalField=focalField, summaryPeriod=summaryPeriod)
  assetsSummaryTable::String = assetsSummary(outassetssummarydf, focalField=focalField, summaryPeriod=summaryPeriod)
  assetsSummaryTableShort::String = assetsSummary(
    outassetssummarydf[[:rownames, :N, :minyear, :maxyear, :total, :mean, :stddev, :p5, :p50, :p95]],
    focalField=focalField, summaryPeriod=summaryPeriod)

  dataReligious::NCCSData=NCCSData(
    data.df[(s::MSymbol->((!ismissing(s)) && s==:ReligionRelated)).(data.df[:categoryname]),:])
  outreligousassetsdf::DataFrame = assetsByYearDF(dataReligious, focalField=focalField, summaryPeriod=summaryPeriod)
  religiousAssetsByYearTable::String =
    assetsByYear(outreligousassetsdf, focalField=focalField, summaryPeriod=summaryPeriod)

  #write the tables
  writeNakedTable(assetsByYearTable, path=outputPath,  outName = "assetsTableByYear.tex")
  writeNakedTable(assetsByYearTableShort, path=outputPath,  outName = "assetsTableByYearShort.tex")
  writeNakedTable(assetsSummaryTable, path=outputPath,  outName = "assetsSummaryTable.tex")
  writeNakedTable(assetsSummaryTableShort, path=outputPath,  outName = "assetsSummaryTableShort.tex")
  writeNakedTable(religiousAssetsByYearTable, path=outputPath,
    outName = "assetsReligiousTable.tex")

  return nothing
end

function getReturnsTables(data::NCCSData; summaryPeriod::UnitRange = RETURNS_SUMMARY_PERIOD,
  focalField::Symbol = :return, sorted::Bool=false, outputPath::String = TABLE_PATH,
  headerName::String = HEADER_NAME, footerName::String = FOOTER_NAME)::Nothing

  #sort if we need to
  (!sorted) && sort!(data.df, [:fisyr, :adjnetassets], (:fisyr, order(:adjnetassets, rev=true)))

  outreturnsyeardf::DataFrame = returnsByYearDF(data, focalField=focalField, summaryPeriod=summaryPeriod)
  returnsByYearTable::String = returnsByYear(outreturnsyeardf, focalField=focalField, summaryPeriod=summaryPeriod)
  returnsByYearTableShort::String = returnsByYear(
    outreturnsyeardf[[:fisyr, :N, :equalmean, :valuemean, :stddev, :skew, :p5, :p50, :p95]],
    focalField=focalField, summaryPeriod=summaryPeriod)


  outreturnssummarydf::DataFrame = returnsSummaryDF(data, focalField=focalField, summaryPeriod=summaryPeriod)
  returnsSummaryTable::String = returnsSummary(outreturnssummarydf, focalField=focalField, summaryPeriod=summaryPeriod)
  returnsSummaryTableShort::String = returnsSummary(
    outreturnssummarydf[[:rownames, :N, :equalmean, :valuemean, :p5, :p50, :p95]],
    focalField=focalField, summaryPeriod=summaryPeriod)

  #write the tables
  writeNakedTable(returnsByYearTable, path=outputPath,  outName = "returnsTableByYear.tex")
  writeNakedTable(returnsByYearTableShort, path=outputPath,  outName = "returnsTableByYearShort.tex")
  writeNakedTable(returnsSummaryTable, path=outputPath,  outName = "returnsSummaryTable.tex")
  writeNakedTable(returnsSummaryTableShort, path=outputPath,  outName = "returnsSummaryTableShort.tex")

  return nothing
end

#NOTE: Consider slicing by type with collors being the years

function makeDescriptiveTables(data::NCCSData)::Nothing

  sort!(data.df, (:fisyr, order(:adjnetassets, rev=true)))

  data.df[:fullsample] = :fullsample #this is a hack for analyzing a full sample

  getReturnsTables(data, sorted=true)
  getAssetsTables(data, sorted=true)
  getTopInstitutions(data, sorted=true)

  deletecols!(data.df, :fullsample)


  return nothing
end
