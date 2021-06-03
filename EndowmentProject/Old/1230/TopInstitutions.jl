const TOP_INSTITUTION_FIELDS = [:name, :fisyr, :netassets, :charitytype,
  :categorynamefiltered]

const TOP_INSTITUTION_METHODS = [:rank, :quantile, :cumulative]

const TOP_INSTITUTIONS_YEAR = 2015
const TOP_INSTITUTIONS_QUANTILE = 0.95
const TOP_INSTITUTIONS_CUMULATIVE = 0.95
const TOP_INSTITUTIONS_RANK = 500
const TOP_INSTITUTIONS_PREFIX = "top_"

processmissing(v::AbstractVector) = collect(skipmissing(v))

#writes out the institutions in the subdataframe
function writeTopInstitutions(subdf::SubDataFrame, label::String;
  outputPath::String = TABLE_PATH,
  topInstitutionFields::Vector{Symbol} = TOP_INSTITUTION_FIELDS)::Nothing

  outPath = "$outputPath\\$(label).csv"
  uCSV.write(outPath, subdf[topInstitutionFields])

  return nothing
end

getTopInstitutionsRank(subdf::SubDataFrame,
  topInstitutionsRank::Int)::SubDataFrame = view(subdf, 1:topInstitutionsRank, :)

function getTopInstitutionsQuantile(subdf::SubDataFrame,
  topInstitutionsQuantile::Float64)
  v::Vector{Float64} = subdf[:netassets]
  qFunc::Function = StatsBase.ecdf(v)
  dfquantile::Vector{Float64} = (qFunc).(subdf[:netassets])

  return view(subdf, dfquantile .> topInstitutionsQuantile, :)
end

function getTopInstitutionsCumulative(subdf::SubDataFrame,
  topInstitutionsCumulative::Float64)

  v::Vector{Float64} = subdf[:netassets]
  cumv::Vector{Float64} = similar(v) #will hold the running total
  cumv[1] = v[1]
  for i::Int ∈ 2:length(v)
    cumv[i] = cumv[i-1] + v[i]
  end
  cumv .= cumv ./ cumv[end] #normalize

  return view(subdf, cumv .< topInstitutionsCumulative, :)
end

#gets the top institutions given a method and SORTED subdataframe
function getTopInstitutions(subdf::SubDataFrame, topMethod::Symbol;
  topInstitutionsYear::Int = TOP_INSTITUTIONS_YEAR,
  topInstitutionsQuantile::Float64 = TOP_INSTITUTIONS_QUANTILE,
  topInstitutionsRank::Int = TOP_INSTITUTIONS_RANK,
  topInstitutionsCumulative::Float64 = TOP_INSTITUTIONS_CUMULATIVE,
  topInstitutionsPrefix::String = TOP_INSTITUTIONS_PREFIX)

  if topMethod == :rank #this assumes we want the top x ranked asset holders
    subdf = getTopInstitutionsRank(subdf, topInstitutionsRank)

  elseif topMethod == :quantile #this assumes we want the top  x percentile asset holders
    subdf = getTopInstitutionsQuantile(subdf, topInstitutionsQuantile)

  elseif topMethod == :cumulative #this assumes we want the top holders which account for x% of assets
    subdf = getTopInstitutionsCumulative(subdf, topInstitutionsCumulative)
  elseif topMethod ≠ :full
      @warn "Top institutions method not detected, assuming full file"
  end

  label::String =
    "$(topInstitutionsPrefix)_$(topInstitutionsYear)_$(topMethod)"
  writeTopInstitutions(subdf, label)

  return nothing
end

#takes in an unsorted data object and gets the top institutions by assets
function getTopInstitutions(data::NCCSData, topMethod::Symbol;
    topInstitutionsYear::Int = TOP_INSTITUTIONS_YEAR)
  sort!(data.df, :netassets, rev=true)
  subdf::SubDataFrame =
    view(data.df, (data.df[:fisyr] .== topInstitutionsYear) .& (
    (!ismissing).(data.df[:netassets])), :)



  getTopInstitutions(subdf, topMethod)

  return nothing
end

#convenience method for running multiple methods at once
function getTopInstitutions(data::NCCSData;
  topMethods::Vector{Symbol} = TOP_INSTITUTION_METHODS,
  topInstitutionsYear::Int = TOP_INSTITUTIONS_YEAR,
  sorted::Bool = false)

  (!sorted) && sort!(data.df, [:fisyr, :adjnetassets], (:fisyr, order(:adjnetassets, rev=true)))
  subdf::SubDataFrame =
    view(data.df, (data.df[:fisyr] .== topInstitutionsYear) .& (
    (!ismissing).(data.df[:netassets])),:)


  #Execute the tasks
  for s::Symbol ∈ topMethods
    getTopInstitutions(subdf, s)
  end

  return nothing
end
