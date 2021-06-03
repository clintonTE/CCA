const TOP_INSTITUTION_FIELDS = [:name, :fisyr, :netassets, :charitytype,
  :categorynamefiltered]

const TOP_INSTITUTION_METHODS = [:rank, :quantile, :cumulative]

const TOP_INSTITUTIONS_YEAR = 2015
const TOP_INSTITUTIONS_QUANTILE = 0.95
const TOP_INSTITUTIONS_CUMULATIVE = 0.95
const TOP_INSTITUTIONS_RANK = 500
const TOP_INSTITUTIONS_PREFIX = "top_"
const WORST_MIN_ASSETS = 100_000_000
const WORST_NUM_YEARS = 10

processmissing(v::AbstractVector) = collect(skipmissing(v))

#writes out the institutions in the subdataframe
function writeTopInstitutions(subdf::SubDataFrame, label::String;
  outputPath::String = TABLE_PATH,
  topInstitutionFields::Vector{Symbol} = TOP_INSTITUTION_FIELDS)::Nothing

  outPath = "$outputPath\\$(label).csv"
  CSV.write(outPath,subdf[topInstitutionFields])

  return nothing
end

getTopInstitutionsRank(subdf::SubDataFrame,
  topInstitutionsRank::Int)::SubDataFrame = view(subdf, 1:topInstitutionsRank, :)

function getTopInstitutionsQuantile(subdf::SubDataFrame,
  topInstitutionsQuantile::Float64)
  v::Vector{Float64} = subdf[:netassets]
  qFunc = StatsBase.ecdf(v)
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

  worstPerformers(data)

  return nothing
end

function worstPerformers(data::NCCSData; worstminassets = WORST_MIN_ASSETS,
  numyears = WORST_NUM_YEARS, lastyear = 2015,
  outputpath::String = TABLE_PATH)::Nothing


  firstyear::Int = lastyear - numyears
  outpath = "$outputpath\\worst_$(firstyear)_$(lastyear)_$(worstminassets ÷ 1_000_000).csv"

  #partition off the part we want
  validfirmsdf = view(data.df, (data.df[:adjnetassets] .> worstminassets) .&
    (data.df[:fisyr] .≤ lastyear) .& (data.df[:fisyr] .> firstyear),
    [:ein, :name, :excess, :return, :t30ret, :sp500ret, :adjnetassets, :monthscovered, :fisyr, :datefiscal])

  #form the consolidated df

  worstdf::DataFrame = DataFrame(ein=validfirmsdf[:ein])
  N::Int = size(worstdf, 1)

  worstdf[:avgwealth] = Vector{MFloat64}(missing, N)
  worstdf[:N] = Vector{MInt}(missing, N)
  worstdf[:monthscovered] = Vector{MInt}(missing, N)
  worstdfindex::Dict = Dict(worstdf[i,:ein] => worstdf[i,:] for i::Int  ∈ 1:N)
  worstdf[:name] = Vector{MString}(missing, N)
  worstdf[:sp500ret] = Vector{MFloat64}(missing, N)
  worstdf[:t30ret] = Vector{MFloat64}(missing, N)
  worstdf[:cagr] = Vector{MFloat64}(missing, N)

  #consolidate the eprformacne data
  for subdf ∈ groupby(validfirmsdf, :ein)
    worstrow = worstdfindex[subdf[1,:ein]]
    worstrow[:name] = subdf[1,:name]

    worstrow[:N] = size(subdf,1)
    worstrow[:monthscovered] = sum(subdf[:monthscovered])

    totalreturn::MFloat64 = prod(1.0 .+ subdf[:return])
    worstrow[:cagr] = totalreturn^(12. / worstrow[:monthscovered]) .- 1

    totalrfr::MFloat64 = prod(1.0 .+ subdf[:t30ret])
    worstrow[:t30ret] = totalrfr^(12. / worstrow[:monthscovered]) .- 1.

    totalsp500::MFloat64 = prod(1.0 .+ subdf[:sp500ret])
    worstrow[:sp500ret] = totalrfr^(12. / worstrow[:monthscovered]) .- 1.

    worstrow[:avgwealth] = mean(subdf[:adjnetassets]) ./ 1_000_000
  end

  #compute excess returns
  worstdf[:excess] = worstdf[:cagr] .- worstdf[:t30ret]
  worstdf = worstdf[completecases(worstdf[[:cagr,:t30ret,:avgwealth,:sp500ret,:excess]]),:]
  worstdf = worstdf[worstdf[:N] .== numyears,:]

  sort!(worstdf, :excess, rev=false)

  worstdf |> CSV.write(outpath)

  return nothing
end
