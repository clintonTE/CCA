
#data type for the investment universe
struct Universe
  df::DataFrame #source dataframe must be properly formatted with a :date column but can include extra observations
  datebegin::Date
  dateend::Date
  returns::Matrix{MFloat64} #should contain the numerical returns data in the DataFrame
  #NOTE for performance reasons, this is M by N!!!

  μ::Vector{Float64} #we can calc this in advance
  N::Int #number of observations in the sample period
  M::Int #number of return series

  missingindices::Vector{Vector{Int}} # a collection of all non-missing indices in returns
end

#creates an investment universe which will be used to generate portfolio returns
function Universe(df::AbstractDataFrame; datebegin::Date = minimum(df[:date]),
  dateend::Date = maximum(df[:date]))::Universe
  subdf::SubDataFrame = view(df,(df[:date] .≥ datebegin) .& (df[:date] .≤ dateend),:)

  local N::Int = size(subdf,1) #get the size in number of observations
  local M::Int = size(subdf, 2) - 1 #exclude date column

  local xcols::Vector{Symbol} = setdiff(names(df), [:date])
  local returns::Matrix{MFloat64} = Matrix{MFloat64}(subdf[xcols])'

  #compute the mean for use later
  local μ::Vector{Float64} = (r::Int->mean(skipmissing(returns[r,:]))).(1:M)

  #map the missing values
  local missingindices::Vector{Vector{Int}} = Vector{Vector{Float64}}()
  for c ∈ 1:N
    push!(missingindices, collect(1:M)[(!ismissing).(returns[:, c])])
  end

  return Universe(df,datebegin, dateend, returns, μ, N, M, missingindices)
end

function Universe(;datapath::String = DATA_PATH, dataname::String=DATA_NAME,
  datebegin::Date = DATE_BEGIN, dateend::Date = DATE_END)

  dateformat::DateFormat = DateFormat(DATE_FORMAT)

  df::DataFrame = CSV.read("$datapath\\$dataname.csv", missingstring="NA")

  #fix the dates
  rename!(df, :date=>:dateold)
  #print(names(df))
  df[:date] = (i->Date("$i", dateformat)).(df[:dateold])
  deletecols!(df, :dateold)

  return Universe(df; datebegin=datebegin, dateend=dateend)
end
