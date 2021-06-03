
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
  missingreturns::Vector{BitVector}
end

#creates an investment universe which will be used to generate portfolio returns
function Universe(df::AbstractDataFrame; datebegin::Date = minimum(df[:date]),
  dateend::Date = maximum(df[:date]))::Universe
  subdf::SubDataFrame = view(df,(df[:date] .≥ datebegin) .& (df[:date] .≤ dateend),:)

  N::Int = size(subdf,1) #get the size in number of observations
  M::Int = size(subdf, 2) - 1 #exclude date column

  xcols::Vector{Symbol} = setdiff(names(df), [:date])
  returns::Matrix{MFloat64} = Matrix{MFloat64}(subdf[xcols])'

  #compute the mean for use later
  μ::Vector{Float64} = (r::Int->mean(skipmissing(returns[r,:]))).(1:M)

  #map the missing values
  missingindices::Vector{Vector{Int}} = Vector{Vector{Float64}}()
  missingreturns::Vector{BitVector} = (n->BitVector(undef, M)).(1:N)
  for c ∈ 1:N
    missingreturns[c] = (!ismissing).(returns[:, c])
    push!(missingindices, collect(1:M)[missingreturns[c]])
  end

  return Universe(df,datebegin, dateend, returns, μ, N, M, missingindices, missingreturns)
end
