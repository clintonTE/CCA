using Revise
using DataFrames, Distributions, LinearAlgebra, Finometrics, Dates

struct Universe
  df::DataFrame #source dataframe must be properly formatted with a :date column but can include extra observations
  datebegin::Date
  dateend::Date
  returns::Matrix{MFloat64} #should contain the numerical returns data in the DataFrame
  #NOTE for performance reasons, this is M by N!!!

  μ::Vector{Float64} #we can calc this in advance
  N::Int #number of observations in the sample period
  M::Int #number of return series

  missingmap::Vector{Vector{Int}} # a collection of all non-missing indices in returns
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
  missingmap::Vector{Vector{Int}} = Vector{Vector{Float64}}()
  for c ∈ 1:N
    push!(missingmap, collect(1:M)[(!ismissing).(returns[:, c])])
  end

  return Universe(df,datebegin, dateend, returns, μ, N, M, missingmap)
end

#lightweight structure for portfolio math
struct Portfolio
  U::Universe
  wₜ::Vector{Float64}
  R::Vector{Float64}
  pμ::Base.RefValue{Float64}

  wmat::Matrix{Float64} #store this so we don't have to reallocate
end

#update the portfolio object
function Portfolio!(P, wₜ)::Portfolio

  P.wₜ .= wₜ
  for c ∈ 1:P.U.N
    for r ∈ P.U.missingmap[c] #iterate over existing values, skipping missing values
      P.wmat[r,c] = P.wₜ[r]
    end

    adjfactor::Float64 = sum(P.wmat[P.U.missingmap[c],c])
    P.wmat[P.U.missingmap[c],c] ./= adjfactor #rescale to account for missing values
  end

  #build the returns (Note: Threads.@threads provides signficant performance gains)
  for c ∈ 1:P.U.N
    P.R[c] = sum(P.wmat[P.U.missingmap[c],c] .* P.U.returns[P.U.missingmap[c],c])
  end

  #update the pointer
  P.pμ[] = mean(P.R)
  return P
end

#construct a brand new portfolio object
function Portfolio(U, wₜ)::Portfolio
  P::Portfolio = Portfolio( #pre-allocate if needed, reduces code reuse
    U,
    Vector{Float64}(undef, U.M),
    Vector{Float64}(undef, U.N),
    Ref(-1.0),
    Matrix{Float64}(undef, U.M, U.N))

  return Portfolio!(P, wₜ)
end



function testFMCov(N=560, M=5000, iter::Int=20)
  cols::Vector{Symbol} = (n->Symbol("X",n)).(1:M)
  crosscor::MFloat64 = 0.25
  stdev::MFloat64 = 0.1
  allmean = .02

  ρ::Matrix{MFloat64} = ones(M,M) .* crosscor
  ρ[diagind(ρ)] .= 1.0 #fix the diagonals
  σ::Vector{MFloat64} = (n->stdev).(1:M)
  Σ::Matrix{MFloat64} = ρ .* stdev^2
  L = cholesky(Σ).L

  #make the test matrix
  testmatrix::Matrix{MFloat64} = rand(Normal(), N, M)

  for r ∈ 1:N
    testmatrix[r,:] .= L * testmatrix[r,:] .+ allmean
  end

  #make missing values
  for c ∈ 1:M
    for r ∈ 1:N
      rand([true, false]) && (testmatrix[r,c] = missing)
    end
  end

  #make the dataframe
  df = DataFrame(testmatrix)
  names!(df, cols)
  df[:date] = Vector{MDate}(missing, N)

  datebegin::Date = Date(2011,11,11)

  df[:date] .= datebegin .+ (Day).(1:N)

  U = Universe(df)
  wₜ = ones(M) ./ M

  P = Portfolio(U,wₜ)

  @time for i ∈ 1:iter
    wₜ .= rand(M)
    Portfolio!(P, wₜ)
  end

  return nothing
end

testFMCov()
