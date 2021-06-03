
#testing
#simulates an artificial portfolio and reweights it iter times
function testPortfolio(;N=560, M=5000, iter=100)
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
