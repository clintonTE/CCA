




function testuniverse(;N=560, M=5000)::Universe
  local cols::Vector{Symbol} = (n->Symbol("X",n)).(1:M)
  local crosscor::MFloat64 = 0.25
  local stdev::MFloat64 = 0.1
  local allmean = .02

  ρ::Matrix{MFloat64} = ones(M,M) .* crosscor
  ρ[diagind(ρ)] .= 1.0 #fix the diagonals
  σ::Vector{MFloat64} = (n->stdev).(1:M)
  Σ::Matrix{MFloat64} = ρ .* stdev^2
  L = cholesky(Σ).L

  #make the test matrix
  local testmatrix::Matrix{MFloat64} = rand(Normal(), N, M)

  for r ∈ 1:N
    testmatrix[r,:] .= L * testmatrix[r,:] .+ allmean
  end

  #make missing values
  for c ∈ 1:M
    for r ∈ 1:N
      (rand(1:10) == 1) && (testmatrix[r,c] = missing)
    end
  end

  #make the dataframe
  df = DataFrame(testmatrix)
  names!(df, cols)
  df[:date] = Vector{MDate}(missing, N)

  local datebegin::Date = Date(2011,11,11)

  df[:date] .= datebegin .+ (Day).(1:N)

  U = Universe(df)

  return U
end

#testing
#simulates an artificial portfolio and reweights it iter times
function testportfolio(;N=560, M=5000, iter=1,
  U::Universe = testuniverse(N=N, M=M))::Portfolio

  wₜ = ones(M) ./ M

  P = Portfolio(U,wₜ)

  local μs::Vector{Float64} = Vector{Float64}(undef, iter)


  for i ∈ 1:iter
    wₜ .= rand(Normal(2,1), M)
    wₜ ./= sum(wₜ)
    Portfolio!(P, wₜ, checkadjustment=true)
    μs[i] = P.pμ[]
  end

  iter>1 && println("Sum of means: $(sum(μs))")

  return P
end

function testcov(;P::Portfolio=testportfolio(iter=1), iter = 1000)::Nothing
  local P1::Portfolio = deepcopy(P)
  local P2::Portfolio = deepcopy(P)
  local N::Int = P.U.N

  covs::Vector{Float64} = Vector{Float64}(undef, iter)
  @time for i::Int ∈ 1:iter
    P1.R .= rand(N)
    P2.R .= rand(N)

    P1.pμ[] = mean(P1.R)
    P2.pμ[] = mean(P2.R)

    covs[i] = cov(P1,P2)
    var(P1)
  end

  print("Sum of covs: $(sum(covs))")

  return nothing
end

function testtrisectedportfolio(;P::Portfolio=testportfolio(iter=1),
  iter = 10, increment = 2, displayiter = 3)

  TP = TrisectedPortfolio(P)
  if displayiter > 0
    display(TP.assignments)
    display(TP.wₛ)
    display(TP.Σ)
  end

  @time for i ∈ 1:iter
    rotatetrisected!(TP, increment)
    if i ≤ displayiter
      display(TP.assignments)
      display(TP.wₛ)
      display(TP.Σ)
    end
  end

  return nothing
end

function testω(;iter::Int = 1)

  errSG::Vector{Float64} = Vector{Float64}(undef, iter)
  errSGP::Vector{Float64} = Vector{Float64}(undef, iter)
  errSGPratio::Vector{Float64} = Vector{Float64}(undef, iter)
  errsum::Vector{Float64} = Vector{Float64}(undef, iter)

  U::Universe=testuniverse()
  G::Portfolio = testportfolio(iter=1, U=U)
  P::Portfolio = initializetestportfolio(G)
  TG::TrisectedPortfolio = TrisectedPortfolio(G)

  P = initializetestportfolio(G)




  @time for i ∈ 1:iter

    initializetestportfolio(G)
    redrawtestportfolio!(P)
    refreshtrisected!(TG)
    #print("var(TG.P: $(var(TG.P)), 1'S1: $(ones(3)' * TG.S * ones(3))\n\n")

    Θ::Parameters = Parameters(
      σ²G = D_σ²G,
      ζ²G = D_ζ²G,
      ζ²P = D_ζ²P,
      SG = var(G),
      SGP = cov(G,P)
    )

    oldS::Matrix{Float64} =  deepcopy(TG.S) #need thsi for checks later

    SGP::Vector{Float64} = cov(TG, P)
    ω::Vector{Float64} = getω(Θ, TG, SGP)
    rescaletrisected!(TG, ω)

    errSG[i] = abs( ω' * oldS * ω - var(TG.P))
    errSGP[i] = abs( ω' * SGP - cov(TG.P,P))
    errSGPratio[i] = abs( cov(TG.P,P) / (ω' * SGP))
    errsum[i] = abs( 1.0 - sum(TG.P.wₜ))

    (iter==1) && print("SG: $(Θ.SG)\n")
    (iter==1) && print("= ω' S ω: $(ω' * oldS * ω)\n")
    (iter==1) && print("= var(TG.P): $(var(TG.P))\n\n")


    (iter==1) && print("SGP: $(Θ.SGP)\n")
    (iter==1) && print("= ω' SGP: $(ω' * SGP)\n")
    (iter==1) && print("= cov(TG.P,P): $(cov(TG.P,P))\n\n")

    (iter==1) && print("1.0 \n")
    (iter==1) && print("= sum(TG.P.wₜ): $(sum(TG.P.wₜ))\n\n")

  end

  if (iter > 1) && (iter<100)
    print("Error SG - w'Sw: $errSG\n")
    print("Error SGP - w'SGP: $errSGP\n")
    print("Error sum 1'wt: $errsum\n")
  elseif iter ≥ 100
    print("Error SG - w'Sw max: $(maximum(errSG)) mean: $(mean(errSG))\n")
    print("Error SGP - w'SGP max: $(maximum(errSGP)) mean: $(mean(errSGP))\n")
    print("Error sum 1'wt max: $(maximum(errsum)) mean: $(mean(errsum))\n")
  end


end

function testposteriors()
  Π::Priors = Priors()
  Θ1::Parameters = Parameters()
  Θ2::Parameters = Parameters()
  F::Posteriors = Posteriors()

  Random.seed!(11)
  F.σ²G!(Θ1, Π)
  F.ζ²G!(Θ1, Π)
  F.ζ²P!(Θ1, Π)
  F.SG!(Θ1, Π)
  F.SGP!(Θ1, Π)

  F = PosteriorsTest()
  Random.seed!(11)
  F.σ²G!(Θ2, Π)
  F.ζ²G!(Θ2, Π)
  F.ζ²P!(Θ2, Π)
  F.SG!(Θ2, Π)
  F.SGP!(Θ2, Π)

  println("Check Θ1: $(println(Θ1))")
  println("Check Θ2: $(println(Θ2))")


end
