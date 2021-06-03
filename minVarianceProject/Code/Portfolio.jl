#primary structure for portfolio math
#some auxillery structures to improve missing value efficiency
struct Portfolio
  U::Universe
  wₜ::Vector{Float64}
  R::Vector{Float64}
  pμ::Base.RefValue{Float64}

  #the below vector contains all valid weights, and does not include invalid weights
  #thus the shape of wvec is a ragged array, with each array of length equal to the
  #number of valid values
  wvec::Vector{Vector{Float64}} #store this so we don't have to reallocate
end

#update the portfolio object
#innermost loop and likely a performance bottleneck, so highly optimized
function Portfolio!(P::Portfolio, wₜ::AbstractVector; normalizeto::Float64 = 1.0,
  checkadjustment::Bool = true)::Portfolio

  P.wₜ .= wₜ

  for c ∈ 1:P.U.N
    for r ∈ P.U.missingindices[c]
      P.wvec[c][r] = wₜ[r]
    end

    adjfactor::Float64 = sum(P.wvec[c]) / normalizeto
    P.wvec[c] ./= adjfactor
    checkadjustment && (adjfactor < .05) && @warn("Extreme adjustment factor $adjfactor -> scaling =
      $(1/adjfactor). Expect unstable results. Consider shifting the test distribution more
      positive or increasing data completeness.")

  end

  #build the returns (Note: Threads.@threads provides signficant performance gains)
  P.R .= 0.0
  for c ∈ 1:P.U.N
    for r ∈ P.U.missingindices[c]
      P.R[c] += P.wvec[c][r] * P.U.returns[r,c]
    end
  end

  #update the pointer
  P.pμ[] = mean(P.R)
  return P
end

#construct a brand new portfolio object
function Portfolio(U::Universe, wₜ::AbstractVector)::Portfolio
  local wvec::Vector{Vector{Float64}} = (n->zeros(U.M)).(1:U.N)

  P::Portfolio = Portfolio( #pre-allocate and call the update function
    U,
    Vector{Float64}(undef, U.M),
    Vector{Float64}(undef, U.N),
    Ref(-1.0),
    wvec)

  return Portfolio!(P, wₜ)
end

#helper method which creates an equi-weighted portoflio
Portfolio(U::Universe)::Portfolio = Portfolio(U, fill(1.0/U.M, U.M))

#creates a new portfolio from another one, re-using some of the references
function Portfolio(sourceP::Portfolio, wₜ::AbstractVector)

  P::Portfolio = Portfolio(sourceP.U,
    Vector{Float64}(undef, sourceP.U.M),
    Vector{Float64}(undef, sourceP.U.N),
    Ref(-1.0),
    sourceP.wvec)


  Portfolio!(P, wₜ)
  return P
end


#draws a random test portfolio using the universe from G
function redrawtestportfolio!(P::Portfolio;numtestassets::Int = 100,
  randomweightdist::T = RANDOM_WEIGHT_DIST,
  weightoftestasset::Function = i::Int -> rand(randomweightdist), #easy to modify if needed
  minadjustment::Float64 = MIN_ADJUSTMENT
  )::Portfolio where T<:Union{Distribution, Vector}

  wₜ::Vector{Float64} = zeros(P.U.M)
  for i ∈ 1:numtestassets
    wₜ[i] = weightoftestasset(i)
  end


  tot::Float64 = sum(wₜ)

  #Don't want extreme weights or a flipped distribution (numerical instability)
  #Value of 0.75 effectively caps the max weight to 200%
  if tot ≤ minadjustment
    return redrawtestportfolio!(P)
  end
  wₜ ./= tot

  shuffle!(wₜ)
  try
    return Portfolio!(P, wₜ, checkadjustment=false)
  catch
    return redrawtestportfolio!(P)
  end
end

#initialize the test portfolio from G
function initializetestportfolio(G::Portfolio)::Portfolio
  P::Portfolio = Portfolio(G, G.wₜ)

  return redrawtestportfolio!(P)
end

function Statistics.cov(P1::Portfolio, P2::Portfolio)::Float64

  return dot(P1.R, P2.R) / P1.U.N - P1.pμ[]*P2.pμ[]
end

function Statistics.var(P1::Portfolio)::Float64
  return cov(P1,P1)
end

function cov!(portfolios::Vector{Portfolio},
  S::AbstractMatrix = Matrix{Float64}(undef, length(portfolios), length(portfolios)))

  #compute the covariance matrix
  for c ∈ 1:length(portfolios)
    for r ∈ 1:c
      S[r,c] = cov(portfolios[r], portfolios[c])
      if r ≠ c
        S[c,r] = S[r,c] #fill in the lower left triangle
      end
    end
  end

  return S
end
