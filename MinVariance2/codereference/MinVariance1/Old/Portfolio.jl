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
  wmat::Matrix{Float64} #stores the weight vectors
end

#update the portfolio object
#innermost loop and likely a performance bottleneck, so highly optimized
function Portfolio!(P, wₜ)::Portfolio

  P.wₜ .= wₜ
  P.wmat .= P.U.missingreturns .* wₜ
  P.wmat ./= sum(P.wmat, dims=1)

  #build the returns (Note: Threads.@threads provides signficant performance gains)
  P.R .= 0.0
  for c ∈ 1:P.U.N
    for r ∈ P.U.missingindices[c]
      P.R[c] += P.wmat[r,c] * P.U.returns[r,c]
    end
  end

  #update the pointer
  P.pμ[] = mean(P.R)
  #print(P.R)
  return P
end

#construct a brand new portfolio object
function Portfolio(U, wₜ)::Portfolio
  #wvecsize::Vector{Int} = (v->length(v)).(U.missingmap)
  #wvec::Vector{Vector{Float64}} = (n->Vector{Float64}(undef, n)).(wvecsize)
  wmat::Matrix{Float64} = Matrix{Float64}(undef, U.M, U.N)

  P::Portfolio = Portfolio( #pre-allocate if needed, reduces code reuse
    U,
    Vector{Float64}(undef, U.M),
    Vector{Float64}(undef, U.N),
    Ref(-1.0),
    #Matrix{Float64}(undef, U.M, U.N))
    wmat)

  return Portfolio!(P, wₜ)
end
