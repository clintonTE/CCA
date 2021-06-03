

function PanelNeweyWest(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}

  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]
  K::Int = length(λᵢ)
  R::Matrix{Float64} = qr(X).R
  Rinv::Matrix{Float64} = R\I #this will save some computational time
  Sₜ::Matrix{Float64} = zeros(K,K) #holds teh central matrix

  #pre-allocate for the spectral matrix
  Σ::Matrix{Float64} = Matrix{Float64}(undef, K, K)

  temp::Matrix{Float64} = Matrix{Float64}(undef, K, K) #pre-allocate working matrix
  RRinv::Matrix{Float64} = BLAS.gemm('N', 'T', Rinv, Rinv) #this is equivelent to [X'X]^-1

  #iterate over all stocks
  for sdf ∈ fs.eachport
    #need to multiply through by the error
    T::Int = size(sdf,1)
    Xₙ::Matrix{Float64} = [ones(T) Matrix{Float64}(sdf[!, coefnames[2:end]])]
    ε::Vector{Float64} = sdf[!, fs.Fret] .- Xₙ * λᵢ
    Xₑ::Matrix{Float64} = Xₙ .* ε

    Sₜ += BLAS.gemm('T','N',1.0/N, Xₑ, Xₑ)

    for v::Int ∈ 1:lags
      #overwrites temp with (1/N)R'R
      BLAS.gemm!('T', 'N', 1.0/N, view(Xₑ, (v+1):T, :),view(Xₑ, 1:(T-v), :), 0.0, temp)
      Sₜ .+= (lags + 1 - v)/(lags+1.) .* (temp .+ temp')
    end

  end
  #Sₜ = Matrix{Float64}(I,K,K)
  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRinv, Sₜ)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(N), RRinvS, RRinv, 0.0, Σ)
  #println("$(diag(Σ .* dofCorrect))")sasa
  return Σ
end

#this version provides an intuitive check
#check #1
function PanelNeweyWestSlow(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}
  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]

  #pre-allocate
  K::Int = size(X,2)
  Sₜ::Matrix{Float64} = zeros(K, K)
  #xxt::Matrix{Float64} = Matrix{Float64}(undef, K,K)
  temp::Matrix{Float64} = Matrix{Float64}(undef, K,K)

  #helper function using the Lars method
  @inline function Rₜ!(v::Int, X::Matrix{Float64}, ε::Vector{Float64};
      T::Int = size(X,1), out::Matrix{Float64} = Matrix{Float64}(undef, K,K))

    out .= 0.0
    for t ∈ (1+v):T
      out .+= X[t,:] * X[t-v,:]' * ε[t] * ε[t-v]
    end

    #out ./= T
    return out
  end

  for sdf ∈ fs.eachport
    T::Int = size(sdf, 1)

    Xₙ::Matrix{Float64} = [ones(T) Matrix{Float64}(sdf[!, coefnames[2:end]])]
    ε::Vector{Float64} = sdf[!, fs.Fret] .- Xₙ * λᵢ
    Sₜ .+= Rₜ!(0, Xₙ, ε, out=temp)


    for v ∈ 1:lags #iterate over lags
      Rₜ!(v, Xₙ, ε, out=temp)
      Sₜ .+= (lags + 1. - v)/(lags + 1.) .* (temp .+ temp')
    end
  end

  Sₜ ./= N

  XXInv::Matrix{Float64} = (X' * X) \ I


  return N * XXInv * Sₜ * XXInv
end

#alg from Correcting for Both Cross-Sectional and
#  Time-Series Dependence in Accounting Research by Ian Gow, Gaizka Ormazabal and Daniel Taylor.
# check #2
function PanelNeweyWestSlow2(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}
  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]
  ε::Vector{Float64} = fs.df[!, fs.Fret] .- X * λᵢ
  ports::typeof(fs.df[!, fs.Fport]) = fs.df[!, fs.Fport]
  #pre-allocate
  K::Int = size(X,2)
  Sₜ::Matrix{Float64} = zeros(K, K)

  #helper function using the Lars method

  for t ∈ 1:N
    Sₜ .+= X[t,:] * X[t,:]' * ε[t]^2
  end

  for v::Int ∈ 1:lags
    for t ∈ (v+1):N
      (ports[t] == ports[t-v]) && (
        Sₜ .+= (lags + 1 - v)/(lags+1.) .* (X[t,:] * X[t-v,:]' .+ X[t-v,:] * X[t,:]') * ε[t] * ε[t-v])
    end
  end
  Sₜ ./= N
  #Sₜ = Matrix{Float64}(I,K,K)

  #the final step is to multiply out the var-covar sandwhich
  Σ::Matrix{Float64} = Matrix{Float64}(undef, K, K)
  R::Matrix{Float64} = qr(X).R
  Rinv::Matrix{Float64} = (R)\Matrix{Float64}(I,K,K) #this will save some computational time
  RRinv::Matrix{Float64} = BLAS.gemm('N', 'T', Rinv, Rinv) #this is equivelent to [X'X]^-1

  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRinv, Sₜ)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(N), RRinvS, RRinv, 0.0, Σ)

  return Σ
end

function FMBartlett(λ::Matrix{Float64},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}

  T::Int = size(λ,1)
  K::Int = length(λᵢ)

  u::Matrix{Float64} = λ .- λᵢ'

  #pre-allocate for the spectral matrix
  Sₜ::Matrix{Float64} = zeros(K,K) #holds teh central matrix
  temp::Matrix{Float64} = Matrix{Float64}(undef, K, K) #pre-allocate working matrix
  Sₜ::Matrix{Float64} += BLAS.gemm('T','N',1.0, u, u)

  for v::Int ∈ 1:lags
    #overwrites temp with (1/N)R'R
    BLAS.gemm!('T', 'N', 1.0, view(u, (v+1):T, :),view(u, 1:(T-v), :), 0.0, temp)
    Sₜ .+= (lags + 1 - v)/(lags+1.) .* (temp .+ temp')
  end

  Sₜ ./= T

  return Sₜ
end

#this version provides an intuitive check
#check #1
function FMBartlettSlow(λ::Matrix{Float64},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}

  #based on Cochrane chap 11/12- Bartlett kernel
  u::Matrix{Float64} = λ .- λᵢ'
  T::Int = size(u,1)

  #pre-allocate
  K::Int = size(u,2)
  Sₜ::Matrix{Float64} = zeros(K, K)
  #xxt::Matrix{Float64} = Matrix{Float64}(undef, K,K)
  temp::Matrix{Float64} = Matrix{Float64}(undef, K,K)

  #helper function using the Lars method
  @inline function Rₜ!(v::Int, u::Matrix{Float64},
      out::Matrix{Float64} = Matrix{Float64}(undef, K,K))

    out .= 0.0
    for t ∈ (1+v):T
      out .+= u[t,:] * u[t-v,:]'
    end

    return out
  end

  Sₜ .+= Rₜ!(0, u, temp)


  for v ∈ 1:lags #iterate over lags
    Rₜ!(v, u, temp)
    Sₜ .+= (lags + 1. - v)/(lags + 1.) .* (temp .+ temp')
  end

  Sₜ ./= T

  return Sₜ
end

#alg from Correcting for Both Cross-Sectional and
#  Time-Series Dependence in Accounting Research by Ian Gow, Gaizka Ormazabal and Daniel Taylor.
# check #2
function FMBartlettML(λ::Matrix{Float64},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}
  T::Int = size(λ,1)

  u::Matrix{Float64} = λ .- λᵢ'

  #pre-allocate
  K::Int = length(λᵢ)
  Sₜ::Matrix{Float64} = zeros(K, K)
  uu::Matrix{Float64} = Matrix{Float64}(undef, K, K)

  #helper function using the Lars method

  for t ∈ 1:T
    Sₜ .+= u[t,:] * u[t,:]'
  end

  for v::Int ∈ 1:lags
    for t ∈ (v+1):T
        #Sₜ .+= (lags + 1 - v)/(lags+1.) .* (u[t,:] * u[t-v,:]' .+ u[t-v,:] * u[t,:]'))
        uu .= u[t,:] * u[t-v,:]'
        Sₜ .+= (lags + 1 - v)/(lags+1.) .* (uu .+ uu')
    end
  end
  Sₜ ./= T
  #Sₜ = Matrix{Float64}(I,K,K)

  return Sₜ
end
