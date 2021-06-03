
function pdinv(M::AbstractMatrix; warnonfail::Bool=false) where Tv

  @assert size(M,1) == size(M,2)

  try
    return cholesky(M)\I
  catch err
    warnonfail && @warn("pdinv cholesky failed ($err), trying other techniques")
  end

  try
    return svd(M)\I
  catch err
    warnonfail && @warn("pdinv svd failed ($err), trying other techniques")
  end

  try
    return pinv(M)
  catch err
    warnonfail && @warn("pinv failed ($err)")
  end

  error("Inversion failed.")
end

#for whatever reason, it seems only qr decomposition methods work for inverting a pd submatrix
#try qr, then fallback on alternate techniques
function pdinv(M::SubArray; warnonfail::Bool = false)
  @assert size(M,1) == size(M,2)
  try
    return qr(M)\I
  catch err
    warnonfail && @warn("pdinv qr failed ($err), trying other techniques")
  end

  return pdinv(M|>typeof(parent(M)); warnonfail)
end

function pdlsq(M::TM, v::TV, ::Type{T}=eltype(M);
    MtM::TMtM = M'*M, warnonfail::Bool=false) where {TM,TV,T, TMtM}
  b = similar(v, size(M,2))
  try
    b .= cholesky(MtM)\(M'*v)
    return b
  catch err
    warnonfail && @warn("lsq cholesky failed ($err), trying other techniques")
  end

  try
    b .= svd(M)\v
    return b
  catch err
    warnonfail && @warn("lsq svd failed ($err), trying other techniques")
  end

  try
    b .= pdinv(MtM; warnonfail)*(M'*v)
    return b
  catch err
    error("pdinv failed ($err), least squares failed!!!")
  end
end


#TODO create ṽ from ω, define the correct zsections
function updateA₀!(Θ::AbstractVolumePartsMCMC{TM,TV,T};
    mcmcdebug::Bool=true,
    Nthreads::Int = Threads.nthreads()) where {TM<:Matrix,TV<:Vector,T}
  @unpack Ψ, Ã, expand, xsection = Θ
  @unpack σ², A₀, G, ω = Ψ
  @unpack μA₀, VA₀inv = Θ.cc.hyper[:A₀]

  dims = (K=size(Ã,1), T=size(Ã,2))

  #first form a matrix for the cumulative product of G
  Ã .= hcat(ones(T, dims.K), cumprod(G,dims=2))

  #compute the precisions
  τₜ = 1 ./ σ²

  #Λₜ = Vector{TM}(undef, dims.T-1)
  #allocate for the cross-sectional precision matrices and means
  Λₜ::Vector{Matrix{T}} = [Matrix{T}(undef, dims.K, dims.K) for t in 1:(dims.T-1)]
  ΛₜĀₜ::Vector{Vector{T}} = [Vector{T}(undef, dims.K) for t in 1:(dims.T-1)]

  #pre-allocate for better threading
  XtX::Vector{Matrix{T}} = [Matrix{T}(undef, dims.K, dims.K) for i ∈ 1:Nthreads]
  Â::Vector{Vector{T}} = [Vector{T}(undef, dims.K) for i ∈ 1:Nthreads]

  #ΛₜĀₜ = Vector{TV}(undef, dims.T-1)

  #cache the following calculation for repeated use
  VA₀invmat = Matrix(VA₀inv, dims.K, dims.K)
  Vinvμ = VA₀inv*fill(μA₀, dims.K)

  vMω(v,ω) = v - ω
  Threads.@threads for t ∈ 1:(dims.T-1)
    n::Int = Threads.threadid()
    #Xₜ = genabsμₖ.(xsection[:Ã,t]', xsection[:LÃ,t]', xsection[:ws,t], xsection[:RLws,t])
    #XtX[n] = Xₜ' * Xₜ
    #mul!(XtX[n], Xₜ', Xₜ)
    
    X̃ₜ = broadcast!(genabsμₖ,
      xsection[:X̃, t], xsection[:Ã,t]', xsection[:LÃ,t]', xsection[:ws,t], xsection[:RLws,t])


    mul!(XtX[n], X̃ₜ', X̃ₜ)
    #BLAS.gemm!('T', 'N', 1.0, X̃ₜ, X̃ₜ, 0.0, XtX[n])
    ṽ = broadcast!(vMω, xsection[:ṽ,t], xsection[:v,t], xsection[:ω,t])
    #ṽ = xsection[:v,t] .- xsection[:ω,t]
    Â[n] = pdlsq(X̃ₜ, ṽ, MtM=XtX[n], warnonfail=mcmcdebug)
    Λₜ[t] .= τₜ[t] .* (XtX[n] .+ VA₀invmat)
    #Āₜ[t] = pdinv(XtX .+ VA₀invmat, XtX*Â .+ Vinvμ, warnonfail=mcmcdebug)
    #Āₜ = pdinv(XtX .+ VA₀invmat, warnonfail=mcmcdebug) *(XtX*Â .+ Vinvμ)
    #ΛₜĀₜ[t] .= Λₜ[t] * Āₜ
    ΛₜĀₜ[t] .=  τₜ[t] .* (XtX[n]*Â[n] .+ Vinvμ)
    #Ā[t] = (XtX + VA₀inv)\(XtX
  end

  #compute the multivariate parameters of the distribution
  #Λ = sum(τₜ .* Λₜ)
  Λ = sum(Λₜ)
  Σ = pdinv(Λ, warnonfail=mcmcdebug)
  if Σ ≈ Σ' #fixes some roundoff error
    Σ .= (Σ .+ Σ') ./ 2
  else
    printmln(Σ)
    error("something is wrong with Σ: !(Σ ≈ Σ')")
  end

  Ā = Σ * sum(ΛₜĀₜ)

  #now draw conditionally
  #=NOTE: from wikipedia, given vectors x₁ and x₂, x₁|x₂~N(μ,Σ) s.t.
  μ = μ₁+ Σ₁₂ (Σ₂₂)⁻¹ (x₂-μ₂)
  Σ = Σ₁₁-Σ₁₂ (Σ₂₂)⁻¹ Σ₂₁
  "This matrix is the Schur complement of Σ22 in Σ. This means that to calculate the conditional
    covariance matrix, one inverts the overall covariance matrix, drops the rows and columns
    corresponding to the variables being conditioned upon, and then inverts back to get the
    conditional covariance matrix." - Wikipedia
  =#

  #WARNING- uncomment below
  for k ∈ shuffle(1:dims.K)
    notk = 1:dims.K .≠ k
    σ²k = Λ[k,k]^-1
    #ΛMkk = pdinv(@view(Σ[notk, notk]), warnonfail=mcmcdebug)
    #ΛMkk = pdinv(Σ[notk, notk], warnonfail=mcmcdebug)

    ΛMkk = cholesky(Σ[notk, notk])\I
    μk = Ā[k] + (Σ[k:k,notk]*ΛMkk*(A₀[notk] .- Ā[notk]))[] #unwraps a 1x1 array
    A₀[k] = rand(Normal(μk, σ²k))
  end

  return A₀

  #=WARNING- below is temporary
  if Σ ≈ Σ' #fixes some roundoff error
    Σ .= (Σ .+ Σ') ./ 2
  else
    printmln(Σ)
    error("something is wrong with Σ: !(Σ ≈ Σ')")
  end

  #???  maybe this is the issue Ā = Σ*ΛĀ
  #@info sum(abs.(Σ .- Σ'))
  d = MvNormal(Ā, Σ)
  return rand(d)=#
end
