

#holds the time-indexed parameters and a pre-allocation for the row-indexed versions
#note- we ahve to be careful about using methods from AbstractVolumeParts since we ahve no V₀

abstract type AbstractVolumePartsMCMC{TM, TV, elTΠ} <: AbstractVolumeParts{TM, TV, elTΠ} end

struct VolumePartsMCMC{
    TΨ<:AbstractVector,
    Tcc<:CapacityChains,
    TXv<:AbstractXYIter,
    TA<:AbstractMatrix,
    #Tv<:AbstractVector,
    Tts<:AbstractVector,
    Texpand<:Dict,
    Txsection<:Dict,
    TM<:AbstractMatrix,
    TV<:AbstractVector,
    elTΠ<:Real} <: AbstractVolumePartsMCMC{TM, TV, elTΠ}

  Ψ::TΨ
  cc::Tcc
  Xv::TXv

  Ã::TA

  ts::Tts
  tsM1::Tts
  expand::Texpand
  xsection::Txsection

  #inner constructor needed to identify elTΠ
  function VolumePartsMCMC(
    Ψ::TΨ, cc::Tcc, Xv::TXv, Ã::TA,
    ts::Tts, tsM1::Tts, expand::Texpand, xsection::Txsection,
    ::Type{TM}, ::Type{TV}, ::Type{T}) where {
      TΨ, Tcc, TXv, TA, Tts, Texpand, Txsection, TM, TV, T}

    (eltype(TΨ) == T) || error("Eltype of Ψ must match T in TM,TV,T")
      #replace the below
    #@assert (Π == reduce(vcat, [vec(A), vec(Ã), vec(G)]))

    return new{TΨ, Tcc, TXv, TA, Tts, Texpand, Txsection,TM, TV, T}(
      Ψ, cc, Xv, Ã, ts, tsM1, expand, xsection)
  end
end



#pre-allocates for the vpexpand
function VolumePartsMCMC(Ψ::CCState, cc::CapacityChains,  Xv::AbstractXYIter,
  ts::Tts, ::Type{TM}, ::Type{TV}, ::Type{T};
  tsM1::Tts=ts .- 1) where {TM, TV, T, Tts<:AbstractVector}

  dims=(K=size(Ψ.G,1), T=size(Ψ.G,2)+1, N=size(Xv.ws,1))

  Ã = rand(T, dims.K, dims.T)

  expand = Dict{Union{Symbol, Tuple{Symbol, Int}},AbstractArray{T}}()
  xsection = Dict{Union{Tuple{Symbol, Int}},AbstractArray{T}}()

  for t ∈ 1:(dims.T-1)

    #these expansions include all values equal or greaer than t
    tpart = ts[ts .≥ t+1]
    tpartM1 = tsM1[tsM1 .≥ t]
    expand[:ws,t] = view(Xv.ws, ts .≥ t+1, :)
    expand[:RLws,t] = view(Xv.RLws, ts .≥ t+1, :)
    expand[:v,t] = view(Xv.v, ts .≥ t+1)

    expand[:G,t] = view(Ψ.G, :, tpartM1)
    expand[:σ²,t] = view(Ψ.σ², tpartM1)
    expand[:ω,t] = view(Ψ.ω, tpartM1)

    expand[:Ã,t] = view(Ã, :, tpart)
    expand[:LÃ,t] = view(Ã, :, tpartM1)

    #these expansions are simple xsections for a particular t
    txsection = ts[ts .== t+1]
    txsectionM1 = tsM1[tsM1 .== t]
    xsection[:ws,t] = view(Xv.ws, ts .== t+1, :)
    xsection[:RLws,t] = view(Xv.RLws, ts .== t+1, :)
    xsection[:v,t] = view(Xv.v, ts .== t+1)
    xsection[:X̃,t] = TM(undef, sum(ts .== t+1), dims.K) #cross-sectional allocations
    xsection[:ṽ,t] = TV(undef, sum(ts .== t+1))

    xsection[:G,t] = view(Ψ.G, :, txsectionM1)
    xsection[:σ²,t] = view(Ψ.σ², txsectionM1)
    xsection[:ω,t] = view(Ψ.ω, txsectionM1)
    xsection[:Ã,t] = view(Ã, :, txsection)
    xsection[:LÃ,t] = view(Ã, :, txsectionM1)
  end

  #convenience alias for the complete case
  expand[:ws] = expand[:ws,1]
  expand[:RLws] = expand[:RLws,1]
  expand[:v] = expand[:v,1]
  #expand[:X̃] = expand[:X̃,1]
  #expand[:ṽ] = expand[:ṽ,1]

  expand[:G] = expand[:G,1]
  expand[:σ²] = expand[:σ²,1]
  expand[:ω] = expand[:ω,1]
  expand[:Ã] = expand[:Ã,1]
  expand[:LÃ] = expand[:LÃ,1]


  Θ = VolumePartsMCMC(Ψ, cc, Xv, Ã, ts, tsM1, expand, xsection, TM, TV, T)
  return Θ
end


#basic constructor from dimension arguments
function VolumePartsMCMC(
  ::Type{TM}=PARAM[:mcmcgpu] ? CUDA.CuMatrix{T} : Matrix{T},
  ::Type{TV}=PARAM[:mcmcgpu] ? CUDA.CuVector{T} : Vector{T},
  ::Type{T}=PARAM[:mcmctype];
    cc::CapacityChains=error("cc is required"),
    Xv::AbstractXYIter=error("Xv is required"),
    ts::Vector{Int}=error("ts is required"),
    initial::Union{NamedTuple, Dict}=error("initial ")
    ) where {T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  Ψ = CCState(cc, initial)
  #=params=(
    A₀ = TV(rand(T, dims.K)),
    G = TM(ones(T, dims.K, dims.T-1)),
    σ² = rand(T, dims.T-1),
    ω = rand(T, dims.T-1)
    )=#


  return VolumePartsMCMC(Ψ, cc, Xv, ts, TM, TV, T)
end

#deepcopy avoids copying the static objects
function Base.deepcopy(Θ::AbstractVolumePartsMCMC{TM, TV, T}) where {TM,TV,T}

  Ψ=deepcopy(Θ.Ψ)
  cc=deepcopy(Θ.cc)
  Xv=Θ.Xv
  ts=Θ.ts
  tsM1=Θ.tsM1

  return VolumePartsMCMC(Ψ, cc, Xv, ts, TM, TV, T, tsM1=tsM1)
end

#=
#the below funciton is now only for testing purposes
function (Θ::VolumePartsIter)(Xv::XYIterLevel, _Xv::AbstractXYIterAlloc{TM,TV,T},
  ::Val{:level}, RegressionType::Val, absμ::Tabsμ = _Xv.absμ_alloc
  ) where {TM<:Matrix,TV<:Vector,T,Tabsμ<:Vector}
  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = initializeÃandfactor(Θ,Xv,_Xv)

  #run the regression for A₁
  try
    A[:,1] .= reg(factored, Xv.v, RegressionType)
  catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage > 1000) ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || error("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])")

    A[:,1] .= reg(factored, Xv.v, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallbackreg fallback val=$(A[:,1])"
  end

  #println("A[:,1]: $(A[:,1])")

  #absμ .= vec(sum(absμₖ,dims=2))
  absμ .= factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end


#if absμₖ is already formed, generating the demeaned predictions should be cheap
function (Θ::VolumePartsIter)(Xv::AbstractXYIter{TM,TV,T}, _Xv::AbstractXYIterAlloc{TM,TV,T},
  absμₖ::Tabsμₖ, ::Val{:leveldemean},
  ::Type{Tabsμ}=typeof(_Xv.absμ_alloc),
  ::Type{elTabsμₖ}=eltype(Tabsμₖ)) where {TM,TV,T,Tabsμₖ,Tabsμ,elTabsμₖ}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))
  absμ::Tabsμ = _Xv.absμ_alloc
  LinearAlgebra.BLAS.gemv!('N',elTabsμₖ(1.0),absμₖ,ones(elTabsμₖ,dims.K),zero(elTabsμₖ),absμ)

  return absμ
end

#primary prediction alogrithm for the demeaned version
function (Θ::VolumePartsIter)(Xv::XYIterDemean, _Xv::AbstractXYIterAlloc{TM,TV,T},
    ::Val{:leveldemean}, RegressionType::Val, absμ::Tabsμ = _Xv.absμ_alloc
    ) where {TM,TV,T,Tabsμ}
  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = initializeÃandfactor(Θ,Xv,_Xv)

  #run the regression for A₁
  try
    A[:,1] .= reg(factored, Xv.ṽ, RegressionType)
  catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage > 1000) ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || error("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])")

    A[:,1] .= reg(factored, Xv.ṽ, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallback val=$(A[:,1])"
  end

  #println("A[:,1]: $(A[:,1])")

  #absμₖ = foldl(Completing(vcat), Map(M->vec(sum(M,1))), xabsμₖ)
  #absμ .= vec(sum(absμₖ,dims=2))
  absμ .= factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end

#this is a cheap version which relies on the already calculated absμₖ
(Θ::VolumePartsIter)(Xv, _Xv, absμₖ, ::Val{:leveldemean}) = (Θ::VolumePartsIter)(
  Xv, _Xv, absμₖ, Val{:level}())



(Θ::VolumePartsIter)(Xv::AbstractXYIter{TM,TV,T},
  _Xv::AbstractXYIterAlloc{TM,TV,T},
    ::Val{:leveldemean}) where {TM<:CuMatrix,TV<:CuVector,T} = error(
    "Do not use CUDA for the demeaning version")


function (Θ::VolumePartsIter)(Xv::XYIterLevel,
  _Xv::AbstractXYIterAlloc{TM,TV,T}, ::Val{:level}, RegressionType::Val,
  absμ::Tabsμ = _Xv.absμ_alloc
  ) where {TM<:CuMatrix,TV<:CuVector,T,Tabsμ<:CuArray}

  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = initializeÃandfactor(Θ,Xv,_Xv)

  #run the regression for A₁
  try
    A[:,1] .= reg(factored, Xv.v, RegressionType)
    #@info "culevel A[:,1]: $(A[:,1])"
    #@info "culevel factored[20,:]: $(factored[20,:])"
    #@info "culevel A[:,20]: $(A[:,20])"
  catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage > 1000) ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || error("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])")

    A[:,1] .= reg(factored, Xv.v, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallback val=$(A[:,1])"
  end

  #println("A[:,1]: $(A[:,1])")

  absμ .= factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end

#this is a cheap version which relies on the already calculated absμₖ
function (Θ::VolumePartsIter)(Xv::XYIterLevel,
    _Xv::AbstractXYIterAlloc{TM,TV,T},
    absμₖ::Tabsμₖ, ::Type{elTabsμₖ}=eltype(Tabsμₖ)
    ) where {TM<:CuMatrix,TV<:CuVector,T, Tabsμₖ<:CuMatrix, elTabsμₖ}

  #first map each coefficient to its appropriate
  dims = (K=size(Θ.A,1), T = size(Θ.A,2))
  absμ::CuVector{elTabsμₖ} = _Xv.absμ_alloc
  CUBLAS.gemv!('N',elTabsμₖ(1.0),absμₖ,CUDA.ones(elTabsμₖ,dims.K),zero(elTabsμₖ),absμ)

  return absμ
end


#computs the predicted values for each strategy k and time t and demeans the result
#via the projection matrix
function projectvolumefromA(Θ::VolumePartsIter, Xv::XYIterDemean,
  _Xv::AbstractXYIterAlloc{TM,TV,T}, absμₖ::Tabsμₖ = _Xv.absμₖ_alloc) where {TM,TV,T,Tabsμₖ}

  @unpack A, LA = Θ.expand
  @unpack ws, RLws, M = Xv
  #first map each coefficient to its appropriate

  absμₖ::Tabsμₖ = _Xv.absμₖ_alloc
  absμₖ .= genabsμₖ.(A', LA', ws, RLws)

  multiply!(M, absμₖ)
  return absμₖ
end

#same as the above without demeaning
function projectvolumefromA(Θ::VolumePartsIter, Xv::XYIterLevel,
  _Xv::AbstractXYIterAlloc{TM,TV,T}, absμₖ::Tabsμₖ = _Xv.absμₖ_alloc) where {TM,TV,T,Tabsμₖ}

  @unpack A, LA = Θ.expand
  @unpack ws, RLws= Xv
  #first map each coefficient to its appropriate

  absμₖ::Tabsμₖ = _Xv.absμₖ_alloc
  broadcast!(genabsμₖ, absμₖ, A', LA', ws, RLws)

  return absμₖ
end

function projectvolumefromA(Θ::VolumePartsIter, Xv::XYIterLevel,
  _Xv::AbstractXYIterAlloc{TM,TV,T}, absμₖ::Tabsμₖ = _Xv.absμₖ_alloc
  ) where {TM<:CuMatrix, TV<:CuVector, T, Tabsμₖ}

  @unpack A, LA = Θ.expand
  @unpack ws, RLws = Xv
  #first map each coefficient to its appropriate

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))

  broadcast!(cudagenabsμₖ, absμₖ, A', LA', ws, RLws)
  #absμₖ .= cudagenabsμₖ.(A', LA', ws, RLws)

  return absμₖ
end


#Updates G based on values of A
#of course, this changes the RHS side of teh regressions, so any
#factored cross-sections need to be recalculated
updateGfromA!(args...;kwargs...) = error("I dropped updateGfromA!. If needed, just replace
  with Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] (though be careful about Ã vs A)")


#uses a base value of A and the G matrix to update all other values of A
function updateAfromGAₜ!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  _Xv::AbstractXYIterAlloc{TM,TV,T}, t::Int,
  boundA::Bool = PARAM[:iterboundA]) where {TM,TV,T}

  local A₀::TV
  local Aₜ::TV = Θ.A[:,t]
  dims = (T=size(Θ.A,2),K=size(Θ.A,1))

  #compute the cumulative growth
  ΠG::TM = _Xv.a_alloc
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G, Xv.Ascale)

  #if boundA #maintain consistency- comprodrow may adjust for small values
  #  Θ.G .= ΠG[:,2:end] ./ ΠG[:,1:(end-1)]
  #end

  A₀ = Aₜ ./ ΠG[:,t] #identify the intial value of assets

  Θ.A .= A₀ .* ΠG #use the growth to compute A
  (vec(Θ.A[:,t]) ≈ Aₜ) || error("vec(Θ.A[:,t]) ($(vec(Θ.A[:,t]))) ≠ Aₜ ($(Aₜ))! ") #quick integrity check
  Θ.A[:,t] .= Aₜ #preserve exact values to partly mitigate roundoff error

  return nothing
end

function updateAfromGAₜ!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  _Xv::AbstractXYIterAlloc{TM,TV,T}, t::Int,
  boundA::Bool = PARAM[:iterboundA]) where {TM<:CuMatrix,TV<:CuVector,T}

  local Aₜ::TV = Θ.A[:,t]
  dims = (T=size(Θ.A,2),K=size(Θ.A,1))

  #compute the cumulative growth
  ΠG::TM = _Xv.a_alloc
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G, Xv.Ascale)
  #ΠG .= hcat(CUDA.ones(T, dims.K, 1), cumprod(Θ.G,dims=2))

  #A₀ = Aₜ ./ ΠG[:,t] #identify the intial value of assets
  #Θ.A .= A₀ .* ΠG
  A = (Aₜ ./ ΠG[:,t]) .* ΠG
  copyto!(Θ.A, A) #use the growth to compute A
  #@assert vec(Θ.A[:,t]) ≈ Aₜ #quick integrity check
  Θ.A[:,t] .= Aₜ #preserve exact values to partly mitigate roundoff error

  return nothing
end


#NOTE: could be used in some future GPU version
cudaGᵢ(Aᵢ,LAᵢ)= Aᵢ/LAᵢ
cudafactorAᵢ(wsᵢ, RLwsᵢ, Gᵢ) = CUDA.abs(wsᵢ - RLwsᵢ/Gᵢ)
cudafactorLAᵢ(wsᵢ, RLwsᵢ, Gᵢ) = CUDA.abs(Gᵢ*wsᵢ - RLwsᵢ)

#takes a volume part object and writes the results to a df
function formresults(panel::AbstractDataFrame, ms::MeasureSpec,
  Θ::VolumePartsIter,
  grossexposures::AbstractMatrix,
  growthrates::AbstractMatrix)

  Fcoefs::Vector{Symbol} = Fcoefficients0(ms)
  Fgcoefs::Vector{Symbol} = (s->Symbol(:G,s)).(Fcoefs)

  tidx::Dict = gettidx(panel, ms)
  tidxinv::Vector{Date} = sort(collect(keys(tidx)))

  if length(Fcoefs) ≠ size(Θ.A,1)
    error("Coeffcients not in expected form. Fcoefs: $Fcoefs")
  end

  results::DataFrame = DataFrame(date=tidxinv)
  resultsmat = hcat(grossexposures,growthrates)
  results = hcat(results, DataFrame(resultsmat, [Fcoefs; Fgcoefs;]))

  return results

end
=#
