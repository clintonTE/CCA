

struct VolumePartsIter{
    TM<:AbstractMatrix,
    TV<:AbstractVector,
    Tg <: AbstractVector,
    Tts<:AbstractVector,
    Texpand<:NamedTuple,
    T<:Real} <: AbstractVolumePartsIter{TM, TV, T}

  #parameters
  A::TM
  Ã::TM #factored A, used in regressions
  G::TM

  g::Tg #vector version by referene of G, should be equivelent

  ts::Tts
  tsM1::Tts
  expand::Texpand

  dims::NamedTuple

  #inner constructor needed to identify elTΠ
  function VolumePartsIter(
    A::TM, Ã::TM, G::TM, g::Tg, ts::Tts, tsM1::Tts, expand::Texpand,
    ::Type{TV} = typeof(similar(A,1))) where {
      TM, TV, Tg, Tts, Texpand}

    return new{TM, TV, Tg, Tts, Texpand, eltype(TM)}(A,Ã,G,g,ts,tsM1, expand)
  end
end

#pre-allocates for the vpexpand
function VolumePartsIter(A::TA,  Ã::TA, G::TA, ts::Tts, tsM1::Tts=ts.-1,#=
  ::Type{TM}, ::Type{TV}, ::Type{T}=eltype(TM)=#) where {
    TA<:AbstractMatrix, TM, TV, T, Tts<:AbstractVector}

  g = vec(G)

  expand = (
    A = view(A, :, ts),
    LA = view(A, :, tsM1),
    Ã = view(Ã, :, ts),
    LÃ = view(Ã, :, tsM1),
    G = view(G, :, tsM1),
    index  = (
      A = to_indices(A, (:, ts)),
      LA = to_indices(A, (:, tsM1)),
      Ã = to_indices(Ã, (:, ts)),
      LÃ = to_indices(Ã, (:, tsM1)),
      G = to_indices(G, (:, tsM1)),
    ))

  Θ = VolumePartsIter(A, Ã, G, g, ts, tsM1, expand)
  return Θ
end


#basic constructor from dimension arguments
function VolumePartsIter(T::Int, K::Int, ts::Vector{Int}, ::Type{Tvp}=PARAM[:itertype],
  ::Type{TM}=PARAM[:itergpu] ? CUDA.CuMatrix{Tvp} : Matrix{Tvp},
  ::Type{TV}=PARAM[:itergpu] ? CUDA.CuVector{Tvp} : Vector{Tvp},
    ) where {Tvp<:Real, TV<:AbstractVector{Tvp}, TM<:AbstractMatrix{Tvp}}

  A = TM(abs.(randn(Tvp, K,T)))
  Ã = TM(abs.(randn(Tvp, K,T)))
  G = TM(ones(Tvp, K, T-1))

  return VolumePartsIter(A, Ã, G, ts, #=TM, TV, Tvp=#)
end

#needed due to the segemented array structure
function Base.deepcopy(Θ::VolumePartsIter)
  A  = deepcopy(Θ.A)
  Ã = deepcopy(Θ.Ã)
  G = deepcopy(Θ.G)


  return VolumePartsIter(A,Ã,G,Θ.ts,Θ.tsM1)
end

#invididual functions for computing volume constriubtions
#genabsμₖ(Aᵢ, LAᵢ, wsᵢ, Rtwsᵢ, RLwsᵢ) = abs(LAᵢ*(Rtwsᵢ - RLwsᵢ)) + abs(Aᵢ*wsᵢ-LAᵢ*Rtwsᵢ)
genabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = abs(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)
genabsμₖ(Aᵢ, LAᵢ, wsᵢ, Rtwsᵢ, RLwsᵢ)  = genabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ)
signabsμₖ(Aᵢ, LAᵢ, wsᵢ, Rtwsᵢ, RLwsᵢ) = sign(Aᵢ*wsᵢ-LAᵢ*Rtwsᵢ) #derivitive

#project the volume at a particular time t
#WARNING- the below only works when sum(Lwsₜ)==sum(wsₜ). If this is not true, need to carry
#an extra vector for the previous weights
function projectt(Aᵢ, LAᵢ, wsₜ, Rtwsₜ, RLwsₜ,  Mₜ)

  return Mₜ * genabsμₖ.(Aᵢ', LAᵢ', wsₜ, Rtwsₜ, RLwsₜ,)
end
function projectt(Aᵢ, LAᵢ, wsₜ, Rtwsₜ, RLwsₜ, ::Nothing)
  return genabsμₖ.(Aᵢ', LAᵢ', wsₜ, Rtwsₜ, RLwsₜ,)
end

#projectt(Aᵢ, LAᵢ, wsₜ, RLwsₜ, ::Nothing) = genabsμₖ.(Aᵢ', LAᵢ', wsₜ, RLwsₜ, 1.0 +sum(RLwsₜ)-sum(wsₜ))

#project the volume for each time t
function projectbyt(A, xsection)
  dims = (K = size(A,1), T = size(A,2))
  @unpack xws, xRtws, xRLws, xM = xsection

  #project the volume at each time t, then concatenate the results
  factored =  reduce(vcat, broadcast(1:(dims.T-1), xws, xRtws,xRLws, xM) do t, xwsₜ, xRtwsₜ, xRLwsₜ, xMₜ
    projectt(A[:,t+1], A[:,t], xwsₜ, xRtwsₜ, xRLwsₜ, xMₜ)
  end)

  return factored
end

#the below is used to get the active volume
genμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = Aᵢ * wsᵢ - LAᵢ * RLwsᵢ
genμₖ(Aᵢ, LAᵢ, wsᵢ, Rtwsᵢ, RLwsᵢ)  = genμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ)
function projectpurchasest(Aᵢ, LAᵢ, wsₜ, Rtwsₜ, RLwsₜ,  Mₜ)

  return Mₜ * genμₖ.(Aᵢ', LAᵢ', wsₜ, Rtwsₜ, RLwsₜ,)
end
function projectpurchasest(Aᵢ, LAᵢ, wsₜ, Rtwsₜ, RLwsₜ, ::Nothing)
  return genμₖ.(Aᵢ', LAᵢ', wsₜ, Rtwsₜ, RLwsₜ,)
end
function projectpurchasesbyt(A, xsection)
  dims = (K = size(A,1), T = size(A,2))
  @unpack xws, xRtws, xRLws, xM = xsection

  #project the volume at each time t, then concatenate the results
  factored =  reduce(vcat, broadcast(1:(dims.T-1), xws, xRtws,xRLws, xM) do t, xwsₜ, xRtwsₜ, xRLwsₜ, xMₜ
    projectpurchasest(A[:,t+1], A[:,t], xwsₜ, xRtwsₜ, xRLwsₜ, xMₜ)
  end)

  return factored
end


#this just runs the regression for the purposes of extracting the coefficients
function coef(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl,
  RegressionType::Val) where {TM, TV, T}

  @unpack G,Ã,A = Θ
  dims = (K = size(A,1), T = size(A,2))
  initializeÃandG!(Θ, Xv)
  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  β = reg(RHS, Xv.ṽ, RegressionType)
  return β
end


#version that performs cross-sectional regressions as needed
function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterControl,
  RegressionType::Val) where {TM,TV,T}
  @unpack G,Ã,A = Θ

  dims = (K = size(A,1), T = size(A,2))

  #set Ã to the cumulative product of G
  initializeÃandG!(Θ, Xv)
  #@info "type of Xv: $(typeof(Xv))"

  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #run the regression for A₁
  β = reg(RHS, Xv.ṽ, RegressionType)
  #below shoudl be equivelent to A[:,1] .= β[1:dims.K]
  #β = (b->ifelse(b>0, b, b-100.)).(β) |> TV

  copyto!(A, 1, β, 1, dims.K)
  #A[:,1] .= parent(β)[1:dims.K]


  #absμ .= vec(sum(absμₖ,dims=2))
  #@eval Main controlbeta = $(β |> Vector)
  absμ = RHS * β


  return absμ
end

function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterControl,
  RegressionType::Val{:none},
  noregressionbase::T = T(PARAM[:iternoregressionbase])) where {TM,TV,T}

  @unpack G,Ã,A = Θ

  dims = (K = size(A,1), T = size(A,2))

  #set Ã to the cumulative product of G
  initializeÃandG!(Θ, Xv)
  #@info "type of Xv: $(typeof(Xv))"

  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #run the regression for A₁

  βW̃ = (size(Xv.W̃,2) > 0) ? reg(Xv.W̃, Xv.ṽ, Val{:cholesky}()) : TV()
  β = [TV(ones(T, dims.K) .* noregressionbase); βW̃]
  #below shoudl be equivelent to A[:,1] .= β[1:dims.K]
  copyto!(A, 1, β, 1, dims.K)
  #A[:,1] .= parent(β)[1:dims.K]


  #absμ .= vec(sum(absμₖ,dims=2))
  #@eval Main controlbeta = $(β |> Vector)
  absμ = RHS * β


  return absμ
end

#we can probably dispatch all of this just from Xv
#NOTE: this is mainly for ZYGOTE
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T},
    Xv::XYIterControl,
    RegressionType::Val) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G)
  Ã = hcat(someones, prodmatG) |> TM

  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #run the regression for A₁
  β = reg(RHS, Xv.ṽ, RegressionType)
  #β = (b->ifelse(b>0, b, b-100.)).(β) |> TV

  #β = ones(size(RHS,2))
  absμ = RHS * β


  return absμ
end

#we can probably dispatch all of this just from Xv
#NOTE: this is mainly for ZYGOTE
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T},
    Xv::XYIterControl,
    RegressionType::Val{:none},
    noregressionbase::T = T(PARAM[:iternoregressionbase])) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G)
  Ã = hcat(someones, prodmatG) |> TM

  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #run the regression for A₁
  βW̃ = (size(Xv.W̃,2) > 0) ? reg(Xv.W̃, Xv.ṽ, Val{:choleskyzygote}()) : TV()
  β = [TV(ones(T, dims.K) .* noregressionbase); βW̃]

  absμ = RHS * β


  return absμ
end

#NOTE: this is mainly for ZYGOTE
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T},
    Xv::XYIterControl,
    RegressionType::Val{:zygotefallbackreg}) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G, limitA=true)
  Ã = hcat(someones, prodmatG) |> TM

  RHS::TM = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #run the regression for A₁
  β = reg(RHS, Xv.ṽ, RegressionType)
  #β = (b->ifelse(b>0, b, b-100.)).(β) |> TV

  absμ = RHS * β


  return absμ
end

function lossforzygote(G::Matrix{T}, Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl,;
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}()) where {TM,TV,T}

  #throw("hello")

  #ε = factored * Aₜ .- Xv.ṽ
  Σε2 = sum((Θ(G, Xv, RegressionType) .- Xv.ṽ).^2) #0.3ms

  return Σε2
end

function ∇lossforzygote(Θ, Xv::AbstractXYIter{TM,TV,T},  GrowthType::Val;
    RegressionType=Val(PARAM[:iterregressiontypezygote]),
    GradType=Val(PARAM[:itergradtype])) where {TM,TV,T}

  #the actual gradient function

  function ∇loss(G::Matrix, GradType::Val{:zygote})
    local grad

    try
      grad = gradient((G)->lossforzygote(G, Θ, Xv,
        RegressionType=RegressionType), G)[1]
    catch err
      (err == InterruptException()) && throw(err)
      errmessage = "$err"
      errmessage = length(errmessage)  > 1000 ? errmessage[1:1000] : errmessage
      @warn "Gradient failed with error $errmessage\nAttempting fallback method"
      try
        grad = gradient((G)->lossforzygote(G, Θ, Xv,
          RegressionType=Val{:zygotefallbackreg}()), G)[1]
      catch err2
        @error "Gradient failed with error $errmessage" exception=(err, catch_backtrace())
        errmessage2 = "$err2"
        errmessage2 = length(errmessage2)  > 1000 ? errmessage2[1:1000] : errmessage2
        @error "Gradient fallback also failed with error $errmessage2" exception=(err2, catch_backtrace())
      end
    end

  end

  function ∇loss(G::Matrix, GradType::Val{:ag})
   return ∇lossag(G, Θ, Xv)
  end


  ∇loss(G::Matrix) = ∇loss(G::Matrix, GradType)


  function ∇vec!(out, x::TV, ::Val{:identity}) where TV <: Vector

    if !(x===Θ.g) #if they are the same array, no need to update
      copyto!(Θ.g, x)
    end

      copyto!(out, ∇loss(Θ.G)) #converts to TV

    return nothing
  end

  function ∇vec!(out, x::TV, ::Val{:identity}) where TV <: CuVector

    copyto!(Θ.g, x)
    G = Θ.G |> Matrix
    copyto!(out, ∇loss(G)) #converts to TV

    return nothing
  end

  function ∇vec!(out, x::TV, ::Val{:log}) where TV <: Vector
    #@info "∇vec! log"
    (x ≡ Θ.g) && throw("x and g should not be the same vector if growthtype≠identity")
    copyto!(Θ.g, exp.(x))

    copyto!(out, ∇loss(Θ.G) .* Θ.G) #adjust by G since d(expy)/dx=expy*dy/dx

    return nothing
  end

  function ∇vec!(out, x::TV, ::Val{:log}) where TV <: CuVector
    (x ≡ Θ.g) && throw("x and g should not be the same vector if growthtype≠identity")
    copyto!(Θ.g, exp.(x))
    G = Θ.G |> Matrix
    copyto!(out, ∇loss(G) .* G)
    #@info "sumabsout: $(sum(abs.(out)))"
  end


  ∇vec!(out, x) = ∇vec!(out, x |> TV, GrowthType) #this version performs the conversion
  ∇vec!(out, x::TV) = ∇vec!(out, x, GrowthType)

  return ∇vec!
end

#Updates G based on values of A
#of course, this changes the RHS side of teh regressions, so any
#factored cross-sections need to be recalculated
updateGfromA!(args...;kwargs...) = throw("I dropped updateGfromA!. If needed, just replace
  with Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] (though be careful about Ã vs A)")


#uses a base value of A and the G matrix to update all other values of A
function updateAfromGAₜ!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},t::Int) where {TM,TV,T}

  local A₀::TV
  local Aₜ::TV = Θ.A[:,t]
  dims = (T=size(Θ.A,2),K=size(Θ.A,1))

  #compute the cumulative growth
  ΠG::TM = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G)


  A₀ = Aₜ ./ ΠG[:,t] #identify the intial value of assets

  Θ.A .= A₀ .* ΠG #use the growth to compute A
  (vec(Θ.A[:,t]) ≈ Aₜ) || throw("vec(Θ.A[:,t]) ($(vec(Θ.A[:,t]))) ≠ Aₜ ($(Aₜ))! ") #quick integrity check
  Θ.A[:,t] .= Aₜ #preserve exact values to partly mitigate roundoff error

  return nothing
end

function updateAfromGAₜ!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  t::Int) where {TM<:CuMatrix,TV<:CuVector,T}

  local Aₜ::TV = Θ.A[:,t]
  dims = (T=size(Θ.A,2),K=size(Θ.A,1))

  #compute the cumulative growth
  ΠG::TM = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G)
  #ΠG .= hcat(CUDA.ones(T, dims.K, 1), cumprod(Θ.G,dims=2))

  #A₀ = Aₜ ./ ΠG[:,t] #identify the intial value of assets
  #Θ.A .= A₀ .* ΠG
  A = (Aₜ ./ ΠG[:,t]) .* ΠG
  copyto!(Θ.A, A) #use the growth to compute A
  #@assert vec(Θ.A[:,t]) ≈ Aₜ #quick integrity check
  Θ.A[:,t] .= Aₜ #preserve exact values to partly mitigate roundoff error

  return nothing
end

#computs the loss, but does so
#given the current values of A irrespective of G
function lossfromA(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControl, RegressionType::Val
    ) where {TM,TV,T}

  @unpack G,Ã,A = Θ
  dims = (K = size(A,1), T = size(A,2))
  RHS::TM = hcat(projectbyt(A, Xv.xsection), Xv.W̃) #NOTE the use of A, not Ã
  β = reg(RHS, Xv.ṽ, RegressionType)
  absμ = RHS * β

  loss = sum((Xv.ṽ .- absμ) .^2 )

  return loss
end



#NOTE: could be used in some future GPU version
cudaGᵢ(Aᵢ,LAᵢ)= Aᵢ/LAᵢ
cudafactorAᵢ(wsᵢ, RLwsᵢ, Gᵢ) = CUDA.abs(wsᵢ - RLwsᵢ/Gᵢ)
cudafactorLAᵢ(wsᵢ, RLwsᵢ, Gᵢ) = CUDA.abs(Gᵢ*wsᵢ - RLwsᵢ)


#takes a volume part object and writes the results to a df
function formresults(panel::AbstractDataFrame, ms::MeasureSpec, Θ::AbstractVolumeParts{TM,TV,T},
    Xv::XYIterControl,  ::Val{:iter}) where {TM,TV,T}

  grossexposures = (Θ.A |> Matrix)' .* grossleverages(ms)'
  growthrates = vcat(zeros(eltype(grossexposures), size(grossexposures,2))', (Θ.G|>Matrix!)')

  Fcoefs::Vector{Symbol} = ms.Fξs
  Fgcoefs::Vector{Symbol} = (s->replace(string(s),r"^Z"=>"G_") |> Symbol).(Fcoefs)
  FWcoefs::Vector{Symbol} = [ms.FW...;]

  if length(Fcoefs) ≠ size(Θ.A,1)
    throw("Coeffcients not in expected form. Fcoefs: $Fcoefs")
  end

  resultsgrossexposures = DataFrame(grossexposures,Fcoefs)
  resultsgrowthrates = DataFrame(growthrates,Fgcoefs)

  tidx::Dict = gettidx(panel, ms)
  tidxinv::Vector{Date} = sort(collect(keys(tidx)))


  results::DataFrame = DataFrame(date=tidxinv)

  results = [results resultsgrossexposures resultsgrowthrates]

  return results

end
