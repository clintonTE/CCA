
#RegressionTypes = (#=glm=Val{:glm}(), =#cholesky=Val{:cholesky}(), qr=Val{:qr}())
reg(X::AbstractMatrix{T},y::Vector{T}, ::Val{:glm}) where T = GLM.coef(lm(X, y))
reg(X::AbstractMatrix{T},y::AbstractArray{T}, ::Val{:cholesky}
  ) where T = cholesky!(Symmetric(X'*X))\ (X' * y)
reg(X::AbstractMatrix{T},y::Vector{T}, ::Val{:qr}) where T = qr(X)\y
reg(X::AbstractMatrix{T}, y::AbstractArray{T}, ::Val{:svd}) where T = svd(X)\y
reg(X::CuMatrix{T}, y::CuArray{T}, ::Val{:cholesky}) where T = hyperreg(X,y) |> CUDA.vec
reg(X::AbstractMatrix{T}, y::AbstractArray{T}, ::Val{:choleskyzygote}) where T = cholesky(X'*X)\ (X' * y)
reg(X::CuMatrix{T}, y::AnyCuVecOrMat{T}, ::Val{:choleskyzygote}) where T = CUDA.cholesky(X'*X)\ (X' * y)
#reg(X::CuMatrix{T}, y::CuVector{T}, ::Val{:none}) where {TV, T} = ones(T, size(X,2)) .* 0.75 |> CuVector{T}
#reg(X::AbstractMatrix{T}, y::AbstractArray{T}, ::Val{:none}) where {TV, T} = ones(T, size(X,2)) .* 0.75
#these methods run with an intercept, which is thrown out and not returned.
reg(X::AbstractMatrix{T}, y::AbstractArray{T},::Val{:intercept}) where T = throw(
  "intercept regressions should not be directly called")
reg(X::AbstractMatrix{T}, y::AbstractArray{T},::Val{:interceptzygote}) where T = throw(
  "intercept regressions should not be directly called")

#=reg(X::AbstractMatrix{T}, y::AbstractArray{T},::Val{:intercept}) where T = reg(
  hcat(ones(T, length(y)), X), y, Val{:cholesky}())[2:end]
reg(X::CuMatrix{T}, y::CuArray{T}, ::Val{:intercept}) where T = reg(
  CUDA.hcat(CUDA.ones(T, length(y)), X), y, Val{:cholesky}())[2:end]
reg(X::AbstractMatrix{T}, y::CuArray{T}, ::Val{:interceptzygote}) where T = reg(
  CUDA.hcat(zygotecuones(T, length(y)), X), y, Val{:choleskyzygote}())[2:end]=#

function reg(X::Matrix{T}, y::Array{T}, ::Val{:fallbackreg}) where T

  #@eval Main X=$X|>deepcopy
  #@eval Main y=$y|>deepcopy
  #= NOTE: Let USVt be the SVD of X, . Note by propertys of svd, Ut*U=Vt*V=I and S  is diagonal
  Then we seek to solve:
    USVt*β = y
    SVt*β = U'y
    β = V*S^-1*U'y = V*(diag(S^-1) .* U' * y)
  =#
  #
  #sX = svd(X)
  #β = sX.V * (sX.S .^-1 .* (sX.U' * y))


  return pinv(X'*X, rtol=sqrt(eps(T)))*(X'*y)
end
#these are just conversion methods for the fallback
reg(X::TX, y::Ty, RegressionType::Val{:fallbackreg}, ::Type{T} = eltype(X) #fallback gpu method
  ) where {TX,Ty,T} = reg(X|>Matrix, y|>Array{T}, RegressionType) |> Ty
reg(X::Matrix{T}, y::Array{T}, RegressionType::Val{:zygotefallbackreg}
  ) where T = reg(X,y,Val{:fallbackreg}())
reg(X::TX, y::Ty, ::Val{:zygotefallbackreg}, ::Type{T} = eltype(X) #fallback gpu method
  ) where {TX,Ty,T} = reg(X|>zygotedevice2host,y|>zygotedevice2host,Val{:fallbackreg}()) |> Ty

#reg(X::AbstractMatrix, y::AbstractArray, ::Val{:qrzygote}) = qr(X)\y

#not sure if there is a potential gpu benefit here- probably would need a very custom kernal
#NOTE: fairly sure these are not needed. If they are needed, just uncomment
#const MIN_MAG = .01
#const MAX_MAG = 3.0
#nonzero(x::T,y) where T = ifelse(isfinite(x) && abs(x) < MIN_MAG, y, x)

#this sets A₁ to 1.0 and all other values to the amount implied from G
#of course, this destroys any current values of A
function initializeÃandG!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T}) where {TM,TV,T}
  Θ.Ã[:,1] .= T(1.0)

  #below is better but it needs to sync with Zygote
  #someones::TM = ones(size(Θ.A,1),1) |> TM

  #copyto!(Θ.Ã, cumprod(hcat(someones, Θ.G),dims=2))

  #below is a bandaid to sync up the initialization with Zygote
  Θ.Ã[:,2:end] .= cumprodbyrow(Θ.G)

  return nothing
end


#this factors X. NOTE assumes we are solving for A₁ AND that A₁ .== 1
function initializeÃandfactor(Θ::AbstractVolumePartsIter{TM, TV,T},
  Xv::AbstractXYIter) where {TM, TV,T}

  @unpack Ã, LÃ = Θ.expand
  @unpack ws, Rtws, RLws= Xv

  initializeÃandG!(Θ, Xv)
  @assert all(Θ.Ã[:,1] .≈ T(1.0)) #always true if Ã is initialized

  factored = genabsμₖ.(Ã', LÃ', ws, Rtws, RLws)

  return factored
end


#this factors X but assumes we are solving for A₁ AND that A₁ .== 1
function initializeÃandfactor(Θ::VolumePartsIter,
  Xv::AbstractXYIter{TM,TV,T}) where {TM<:CuMatrix, TV<:CuVector, T}

  @unpack Ã, LÃ = Θ.expand
  @unpack ws, Rtws, RLws= Xv

  initializeÃandG!(Θ, Xv)
  @assert all(Θ.Ã[:,1] .≈ T(1.0)) #always true if Ã is initialized
  #first map each coefficient to its appropriate

  factored = genabsμₖ.(Ã', LÃ', ws, Rtws, RLws)

  return factored
end


#updates A₁ for a single pass by running the appropriate regression
#returns the loss function
#WARNING- we do not update all of A here after running the regression
function updateA₁!(loss, Θ::VolumePartsIter,
  Xv::AbstractXYIter{TM, TV,T},
  RegressionType::Val) where {TM, TV,T}

  local λ
  #try
    λ = loss(Θ, Xv, RegressionType)
  #=catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage) > 1000 ? errmessage[1:1000] : errmessage
    WARNING- I believe the below code is obsolete, so it is commented out
    @warn "Regression failed with error $errmessage\nAttempting fallback method"
    λ = loss(Θ, Xv, Val{:fallbackreg}())
  #@assert typeof(λ) <: Float64
  #@info "update: λ=$λ"
end=#
  return λ
end


  #set the initial value for g from the growth type
  #called for the inital vector at the beginning of optimization
function initializeΘgvector(g::TV, ::Type{T}, GrowthType::Val, G::Union{TV,Nothing}=nothing;
  initializeto1::Bool=true) where {TV, T}
  _initializeΘg(g::AbstractVector, ::Val{:identity}) = g |> Vector |> deepcopy
  _initializeΘg(g::AbstractVector, ::Val{:log}) = g |> Vector |> deepcopy .|> log
  function _initializeΘg(g::AbstractVector, ::Val{:twosidedlog},
    lowerbound::T=PARAM[:optimlowerbound], twosidedloglower::T=PARAM[:itertwosidedloglower],
    upperbound::T=PARAM[:optimupperbound], twosidedlogupper::T=PARAM[:itertwosidedlogupper])
    gout = g |> Vector |> deepcopy
    Δbound = upperbound
    Δtwosidedlog = twosidedlogupper - twosidedloglower
    gout .= -(twosidedloglower)*Δbound/Δtwosidedlog #solved to set G=1
    gout
  end

  if initializeto1
    g .= T(1.0)
  end
  gidentity = g |> Vector |> deepcopy #use the pretransformation version for verification later
  gout = _initializeΘg(g, GrowthType)
  #g .= gout |> TV
  copyto!(g, gout)
  @assert all(conditiongrowth(gout, GrowthType) |> Vector .≈ gidentity)
  return gout
end


conditiongrowth(g, ::Val{:identity})= g
conditiongrowth(g, ::Val{:log}) = exp.(g)
conditiongrowth(g::CuVector, ::Val{:log}) = CUDA.exp.(g)
function conditiongrowth(g::AbstractVector{T}, ::Val{:twosidedlog},
  lowerbound::T=PARAM[:optimlowerbound], twosidedloglower::T=PARAM[:itertwosidedloglower],
  upperbound::T=PARAM[:optimupperbound], twosidedlogupper::T=PARAM[:itertwosidedlogupper]
  ) where T

  (upperbound + lowerbound ≈ 0.0) || throw("bounds of twosided log must be symmetric about 0")
  Δbound = upperbound
  Δtwosidedlog = twosidedlogupper - twosidedloglower
  function maptolog(x) #maps the input to the two sided log of g e.g. [-10,10]->±exp([-2.5,2.5])
    normalizedabsx::T = abs(x)/Δbound
    mappedabsx::T = exp(normalizedabsx*Δtwosidedlog + twosidedloglower)
    return (signsans0(x)) * mappedabsx
  end

  return maptolog.(g)
end

function updateA₁!(x::TV, loss, Θ::VolumePartsIter,
  Xv::AbstractXYIter{TM, TV,T},
  RegressionType::Val, ::Val{:identity}
  ) where {TM, TV,T}

  #no need to copy if they are equivelent
  (Θ.g===x) && return throw("I disallowed Θ.g===x due to race error potential
    (It shouldn't make a big difference anyway)")
  #@info "updateA₁!: growth type identity"
  copyto!(Θ.g,x)
  return updateA₁!(loss,Θ,Xv,RegressionType)
end


function updateA₁!(x::TV, loss, Θ::VolumePartsIter,
  Xv::AbstractXYIter{TM, TV,T},
  RegressionType::Val, GrowthType::Val;
  ) where {TM, TV,T}

  (Θ.g===x) && throw("x should be a different vector if solving for GrowthType $GrowthType")
  #no need to copy if they are equivelent
  copyto!(Θ.g, conditiongrowth(x, GrowthType))
  #@info "updateA₁!: growth type identity"
  return updateA₁!(loss,Θ,Xv,RegressionType)
end

#handles the CUDA case by converting the supplied vector to a CUDA vector
updateA₁!(x::Vector, loss, Θ::VolumePartsIter,
  Xv::AbstractXYIter{TM, TV,T},
  RegressionType::Val, GrowthType::Val;
  ) where {TM<:CuMatrix, TV<:CuVector,T} = updateA₁!(
  x |> TV, loss,Θ,Xv,RegressionType,GrowthType)
