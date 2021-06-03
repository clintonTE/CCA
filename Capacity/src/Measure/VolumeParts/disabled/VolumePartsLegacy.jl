

#define the critical parts of the interface
#note- we may need another layer if we make a version without an Intercept
abstract type AbstractVolumeParts{
    TM<:AbstractMatrix, TV<:AbstractVector, T<:Real} <: AbstractVector{T} end

#convenience assignment operator
function updatexfromy!(dest::T, source::T) where T<:AbstractVolumeParts
  dest .= source
  return nothing
end

#this is the predict function
function (Θ::AbstractVolumeParts)(X::NamedTuple, ::AbstractVector, ::Val{:level})
  @unpack ws, RLws = X
  #first map each coefficient to its appropriate
  A = view(Θ.A, Θ.expand.Aidx...)
  LA = view(Θ.A, Θ.expand.LAidx...)
  V₀ = view((exp).(Θ.V₀), Θ.expand.V₀idx)

  absμₖ = abs.(A .* ws' .- LA .* RLws')
  absμ = vec(sum(absμₖ,dims=1)) .+ V₀

  return absμ
end

#this is the predict function
function (Θ::AbstractVolumeParts)(X::NamedTuple, ::AbstractVector, ::Val{:levellog})
  @unpack ws, RLws = X
  #first map each coefficient to its appropriate
  A = view(exp.(Θ.A), Θ.expand.Aidx...)
  LA = view(exp.(Θ.A), Θ.expand.LAidx...)
  V₀ = view((exp).(Θ.V₀), Θ.expand.V₀idx)

  absμₖ = abs.(A .* ws' .- LA .* RLws')
  absμ = vec(sum(absμₖ,dims=1)) .+ V₀

  return absμ
end

genabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = abs(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)
signabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = sign(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)
cudagenabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = CUDA.abs(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ) #one of these should go away
cudaabsμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = CUDA.abs(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)

#makes volumeparts callable for a level estimate
function (Θ::AbstractVolumeParts)(X::NamedTuple, ::CuVector{T}, ::Val{:level}) where T
  @unpack ws, RLws = X
  #re-creating the matrices gets around
  #https://github.com/JuliaGPU/CUDA.jl/issues/684
  #and https://github.com/FluxML/Zygote.jl/issues/600

  #REMINDER- A, LA, and V₀ are in logs here
  A = CuMatrix{T}(Θ.A[Θ.expand.Aidx...])
  LA = CuMatrix{T}(Θ.A[Θ.expand.LAidx...])
  V₀ = CUDA.exp.(Flux.CuArrays.CuVector{T}((Θ.V₀)[Θ.expand.V₀idx]))

  #absμₖ = CUDA.abs.(A .* ws' .- LA .* RLws')
  absμₖ = cudaabsμₖ.(A, LA, ws', RLws')
  absμ = CUDA.vec(CUDA.sum(absμₖ,dims=1))

  return absμ .+ V₀
end


#makes volumeparts callable for a level estimate
function (Θ::AbstractVolumeParts)(X::NamedTuple, ::CuVector{T}, ::Val{:levellog}) where T
  @unpack ws, RLws = X
  #re-creating the matrices gets around
  #https://github.com/JuliaGPU/CUDA.jl/issues/684
  #and https://github.com/FluxML/Zygote.jl/issues/600

  #REMINDER- A, LA, and V₀ are in logs here
  A = CUDA.exp.(CuMatrix{T}(Θ.A[Θ.expand.Aidx...]))
  LA = CUDA.exp.(CuMatrix{T}(Θ.A[Θ.expand.LAidx...]))
  V₀ = CUDA.exp.(CuVector{T}((Θ.V₀)[Θ.expand.V₀idx]))

  #absμₖ = CUDA.abs.(A .* ws' .- LA .* RLws')
  absμₖ = cudaabsμₖ.(A, LA, ws', RLws')
  absμ = CUDA.vec(CUDA.sum(absμₖ,dims=1)) .+ V₀

  return absμ
end

formresults(panel::AbstractDataFrame, ms::MeasureSpec, Θ::AbstractVolumeParts,
  ::Val{:levellog}) = (
  formresults(panel, ms, Θ, exp.(Θ.A)' .* grossleverages(ms)'))

formresults(panel::AbstractDataFrame, ms::MeasureSpec, Θ::AbstractVolumeParts,
  ::Val{:level}) = (
  formresults(panel, ms, Θ, (Θ.A |> Matrix!)' .* grossleverages(ms)'))

function formresults(panel::AbstractDataFrame, ms::MeasureSpec, Θ::AbstractVolumeParts,
    ::Val{:iter})

  Amat = (Θ.A |> Matrix!)' .* grossleverages(ms)'
  Gmat = vcat(zeros(eltype(Amat), size(Amat,2))', (Θ.G|>Matrix!)')
  return formresults(panel, ms, Θ, Amat, Gmat)
end


###The below provides 1d indexing of the parameters
Base.length(Θ::AbstractVolumeParts) = length(Θ.Π)
Base.size(Θ::AbstractVolumeParts) = (Base.length(Θ),)
Base.IndexStyle(::Type{<:AbstractVolumeParts}) = IndexLinear()

Base.getindex(Θ::AbstractVolumeParts, inds) = Θ.Π[inds]
Base.setindex!(Θ::AbstractVolumeParts, x::Real, i::Int) = (Θ.Π[i] = x)


∇rms(Θ::AbstractVolumeParts, ∇Θ) = sqrt(mean((vec∇(Θ, ∇Θ)) .^ 2))


####Expanders

levelexpandlsq(A,V₀, ts) = (
  Aidx = to_indices(A, (:, ts)),
  LAidx = to_indices(A, (:, ts .- 1)),
  V₀idx = ts .- 1)

levelexpandmle(A, V₀, ts) = (
  Aidx = to_indices(A, (:, ts)),
  LAidx = to_indices(A, (:, ts .- 1)),
  V₀idx = ts.-1,
  τidx = ts.-1)


volumepartsexpanders = Dict{Symbol, Function}(
  :levelmle => levelexpandmle,
  :levellsq => levelexpandlsq,
)


cudasignμₖ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = CUDA.sign(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)
cudanoop(x) = x
