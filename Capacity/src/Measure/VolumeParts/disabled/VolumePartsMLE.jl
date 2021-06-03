

#holds the time-indexed parameters and a pre-allocation for the row-indexed versions
struct VolumePartsMLE{
    TΠ<:AbstractVector,
    TA<:AbstractMatrix,
    TV₀<:AbstractVector,
    Tτ<:AbstractVector,
    Tts<:AbstractVector,
    Texpand<:NamedTuple,
    elTΠ<:Real} <: AbstractVolumeParts{TΠ, TA, Tts, Texpand, elTΠ}

  Π::TΠ

  A::TA
  V₀::TV₀
  τ::Tτ

  #not to be accessed directly- indices for the expand function
  ts::Tts
  expand::Texpand

  #inner constructor needed to identify elTΠ
  function VolumePartsMLE(Π::TΠ, A::TA, V₀::TV₀, τ::Tτ, ts::Tts, expand::Texpand) where {
      TΠ, TA, TV₀, Tτ, Tts, Texpand}

    @assert (Π == reduce(vcat, [vec(A), vec(V₀), vec(τ)]))

    return new{TΠ, TA, TV₀, Tτ, Tts, Texpand, eltype(TΠ)}(Π,A,V₀,τ, ts,expand)
  end
end


#pre-allocates for the vpexpand
function VolumePartsMLE(A::TA, V₀::TV₀, τ::Tτ, ts::Tts, expand::Texpand) where {
    TA<:AbstractMatrix,
    TV₀<:AbstractVector,
    Tτ<:AbstractVector,
    Tts<:AbstractVector,
    Texpand<:NamedTuple}

  #create vector versions for the parameter vector
  vecA = vec(A)
  vecV₀ = vec(V₀)
  vecτ = vec(τ)

  #this will hold the parameter data
  Π = reduce(vcat, [vecA, vecV₀, vecτ])
  start=1
  linkedA = less_unsafe_wrap(TA, source=Π, start=start, dims=size(A))
  start+=length(vecA)
  linkedV₀ = less_unsafe_wrap(TV₀, source=Π, start=start, dims=size(V₀))
  start+=length(vecV₀)
  linkedτ = less_unsafe_wrap(Tτ, source=Π, start=start, dims=size(τ))
  start+=length(vecτ)
  @assert start - 1 == length(Π)

  #final checks on the linked dataframes
  @assert sum(linkedA .== A) + sum(linkedV₀ .== V₀) + sum(linkedτ .== τ)== length(Π)

  Θ = VolumePartsMLE(Π, linkedA, linkedV₀, linkedτ, ts, expand)

  return Θ
end


#basic constructor from dimension arguments
function VolumePartsMLE(T::Int, K::Int, ts::Vector{Int}, ::Type{Tvp}=PARAM[:fluxtype];
    expander::Function=volumepartsexpanders[PARAM[:fluxmodel]]) where {Tvp<:Real}

  #note- includes a work-around until Flux.gpu works
  A = abs.(randn(Tvp, K,T))
  V₀ = abs.(randn(Tvp, T-1))
  τ = abs.(randn(Tvp, T-1))

  return VolumePartsMLE(A, V₀, τ, ts, expander(A, V₀, ts))
end


#creates new VolumeParts with copies of the trainable portion and pass-by-reference
# of the non-trainable indices
function copytrainable(Θ::VolumePartsMLE)
  #copy trainable portions by value
  A  = deepcopy(Θ.A)
  V₀ = deepcopy(Θ.V₀)
  τ = deepcopy(Θ.τ)

  #copy the non-trainable portions by reference
  return VolumePartsMLE(A, V₀, τ, Θ.ts, Θ.expand)
end

#needed due to the segemented array structure
function Base.deepcopy(Θ::VolumePartsMLE)
  A  = deepcopy(Θ.A)
  V₀  = deepcopy(Θ.V₀)
  τ  = deepcopy(Θ.τ)
  ts = deepcopy(Θ.ts)
  expand = deepcopy(Θ.expand)

  return VolumePartsMLE(A,V₀,τ,ts,expand)
end


Flux.@functor VolumePartsMLE
Flux.params(Θ::VolumePartsMLE) = Flux.params(Θ.A, Θ.V₀,Θ.τ,)

function testvolumepartsindicesmle()
  Θ = VolumePartsMLE(200,5, rand(2:200, 10^5))
  Θcopy = deepcopy(Θ)

  #first test getidx
  (Θ.A == reshape(Θ[1:1000], :, 200)) || error("inconsistent getidx dimensions Θ.A")
  (Θ.V₀ == Θ[1001:1199]) || error("inconsistent getidx dimensions Θ.V₀")
  (Θ.τ == Θ[1200:1398]) || error("inconsistent getidx dimensions Θ.τ")

  #now test setidx
  testvals = rand(PARAM[:fluxtype], 1398)
  testidx = Random.randperm(1398)

  @assert allunique(testidx)
  Θ[testidx] .= testvals

  testcombined = hcat(testidx, testvals)
  testcombined = testcombined[sortperm(testcombined[:,1]),:]
  orderedtestvals = vec(testcombined[:,2])
  display(zip(Θ.Π, orderedtestvals))
  @assert sum(orderedtestvals .== Θ.Π) == length(Θ.Π)

  vec(Θcopy.A) .= orderedtestvals[1:1000]
  vec(Θcopy.V₀) .= orderedtestvals[1001:1199]
  vec(Θcopy.τ) .= orderedtestvals[1200:1398]

  (Θcopy.A == Θ.A) || error("inconsistent setidx! dimensions Θ.A")
  (Θcopy.V₀ == Θ.V₀) || error("inconsistent setidx! dimensions Θ.V₀")
  (Θcopy.τ == Θ.τ) || error("inconsistent setidx! dimensions Θ.τ")


  @info "test vol parts indices completed successfully"
end

#takes a volume part object and writes the results to a df

function formresults(panel::AbstractDataFrame, ms::MeasureSpec, Θ::VolumePartsMLE,
    grossexposures::AbstractMatrix)
  Fcoefs::Vector{Symbol} = Fcoefficients(ms)

  tidx::Dict = gettidx(panel, ms)
  tidxinv::Vector{Date} = sort(collect(keys(tidx)))

  if (length(Fcoefs) ≠ size(Θ.A,1) + 1) || Fcoefs[end] ≠ ms.Fintercept
    error("Coeffcients not in expected form. Fcoefs: $Fcoefs")
  end

  results::DataFrame = DataFrame(date=tidxinv)

  #delog the appropriate series
  V₀ = (exp).(Θ.V₀)
  τ = (exp).(Θ.τ)

  results = hcat(results, DataFrame(hcat(grossexposures, [missing; V₀...;], [missing; τ...;]), [Fcoefs; :tau]))

  return results

end



###########∇VolumeParts
#this struct/function computes the Jacobian
vec∇(Θ::VolumePartsMLE, ∇Θ) = reduce(vcat, [vec(∇Θ[Θ.A]),vec(∇Θ[Θ.V₀]),vec(∇Θ[Θ.τ])])

function test∇VolumeParts!(loss, X::NamedTuple, y, Θ::VolumePartsMLE)
  @info "test ∇VolumeParts!"
  Π = Flux.params(Θ) #formatted parameters
  @info "length: $(length(Π))"
  gs = gradient(()->loss(X,y),Π)

  rms(v) = sqrt(mean(v.^2))
  vrms(m) = rms(vec(m))

  @info "test gradrms: Θ.A-$(vrms(gs[Θ.A])); Θ.V₀-$(vrms(gs[Θ.V₀])); Θ.τ-$(vrms(gs[Θ.τ]));"

  ∇Θ! = ∇VolumeParts!(loss, X, y, Θ)

  #out .= reduce(vcat, [∇Θ[Θ.A], vec(∇Θ[Θ.V₀])])#vec∇(Θ, ∇Θ)
  out = similar(Θ.Π)
  ∇Θ!(out, Θ)
  NvecA::Int = length(vec(Θ.A))
  NvecV₀::Int = length(vec(Θ.V₀))
  Nvecτ::Int = length(vec(Θ.τ))
  @info "test ∇Θ!: Θ.A-$(rms(out[1:NvecA]));" *
    " Θ.V₀-$(vrms(out[(NvecA+1):(NvecA + NvecV₀)]));" *
    " Θ.τ-$(vrms(out[(NvecA+NvecV₀+1):(NvecA+NvecV₀+Nvecτ)]));"
end
