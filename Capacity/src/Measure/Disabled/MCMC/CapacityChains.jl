using Revise, MCMCChains, LinearAlgebra, UnPack


####################

#holds the main state vector
#unfinished and untested
abstract type AbstractCCState{T} <: AbstractVector{T} end
struct CCState{TΠ<:AbstractVector, Tparams<:Dict, Tparamsidx, T} <: AbstractCCState{T}
  Π::TΠ
  params::Tparams
  paramsidx::Tparamsidx

  CCState(Π::TΠ, params::Tparams, paramsidx::Tparamsidx) where {
    TΠ, Tparams, Tparamsidx}=new{TΠ, Tparams, Tparamsidx, eltype(Π)}(Π, params, paramsidx)
end


#forms a state vector given a parameter index
#main constructor, but internally facing only
#the paramsidx object should map all params to indices in a vector
function CCState(paramsidx::Dict, ::Type{T}) where T
  Nparams = sum(length.(values(paramsidx)))
  Π = rand(T, Nparams)

  #build high-performance array wraps for indexing into the state vector
  params::Dict{Symbol, Union{Ref,AbstractArray{T}}} = Dict()
  for (n,p) ∈ pairs(paramsidx)
    Tp = typeof(p)
    if Tp <: Int
      params[n] = Ref(Π, p)
    elseif Tp <: AbstractVector
      params[n] = less_unsafe_wrap(Vector{T}, source=Π, start=minimum(p), dims=size(p))
    elseif Tp <: AbstractMatrix
      params[n] = less_unsafe_wrap(Matrix{T}, source=Π, start=minimum(p), dims=size(p))
    else
      @assert false
    end
  end

  #integrity checks
  length(Π)==sum(length.(values(params)))
  coversindex = falses(length(Π)) #keep track of which indices params references
  for (n,p) ∈ pairs(paramsidx)
    @assert size(p) ≡ size(params[n]) #check dimensions
    indicesused = [p...;]
    @assert all(reshape(Π[indicesused], size(params[n])) .≡ params[n]) #checks for substantial identity
    @assert all((!).(coversindex[indicesused])) #makes sure we have no overlaps
    coversindex[indicesused] .= true
  end
  @assert all(coversindex) #makes sure params references all indices of Π

  return CCState(Π,params,paramsidx)
end

###The below provides 1d indexing of the parameters
Base.length(Ψ::AbstractCCState) = length(Ψ.Π)
Base.size(Ψ::AbstractCCState) = (Base.length(Ψ),)
Base.IndexStyle(::Type{<:AbstractCCState}) = IndexLinear()

Base.getindex(Ψ::AbstractCCState, inds) = Ψ.Π[inds]
Base.setindex!(Ψ::AbstractCCState, x::Real, i::Int) = (Ψ.Π[i] = x)

#shortcut for directly indexing into the capacity state
function Base.getproperty(Ψ::TΨ, Fparam::Symbol) where TΨ <: AbstractCCState
  hasfield(TΨ, Fparam) && return getfield(Ψ, Fparam)
  Ψ.params[Fparam]
end

function Base.show(io::IO, ::MIME"text/plain", Ψ::AbstractCCState{T}) where T
  for k ∈ keys(Ψ.params)
    println(io, "Showing contents of $k: ")
    printmln(Ψ.params[k])
  end
end

function Base.deepcopy(Ψ::AbstractCCState{T}) where T
  Ψcpy = CCState(Ψ.paramsidx, T)

  #integrity checks
  Ψcpy .= Ψ
  @assert all(Ψcpy.Π .≡ Ψ.Π)
  for s ∈ keys(Ψ.params)
    @assert all(Ψcpy.params[s] .≡ Ψ.params[s])
  end

  return Ψcpy
end


#contains the MCMC chain
mutable struct CapacityChains{
    Tchain<:Chains, Tvalue, Tparams<:Dict, Tparamsidx<:Dict,
      Thyper<:Dict, Tidxiter, T}

  chain::Tchain #holds an MCMC chain object
  value::Tvalue
  params::Tparams
  paramsidx::Tparamsidx
  hyper::Thyper

  i::Int #the iteration counter
  previdx::Int
  idxiter::Tidxiter #index iterator, relevant during burn-in

  Nburn::Int
  burncomplete::Bool
  function CapacityChains(;
      chain::Tchain,
      params::Tparams,
      paramsidx::Tparamsidx,
      hyper::Thyper,
      Nburn::Int) where {Tchain, Tparams,Tparamsidx,Thyper}

    value = chain.value
    previdx=0
    idxiter=Iterators.cycle(1:(size(chain.value,1))) |> Iterators.Stateful

    burncomplete=Nburn≤0

    return new{Tchain, typeof(value), Tparams, Tparamsidx, Thyper, typeof(idxiter), eltype(chain.value)}(
      chain, value, params, paramsidx, hyper, 1, previdx, idxiter, Nburn, burncomplete)
  end
end

#main constructor for the MCMC chain
function CapacityChains(::Type{T};
  Nburn::Int = error("Nburn is required"),
  Ndraws::Int = error("Ndraws is required"),
  paramstemplate::Union{NamedTuple,Dict} = error("paramstemplate is required"),
  hyper::Union{Dict,NamedTuple}) where T <: Real

  #total number of parameters in the chain
  Nparams::Int = sum([length(p) for p ∈ values(paramstemplate)])
  paramnames::Vector{String} = Vector{String}(undef, Nparams)
  col::Int = 1 #the column in the values matrix
  paramsidx::Dict{Symbol, Union{Int, Vector{Int}, Matrix{Int}}} = Dict()
  for (i,(n,p)) ∈ enumerate(pairs(paramstemplate))
    nparams = length(p)
    @assert nparams == prod(size(p))

    if typeof(p) <: Real #the parameter is a single value
      paramnames[col] = n|>string
      paramsidx[n] = col
      col += 1

    #a corresponding  vector of column indices for a parameter vector
    elseif typeof(p) <: AbstractVector
      paramsidx[n] = similar(p, Int)
      colend = col+length(p)-1 #the end of the range
      paramsidx[n] .= col:colend

      #assign a string name to each parameter
      paramnames[col:colend] .= (coord->"$n[$coord]").(col:colend)

      col = col + length(p)
      @assert col == colend + 1
    elseif typeof(p) <: AbstractMatrix
      paramsidx[n] = similar(p, Int)
      colend = col+length(p)-1 #the end of the range
      paramsidx[n] .= reshape(col:colend, size(p)) #creates the indices

      #assign a string name to each parameter. The string includes name and index e.g. a[2,3]
      paramnames[col:colend] .= (coord->"$n[$(coord[1]),$(coord[2])]").(
        Iterators.product(1:size(p,1), 1:size(p,2)) |> collect |> vec)
      col = col + length(p)
      @assert col == colend + 1
    else
      error("Unrecognized type of parameter $n: $(typeof(p))")
    end
  end
  @assert Nparams == col - 1

  valuetemplate = rand(T, Ndraws, Nparams, 1)
  chain = Chains(valuetemplate, paramnames)

  #the below views are carefully set to allow linear indexing
  params = Dict(n => view(chain.value,:,vec([p...;]),1) for (n,p) ∈ pairs(paramsidx))

  cleanedhyper::Dict = typeof(hyper) <: Dict ? hyper : (hyper |> Dict)

  return CapacityChains(hyper=cleanedhyper; chain, params, paramsidx, Nburn)
end

#gets the index within the value array
function Base.getproperty(cc::CapacityChains, sym::Symbol)
  (sym === :idx) && return peek(cc.idxiter)
  return getfield(cc, sym)
end

#increments the state of the chain
function increment!(cc::CapacityChains)
  if cc.i == cc.Nburn #if we are at the end of the burning period, reset
    @assert !cc.burncomplete
    cc.burncomplete=true
    resetidx!(cc)
  else
    cc.previdx = cc.idx
    popfirst!(cc.idxiter)
    (cc.burncomplete) && (cc.idx < cc.previdx) && error("Attempted to increment chain past chain length")
  end

  cc.i+=1
  return nothing
end

#updates the chain based on the state
function updateandincrement!(cc::CapacityChains, Ψ::AbstractCCState; verify::Bool=true)

  cc.value[cc.idx, :, 1] .= Ψ.Π
  increment!(cc)

  @unpack previdx = cc #check the integrity of the assignment
  if verify
    for n ∈ keys(cc.params)
      (all(cc.params[n][previdx, :, 1] .== [getproperty(Ψ, n)...;])) || error(
        "!all(cc.params[n][previdx, :, 1] .== [getproperty(Ψ, n)...;])
        n=$n, previdx=$previdx
        cc.params[n][previdx, :, 1]=$(cc.params[n][previdx, :, 1])
        [getproperty(Ψ, n)...;]=$([getproperty(Ψ, n)...;])
        "
      )
    end
  end

  return nothing
end

#a couple of convenience methods
Base.length(cc::CapacityChains) = size(cc.value,1)
Base.size(cc::CapacityChains) = size(cc.value)[1:2]
Base.size(cc::CapacityChains, i::Int) = size(cc.value)[i]

#resets the circular indexing to put values in chronological order
function resetidx!(cc::CapacityChains)
  @unpack value = cc
  value .= @view value[[(cc.idx+1):end...;1:cc.idx...;], :, :]

  #reset the index
  while cc.idx ≠ 1
    popfirst!(cc.idxiter)
  end
  cc.previdx = length(cc)

  return nothing
end

CCState(cc::CapacityChains, ::Type{T}=eltype(cc.value)) where T = CCState(cc.paramsidx, T)

#Base.Broadcast.broadcastable(yp::YearPeriod) = Ref(yp)
#helper function that assigns an initial value
function CCState(cc::CapacityChains, initial::Union{NamedTuple,Dict},
  ::Type{T}=eltype(cc.value)) where T

  Ψ = CCState(cc, T)

  for (n,p) ∈ pairs(initial)
    Ψ.params[n] .= initial[n]
  end

  #integrity check
  for (n,p) ∈ pairs(Ψ.params)
    @assert all(Ψ.params[n] .== initial[n])
  end

  return Ψ
end
