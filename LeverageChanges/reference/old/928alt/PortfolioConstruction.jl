


#this will hold all attributes relating to portfolio construction
struct PortfolioConstruction{T}
  Fval::Symbol #the field that holds the values
  Fnval::Symbol #this field will hold a count in the final portfolio
  Fquantile::Symbol #the field that holds the quantiles
  Fgroup::Symbol #the field that holds the group values
  Fweight::Symbol#the field that holds the weights

  valtype::Type
  breakpoints::T #holds a tuple of breakpints
  groupweights::T #holds a tuple of weighting coefficents (after market cap normalization)

  Ngroups::Int #total number of groups
end

#main constructor
#most items are pulled from constants
function PortfolioConstruction(::Type{T}, Fval::Symbol;
  factorbreakpoints::Dict = FACTOR_BREAKPOINTS,
  factorgroupweights::Dict = FACTOR_GROUP_WEIGHTS,
  factorvaluetypes::Dict = FACTOR_VALUE_TYPES) where T

  Fnval = Symbol("N_", Fval)
  Fquantile = Symbol(:Q, Fval)
  Fgroup = Symbol(:G, Fval)
  Fweight = Symbol(:W, Fval)

  valtype::Type = factorvaluetypes[Fval]
  breakpoints::T = factorbreakpoints[Fval]
  groupweights::T = factorgroupweights[Fval]

  return PortfolioConstruction{T}(Fval, Fnval, Fquantile, Fgroup,
    Fweight, valtype, breakpoints, groupweights, length(breakpoints))
end

#constructor which does not require a type
function PortfolioConstruction(Fval::Symbol;
  factorngroups::Dict = FACTOR_N_GROUPS)

  T::Type = NTuple{factorngroups[Fval], Float64}

  return PortfolioConstruction(T, Fval)
end

#define eltype for portfolio construction
Base.eltype(::Type{PortfolioConstruction{T where T}}) = T
Base.eltype(pc::PortfolioConstruction) = eltype(typeof(pc))

#########################
# Portfolio Construction
##########################

#create a container structure for all of the portfolio constructions
struct PortfolioConstructions
  pcs::Vector{PortfolioConstruction}
  idx::Dict{Symbol, PortfolioConstruction}
  N::Int
end

#main constructor
function PortfolioConstructions(factorfields::Vector{Symbol} = FACTOR_FIELDS)
  pcs::Vector{PortfolioConstruction} = (PortfolioConstruction).(factorfields)
  idx::Dict = Dict(pc.Fval=>pc for pc ∈ pcs)

  return PortfolioConstructions(pcs, idx, length(pcs))
end


#make indexing work
Base.getindex(pcs::PortfolioConstructions, syms::Symbol...) = (s::Symbol->
  pcs.idx[s]).(syms)
Base.getindex(pcs::PortfolioConstructions, sym::Symbol) = pcs.idx[sym]
Base.getindex(pcs::PortfolioConstructions, indices::Int...) = pcs.pcs[indices]
Base.getindex(pcs::PortfolioConstructions, indice::Int) = pcs.pcs[indice]

#make length work
Base.length(pcs::PortfolioConstructions)::Int = length(pcs.pcs)

#make iteration work
function Base.iterate(pcs::PortfolioConstructions, state::Int = 1)
  (state > pcs.N) && return nothing
  return (pcs[state], state + 1)
end


Fvals(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fval).(pcs)
Fquantiles(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fquantile).(pcs)
Fgroups(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fgroup).(pcs)
Fweights(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fweight).(pcs)

function computesorts(pcs::PortfolioConstructions)
  Nsorts::Int = prod((pc->pc.Ngroups).(pcs))

  sorts::Vector{NTuple{pcs.N, Int}} = Vector{NTuple{pcs.N, Int}}(undef, Nsorts)
  sortlabels::Vector{Vector{Int}} = (pc->collect(1:pc.Ngroups)).(pcs)
  println(sortlabels)

  ctr::Int = 0
  for p ∈ Iterators.product(sortlabels...)
    ctr+=1
    sorts[ctr] = NTuple{pcs.N, Int}(p)
  end

  return sorts
end
