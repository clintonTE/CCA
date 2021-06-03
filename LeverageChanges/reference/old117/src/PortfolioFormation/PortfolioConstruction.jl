
#############################Move the below elsewhere

struct PortfolioConstruction{T}
  Fval::Symbol #the field that holds the values
  Fcontrol::NSymbol #holds the control if one is used
  Fweighton::NSymbol #holds the variable that can be weighted against (e.g. market cap)

  Fgroup::Symbol #the field that holds the group values
  Fn::Symbol

  valtype::Type
  breakpoints::T #holds a tuple of breakpints
  groupweights::T #holds a tuple of weighting coefficents (after market cap normalization)

  Ngroups::Int
  name::Symbol
end

#main constructor
#most items are pulled from constants
function PortfolioConstruction(::Type{T}, Fval::Symbol;
  name::Symbol = Fval,
  Fcontrol::NSymbol=nothing,
  Fweighton::NSymbol = nothing,
  breakpoints::T = nothing,
  groupweights::T = nothing,
  valtype::Type = nothing) where T

  Fgroup = Symbol("G_", name)
  Fnval = Symbol("N_", name)

  return PortfolioConstruction{T}(Fval, Fcontrol, Fweighton,
    Fgroup, Fnval,
    valtype, breakpoints, groupweights,
    length(breakpoints), name)
end



#########################
const CQuartileValType = NTuple{4,Float64}
const CQuartileConstruction = PortfolioConstruction{CQuartileValType}

const QUARTILE_BREAKPOINTS = (0.25,0.5,0.75,1.0)
const QUARTILE_GROUP_WEIGHTS = (-1.0,0.0,0.0,1.0)

#simplified constructor to make quartiles
function CQuartileConstruction(Fval::Symbol, Fcontrol::Symbol;
  Fweighton::NSymbol=nothing, name::Symbol = Fval)

  return PortfolioConstruction(CQuartileValType, Fval, Fcontrol=Fcontrol, Fweighton=Fweighton,
    breakpoints=QUARTILE_BREAKPOINTS,
    groupweights=QUARTILE_GROUP_WEIGHTS,
    valtype=Float64, name=name)
end

#define eltype for portfolio construction
Base.eltype(::Type{PortfolioConstruction{T where T}}) = T
Base.eltype(pc::PortfolioConstruction) = eltype(typeof(pc))

#########################
# Portfolio Construction
##########################

#create a container structure for all of the portfolio constructions
struct PortfolioConstructions
  pcs::Vector{PortfolioConstruction{<:NTuple}}
  idx::Dict{Symbol, PortfolioConstruction}
  N::Int
end

#main constructor
function PortfolioConstructions(::Type{T},
  Fvals::Vector{Symbol};
  Fcontrols::Union{NSymbol, Vector{NSymbol}} = nothing,
  Fweightons::Union{NSymbol, Vector{NSymbol}} = nothing,
  breakpoints::T = nothing,
  groupweights::T = nothing,
  valtypes::Union{Type, Vector{Type}} = nothing,
  namevec::Vector{Symbol} = deepcopy(Fvals)) where T

  #kinda ugly, but its a compact way to define the portfolio constructions
  #nameing is a bit incosnsitent to avoid ambigueties
  pcs::Vector{PortfolioConstruction} = (
    (Fval,Fcontrol,Fweighton,valtype,name) ->
      PortfolioConstruction(T, Fval,
        Fcontrol=Fcontrol, Fweighton=Fweighton, breakpoints=breakpoints,
        groupweights=groupweights, valtype=valtype, name=name)).(
    Fvals, Fcontrols, Fweightons, valtypes, namevec)

  idx::Dict = Dict(pc.name=>pc for pc ∈ pcs)

  return PortfolioConstructions(pcs, idx, length(pcs))
end

const CQuartileConstructions = PortfolioConstructions
function CQuartileConstructions(Fvals::Vector{Symbol}, Fcontrols::Union{Symbol, Vector{Symbol}};
  Fweightons::Union{NSymbol, Vector{NSymbol}}=nothing,
  namevec::Vector{Symbol} = deepcopy(Fvals))

  PortfolioConstructions(CQuartileValType, Fvals,
    Fcontrols=Fcontrols,
    Fweightons=Fweightons,
    breakpoints=QUARTILE_BREAKPOINTS,
    groupweights=QUARTILE_GROUP_WEIGHTS,
    valtypes=Float64, namevec=namevec)
end

function CQuartileConstructions(d::Dict)
  Fvals::Vector{Symbol} = d[:Fvals]
  Fcontrols::Union{Symbol, Vector{Symbol}} = d[:Fcontrols]
  Fweightons::Union{NSymbol, Vector{NSymbol}}= haskey(d, :Fweightons) ? d[:Fweightons] : nothing
  namevec::Vector{Symbol} = haskey(d, :namevec) ? d[:namevec] : nothing

  PortfolioConstructions(CQuartileValType, Fvals,
    Fcontrols=Fcontrols,
    Fweightons=Fweightons,
    breakpoints=QUARTILE_BREAKPOINTS,
    groupweights=QUARTILE_GROUP_WEIGHTS,
    valtypes=Float64,
    namevec=namevec)
end

function Base.merge(pcs1::PortfolioConstructions, pcs2::PortfolioConstructions)
  pcs::Vector{PortfolioConstruction} = [pcs1.pcs; pcs2.pcs;]
  idx = merge(pcs1.idx, pcs2.idx)
  N::Int = pcs1.N + pcs2.N

  return PortfolioConstructions(pcs, idx, N)
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
Fweightons(pcs::PortfolioConstructions)::Vector{NSymbol} = (
  pc::PortfolioConstruction->pc.Fweighton).(pcs)
Fcontrols(pcs::PortfolioConstructions)::Vector{NSymbol} = (
  pc::PortfolioConstruction->pc.Fcontrol).(pcs)
Fgroups(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fgroup).(pcs)
Fweights(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fweight).(pcs)
Fns(pcs::PortfolioConstructions)::Vector{Symbol} = (
  pc::PortfolioConstruction->pc.Fn).(pcs)

function IdentifyPortfolios(pcs::PortfolioConstructions)
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
