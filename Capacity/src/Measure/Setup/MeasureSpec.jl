
abstract type AbstractMeasureSpec end

#holds the parameters for forming the weights
#not meant to hold any meaningful data- just parameters
struct MeasureSpec{TFvol<:NSymbol, TFDvol<:NSymbol} <: AbstractMeasureSpec
  Fξs::Vector{Symbol}
  ξs::Vector{Characteristic}
  ξidx::Dict{Symbol, Characteristic}
  Fvol::TFvol
  FDvol::TFDvol

  FW::Vector{Symbol} #these are global controls (use case: intercept, controls by k or n)
  FWbyt::Vector{Symbol} #these are controls grouped by t only (use case: time fixed effects)
  FWbykt::Dict{Symbol, Vector{Symbol}} #these are controls grouped by t AND k



  K::Int
end


Fweights(ms::MeasureSpec, s::Symbol) = (ξ->ξ.X[s]).(ms.ξs)
Fgroupedcontrols(ms::MeasureSpec) = [ms.FWbyt...;(Fξ->ms.FWbykt[Fξ]).(ms.Fξs)...;]
Fcontrols(ms::MeasureSpec) = [ms.FW...; Fgroupedcontrols(ms)...;]
grossleverages(ms::MeasureSpec) = (ξ->ξ.ws.gross).(ms.ξs)

#creates the control names as needed
function getcontrolnames(Fξs::Vector{Symbol};
    controlsprefix::Symbol = PARAM[:controlsprefix],
    controlsglobal::Vector{Symbol} = PARAM[:controlsglobal],
    controlsbykt::Vector{Symbol} = PARAM[:controlsbykt],
    controlsbyt::Vector{Symbol} = PARAM[:controlsbyt])


  FW = (s->Symbol(controlsprefix,s)).(controlsglobal)
  FWbyt = (s->Symbol(controlsprefix,s)).(controlsbyt)
  FWbykt::Dict{Symbol,Vector{Symbol}} = Dict{Symbol,Vector{Symbol}}()
  for Fξ ∈ Fξs
    FWbykt[Fξ] = [Symbol(controlsprefix, s, Fξ) for s ∈ controlsbykt]
  end


  ntout = (;FW=FW, FWbyt=FWbyt, FWbykt=FWbykt)
  #println("ntout: $(ntout), typeof: $(typeof(ntout))")
  return ntout
end

function MeasureSpec(ξstuple::Characteristic...;
  Fvol::NSymbol = PARAM[:Fvol],
  FDvol::NSymbol = nothing, #set the keyword to Symbol(:D, Fvol) to enable
  ξs::Vector{Characteristic} = collect(ξstuple),
  Fξs::Vector{Symbol} = (ξ->ξ.Fξ).(ξs),
  controlnames::NamedTuple=getcontrolnames(Fξs)
  )

  ξidx::Dict{Symbol, Characteristic} = Dict(zip(Fξs, ξs))

  #permnodate::Vector{Tuple{Int, Date}} = (Tuple).(zip(panel.permno, panel.date))
  #println("controlnames: $(controlnames), typeof: $(typeof(controlnames))")

  return MeasureSpec(Fξs, ξs, ξidx, Fvol, FDvol,
    controlnames.FW, controlnames.FWbyt, controlnames.FWbykt, length(Fξs))
end

#set up a get index notation
Base.getindex(ms::MeasureSpec, Fξ::Symbol) = ms.ξidx[Fξ]
Base.getindex(ms::MeasureSpec, Fξs::Symbol...) = (Fξ->ms.idx[Fξ]).(Fξs)
Base.getindex(ms::MeasureSpec, i::Int) = ms.ξiBase.dx[i]
Base.getindex(ms::MeasureSpec, is::Int...) = ms.ξidx[is]

Base.length(ms::MeasureSpec) = length(Fξs)

#make iteration work
function Base.iterate(ms::MeasureSpec, state::Int = 1)
  (state > ms.K) && return nothing
  return (ms[state], state + 1)
end

#some convenience methods
getT(panel::AbstractDataFrame, ms::MeasureSpec) = length(unique(panel.date))
getN(panel::AbstractDataFrame, ms::MeasureSpec) = length(unique(panel.permno))
getdims(panel::AbstractDataFrame, ms::MeasureSpec) = (
  T=getT(panel, ms), N = getN(panel, ms), K = ms.K)

function gettidx(panel::AbstractDataFrame, ms::MeasureSpec)
  date::Vector{Date} = sort(unique(panel.date))
  tidx::Dict{Date,Int} = Dict(date[t]=>t for t ∈ 1:length(date))

  return tidx
end



#=makes sure the MeasureSpec matches the ordering in the dataframe
function verifymeasureintegrity(panel::DataFrame, ms::MeasureSpec)
  Nrows::Int = size(panel, 1)

  permno::Vector{Int} = (t->t[1]).(ms.permnodate)
  date::Vector{Date} = (t->t[2]).(ms.permnodate)

  (sum(permno .== panel.permno) == Nrows) || error("permno integrity check failed")
  (sum(date .== panel.date) == Nrows) || error("date integrity check failed")
end=#
