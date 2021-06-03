
struct Characteristic <: AbstractCharacteristic
  ws::WeightSpec
  Fξ::Symbol
  X::Dict{Symbol, Symbol}
end


#builds a characteristic object, which includes cross-sectional weights and
#weight coefficients
function Characteristic(ws::WeightSpec, Fξ::Symbol)

  Fw = Symbol(:w, ws.weightprefix, Fξ)

  #a couple of helper funtions for consistent renaming
  lagname(s::Symbol) = Symbol(:L, s)
  retname(s::Symbol) = Symbol(:R, s)
  rettname(s::Symbol) = Symbol(:Rt, s)

  #WARNING: HARDCODE ALERT!!!
  X = Dict{Symbol, Symbol}()
  X[:Fw] = Fw
  X[:FLw] = lagname(X[:Fw])
  X[:FLLw] = lagname(X[:FLw])
  X[:FLLLw] = lagname(X[:FLLw])
  #not used for now
  X[:FRLw] = retname(X[:FLw])
  X[:FRLLw] = retname(X[:FLLw])

  X[:FRtw] = rettname(X[:Fw])
  return Characteristic(ws, Fξ, X)
end

weightfields(ξ::Characteristic) = collect(values(ξ.X))
Characteristic(ξ::Characteristic,  Fξ::Symbol, X::Dict
  ) = Characteristic(ξ.ws, Fξ, X)

demeanξ(ξraw::AbstractVector{<:Union{Real, Missing}}) = ξraw .- mean(skipmissing(ξraw))
#this gives the quantile fo each value in an array
#would have preferred the ecdf function, but it doesn't play nice with
#containers that can hold missing values, even if none are actually held
#=function quantiles(ξraw::AbstractVector{<:Union{Real,Missing}})
  (length(ξraw) .- competerank(ξraw,rev=true) .+ 1)./length(ξraw)
end=#

function quantiles(ξraw::AbstractVector{<:Union{Real,Missing}})
  (length(ξraw) .- competerank(ξraw,rev=true) .+ 1)./length(ξraw)
end

#alternate procedure
#=function quantiles(ξraw::AbstractVector{<:Union{Real,Missing}})
  F = ecdf(ξraw |> Vector{Float64})
  return map(F, ξraw)
end=#


#conditions each slice of the panel for a particular characteristic
function conditionξ!(F::Function, panel::AbstractDataFrame, Fξraw::Symbol, Fξ::Symbol)::Nothing

  (sum((!ismissing).(panel[!, Fξ])) > 0) && error("$Fξ in panel already contains data")
  spanels::GroupedDataFrame = groupby(view(panel, (!ismissing).(panel[!,Fξraw]), :), :date)
  Threads.@threads for i ∈ 1:length(spanels)
      spanel::SubDataFrame = spanels[i]
      Fξsubcol::SubArray = spanel[!, Fξ]
      Fξsubcol .= F(spanel[!, Fξraw])
  end

  return nothing
end

#=standard conditioning of a characteristic column
function conditionξ!(panel::AbstractDataFrame, Fξraw::Symbol;
  prefix::Symbol = PARAM[:charprefix],
  Fξ::Symbol = Symbol(prefix, Fξraw))::Nothing

  #allocate the data
  ξtype::Type = Union{eltype(panel[!,Fξraw]), Missing}
  panel[!, Fξ] = missings(ξtype, size(panel,1))

  conditionξ!(demeanξ, panel, Fξraw, Fξ)
  return nothing
end

#apply a quantile transformation before demeaning
function conditionξquantile!(panel::AbstractDataFrame, Fξraw::Symbol;
  prefix::Symbol = PARAM[:charprefix],
  Fξ::Symbol = Symbol(prefix, Fξraw))::Nothing

  #allocate the data
  panel[!, Fξ] = Vector{Union{eltype(panel[!,Fξraw])}}(undef, size(panel,1))

  conditionξ!(demeanξ ∘ quantiles, panel, Fξraw, Fξ)
  return nothing
end=#
