struct FactorSpec
    df::SubDataFrame #the data used
    Fdate::Symbol #date or period field
    Ffactors::Vector{Symbol} #factor values
    Ftest::Vector{Symbol} #test portfolios
    Θ::Dict{Symbol, Any} #will hold other parameters relating to FM
end

#a simple constructor
function FactorSpec(df::D, Fdate::Symbol, Ffactors::Vector{Symbol},
    Ftest::Vector{Symbol}; bounds::Union{AbstractRange, Nothing} = nothing,
  Θ::Dict{Symbol, Any} = Dict{Symbol, Any}()
  ) where {D<:AbstractDataFrame}

  local sdf::SubDataFrame
  local dfnames::Vector{Symbol} = names(df)

  #make sure the dataframe has everything we want
  (Fdate ∉ dfnames) && error("$Fdate not found in dataframe!")
  (s->(s ∉ dfnames) && error("Factor $s not found in dataframe!")).(Ffactors)
  (s->(s ∉ Ftest) && error("Test portfolio $s not found in dataframe!")).(Ftest)

  if isnothing(bounds)
    sdf = view(df, :, :)
  else
    sdf = view(df, (d::MDate-> (!ismissing(d)) &&
      (d ≥ minimum(bounds)) && (d≥maximum(bounds))).(df[!,Fdate]), :)
  end

  sdf = view(sdf, completecases(sdf), [Fdate; Ffactors; Ftest])


  return FactorSpec(sdf, Fdate, Ffactors, Ftest, Θ)
end

#creates the formula
function rhs(fs::FactorSpec; intercept::Bool = true)
  local rhs::Union{Symbol, Expr}
  local rhsstring::String

  rhsstring = join((string).(fs.Ffactors), " + ")
  (!intercept) && (rhsstring = "$rhsstring + 0")

  rhs = Meta.parse(rhsstring)

  return rhs
end

function rhsnames(fs::FactorSpec; intercept::Bool = true)

  local nms::Vector{Symbol}

  if intercept
    nms = [:intercept; fs.Ffactors;]
  else
    nms = fs.Ffactors
  end

  return fs.Ffactors
end

#run time series regressions with an arbitrary aggregation function
function timeseriesregressions(fs::FactorSpec, ::Type{T} = FMLM;
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel = false)::Vector{T} where T

  Nreg::Int = length(fs.Ftest)
  specs::FMSpecs = FMSpecs(Nreg, T, aggfunc=aggfunc)

  for p ∈ fs.Ftest
    push!(specs, yspec = p, xspec = rhs(fs), xnames = rhsnames(fs))
  end


  #println("Nreg: $Nreg, num specs: $(length(specs))")

  computeFMLMresults!(fs.df, specs, parallel=parallel,
    containsmissings=false)

  return specs.results
end
