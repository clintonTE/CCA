#=struct TestPort
  Fret::Symbol
  Ffactors::Vector{Symbol}
end=#

const BLAS_MULTITHREAD = 6

struct FactorSpecification
    df::SubDataFrame #the data used
    Fret::Symbol
    Fdate::Symbol #date or period field
    Fport::Symbol

    #these are considered factors, generally used in the time series regressions
    Ffactors::Union{Nothing,Vector{Symbol}}

    #These are chracteristics, typically used for the cross-sectional regressions
    Fcharacteristics::Union{Nothing,Vector{Symbol}}

    Fcontrols::Vector{Symbol} #holds controls only

    #accessor variables
    eachport::GroupedDataFrame
    eachdate::GroupedDataFrame
    portidx::Dict{Int, Int}
    dateidx::Dict{Date, Int}

    name::NSymbol
end

#NOTE: Brackets index the test portfolios
Base.getindex(fs::FactorSpecification, sym::Symbol) = fs.eachport[portidx[sym]]
Base.getindex(fs::FactorSpecification, dt::Date) = fs.eachdate[dateidx[dt]]

#the below methods could probably be improved
Base.getindex(fs::FactorSpecification, syms::Symbol...) = view(
  fs.df, (s::Symbol-> s ∈ syms).(fs.df[!,fs.Fport]), :)
Base.getindex(fs::FactorSpecification, dts::Date...) = view(
  fs.df, (d::Date-> s ∈ dts).(fs.df[!,fs.Fdate]), :)


@inline function notsmallgroups(df::AbstractDataFrame, Fgroup::Symbol;
  minrowspergroup::Int = error("minrowspergroup must be assigned"))::Vector{Bool}

  dfcount::DataFrame = DataFrame(group = df[!,Fgroup])
  dfcount.count = Vector{MInt}(undef, size(dfcount,1))

  for sdfcount ∈ groupby(dfcount, :group)
    sdfcount.count .= size(sdfcount,1)
  end

  return dfcount.count .≥ minrowspergroup
end

#main constructor
function FactorSpecification(df::AbstractDataFrame, Fret::Symbol, Fdate::Symbol, Fport::Symbol;
    Ffactors::Union{Nothing, Vector{Symbol}}=nothing,
    Fcharacteristics::Union{Nothing, Vector{Symbol}}=nothing,
    bounds::Union{AbstractRange, Nothing} = nothing,
    Fcontrols = Vector{Symbol}(), name::NSymbol = nothing,
    minrowspergroup::Int = MIN_ROWS_PER_GROUP,
    skipcrosssectiondfcheck::Bool = false,
    skiptimeseriesdfcheck::Bool = false
    )

  local sdf::SubDataFrame
  local dfnames::Vector{Symbol} = names(df)

  #make sure the dataframe has everything we want
  (Fret ∈ dfnames) || error("return field $Fret not found in dataframe!")
  (Fdate ∈ dfnames) || error("datefield $Fdate not found in dataframe!")
  (Fport ∈ dfnames) || error("port field $Fport not found in dataframe!")
  (!isnothing(Ffactors)) && (
    (s->(s ∈ dfnames) || error("Factor $s not found in dataframe!")).(Ffactors))
  (!isnothing(Ffactors)) && (
    (s->(s ∈ dfnames) || error("Factor $s not found in dataframe!")).(Fcharacteristics))
  (s->(s ∈ dfnames) || error("Control $s not found in dataframe!")).(Fcontrols)
  #println("num missing = $(size(df[(!).(completecases(df[!,[Fret,Fdate,Fport]])),:]))")
  (sum((!).(completecases(df[!,[Fret,Fdate,Fport]]))) == 0) || error("missing values found in identifiers")

  if isnothing(bounds) && typeof(df) == DataFrame
    sdf = view(df, :, :)
  else
    sdf = view(df, (d::MDate->
      (d ≥ minimum(bounds)) && (d≥maximum(bounds))).(df[!,Fdate]), :)
  end

  allfields = filter(s->!isnothing(s), [Fret; Fdate; Fport; Ffactors; Fcharacteristics; Fcontrols])
  sdf = view(sdf, completecases(sdf[!, allfields]), allfields)

  skiptimeseriesdfcheck || (sdf = view(sdf, notsmallgroups(sdf, Fport, minrowspergroup=minrowspergroup), :))
  skipcrosssectiondfcheck || (sdf = view(sdf, notsmallgroups(sdf, Fdate, minrowspergroup=minrowspergroup), :))


  (size(sdf,1)==0) && error("no rows left with reasonable sample size")

  #build subarrays of the most comon slices we will take for analysis
  eachport = groupby(sdf, Fport)
  eachdate = groupby(sdf, Fdate)
  portidx::Dict = Dict(eachport[i][1, Fport] => i for i ∈ 1:length(eachport))
  dateidx::Dict = Dict(eachdate[i][1, Fdate] => i for i ∈ 1:length(eachdate))

  return FactorSpecification(sdf, Fret, Fdate, Fport, Ffactors, Fcharacteristics, Fcontrols,
    eachport, eachdate, portidx, dateidx, name)
end


function rhs(fs::FactorSpecification, focalvars::Vector{Symbol}; intercept::Bool = true)
  local rhs::Union{Symbol, Expr}
  local rhsstring::String

  rhsstring = join((string).([focalvars; fs.Fcontrols;]), " + ")
  (!intercept) && (rhsstring = "$rhsstring + 0")

  rhs = Meta.parse(rhsstring)

  return rhs
end

function rhsnames(fs::FactorSpecification, focalvars::Vector{Symbol}; intercept::Bool = true)
  local nms::Vector{Symbol}

  if intercept
    nms = [:intercept; focalvars; fs.Fcontrols;]
  else
    nms = [focalvars; fs.Fcontrols;]
  end

  return nms
end

#creates the formulas
factorrhs(fs::FactorSpecification; intercept::Bool = true)::FMExpr = rhs(
  fs, fs.Ffactors, intercept=intercept)
factorrhsnames(fs::FactorSpecification; intercept::Bool = true)::Vector{Symbol} = rhsnames(
  fs, fs.Ffactors, intercept=intercept)
characteristicrhs(fs::FactorSpecification; intercept::Bool = true)::FMExpr = rhs(
  fs, fs.Fcharacteristics, intercept=intercept)
characteristicrhsnames(fs::FactorSpecification; intercept::Bool = true)::Vector{Symbol} = rhsnames(
  fs, fs.Fcharacteristics, intercept=intercept)

#run time series regressions with an arbitrary aggregation function
function timeseriesregressions(fs::FactorSpecification, ::Type{T} = FMLM;
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel::Bool = false)::Vector{T} where T

  Nreg::Int = length(fs.eachport)
  specs::FMSpecs = FMSpecs(Nreg, T, aggfunc=aggfunc)

  timeseriesrhs = factorrhs()
  timeseriesrhsnames = factorrhsnames()

  for i ∈ 1:length(fs.eachport)
    push!(specs, yspec = fs.Fret, xspec = timeseriesrhs, xnames = timeseriesrhsnames)
  end

  if parallel
    BLAS.set_num_threads(1)
  else
    BLAS.set_num_threads(BLAS_MULTITHREAD)
  end

  computeFMLMresults!(fs.eachport, specs, parallel=parallel, containsmissings=false)

  return specs.results
end

#run time series regressions with an arbitrary aggregation function
function crosssectionalregressions(fs::FactorSpecification, ::Type{T} = FMLM;
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel::Bool = false)::Vector{T} where T

  Nreg::Int = length(fs.eachport)
  specs::FMSpecs = FMSpecs(Nreg, T, aggfunc=aggfunc)

  crosssectionalrhs = characteristicrhs(fs)
  crosssectionalrhsnames = characteristicrhsnames(fs)

  for i ∈ 1:length(fs.eachdate)
    push!(specs, yspec = fs.Fret, xspec = crosssectionalrhs, xnames = crosssectionalrhsnames)
  end

  if parallel
    BLAS.set_num_threads(1)
  else
    BLAS.set_num_threads(BLAS_MULTITHREAD)
  end
  computeFMLMresults!(fs.eachdate, specs, parallel=parallel, containsmissings=false)

  return specs.results
end

struct FactorRegression
  λ::DataFrame
  λₜ::Vector{Float64}
  σₜ::Vector{Float64}
  t::Int
end



function characteristicFM(fs::FactorSpecification; parallel::Bool = false)::FactorRegression
  local t::Int
  local mλ::Matrix{Float64}
  local λ::DataFrame
  local λₜ::Vector{Float64}
  local σₜ::Vector{Float64}

  ##STEP 1: Time Series Regressions
  ########################
  #labels::Vector{Symbol} = [:intercept; fmlm.Ffactors]

  #make the time-series regression aggregation function
  T::Type = Vector{Float64}
  @inline aggfunc(fmlm::FMLM) = fmlm.β

  results = crosssectionalregressions(fs, T, aggfunc=aggfunc, parallel=parallel)
  coefnames::Vector{Symbol} = characteristicrhsnames(fs)

  #now build the factor regression object
  try
    t = length(results)
    mλ = vcat((transpose).(results)...) #create rows and make result matrix
    λ = DataFrame(mλ, coefnames)
    λₜ = (c->mean(c)).(eachcol(λ)) ./ t
    σₜ = (c->std(c)).(eachcol(λ)) ./ t^0.5
  catch err
    println(coefnames)
    println(results)
    println(fs.df[1:5,:])
    error(err)
  end

  return FactorRegression(λ, λₜ, σₜ, t)
end
