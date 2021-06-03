#=struct TestPort
  Fret::Symbol
  Ffactors::Vector{Symbol}
end=#



const BLAS_MULTITHREAD = 1

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
    portidx::Dict{Union{Int,Symbol}, Int}
    dateidx::Dict{Union{Date,Int}, Int}

    name::Symbol
end

#NOTE: Brackets index the test portfolios
Base.getindex(fs::FactorSpecification, sym::Symbol) = fs.eachport[portidx[sym]]
Base.getindex(fs::FactorSpecification, dt::Union{Date,Int}) = fs.eachdate[dateidx[dt]]

#the below methods could probably be improved
Base.getindex(fs::FactorSpecification, syms::Symbol...) = view(
  fs.df, (s::Symbol-> s ∈ syms).(fs.df[!,fs.Fport]), :)
Base.getindex(fs::FactorSpecification, dts::Union{Date,Int}...) = view(
  fs.df, (d::Union{Date,Int}-> s ∈ dts).(fs.df[!,fs.Fdate]), :)


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
    Fcontrols = Vector{Symbol}(),
    name::Symbol = nothing,
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
  (!isnothing(Fcharacteristics)) && (
    (s->(s ∈ dfnames) || error("Factor $s not found in dataframe!")).(Fcharacteristics))
  (s->(s ∈ dfnames) || error("Control $s not found in dataframe!")).(Fcontrols)
  #println("num missing = $(size(df[(!).(completecases(df[!,[Fret,Fdate,Fport]])),:]))")
  #(sum((!).(completecases(df[!,[Fdate,Fport]]))) == 0) || error("missing values found in identifiers")
  issorted(df, [Fport, Fdate]) || error("DataFrame must be sorted by $Fport and $Fdate")

  if isnothing(bounds) && typeof(df) == DataFrame
    sdf = view(df, :, :)
  elseif !isnothing(bounds)
    sdf = view(df, (d::Any ->
      (d ≥ minimum(bounds)) && (d≥maximum(bounds))).(df[!,Fdate]), :)
  else
    sdf = df
  end

  allfields::Vector{Symbol} = filter(s->!isnothing(s), [Fret; Fdate; Fport; Ffactors; Fcharacteristics; Fcontrols])
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


function rhs(fs::FactorSpecification, focalvars::Union{Vector{Symbol},Nothing};
    intercept::Bool = true)
  local rhs::Union{Symbol, Expr, Nothing}
  local rhsstring::String

  joinvec::Vector{NSymbol} = [focalvars; fs.Fcontrols;]
  joinvec = joinvec[(!isnothing).(joinvec)]

  rhsstring = join((string).(joinvec), " + ")
  (!intercept) && (rhsstring = "$rhsstring + 0")

  rhs = Meta.parse(rhsstring)

  return rhs
end

function rhsnames(fs::FactorSpecification, focalvars::Union{Vector{Symbol},Nothing};
    intercept::Bool = true)::Vector{Symbol}

  local nms::Vector{Union{Symbol,Nothing}}

  if intercept
    nms = [:intercept; focalvars; fs.Fcontrols;]
  else
    nms = [focalvars; fs.Fcontrols;]
  end

  return nms[(!isnothing).(nms)]
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
function timeseriesregressions(fs::FactorSpecification,
  ::Type{T} = FMLM,
  ::Type{M} = Matrix{Float64},
  ::Type{V} = Vector{Float64};
  timeseriesdep::Symbol=fs.Fret,
  timeseriesrhs = factorrhs(fs),
  timeseriesrhsnames = factorrhsnames(fs),
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel::Bool = false,
  qrtype::Type = Matrix{Float64})::Vector{T} where {
    T, M<:AbstractMatrix, V<:AbstractVector}

  Nreg::Int = length(fs.eachport)
  specs::FMSpecs = FMSpecs(Nreg, T, aggfunc=aggfunc)

  timeseriesrhs = factorrhs(fs)
  timeseriesrhsnames = factorrhsnames(fs)

  for i ∈ 1:length(fs.eachport)
    push!(specs, yspec = fs.Fret, xspec = timeseriesrhs, xnames = timeseriesrhsnames)
  end

  if parallel
    BLAS.set_num_threads(1)
  else
    BLAS.set_num_threads(BLAS_MULTITHREAD)
  end

  computeFMLMresults!(fs.eachport, specs, M, V,
    parallel=parallel, containsmissings=false, qrtype = qrtype)

  return specs.results
end

#run time series regressions with an arbitrary aggregation function
function crosssectionalregressions(fs::FactorSpecification, ::Type{T} = FMLM,
    ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel::Bool = false, qrtype::Type = Matrix{Float64})::Vector{T} where {
    T, M<:AbstractMatrix, V<:AbstractVector}

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
  computeFMLMresults!(fs.eachdate, specs, M, V,
    parallel=parallel, containsmissings=false,
    qrtype=qrtype)

  #alldates::Vector{Date} = (df::SubDataFrame->df[1,fs.Fdate]).(fs.eachdate)

  #dateidx::Dict = Dict(zip(alldates, specs.results))

  return specs.results
end

struct FactorRegression
  λ::Union{Nothing, Matrix{Float64}}
  λᵢ::Vector{Float64}
  σᵢ::Vector{Float64}
  t::NInt
  name::Symbol
  coefidx::Dict{Symbol, Int}
  otherstats::Dict{Symbol, Any}
end

FactorRegression( λ::Union{Nothing, Matrix{Float64}}, λᵢ::Vector{Float64},
  σᵢ::Vector{Float64}, t::NInt, name::Symbol,
  coefidx::Dict{Symbol, Int}) = FactorRegression(λ, λᵢ, σᵢ, t, name, coefidx, Dict{Symbol, Any}())


function manyFF(fs::FactorSpecification,
  ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
  parallel::Bool = false,
  errorfunction::Function = neweywestΣfunc(2),
  qrtype::Type = Matrix{Float64})::Vector{FactorRegression} where {
    M<:AbstractMatrix, V<:AbstractVector}

  local frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #make the time-series regression aggregation function
  T::Type = FMLM
  @inline aggfunc(fmlm::FMLM) = fmlm

  results = timeseriesregressions(fs, T, M, V,
    aggfunc=aggfunc, parallel=parallel, qrtype=qrtype)
  coefnames::Vector{Symbol} = factorrhsnames(fs)
  coefidx = Dict{Symbol, Int}(zip(coefnames, 1:(length(coefnames))))

  #now build the factor regression objects, one for each regression
  for (i,port::SubDataFrame) ∈ enumerate(fs.eachport)
    local λᵢ::Vector{Float64} = results[i].β
    local σᵢ::Vector{Float64} = diag(errorfunction(results[i])).^0.5
    local t::Int = size(port,1)
    local name::Symbol = Symbol("$(fs.name)_col$(port[1,fs.Fport])")
    push!(frs, FactorRegression(nothing, λᵢ, σᵢ, t, name, coefidx))
  end

  return frs
end

function characteristicFM(fs::FactorSpecification,
    ::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    parallel::Bool = false, lags::Int=2,
    qrtype::Type = Matrix{Float64})::FactorRegression where {
      M<:AbstractMatrix, V<:AbstractVector}

  local t::Int
  local λ::Matrix{Float64}
  local λᵢ::Vector{Float64}
  local σᵢ::Vector{Float64}
  local name::Symbol
  local coefidx::Dict{Symbol, Int}


  #make the time-series regression aggregation function
  T::Type = Vector{Float64}
  @inline aggfunc(fmlm::FMLM) = fmlm.β

  results = crosssectionalregressions(fs, T, M, V,
    aggfunc=aggfunc, parallel=parallel, qrtype=qrtype)
  coefnames::Vector{Symbol} = characteristicrhsnames(fs)


  coefidx = Dict{Symbol, Int}(zip(coefnames, 1:(length(coefnames))))

  #now build the factor regression object
  λ = vcat((transpose).(results)...) #create rows and make result matrix
  origrows::Int = size(λ,1)
  λ = λ[(i->isfinite(sum(λ[i,:]))).(1:size(λ,1)), :] #checks for na/missing and drops such rows
  newrows::Int = size(λ,1)
  (origrows ≠ newrows) && (println("$(fs.name): removed " *
    "$(origrows - newrows) rows during DQ. $(newrows) remain."))
  t = size(λ,1)
  λᵢ = (c->mean(c)).(eachcol(λ))
  #println("coefnames: $coefnames")
  #println("λᵢ: $λᵢ\n")
  Σ::Matrix{Float64} = FMBartlett(λ, λᵢ, lags)
  σᵢ = diag(Σ).^ 0.5 ./ t^0.5

  #=println("\nFastNW: $σᵢ " *
    #"\nSlowNW: $(diag(FMBartlettSlow(λ, λᵢ, lags)).^0.5 ./ t^0.5)" *
    #"\nFMNeweyWestML: $(diag(FMBartlettML(λ, λᵢ, lags)).^0.5 ./ t^0.5)" *
    "\nNormalSE: $((c->std(c)).(eachcol(λ)) ./ t^0.5)" *
    "\nλᵢ: $λᵢ"
    )=#


  name = fs.name


  return FactorRegression(λ, λᵢ, σᵢ, t, name, coefidx)
end


function panelFM(fs::FactorSpecification,
    ::Type{M}=Matrix{Float64}, ::Type{V}=Vector{Float64};
    errorfunctions::Vector{Function}=[homoskedasticΣ!],
    withinsym::NSymbol = nothing,
    clustersyms::Union{NSymbol, Vector{Symbol}} = [fs.Fdate, fs.Fport],
    qrtype::Type = Matrix{Float64})::Vector{FactorRegression} where {
      M<:AbstractMatrix, V<:AbstractVector}

  local frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #make the time-series regression aggregation function
  T::Type = Vector{Float64}
  @inline aggfunc(fmlm::FMLM) = fmlm

  rhs::FMExpr = characteristicrhs(fs)
  coefnames::Vector{Symbol} = characteristicrhsnames(fs)
  fmlm::FMLM = FMLM(fs.df, rhs, fs.Fret, M, V,
    withinsym=withinsym, clustersyms=clustersyms, Xnames=coefnames,
    containsmissings=false, qrtype=qrtype)
  coefidx = Dict{Symbol, Int}(zip(coefnames, 1:(length(coefnames))))
  local name = fs.name

  local λᵢ::Vector{Float64} = fmlm.β

  for (j,f) ∈ enumerate(errorfunctions)

    local σᵢ::Vector{Float64}

    #print("j=$j time: ")
    Σ::Matrix{Float64} = f(fmlm)
    σᵢ = diag(Σ).^ 0.5

    push!(frs, FactorRegression(nothing, λᵢ, σᵢ, nothing, name, coefidx))
  end



  return frs
end

#run time series regressions with an arbitrary aggregation function
function BJSFFrow(fs::FactorSpecification, df::AbstractDataFrame,
  ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
  parallel::Bool = false,
  errorfunction::Function = neweywestΣfunc(2),
  Ffactors::Vector{Symbol} = [fs.Fret; fs.Ffactors],
  qrtype::Type = Matrix{Float64},
  regressiongroupname::Symbol = df.portname[1])::Vector{FactorRegression} where {
    M<:AbstractMatrix, V<:AbstractVector}

  #this will hold the output
  local frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #iterate through each dependent variable
  for (i, Ysym) ∈ enumerate(Ffactors)
    focals::Vector{Symbol} = setdiff(Ffactors, [Ysym]) #The Y-var is excluded here
    Xexpr = rhs(fs, focals)
    Xnames = rhsnames(fs, focals)

    #run the regression
    fmlm = FMLM(df, Xexpr, Ysym, M, V, Xnames = Xnames, containsmissings=false, qrtype=qrtype)

    #then make the FactorRegression object
    λᵢ::Vector{Float64} = fmlm.β
    σᵢ::Vector{Float64} = diag(errorfunction(fmlm)).^0.5
    t::Int = size(fs.df,1)
    name::Symbol = Symbol("$(regressiongroupname)_col$(Ysym)")
    coefidx::Dict = Dict{Symbol, Int}(zip(Xnames, 1:(length(Xnames))))
    otherstats::Dict = Dict(:R²=>R²(fmlm, adjusted=false))

    push!(frs, FactorRegression(nothing, λᵢ, σᵢ, t, name, coefidx, otherstats))
  end

  return frs
end

function BJSFF(fs::FactorSpecification,
  ::Type{M} = Matrix{Float64}, ::Type{V} = Vector{Float64};
  parallel::Bool = false,
  errorfunction::Function = neweywestΣfunc(2),
  qrtype::Type = Matrix{Float64})::Dict{Symbol, FactorRegression} where {
    M<:AbstractMatrix, V<:AbstractVector}

  #this will hold the output
  results::Dict = Dict{Symbol, FactorRegression}()

  #first build the baseline
  baselinedf::DataFrame = unique(fs.df[!, [:date; fs.Ffactors;]])
  baselinefrs::Vector{FactorRegression} = BJSFFrow(fs, baselinedf, Ffactors=fs.Ffactors, regressiongroupname=:baseline)
  for fr ∈ baselinefrs
    results[fr.name] = fr
  end

  #then build the rest
  for sdf ∈ fs.eachport
    frs::Vector{FactorRegression} = BJSFFrow(fs, sdf)
    for fr ∈ frs
      results[fr.name] = fr
    end
  end

  return results
end
