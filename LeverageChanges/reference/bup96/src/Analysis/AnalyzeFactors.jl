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
    dateidx::Dict{Date, Int}

    name::Symbol
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
    Fcontrols = Vector{Symbol}(),
    name::NSymbol = nothing,
    minrowspergroup::Int = MIN_ROWS_PER_GROUP,
    skipcrosssectiondfcheck::Bool = false,
    skiptimeseriesdfcheck::Bool = false
    )

  local sdf::SubDataFrame
  local dfnames::Vector{Symbol} = names(df)

  #make sure the dataframe has everything we want
  isnothing(name) && error("name is a required key-word field!")
  (Fret ∈ dfnames) || error("return field $Fret not found in dataframe!")
  (Fdate ∈ dfnames) || error("datefield $Fdate not found in dataframe!")
  (Fport ∈ dfnames) || error("port field $Fport not found in dataframe!")
  (!isnothing(Ffactors)) && (
    (s->(s ∈ dfnames) || error("Factor $s not found in dataframe!")).(Ffactors))
  (!isnothing(Fcharacteristics)) && (
    (s->(s ∈ dfnames) || error("Factor $s not found in dataframe!")).(Fcharacteristics))
  (s->(s ∈ dfnames) || error("Control $s not found in dataframe!")).(Fcontrols)
  #println("num missing = $(size(df[(!).(completecases(df[!,[Fret,Fdate,Fport]])),:]))")
  (sum((!).(completecases(df[!,[Fret,Fdate,Fport]]))) == 0) || error("missing values found in identifiers")
  issorted(df, [Fport, Fdate]) || error("DataFrame must be sorted by $Fport and $Fdate")

  if isnothing(bounds) && typeof(df) == DataFrame
    sdf = view(df, :, :)
  elseif !isnothing(bounds)
    sdf = view(df, (d::MDate->
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
function timeseriesregressions(fs::FactorSpecification, ::Type{T} = FMLM;
  aggfunc::Function = (f(m::FMLM)::T=m),
  parallel::Bool = false)::Vector{T} where T

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

  #alldates::Vector{Date} = (df::SubDataFrame->df[1,fs.Fdate]).(fs.eachdate)

  #dateidx::Dict = Dict(zip(alldates, specs.results))

  return specs.results
end

struct FactorRegression
  λ::Union{Nothing, Matrix{Float64}}
  λᵢ::Vector{Float64}
  σᵢ::Vector{Float64}
  t::Int
  name::Symbol
  coefidx::Dict{Symbol, Int}
end

function manyFF(fs::FactorSpecification;
  parallel::Bool = false,
  errorfunction::Function = getNeweyWestFunc(2))::Vector{FactorRegression}

  local frs::Vector{FactorRegression} = Vector{FactorRegression}()

  #make the time-series regression aggregation function
  T::Type = FMLM
  @inline aggfunc(fmlm::FMLM) = fmlm

  results = timeseriesregressions(fs, T, aggfunc=aggfunc, parallel=parallel)
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

function returnfactorFF(fs::FactorSpecification; parallel::Bool = false)::FactorRegression
  local t::Int
  local λ::Matrix{Float64}
  local λᵢ::Vector{Float64}
  local σᵢ::Vector{Float64}
  local name::Symbol
  local coefidx::Dict{Symbol, Int}

  #make the time-series regression aggregation function
  T::Type = Vector{Float64}
  @inline aggfunc(fmlm::FMLM) = fmlm.β

  results = timeseriesregressions(fs, T, aggfunc=aggfunc, parallel=parallel)
  coefnames::Vector{Symbol} = factorrhsnames(fs)
  coefidx = Dict{Symbol, Int}(zip(coefnames, 1:(length(coefnames))))

  #now build the factor regression object
  mλ = vcat((transpose).(results)...) #create rows and make result matrix
  origrows::Int = size(mλ,1)
  mλ = mλ[(i->isfinite(sum(mλ[i,:]))).(1:size(mλ,1)), :] #checks for na/missing and drops such rows
  newrows::Int = size(mλ,1)
  (origrows ≠ newrows) && (println("$(fs.name): removed " *
    "$(origrows - newrows) rows during DQ. $(newrows) remain."))
  #λ = DataFrame(mλ, coefnames)
  t = size(λ,1)
  λᵢ = (c->mean(c)).(eachcol(λ))
  σᵢ = (c->std(c)).(eachcol(λ)) ./ t^0.5
  name = fs.name

  return FactorRegression(mλ, λᵢ, σᵢ, t, name, coefidx)
end

#NOTE: The problem is the garbage data is in the df. Clean the df, clean the SEs.
function FMNeweyWest(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}

  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]
  K::Int = length(λᵢ)
  R::Matrix{Float64} = qr(X).R
  Rinv::Matrix{Float64} = (R)\Matrix{Float64}(I,K,K) #this will save some computational time
  Sₜ::Matrix{Float64} = zeros(K,K) #holds teh central matrix

  #pre-allocate for the spectral matrix
  Σ::Matrix{Float64} = Matrix{Float64}(undef, K, K)

  temp::Matrix{Float64} = Matrix{Float64}(undef, K, K) #pre-allocate working matrix
  RRinv::Matrix{Float64} = BLAS.gemm('N', 'T', Rinv, Rinv) #this is equivelent to [X'X]^-1

  #iterate over all stocks
  for sdf ∈ fs.eachport
    #need to multiply through by the error
    T::Int = size(sdf,1)
    Xₙ::Matrix{Float64} = [ones(T) Matrix{Float64}(sdf[!, coefnames[2:end]])]
    ε::Vector{Float64} = sdf[!, fs.Fret] .- Xₙ * λᵢ
    Xₑ::Matrix{Float64} = Xₙ .* ε

    Sₜ::Matrix{Float64} += BLAS.gemm('T','N',1.0/N, Xₑ, Xₑ)

    for v::Int ∈ 1:lags
      #overwrites temp with (1/N)R'R
      BLAS.gemm!('T', 'N', 1.0/N, view(Xₑ, (v+1):T, :),view(Xₑ, 1:(T-v), :), 0.0, temp)
      Sₜ .+= (lags + 1 - v)/(lags+1.) .* (temp .+ temp')
    end

  end
  #Sₜ = Matrix{Float64}(I,K,K)
  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRinv, Sₜ)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(N), RRinvS, RRinv, 0.0, Σ)
  #println("$(diag(Σ .* dofCorrect))")sasa
  return Σ
end

#this version provides an intuitive check
#check #1
function FMNeweyWestSlow(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}
  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]

  #pre-allocate
  K::Int = size(X,2)
  Sₜ::Matrix{Float64} = zeros(K, K)
  #xxt::Matrix{Float64} = Matrix{Float64}(undef, K,K)
  temp::Matrix{Float64} = Matrix{Float64}(undef, K,K)

  #helper function using the Lars method
  @inline function Rₜ!(v::Int, X::Matrix{Float64}, ε::Vector{Float64};
      T::Int = size(X,1), out::Matrix{Float64} = Matrix{Float64}(undef, K,K))

    out .= 0.0
    for t ∈ (1+v):T
      out .+= X[t,:] * X[t-v,:]' * ε[t] * ε[t-v]
    end

    #out ./= T
    return out
  end

  for sdf ∈ fs.eachport
    T::Int = size(sdf, 1)

    Xₙ::Matrix{Float64} = [ones(T) Matrix{Float64}(sdf[!, coefnames[2:end]])]
    ε::Vector{Float64} = sdf[!, fs.Fret] .- Xₙ * λᵢ
    Sₜ .+= Rₜ!(0, Xₙ, ε, out=temp)


    for v ∈ 1:lags #iterate over lags
      Rₜ!(v, Xₙ, ε, out=temp)
      Sₜ .+= (lags + 1. - v)/(lags + 1.) .* (temp .+ temp')
    end
  end

  Sₜ ./= N

  XXInv::Matrix{Float64} = (X' * X) \ Matrix{Float64}(I,K,K)


  return N * XXInv * Sₜ * XXInv
end

#alg from Correcting for Both Cross-Sectional and
#  Time-Series Dependence in Accounting Research by Ian Gow, Gaizka Ormazabal and Daniel Taylor.
# check #2
function FMNeweyWestML(fs::FactorSpecification,
  coefnames::Vector{Symbol},
  λᵢ::Vector{Float64},
  lags::Int)::Matrix{Float64}
  N::Int = size(fs.df,1)

  X::Matrix{Float64} = [ones(N) Matrix{Float64}(fs.df[!, coefnames[2:end]])]
  ε::Vector{Float64} = fs.df[!, fs.Fret] .- X * λᵢ
  ports::typeof(fs.df[!, fs.Fport]) = fs.df[!, fs.Fport]
  #pre-allocate
  K::Int = size(X,2)
  Sₜ::Matrix{Float64} = zeros(K, K)

  #helper function using the Lars method

  for t ∈ 1:N
    Sₜ .+= X[t,:] * X[t,:]' * ε[t]^2
  end

  for v::Int ∈ 1:lags
    for t ∈ (v+1):N
      (ports[t] == ports[t-v]) && (
        Sₜ .+= (lags + 1 - v)/(lags+1.) .* (X[t,:] * X[t-v,:]' .+ X[t-v,:] * X[t,:]') * ε[t] * ε[t-v])
    end
  end
  Sₜ ./= N
  #Sₜ = Matrix{Float64}(I,K,K)

  #the final step is to multiply out the var-covar sandwhich
  Σ::Matrix{Float64} = Matrix{Float64}(undef, K, K)
  R::Matrix{Float64} = qr(X).R
  Rinv::Matrix{Float64} = (R)\Matrix{Float64}(I,K,K) #this will save some computational time
  RRinv::Matrix{Float64} = BLAS.gemm('N', 'T', Rinv, Rinv) #this is equivelent to [X'X]^-1

  #this is [X'X]^-1S=[R'R]^-1S
  RRinvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRinv, Sₜ)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(N), RRinvS, RRinv, 0.0, Σ)

  return Σ
end

function characteristicFM(fs::FactorSpecification; parallel::Bool = false, lags=2)::FactorRegression
  local t::Int
  local λ::Matrix{Float64}
  local λᵢ::Vector{Float64}
  local σᵢ::Vector{Float64}
  local name::Symbol
  local coefidx::Dict{Symbol, Int}


  #make the time-series regression aggregation function
  T::Type = Vector{Float64}
  @inline aggfunc(fmlm::FMLM) = fmlm.β

  results = crosssectionalregressions(fs, T, aggfunc=aggfunc, parallel=parallel)
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

  Σ::Matrix{Float64} = FMNeweyWest(fs, coefnames, λᵢ, lags)
  σᵢ = diag(Σ).^ 0.5

  println("\nFastNW: $σᵢ " *
    #"\nSlowNW: $(diag(FMNeweyWestSlow(fs, coefnames, λᵢ, lags)).^0.5)" *
    #"\nFMNeweyWestML: $(diag(FMNeweyWestML(fs, coefnames, λᵢ, lags)).^0.5)" *
    "\nNormalSE: $((c->std(c)).(eachcol(λ)) ./ t^0.5)" *
    "\nλᵢ: $λᵢ"
    )

  #=if length(λᵢ) > 8
    (vmax, ind) = findmax(λ[:,end-2])
    println("coefname: $(coefnames[end-2])")
    println("Max val = $vmax")
    println(fs.eachdate[ind])
    for sdf ∈ fs.eachdate
      println("$(sdf[1,:date]): $(size(sdf))")
    end
    error("Debug info printed.")
  end=#

  name = fs.name


  return FactorRegression(λ, λᵢ, σᵢ, t, name, coefidx)
end

#=
86-50-31
struct FactorRegressions
  frs::Vector{FactorRegression}
  nameidx::Dict
  fssidx::Dict{}
end

function FactorRegressions=#
