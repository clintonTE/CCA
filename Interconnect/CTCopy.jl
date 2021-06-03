module CTCopy

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


##################Dependencies

using  DataFrames, Distributions, StatsBase, GLM, CategoricalArrays,
  Dates, NLopt, ForwardDiff, Formatting, LinearAlgebra


######################Methods####################
export CTLM, #Regression methods
  CT2SLS,
  project!,
  getResid!,
  pullModelMatrix,
  getCoeff!,
  getHomoskedΣ!,
  getModWhiteΣ!,
  getWhiteΣ!,
  getNeweyWest!,
  getModWhiteΣSlow,
  getHomoskedΣSlow,
  getWhiteΣSlow,
  getNeweyWestSlow,
  dropNullsFromDF,
  CTExpr,
  CTSym,
  CTQR,
  get1stStage,
  getTerm,
  getR,
  getR²,
  getNeweyWestFunc,
  getClustered!,

  texTable, #IO Mthods
  writeTables2File,
  array2String,
  num2Str,
  vec2String,

  skewnessStat, #CTstat methods
  kurtosisStat,
  JBStat,
  UnivariateTest,
  autoρ,
  ljungBox,
  test,
  pAutoρ,
  ARMA,
  ARMAFillYHat!,
  ARMALogLike,
  ARMAEstimate,
  generalizeARMA!,
  narrowARMA!,
  lagDF!,

  CTΦInv,#CTNumerical
  CTΦ,

  VAR, #CTVAR
  VAREstimate,
  VARPrediction,
  varΣ,
  propagateImpulse




##################Custom types
###NOTE: 0.7 Compat Hacks, delte on upgrade
#const Nothing = Void
#const LinearAlgebra = Base.LinAlg #NOTE:Compat Hack
#const AbstractRange = Range #NOTE compat hack

###CTReg Types
const CTExpr = Union{Symbol,Expr,Nothing}
const CTSym = Union{Symbol,Nothing}
const CTData = Union{Float64,Int,Date,Symbol, CategoricalValue, Missing, Nothing}#,
abstract type CTModel end

###VAR Types
const Structure = Union{Nothing, Vector{Float64}}
const CoefficientErrors = Union{Nothing, Vector{Float64}}
const RegressionSS = Union{Nothing, Float64}
const RegressionN = Union{Nothing, Int}

################Constants

#shared cosntants
const DEFAULT_STAR_LEGEND = "*p<0.1, **p<0.05, ***p<0.01"
const DEBUG_CTMOD = false
const DEFAULT_DECIMALS = 4

#Stat constants
const MAX_TIME = 30.
const F_TOL_ABS = 10. ^ -12.
const F_TOL_REL = 10. ^ -10.
const OPT_BOUND = 2.

#numerical constants
const sqrt2Inv = 1.0/(2.0^.5)
const sqrt2PiInv = 1.0/(2.0*π)^0.5
#  MFloat64, MInt, MBool, MSymbol}


###################Files############

include("CTReg.jl")
include("CTIO.jl")
include("CTStat.jl")
include("CTVAR.jl")
include("CTNumerical.jl")

end
