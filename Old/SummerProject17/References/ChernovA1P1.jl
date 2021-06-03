module  ChernovA1P1

#this is useful if we have mdoules
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end
#=
weave(Pkg.dir("$(pwd())\\ChernovA1P1.jl"),
  informat="script",
  out_path = "$(pwd())\\ChernovA1P1.html",
  doctype = "md2html")
=#
using DataFrames, Distributions, StatsBase, GZip, JLD, Gadfly, CTModCopy, JuMP, NLopt


####################Working constants
const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "m/d/y"
const TBILL_NAME = "t-bill"

#useful type
const CTSym = Union{Symbol,Void}
const PlotContainer = Union{Plot,Gadfly.Compose.Context}
const YEAR_RANGES = [1954:2016,1954:1975, 1976:1981, 1982:2006, 2007:2016]



#a simple routine for making an easy to read dataframe binary
function preProcessTBills()::Void
  dateFormat::String = DATE_FORMAT_STR
  tBillDF::DataFrame = readtable("$DATA_PATH\\$TBILL_NAME.csv")

  #fix the dates
  tBillDF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(tBillDF[:,:DATE])::DataVector{Date}
  sort!(tBillDF, cols = [:DATE])
  tBillDF[:,:TB3MS] .*= .01
  tBillDF[:,:rates] ./=12.0 #(1.0+tBillDF[:,:TB3MS]).^(1.0/12.0) .- 1.0

  stream::IOStream = open("$DATA_PATH\\$TBILL_NAME.jls", "w")
  serialize(stream, tBillDF)
  close(stream)

  return nothing
end

#this function takes in a year range and a list of dates and returns a boolean vector
# which is true iff a given date is inside the range.
getDateIndFromYears(dates::DataVector{Date}, yearRng::UnitRange{Int})::Vector{Bool} =
  (d::Date->Dates.year(d) ∈ yearRng).(dates)

#this type is designed to hold the specifications for a single iteration of the problem
mutable struct Problem
  yearRange::UnitRange
  rates::Vector{Float64}
  T::Int
end

#this constructor takes in a year range and creates a problem specification
function Problem(df::DataFrame, yearRange::UnitRange)::Problem

  #create a simple vector from the data frame given the year range
  rates::Vector{Float64} = Vector{Float64}(df[getDateIndFromYears(df[:,:DATE], yearRange),:rates])
  T::Int = length(rates)
  return Problem(yearRange, rates, T)

end


function estimateVasicek(problem::Problem)

  #Define the objective function
  T::Int = problem.T
  x::Vector{Float64} = problem.rates

  #u::Float64 = mean(x)
  #v::Float64 = var(x, corrected=false)
  function obj(u, b, v)
    #println("u=$u, b=$b, v=$v")
    l = -(x[1]-u)*(x[1]-u)*(1.0-b*b)
    uTimes1Minusb = u*(1.0-b)

    for t::Int ∈ 2:T
      val = (x[t] - uTimes1Minusb - b*x[t-1])
      l += -val*val
    end

    l /= 2.0*v
    l += 0.5*(log(1.0-b^2.0))-0.5*T*log(v)

    return l
  end

  #println(obj(0.01, 0.5, 0.0001))

  #build the optimizer
  #Algorithms: :LD_LBFGS (fails), :LN_COBYLA (bad), :LN_BOBYQA (bad)
  #  :LD_SLSQP (fails), :LD_TNEWTON (fails), :LD_VAR2/:LD_VAR1 (fails),
  #  :LD_MMA (bad), :GN_ESCH (hang), :GN_DIRECT (questionable), GN_DIRECT_L (questionable),
  # GN_ISRES (bad), :LN_NELDERMEAD(bad), :LN_SBPLX(bad), :LD_LBFGS (fails), :GD_STOGO (fails)
  #GN_CRS2_LM (OK)

  mod = Model(solver=NLoptSolver(algorithm=:GN_CRS2_LM))

  @variable(mod, 10.0^-6.0 <= u <= .2)
  @variable(mod, 10.0^-8.0 <= b <= 1.0-10.0^-8.0)
  @variable(mod, 10.0^-12.0 <= v <= .1)

  setvalue(u, mean(x))
  setvalue(b, .99)
  setvalue(v, var(x))

  JuMP.register(mod, :obj, 3, obj, autodiff=true)
  @NLobjective(mod, Max, obj(u,b,v))

  #println("Running Solve: ")
  output = solve(mod)

  #println("For range $(problem.yearRange), estimates are: ")
  #println("μ = $(getvalue(u)), b = $(getvalue(b)),  σ² = $(getvalue(v)), σ= $(getvalue(v)^0.5)")
  #println("$(getvalue(u)), $(getvalue(b)),  $(getvalue(v)), $(getvalue(v)^0.5)",
  #  ", $(obj(getvalue(u),getvalue(b),getvalue(v))), $(mean(problem.rates)), $(std(problem.rates; corrected=false))")

end

#"this is a script function to hold non-reusable code"
function problemScript(;refreshData::Bool = true, yearRanges::Vector{UnitRange{Int}} = YEAR_RANGES)::Void

  #if we need to, write a new file
  if !isfile("$DATA_PATH\\$TBILL_NAME.jls") || refreshData
    preProcessTBills()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$TBILL_NAME.jls")
  tBillDF::DataFrame = deserialize(stream)
  close(stream)

  #build the problem specifications
  problemSpecs::Vector{Problem} = ((r::UnitRange)->Problem(tBillDF, r)).(yearRanges)

  #println("μ, ϕ, σ², σ, Obj (Log Likelihood), x̄ (unconditional), σ (unconditional)")
  estimateVasicek.(problemSpecs)

  return nothing
end

@time problemScript(refreshData = false)

end
