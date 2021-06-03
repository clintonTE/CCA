module Production
##################Dependencies####################
using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end

using Finometrics #my personal module
using DataFrames, Distributions, GZip, ForwardDiff, Dates, LinearAlgebra, Statistics,
   Formatting, StatsBase, DataStructures, Distributed,
  Serialization,  CSV, Random, Plots, UnicodePlots, SIMD
  #, Measures, Gadfly, Fontconfig, Cairo NOTE: uncomment for graphs

############################Exported functions
export Firm, Problem, Simulation, graphscript!

#############################Constants
const DEF_STATES_PER_VAR = 10

###Firm Specific DefaultParameters
const DEF_τ = [0.0, 0.0]
const DEF_δ = 0.3
const DEF_Θ = (0.8, 0.2)
const DEF_α = 1.05
const DEF_β = 0.98
const DEF_γ = Float64(1/3)
const DEF_ζ = (1.0, 0.0)
const DEF_KSTATES = collect(range(0.0, stop=5.0, length=DEF_STATES_PER_VAR))
const DEF_YSTATES = collect(range(0.0, stop=5.0, length=DEF_STATES_PER_VAR))

const TOL_CONSTRAINT = 1.0e-10
const N_TOL = 1.0e-10

#const VALID_SIMD_WIDTHS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]


#Other problem default parameters
const DEF_p = [0.81 0.68; 0.19 0.32]
const DEF_M = 2
const DEF_D_TOL = 10e-6
const DEF_MAX_ITER = 1000
const DEF_N_SIMS = 10^5


const OUT_NAME = "test"
const OUT_PAPER = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\productionproject\\sub"
const OUT_PATH = "output"

const MAX_WORKERS = Base.Sys.CPU_THREADS ÷ 2 # gets num of physical cores
const PENALTY = -10_000.0


##########################code files##########################
include("code\\Numerical.jl")
include("code\\Firm.jl")
include("code\\Problem.jl")
include("code\\Simulation.jl")
include("code\\Graphscripts.jl")

end
