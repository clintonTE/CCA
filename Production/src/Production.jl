module Production

##################Dependencies####################
using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

#=if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end=#

using Finometrics #my personal module
using DataFrames, Distributions, Dates, LinearAlgebra, Statistics,
   Formatting, StatsBase, Serialization,  CSV, Random, Plots
  #, Measures, Gadfly, Fontconfig, Cairo NOTE: uncomment for graphs

############################Exported functions
export Firm, Problem, Simulation, graphscript!

#############################Constants

const TOL_CONSTRAINT = 1.0e-10
const N_TOL = 1.0e-10

const OUT_NAME = "test"
const OUT_PAPER = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\productionproject\\sub"
const OUT_PATH = "output"

const PENALTY = -10_000.0


##########################code files##########################
include("Numerical.jl")
include("Firm.jl")
include("Problem.jl")
include("Simulation.jl")
include("Graphscripts.jl")

end # module
