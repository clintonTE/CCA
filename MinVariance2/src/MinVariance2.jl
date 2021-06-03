using Revise #load this first
module MinVariance2
using Revise


#########################Packages#########################
using Juno, UnPack, Atom
using LinearAlgebra, Distributions, CSV, DataFrames, Statistics
using StatsBase, Random, CodecZlib, CodecZstd, XLSX, Dates

using Finometrics

#########################Globals#########################
#should only use these if signficant advantage to a parameter

#########################Constants#########################
const PARAMETER_FILE = "MVParameters"
const PARAM = Dict{Symbol,Any}()

#########################Files#########################
include("Parameters.jl")
include("Utilities.jl")

end # module
