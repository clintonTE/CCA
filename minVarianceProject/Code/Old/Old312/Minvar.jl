module Minvar
##################Dependencies####################
using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end

using DataFrames, Distributions, GZip, ForwardDiff, Dates, LinearAlgebra, Statistics,
   Measures, Formatting, StatsBase, Distributions, DataStructures, Finometrics, Distributed,
  Serialization,  CSV, Random, Distributed#, Gadfly, Fontconfig, Cairo NOTE: uncomment for graphs


#########################Exported functions##########
export testportfolio,
  testcov,
  testtrisectedportfolio,
  matlab2julia,
  mathematica2julia,
  testω,
  testposteriors,
  estimateminvariance,
  initializepar, #need this to set up the parallelization
  cleanup

################Priors and initial vlaues############
#Priors
const D_θG = 0.16^2/20 #current vix converted to daily
const D_δ²G = .76/20*2 #reasonably diffuse prior from volatiltiy of the vix multiplied by 2
const D_αG = 1.01 #Pick a diffuse prior, but want the expectation to exist
const D_βG = 7.0e-5 / 2.0 * (D_αG - 1.0) #derive from vix of vix (See documentation)
const D_αP = D_αG
const D_βP = D_βG * 2.0 * (D_αP - 1.0) #expect this is more volatile.
const D_θSG = D_θG
const D_δ²SG = D_δ²G*2.0 #more uncertain
const D_θSGP = D_θG
const D_δ²SGP = D_δ²G*4.0 #most uncertain


#Initial Parameter Values
const D_σ²G=.001
const D_ζ²G= .004
const D_ζ²P = .009
const D_SG = .002
const D_SGP = .003

##############Other constants and methedological parameters############
#const TEX_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\Endowment Project"
#const GRAPH_PATH = TEX_PATH * "\\sub\\fig"
#const WORKING_PATH = pwd() * "\\working"

const OUT_NAME = "test"
const OUT_PATH = "output"

const N_BURN = 100
const N_SAMPLES = 100

const DATA_NAME = "sampledata"
const DATA_PATH = pwd() * "\\data"

const DATE_BEGIN = Date(2015,1,1)
const DATE_END =  Date(2016,12,31)
const DATE_FORMAT = "yyyymmdd"

const MAX_WORKERS = Base.Sys.CPU_THREADS ÷ 2 # gets num of physical cores

const NUM_TEST_ASSETS = 100 #this shouldn't make a big difference in the results
const MIN_ADJUSTMENT = 0.75 #smaller values->more extreme weights allowed for test portfolio

#without this throws an error on each run if already defined
if !(@isdefined RANDOM_WEIGHT_DIST)
  const RANDOM_WEIGHT_DIST = TriangularDist(-1.0,1.5,1/NUM_TEST_ASSETS)
end

#############File linkage###############################
include("code\\MinvarPar.jl")
include("code\\Universe.jl")
include("code\\Portfolio.jl")
include("code\\Parameters.jl")
include("code\\Posteriors.jl")
include("code\\TrisectedPortfolio.jl")
include("code\\X2Julia.jl")
include("code\\Estimation.jl")
include("code\\Testsandoneoffs.jl")
end
