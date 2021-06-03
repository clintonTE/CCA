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
  Serialization,  CSV, Random, Distributed, Gadfly, Fontconfig, Cairo# NOTE: uncomment for graphs


#########################Exported functions##########
export testportfolio,
  testcov,
  testtrisectedportfolio,
  matlab2julia,
  mathematica2juliaEq,
  mathematica2juliaW,
  testω,
  testposteriors,
  estimateminvariance,
  initializepar, #need this to set up the parallelization
  cleanup,
  makeconvergencegraphs,
  Estimation


################Priors and initial vlaues############
###Infromative Priors
const D_θG = 0.16^2/20 #current vix converted to daily
const D_δ²G = .76/20*2 #reasonably diffuse prior from volatiltiy of the vix multiplied by 2
const D_θSG = D_θG
const D_δ²SG = D_δ²G*2.0 #more uncertain
const D_θSGP = D_θG
const D_δ²SGP = D_δ²G*4.0 #most uncertain

const D_zG = .76/20/510*2*.01
const D_zP = D_zG * 2
const D_zGP = (.76/20/510)^0.5*(D_zG/2)/2 #std dev * std dev G * std dev P /2
const D_Ψ = [D_zG D_zGP; D_zGP D_zP]
const D_ν = 3.01

const U_θG = 0.0001 #current vix converted to daily
const U_δ²G = 0.0001 #reasonably diffuse prior from volatiltiy of the vix multiplied by 2
const U_θSG = 0.0001
const U_δ²SG = 0.0001 #more uncertain

const U_zG = 0.0001
const U_zP = 0.0001
const U_zGP = 0.0000 #std dev * std dev G * std dev P /2
const U_Ψ = [U_zG U_zGP; U_zGP U_zP]
const U_ν = 3.01



#Initial Parameter Values
const D_σ²G= 0.0001
const D_ζ²G= 0.0001
const D_ζ²P = 0.0001
const D_ζGP = 0.0000
const D_SG = 0.0001
const D_SGP = .0001

const SMALL = 1.0e-20

###Legacy, not used right now
const U_θSGP = 0.0001
const U_δ²SGP = 0.0001 #most uncertain

##############Other constants and methedological parameters############
#const TEX_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\Endowment Project"
#const GRAPH_PATH = TEX_PATH * "\\sub\\fig"
#const WORKING_PATH = pwd() * "\\working"

const OUT_NAME = "test"
const OUT_PATH = "output"
#const OUT_GRAPH_PATH = "C:\\Users\\Clinton Tepper\\Dropbox\\Apps\\Overleaf\\minvar\\sub"
const OUT_GRAPH_PATH = "output"
const OUT_ESTIMATE_NAME = "lastestimate"
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}}

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
const MIN_BURN = 1_000 #max used in convergence graphs
const CONVERGENCE_WINDOWS = [1000,10_000,100_000,1_000_000]
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
include("code\\Output.jl")
include("code\\Testsandoneoffs.jl")


end
