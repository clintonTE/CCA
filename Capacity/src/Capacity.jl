using Revise
module Capacity
using Revise

using Atom, Juno, UnPack
using Optim, Pkg, SparseArrays
using FixedEffectModels
using CSV, DataFrames, Serialization, Dates, Statistics, StatsBase, BenchmarkTools
using LinearAlgebra, Distributions, Distributed, VegaLite, Random
using CodecZlib, CodecZstd, GLM, XLSX, InteractiveUtils, Base.Iterators
using PGFPlotsX, Colors, LaTeXStrings
#, BlackBoxOptim
#using MCMCChains
using CUDA
using Transducers, BangBang, LineSearches, LoopVectorization, Zygote


using Finometrics

#https://github.com/clintonTE/Finometrics#master

#TODO: Find a way to test the Jacobian from Zygote
#maybe see if there is a way to do a hand calculation
using CUDA: CuVector, CuMatrix, CuArray
CUDA.allowscalar(false)


import Base: Symbol, show, length, iterate, broadcastable, *

abstract type AbstractCharacteristic end

const RANDOM_SEED_VALUE = 1111
Random.seed!(RANDOM_SEED_VALUE)

#########################CRITICAL MACROS####################
#NOTE: This macro allows for the parallelization routines to be shut down
#by flicking the above PARALLEL switch
#BROKEN!!!
#=macro mpar(expr)
  quote
    if PARALLEL[]
      :($(Threads.@threads($expr)))
    else
      :($($expr))
    end
  end
end


macro mspawn(expr)
  quote
    if PARALLEL[]
      :($(Threads.@spawn($expr)))
    else
      :($($expr))
    end
  end
end=#


###########################CONSTANTS##########################
const PARAMETER_FILE = pwd() * "\\Parameters.xlsx"
const PARALLEL = Ref{Bool}(true) #perhaps we can keep this the only global state var
const PARAM = Dict{Symbol,Any}()
const LARGE_VAL = prevfloat(Inf)^0.33
const LARGE_VAL32 = (prevfloat(Inf32)^0.5f0)/10f0
const SMALL_VAL = 1e-6
const SMALL_VAL32 = 1f-6


IN_SMALL_CSV_STREAM(p::String) = open(p)
IN_SMALL_CSV_STREAM(F::Function, p::String) = open(F, p)
const SMALL_CSV_EXTENSION = "csv"

IN_CSV_STREAM(p::String) = read(p) |> x->transcode(ZstdDecompressor,x)
IN_CSV_STREAM(F::Function, p::String) = open(F, ZstdDecompressorStream, p)
const CSV_EXTENSION = "csv.zstd"

function IN_BIN_STREAM(p::String)
  local obj
  open(ZstdDecompressorStream, p) do io
      obj = deserialize(io)
  end

  return obj
end


const BIN_EXTENSION = "bin.zstd"
function OUT_BIN_STREAM(p::String, obj, level=1)

  io = ZstdCompressorStream(open(p, "w"), level=level)
  try
    serialize(io, obj)
  finally
    close(io)
  end
end

#####abstract types
abstract type AbstractIterState{T<:Real} end
abstract type AbstractVolumeParts{TM<:AbstractMatrix, TV<:AbstractVector, T<:Real} end

#note- we ahve to be careful about using methods from AbstractVolumeParts since we ahve no V₀
abstract type AbstractVolumePartsIter{TM,TV,T} <: AbstractVolumeParts{TM,TV,T} end

abstract type AbstractAnnhilator{T<:Real} end

include("Utilities\\GeneralUtilities.jl")
#include("Measure\\DistributionUtilities.jl")  #disabled due to flux incompatibility
include("Utilities\\TestUtilities.jl")
include("Utilities\\Annihilator.jl")
include("Utilities\\HyperReg.jl")
include("Utilities\\AR1.jl")

include("Preliminary\\PrelimBeta.jl")
include("Preliminary\\PrelimFunds.jl")
#include("Preliminary\\PrelimPreprocessCRSP.jl")
include("Data\\TrailingCRSPMeasure.jl")
include("Data\\Inputs\\Parameters.jl")
include("Data\\Inputs\\IBES.jl")
include("Data\\Inputs\\CCM.jl")
include("Data\\Inputs\\CRSP.jl")
include("Data\\Inputs\\COMP.jl")
include("Data\\Inputs\\FF.jl")

include("Data\\MergeCRSPCOMP.jl")
include("Data\\CRSPweekly.jl")
include("Data\\CRSPmonthly.jl")
include("Measure\\Setup\\WeightSpec.jl")
include("Measure\\Setup\\Characteristics.jl")
include("Measure\\Setup\\MeasureSpec.jl")
include("Measure\\Setup\\Standardization.jl")
include("Measure\\Verification\\Returns.jl")
include("Measure\\Verification\\Funds.jl")
#include("VolumeParts\\VolumeParts.jl")

#include("VolumeParts\\VolumePartsLSq.jl")
include("Measure\\VolumeParts\\XYIter.jl")
include("Measure\\VolumeParts\\TestXv.jl")
include("Measure\\VolumeParts\\VolumePartsIter.jl")

#include("VolumeParts\\VolumePartsFlux.jl") #disabled due to flux incompatibility
#NOTE: below may be resurrected at some point
#include("Measure\\MCMC\\CapacityChains.jl")
#include("VolumeParts\\VolumePartsMCMC.jl")

include("Measure\\Iter\\IterModels.jl")
include("Measure\\VolumeParts\\VolumePartsIterLevelDemean.jl") #mostly for testing purposes
include("Measure\\VolumeParts\\VolumePartsIterTest.jl")
include("Measure\\Iter\\IterEstimate.jl")
include("Measure\\Iter\\IterPanelAnalysis.jl")
include("Measure\\Iter\\Iter.jl")

#include("Measure\\IterSolveFP.jl")
include("Measure\\Iter\\IterSolveOptim.jl")
include("Measure\\Iter\\IterBootstrap.jl")
#include("Measure\\Iter\\IterSolveBB.jl")

#include("Measure\\MCMC\\Posteriors.jl") #NOTE: these may be resurrected for use later
#include("Measure\\MCMC\\MCMCEstimate.jl")
#include("Measure\\MCMC\\MCMCModels.jl")
#include("Measure\\MCMC\\TestMCMC.jl")

include("Measure\\Verification\\Comomentum.jl")
#include("Measure\\Verification\\Comomentumbb.jl")

include("Measure\\Setup\\Measures.jl")
include("Data\\MomentumStrategies.jl")
include("Measure\\Setup\\Controls.jl")
include("Data\\FormPanel.jl")
include("Measure\\ComputeMeasure.jl")

include("Analysis\\Analyze.jl")
include("Analysis\\LeadLag.jl")
include("Analysis\\LPFilter.jl")
include("Analysis\\MovingBlockBootstrap.jl")
include("Analysis\\Verify.jl")
include("Analysis\\IO\\Tables.jl")
include("Analysis\\IO\\Graphs.jl")
include("Analysis\\IO\\LLTables.jl")
include("Analysis\\IO\\IO.jl")



REGRESSION_METHOD = Dict{Symbol, Any}(
#  :flux=>formfluxregression,
#  :opt=>formoptregression,
  :iter=>formiterregression,
  #:mcmc=>formmcmcregression
)


end # module
