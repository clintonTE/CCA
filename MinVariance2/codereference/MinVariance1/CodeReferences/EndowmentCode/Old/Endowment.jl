module Endowment

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end

#Use this to print
#=
using Weave
codeName = "train"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#
##############0.7/1.0 compatability hacks
#NOTE: Delete these upon 0.7


#############Dependencies####################
using DataFrames, Distributions, CSV, GZip, ForwardDiff,
  Gadfly, NLopt, Measures, Formatting, HTTP,
  CodecZlib

using CT

##################Exports########################
export initializePar, ###Acquire NCCSData
cleanup,
getNCCSFiles,
getNCCSFilesIndex,
preProcessData,
validateSingleTypeCSV,
readSingleTypeCSV,
readMetaIndex,
formMetaData,
MetaTuple,
MetaDictionary,
computeEndowmentReturns



#############################Shared Constants###########
const OUTPUT_PATH = pwd() * "\\output"
const WORKING_PATH = pwd() * "\\working"
const DATA_PATH = pwd() * "\\data"
const NCCS_PATH = DATA_PATH * "\\NCCS"
const DOCUMENTATION_PATH = DATA_PATH * "\\documentation"

const METADATA_NAME = "NCCSMetadata"
const META_INDEX_NAME = "NCCSMetaIndex"
const META_PATH = DOCUMENTATION_PATH * "\\raw meta"
const SUB_TABLES_PATH = META_PATH * "\\subtables"
const FILES_INDEX_NAME = "NCCSFiles"
const KEEPER_NAME = "Keeper List"
const Nothing = Void #NOTE: 0.7 compatabiltiy hack, delete on ugrade
const VAR_TYPES = Dict(:Numeric=>Float64, :Character=>String)


const Conformed = Symbol
const MConformed = Union{Symbol, Missing}

#default compressor algorithms
I_CSV_STREAMER(p::String) = GzipDecompressorStream(open(p))
I_JLS_STREAMER(p::String) = GzipDecompressorStream(open(p))
O_JLS_STREAMER(p::String) = GzipCompressorStream(open(p, "w"))

include("EndowmentCode\\EndowmentPar.jl")
include("EndowmentCode\\MetaTuples.jl")
include("EndowmentCode\\MetaDictionaries.jl")
include("EndowmentCode\\EndowmentGetNCCS.jl")
include("EndowmentCode\\EndowmentData.jl")





end
