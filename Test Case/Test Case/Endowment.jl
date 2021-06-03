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
  Gadfly, NLopt, Measures, Formatting, HTTP, TexTables,
  CodecZlib, StatsBase, Distributions, DataStructures

#using CT

##################Exports########################
export initializePar, ###Acquire NCCSData
cleanup,
downloadNCCSCSVs,
preProcessData,
validateSingleTypeCSV,
readSingleTypeCSV,
readMetaIndex,
formMetaData,
MetaTuple,
MetaDictionary,
NCCSFile,
NCCS_PATH,

MBool,
MInt,
MFloat64,
MSymbol,
MString



#############################Shared Constants###########
const OUTPUT_PATH = pwd() * "\\output"
const WORKING_PATH = pwd() * "\\working"
const DATA_PATH = pwd() * "\\data"
const NCCS_PATH = DATA_PATH * "\\NCCS"
const DOCUMENTATION_PATH = DATA_PATH * "\\documentation"

const HEADER_NAME = "header.tex"
const FOOTER_NAME = "footer.tex"
const METADATA_NAME = "NCCSMetadata"
const META_INDEX_NAME = "NCCSMetaIndex"
const META_PATH = DOCUMENTATION_PATH * "\\raw meta"
const SUB_TABLES_PATH = META_PATH * "\\subtables"
const FILES_INDEX_NAME = "NCCSFiles"
const KEEPER_NAME = "Keeper List"
const Nothing = Void #NOTE: 0.7 compatabiltiy hack, delete on ugrade
const VAR_TYPES = Dict(:Numeric=>Float64, :Character=>String)

const MBool = Union{Bool, Missing}
const MInt = Union{Int, Missing}
const MFloat64 = Union{Float64, Missing}
const MSymbol = Union{Symbol, Missing}
const MString = Union{String, Missing}

#default compressor algorithms
I_CSV_STREAMER(p::String) = GzipDecompressorStream(open(p))
I_JLS_STREAMER(p::String) = GzipDecompressorStream(open(p))
O_JLS_STREAMER(p::String) = GzipCompressorStream(open(p, "w"))

include("EndowmentCode\\EndowmentPar.jl") #holds supporting functions for parallelization
include("EndowmentCode\\MetaTuples.jl") #holds methods related to field-level metadata
include("EndowmentCode\\MetaDictionaries.jl") #holds methods related to file level metadata
include("EndowmentCode\\Metadata.jl") #holds all other metadata information
include("EndowmentCode\\DownloadNCCSCSVs.jl") #holds methods related to downloading the data
include("EndowmentCode\\NCCSFiles.jl") #holds methods related to the data files
include("EndowmentCode\\EndowmentData.jl") #holds high-level scripts and other miscellaneous





end
