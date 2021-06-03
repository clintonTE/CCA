using Revise
module QFactors
using Revise

using  CSV, DataFrames, Serialization, Dates, Statistics,
  Distributed, Finometrics, Measures,
  Formatting, CodecLz4, CodecZlib

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

############################CONSTANTS

const OUT_PATH = pwd() * "\\output"
const DATA_PATH = pwd() * "\\data"


const CRSP_NAME = "CRSP-M"
const CRSP_DATE_FORMAT = dateformat"yyyymmdd"

const COMP_NAME = "COMP"
const COMP_A_NAME = "COMP-A"
const COMP_Q_NAME = "COMP-Q"
const COMP_DATE_FORMAT = dateformat"yyyymmdd"

###########################IO Compression/Decompression routines
IN_CSV_STREAM(p::String) = LZ4DecompressorStream(open(p))
IN_CSV_STREAM(F::Function, p::String) = open(F, LZ4DecompressorStream, p)

IN_JLS_STREAM(p::String) = LZ4DecompressorStream(open(p))
IN_JLS_STREAM(F::Function, p::String) = open(F, LZ4DecompressorStream, p)

OUT_JLS_STREAM(p::String) = LZ4CompressorStream(open(p, "w"))
OUT_JLS_STREAM(F::Function, p::String) = open(F, LZ4CompressorStream, p, "w")


##############################EXPORTS




############################FILES
include("QUtilities.jl")
include("PreprocessCRSP.jl")
include("PreprocessCOMP.jl")



end
