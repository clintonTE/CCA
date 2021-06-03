using Revise
module QFactors
using Revise

using  CSV, DataFrames, Serialization, Dates, Statistics,
  Distributed, Finometrics, Measures, StatsBase,
  Formatting, CodecLz4, CodecZlib,  Distributions, LinearAlgebra

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


import Base.eltype, Base.getindex, Base.iterate, Base.length

#necessary due to exporting conflict
describe(df::DataFrame) = DataFrames.describe(df)

############################CONSTANTS

const OVERLEAF_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\QFactors"
const TABLE_PATH = "$OVERLEAF_PATH\\sub\\tab"

const OUT_PATH = pwd() * "\\output"
const DATA_PATH = pwd() * "\\data"

const CRSP_TYPE = :monthly
const CRSP_SUFFIX = CRSP_TYPE == :monthly ? 'M' : 'D'
const CRSP_NAME = "CRSP-$CRSP_SUFFIX"
const CRSP_DATE_FORMAT = dateformat"yyyymmdd"

const COMP_NAME = "COMP"
const COMP_A_NAME = "COMP-A"
const COMP_Q_NAME = "COMP-Q"
const COMP_DATE_FORMAT = dateformat"yyyymmdd"

const CCM_NAME = "CCM"
const CCM_DATE_FORMAT = dateformat"yyyymmdd"

const LAST_DATE_STRING = "20191231"
const LAST_DATE = Date(LAST_DATE_STRING, CCM_DATE_FORMAT)

const UNIV_NAME = "$(CRSP_NAME)_$(COMP_NAME)_UNIV"

const FACT_NAME = "$(CRSP_NAME)_$(COMP_NAME)_FACT"

const TEST_OUTPUT = true

const TIME_2_STALE = Month(6)
const MIN_DATE = Date(1972,1,1)
const MAX_DATE = Date(2018,12,31)

const VWCRSP_NAME = "VWCRSP-D"
const RFCRSP_NAME = "RFCRSP-M"
const IDXCRSP_PATH = "$DATA_PATH\\crspidx"
const IDXCRSP_DATE_FORMAT = CRSP_DATE_FORMAT
const IDX_DATE_COL = :caldt

const FF_PORT_NAME = "FF-25-SIZE-BM-D"
const FF_PATH = "$DATA_PATH\\french"
const FF_PORT_RENAME =   Dict(:smalllobm=>:me1bm1, :smallhibm=>:me1bm5,
  :biglobm=>:me5bm1, :bighibm=>:me5bm5)
const FF_PORT_VALS = [Symbol("me$(i)bm$(j)") for i ∈ 1:5 for j ∈ 1:5]
const FF_DATE_FORMAT = CRSP_DATE_FORMAT



#############################Multi-threading

#available lock if needed- use lock(spinlock) unlock(spinlock) etc
const spinlock = Threads.SpinLock()

const PARALLEL = Ref{Bool}(true)

#NOTE: This macro allows for the parallelization routines to be shut down
#by flicking the above PARALLEL switch
macro mpar(expr)
  quote
    if PARALLEL[]
      :($(Threads.@threads($expr)))
    else
      :($($expr))
    end
  end
end


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
include("FormFactors\\PreprocessCRSP.jl")
include("FormFactors\\PreprocessCOMP.jl")
include("FormFactors\\MergeCOMPCCM.jl")
include("FormFactors\\MergeCRSPCOMP.jl")
include("FormFactors\\PortfolioConstruction.jl")
include("FormFactors\\FormFactors.jl")
include("FormFactors\\MergeFACTOthers.jl")
include("AnalyzeFactors.jl")
include("AnalyzeQFactors.jl")





end
