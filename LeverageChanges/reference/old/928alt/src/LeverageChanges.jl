using Revise
module LeverageChanges
using Revise

using  CSV, DataFrames, Serialization, Dates, Statistics,
  Distributed, Finometrics, Measures, StatsBase,
  Formatting, CodecLz4, CodecZlib,  Distributions, LinearAlgebra

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


import Base: eltype, getindex, iterate, length, isfinite, merge

#necessary due to exporting conflict
describe(df::DataFrame) = DataFrames.describe(df)


#############################CONSTANTS#####################
TEST_OUTPUT = true

const OUT_PATH = pwd() * "\\output"
const DATA_PATH = pwd() * "\\data"
const ARCHIVE_PATH = pwd() * "\\working"

const CRSP_TYPE = :monthly
const CRSP_DATE_FORMAT = dateformat"yyyymmdd"

const CRSPM_NAME = "CRSP-M"
const CRSPM_PARTS = ["CRSP-M"]

const CRSPD_NAME = "CRSP-D"
const CRSPD_PARTS = ["CRSP-D-II", "CRSP-D-I"]

const COMP_NAME = "COMP"
const COMP_A_NAME = "COMP-A"
const COMP_Q_NAME = "COMP-Q"
const COMP_DATE_FORMAT = dateformat"yyyymmdd"

const CCM_NAME = "CCM"
const CCM_DATE_FORMAT = dateformat"yyyymmdd"
const LAST_DATE_STRING = "21001231"
const LAST_DATE = Date(LAST_DATE_STRING, CCM_DATE_FORMAT)

const MONTHS_2_STALE_LONG = 18
const MONTHS_2_STALE_SHORT = 6
const WINSORIZE_PROP = 0.001

const UNIV_NAME = "UNIV"

const LEV_NAME = "LEV"

const PORT_NAME = "PORT"
const PANEL_NAME = "PANEL"

const VALIDATE_MERGED = true

const USE_QUARTERLY_DATES = true

const MIN_ROWS_PER_GROUP = 14


##constants related to outside data
const VWCRSP_NAME = "VWCRSP-D"
const RFCRSP_NAME = "RFCRSP-M"
const RETAINED_VWCRSP_COLS = [:date, :vwindd, :ewindd]
const RETAINED_RFCRSP_COLS = [:date, :t30ind]
const IDXCRSP_PATH = "$DATA_PATH\\crspidx"
const IDXCRSP_DATE_FORMAT = CRSP_DATE_FORMAT
const IDX_DATE_COL = :caldt

const FF5_NAME = "FF5"
const FF3M_NAME = "FF3M"
const FF_PATH = "$DATA_PATH\\french"
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


############################FILES
include("PanelFormation\\PreprocessCRSP.jl")
include("PanelFormation\\PreprocessCOMP.jl")
include("PanelFormation\\MergeCOMPCCM.jl")
include("PanelFormation\\MergeCRSPCOMP.jl")
include("Utilities.jl")
include("PanelFormation\\FormDataSeries.jl")
include("PortfolioFormation\\PortfolioConstruction.jl")
include("PortfolioFormation\\ConstructPortfolios.jl")
include("PanelFormation\\MergeOthers.jl")
include("Analysis\\AnalyzeFactors.jl")
include("Analysis\\Table1.jl")
include("Analysis\\Replicate.jl")

end # module
