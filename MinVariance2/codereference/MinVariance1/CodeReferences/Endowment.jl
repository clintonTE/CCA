#using Revise
module Endowment
using Revise

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

if "$(pwd())\\data" ∉ LOAD_PATH
  push!(LOAD_PATH,"$(pwd())\\data")
end


#TODO (maybe): Refactor for metadata of processed data fields
#Revert CSV from master (test once this is done)
#Use this to print
#=
using Weave
codeName = "train"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

#############Dependencies####################
using DataFrames, Distributions, GZip, ForwardDiff, Dates,
   NLopt, Measures, Formatting, HTTP, TexTables,
  CodecZlib, StatsBase, Distributions, DataStructures, Finometrics, Distributed,
  Serialization, uCSV, CSV, Gadfly, Fontconfig, Cairo#FileIO, CSVFiles

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
downloadIRSPDFs,
processNCCS,
NCCSData,
constructPrimaryOutput,
makeDescriptive

#############################Shared Constants###########
const TEX_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\Endowment Project"
const GRAPH_PATH = TEX_PATH * "\\sub\\fig"
const WORKING_PATH = pwd() * "\\working"
const DATA_PATH = pwd() * "\\data"
const NCCS_PATH = DATA_PATH * "\\NCCS"
const DOCUMENTATION_PATH = DATA_PATH * "\\documentation"
const IRS_PATH = DOCUMENTATION_PATH * "\\irs"
const FRED_PATH = DATA_PATH * "\\FRED"
const WRDS_PATH = DATA_PATH * "\\WRDS"
const WRDS_DATE_FORMAT = "yyyy-mm-dd"
const FISCAL_DATE_FORMAT = "yyyymm"
const IRS_SOI_PATH = DATA_PATH * "\\IRSSOI"

const INDUSTRY_CROSS_YEAR = 2013
const COMPUSTAT_NAME = "$INDUSTRY_CROSS_YEAR-COMPUSTAT"
const SOI_NAME = "$INDUSTRY_CROSS_YEAR-SOI"
const INDUSTRY_COMPARISON_SOURCE = #=:wrds =# :soiwrds

const IRS_INDEX_NAME = "IRS List"
const HEADER_NAME = "header.tex"
const FOOTER_NAME = "footer.tex"
const METADATA_NAME = "NCCSMetadata"
const META_INDEX_NAME = "NCCSMetaIndex"
const META_PATH = DOCUMENTATION_PATH * "\\raw meta"
const SUB_TABLES_PATH = META_PATH * "\\subtables"
const FILES_INDEX_NAME = "NCCSFiles"
const KEEPER_NAME = "Keeper List"
const SUMMARY_NAME = "summaryTable"
const PROCESSED_DATA_NAME = "endowmentData"
const CONSOLIDATED_NAME = "summaryTableConsolidated"
const PRIMARY_OUTPUT_NAME = "primaryOutput"
#const Nothing = Void #NOTE: 0.7 compatabiltiy hack, delete on ugrade
const VAR_TYPES = Dict(:Numeric=>Float64, :Character=>String)
const FILE_TYPES_NAME_TO_CODE = Dict("Public Charity" => :pc,
    "Private Foundation" => :pf, "Other 501c" => :co)
const FRED_INFLATION_FILE = "inflation"
const FRED_DATE_FORMAT = "m/d/yyyy"
const TABLE_PATH = TEX_PATH * "\\sub\\tab"

const MAX_SERIALIZATION_WORKERS = 2
const MAX_MAPPING_WORKERS = 2

const MBool = Union{Bool, Missing}
const MInt = Union{Int, Missing}
const MFloat64 = Union{Float64, Missing}
const MSymbol = Union{Symbol, Missing}
const MString = Union{String, Missing}
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs

const REVENUE_CUTOFF = 1. #require at least one year with a higher value-
# set to a higher value than the lwoer bound for an effect

#reqquire all asset levels above this mark (set to $10_000_000 for testing, 50,000 breadth)
const REVENUE_LOWER_BOUND = 1_000_000.
const WEALTH_LOWER_BOUND = REVENUE_LOWER_BOUND/2. #half the revenue bound ensures reasonable net assets
const ENDOWMENT_RETURNS_NAME = "endowmentReturns"
const LARGE_DIFFERENCE_CUTOFF = .1
const MAX_RETURN_CUTOFF = 2.
const MIN_RETURN_CUTOFF = -0.9
const MIN_POINTS_PER_CATEGORY = 10
const BASE_YEAR = 2015

const MISSING_STRINGS = ["NUL", "NULL", "NIL", "missing", "None"]
const VERIFY_FIELDS = [:name, :ein, :fisyr]

const MIN_POINTS_FOR_BETA = 5

include("EndowmentCode\\CodeDictionaries.jl")

#default compressor algorithms
#I_CSV_STREAMER(p::String) = GzipDecompressorStream(open(p))
I_CSV_STREAMER(p::String) = gzopen(p)
#I_CSV_STREAMER(p::String) = Stream(format"CSV", GzipDecompressorStream(open(p)))
#I_CSV_STREAMER(p::String) = IOBuffer(read(gzopen(p)))
I_JLS_STREAMER(p::String) = GzipDecompressorStream(open(p))
O_JLS_STREAMER(p::String) = GzipCompressorStream(open(p, "w"))

#csvopen(iStream::Any)::DataFrame = load(Stream(format"CSV", iStream)) |> DataFrame
#csvsave(oStream::Any, df::T where T<:AbstractDataFrame) = save(Stream(format"CSV", oStream),df)

include("EndowmentCode\\EndowmentPar.jl") #holds supporting functions for parallelization
include("EndowmentCode\\MetaTuples.jl") #holds methods related to field-level metadata
include("EndowmentCode\\MetaDictionaries.jl") #holds methods related to file level metadata
include("EndowmentCode\\Metadata.jl") #holds all other metadata information
include("EndowmentCode\\DownloadNCCSCSVs.jl") #holds methods related to downloading the data
include("EndowmentCode\\DownloadIRSPDFs.jl") #holds methods related to downlaoding IRS files (no dependencies)
include("EndowmentCode\\NCCSFiles.jl") #holds methods related to the data files
include("EndowmentCode\\NCCSSummary.jl") #Gets summary information on the files
include("EndowmentCode\\MappingFunctions.jl")
include("EndowmentCode\\Mapping.jl")
include("EndowmentCode\\WRDSData.jl") #contains scripts relating to acquisition of wrds data
include("EndowmentCode\\EndowmentData.jl") #holds high-level scripts and other miscellaneous
include("EndowmentCode\\ComputeNCCS.jl")
include("EndowmentCode\\TopInstitutions.jl")
include("EndowmentCode\\AnalyzeBeta.jl")
include("EndowmentCode\\DescriptiveTables.jl")
include("EndowmentCode\\IndustryFigures.jl")
include("EndowmentCode\\DescriptiveFigures.jl")




end
