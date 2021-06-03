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
using DataFrames, Distributions, GZip, Dates,
   Measures, Formatting, HTTP, ODBC, Statistics,CodecZlib,
  StatsBase, Distributions, Finometrics, Distributed,
  Serialization, uCSV, CSV, #, StatsPlots,#, PlotlyJS#, UnicodePlots
   Gadfly, Cairo#, Fontconfig#FileIO, CSVFiles

Gadfly.push_theme(:default)
#plotlyjs()

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
reportandanalyze,
getfielddf

#############################Shared Constants###########
const TEX_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\Endowment Project"
const GRAPH_PATH = TEX_PATH * "\\sub\\fig"
const GRAPH_EXT = "pdf"
const GRAPH_SUFFIX = "v2"
const GRAPH_SIZE = (1000,750)
const WORKING_PATH = pwd() * "\\working"
const DATA_PATH = pwd() * "\\data"
const NCCS_PATH = DATA_PATH * "\\NCCS"
const DOCUMENTATION_PATH = DATA_PATH * "\\documentation"
const IRS_PATH = DOCUMENTATION_PATH * "\\irs"
const FRED_PATH = DATA_PATH * "\\FRED"

const WRDS_PATH = DATA_PATH * "\\WRDS"
const WRDS_DATE_FORMAT = "yyyy-mm-dd"
const REFRESH_WRDS = false

const FISCAL_DATE_FORMAT = "yyyymm"
const IRS_SOI_PATH = DATA_PATH * "\\IRSSOI"

const FF_PATH = DATA_PATH * "\\FF"
const FF_NAME = "FF"
const FF_DATE_FORMAT = "yyyymm"

const OPEN990_PATH = DATA_PATH * "\\open990"
const OPEN990_FIELD_PATH = OPEN990_PATH * "\\fieldfiles"
const OPEN990_DATE_FORMAT = "yyyymm"
const OPEN990_FIELD_INDEX_NAME = "fieldfileindex"
const OPEN990_FIELD_DATA_NAME = "fieldfile"
#const OPEN990_SUBMITTED__ON_DATE_FORMAT = "m/d/yyyy"

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
const VAR_TYPES = Dict(:Numeric=>Float64, :Character=>String)
const FILE_TYPES_NAME_TO_CODE = Dict("Public Charity" => :pc,
    "Private Foundation" => :pf, "Other 501c" => :co)
const FRED_INFLATION_FILE = "inflation"
const FRED_DATE_FORMAT = "m/d/yyyy"
const TABLE_PATH = TEX_PATH * "\\sub\\tab"

const DECIMALS = 2
const DECIMALS_RETURN = 2

const MAX_SERIALIZATION_WORKERS = 6
const MAX_MAPPING_WORKERS = 6

const MBool = Union{Bool, Missing}
const MInt = Union{Int, Missing}
const MFloat64 = Union{Float64, Missing}
const MSymbol = Union{Symbol, Missing}
const MString = Union{String, Missing}

#const PlotContainer = T where T <: Plots.Plot #holds graphs
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs

const REVENUE_CUTOFF = 1. #require at least one year with a higher value-
# set to a higher value than the lwoer bound for an effect

#reqquire all asset levels above this mark (set to $10_000_000 for testing, 50,000 breadth)
#NOTE: Recently been setting the bound at $1_000_000
const REVENUE_LOWER_BOUND = 500_000.0
const WEALTH_LOWER_BOUND = 1_000_000.0

# NOTE: Below limit leads to a hard drop for records below this mark.
# Roughly deflation adjusted 1988(EOY-126.1)-> 2015 (EOY-236.525)
const PREFILTER_REVENUE_LOWER_BOUND = REVENUE_LOWER_BOUND * 0.5
const PREFILTER_WEALTH_LOWER_BOUND = WEALTH_LOWER_BOUND * 0.5

const ENDOWMENT_RETURNS_NAME = "endowmentReturns"
const LARGE_DIFFERENCE_CUTOFF = .2
const MAX_RETURN_CUTOFF = 1.
const MIN_RETURN_CUTOFF = -0.9
const MIN_POINTS_PER_CATEGORY = 50
const BASE_YEAR = 2015

#technically this limit is on gross receipts s.t. gross receipts > revenue
#currently warns if breached
const LOWER_EZ_REVENUE_LIMIT = 50_000.0

const PERSISTENCE_PERIODS = [1,2,3,4]
const PERSISTENCE_YEARS = [2015]
const PERSISTENCE_CUTOFFS = [0.3333, 0.6667, 1.0]
const LBENCHMARKS = [nothing, :bmsp500, :bmff3, :bmff5]
const TOO_SMALL_TO_LOG = 1e-3 #anything smaller than this is more likely to be a decimal error

const CROSSTAB_CUTOFFS = [0.5, 1.0]
const CROSSTAB_YEAR_RANGE = 2009:2015

const MISSING_STRINGS = ["NUL", "NULL", "NIL", "missing", "None"]
const VERIFY_FIELDS = [:name, :ein, :fisyr]
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]

const MIN_POINTS_PLUS_K_FOR_BETA = 1
const WINSORIZE_PROP = 0.025

include("EndowmentCode\\CodeDictionaries.jl")

#default compressor algorithms
#I_CSV_STREAMER(p::String) = GzipDecompressorStream(open(p))
#I_CSV_STREAMER(p::String) = gzopen(p)
I_CSV_STREAMER(p::String) = open(GzipDecompressorStream, p)
I_CSV_STREAMER(F::Function, p::String) = open(F, GzipDecompressorStream, p)
#I_CSV_STREAMER(F::Function, p::String) = gzopen(F, p, "r")
#I_CSV_STREAMER(p::String) = Stream(format"CSV", GzipDecompressorStream(open(p)))
#I_CSV_STREAMER(p::String) = IOBuffer(read(gzopen(p)))
I_JLS_STREAMER(p::String) = open(GzipDecompressorStream, p)
I_JLS_STREAMER(F::Function, p::String) = open(F, GzipDecompressorStream, p)
O_JLS_STREAMER(p::String) = GzipCompressorStream(open(p, "w"))
O_JLS_STREAMER(F::Function, p::String) = open(F, GzipCompressorStream, p, "w")

#csvopen(iStream::Any)::DataFrame = load(Stream(format"CSV", iStream)) |> DataFrame
#csvsave(oStream::Any, df::T where T<:AbstractDataFrame) = save(Stream(format"CSV", oStream),df)

include("EndowmentCode\\EndowmentPar.jl") #holds supporting functions for parallelization
include("EndowmentCode\\DownloadAndValidate\\MetaTuples.jl") #holds methods related to field-level metadata
include("EndowmentCode\\DownloadAndValidate\\MetaDictionaries.jl") #holds methods related to file level metadata
include("EndowmentCode\\DownloadAndValidate\\Metadata.jl") #holds all other metadata information
include("EndowmentCode\\DownloadAndValidate\\DownloadNCCSCSVs.jl") #holds methods related to downloading the data
include("EndowmentCode\\DownloadAndValidate\\DownloadIRSPDFs.jl") #holds methods related to downlaoding IRS files
include("EndowmentCode\\DownloadAndValidate\\NCCSFiles.jl") #holds methods related to the data files
include("EndowmentCode\\DownloadAndValidate\\NCCSSummary.jl") #Gets summary information on the files
include("EndowmentCode\\PreprocessAndMap\\MappingFunctions.jl")
include("EndowmentCode\\PreprocessAndMap\\Mapping.jl")
include("EndowmentCode\\ReportAndAnalyze\\WRDSData.jl") #contains scripts relating to acquisition of wrds data
include("EndowmentCode\\DownloadAndValidate\\EndowmentData.jl") #holds high-level scripts and other miscellaneous

include("EndowmentCode\\ReportAndAnalyze\\FinancialData.jl")
include("EndowmentCode\\ReportAndAnalyze\\Open990.jl")
include("EndowmentCode\\ReportAndAnalyze\\ComputeNCCS.jl")
include("EndowmentCode\\ReportAndAnalyze\\TopInstitutions.jl")
include("EndowmentCode\\ReportAndAnalyze\\AnalyzeBeta.jl")
include("EndowmentCode\\ReportAndAnalyze\\Persistence.jl")
include("EndowmentCode\\ReportAndAnalyze\\DescriptiveTables.jl")
include("EndowmentCode\\ReportAndAnalyze\\IndustryFigures.jl")
include("EndowmentCode\\ReportAndAnalyze\\DescriptiveFigures.jl")
include("EndowmentCode\\ReportAndAnalyze\\CrossTab.jl")
include("EndowmentCode\\ReportAndAnalyze\\SensitivitySpecifications.jl")
include("EndowmentCode\\ReportAndAnalyze\\SensitivityRegressions.jl")
include("EndowmentCode\\ReportAndAnalyze\\ExtractTails.jl")
include("EndowmentCode\\ReportAndAnalyze\\ReportAndAnalyze.jl")





end
