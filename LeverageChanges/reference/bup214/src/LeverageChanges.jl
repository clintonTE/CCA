using Revise
module LeverageChanges
using Revise

using  CSV, DataFrames, Serialization, Dates, Statistics,
  Distributed, Finometrics, StatsBase,
  CodecZstd, CodecZlib,  Distributions, LinearAlgebra,
  CuArrays,#Plots,
  #Gadfly, Cairo, Fontconfig, Compose
  VegaLite

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end


import Base: eltype, getindex, iterate, length, isfinite, merge

#necessary due to exporting conflict
describe(df::DataFrame) = DataFrames.describe(df)


#############################CONSTANTS#####################
const IVO_REP = true
TEST_OUTPUT = true

const OUT_PATH = pwd() * "\\output"
const DATA_PATH = pwd() * "\\data"
const ARCHIVE_PATH = pwd() * "\\working"
const TEX_PATH = "C:\\Users\\Clinton\\Dropbox\\Apps\\Overleaf\\LeverageChangesRep"
const TABLE_PATH = "$TEX_PATH\\sub\\tab"
const FIGURE_PATH = "$TEX_PATH\\sub\\fig"
const DEFAULT_DECIMALS = 2

const CRSP_DATE_FORMAT = dateformat"yyyymmdd"

const CRSPM_NAME = "CRSP-M"
#const CRSPM_PARTS = ["CRSP-M"]

const CRSPD_NAME = "CRSP-D"
#const CRSPD_PARTS = ["CRSP-D-II", "CRSP-D-I"]

const CRSPD_EXPANDED_NAME = "CRSP-DEX"
#const CRSPD_EXPANDED_PARTS = ["CRSP-D-II", "CRSP-D-I", "CRSP-D-0"]

const COMP_NAME = "COMP"
const COMP_A_NAME = "COMP-A"
const COMP_Q_NAME = "COMP-Q"
const COMP_DATE_FORMAT = dateformat"yyyymmdd"

const CCM_NAME = "CCM"
const CCM_DATE_FORMAT = dateformat"yyyymmdd"
const LAST_DATE_STRING = "21001231"
const LAST_DATE = Date(LAST_DATE_STRING, CCM_DATE_FORMAT)

const MONTHS_2_STALE_LONG = 18
const MONTHS_2_STALE_SHORT = 6 #maybe mess with this?
const WINSORIZE_PROP = 0.001

const UNIV_NAME = "UNIV"
const LEV_NAME = "LEV"
const PORT_NAME = "PORT"
const PANEL_NAME = "PANEL"

const VALIDATE_MERGED = true

const USE_QUARTERLY_DATES = (!IVO_REP) #NOTE


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

const EVENT_NAME = "event"
const EVENT_PATH = "$DATA_PATH\\event"
const CRSP_DISTRIBUTIONS_NAME  = "crspdistributions"

const CUMEX_NAME = "cumex"

const TAXES_PATH = "$DATA_PATH\\taxes"
const TAXES_NAME = "histtaxes-y"

#currently a bug prevents threaded=true from allowing a mutable dataframe
const CSV_THREADED = true

#Missing vlaues
const CRSP_MISSINGS = [
  "A", "B", "C", "D", "E", "S", "T", "P",
  -55.0, -66.0, -77.0, -88.0, -99.0,
  -55, -66, -77, -88, -99,
  "-55.0", "-66.0", "-77.0", "-88.0", "-99.0"]

const CRSP_PARSED_MISSINGS = (-55.0, -66.0, -77.0, -88.0, -99.0)


###Constants related to portfolio speicifications and sorts

#=struct LeverageSpec
  Fvals::Vector{Symbol}
  Fcontrols::Symbol
  Fweightons::NSymbol
  Fnames::Vector{Symbol}
  Fgroups::Vector{Symbol}
end

function LeverageSpec(;
  Fvals::Vector{Symbol} = error("Fvals is required"),
  Fcontrols::Symbol = error("Fcontrols is required"),
  Fweightons::NSymbol = nothing,
  Fnames::Vector{Symbol} = error("Fnames is required"),
  Fgroups::Vector{Symbol} = error("Fgroups is required")
  ) = LeverageSpec(Fvals, Fcontrols, Fweightons, Fnames, Fgroups)

CQUART_FVALS
const CQUART = LeverageSpec(
 Fvals = [:mknegcash, :bknegcash, :mkflev, :bkflev, :mkliab, :bkliab],
 Fcontrols = :LLret12m,
 Fnames = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS),
 Fgroups = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES))

const CQUART_ALT = LeverageSpec(
 Fvals = CQUART.Fvals,
 Fcontrols = :Lmkat,
 Fnames = (s::Symbol->Symbol(s,"_ac")).(CQUART_FNAMES),
 Fgroups = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_ALT))

const CQUART_D =#

###Main specification
const FRET = :ret
const CQUART_FVALS = [:mknegcash, :bknegcash, :mkflev, :bkflev, :mkliab, :bkliab]
const CQUART_FCONTROLS = :LLret12m
const CQUART_FWEIGHTONS = nothing
const CQUART_FNAMES = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS)
const CQUART_FGROUPS = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES)

const CQUART_FCONTROLS_ALT = :Lmkat
const CQUART_FNAMES_ALT = (s::Symbol->Symbol(s,"_ac")).(CQUART_FNAMES)
const CQUART_FGROUPS_ALT = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_ALT)

###differenced specification
const CQUART_FVALS_D = [:Dmknegcash, :Dbknegcash, :Dmkflev, :Dbkflev, :Dmkliab, :Dbkliab]
const CQUART_FCONTROLS_D = :LLret12m
const CQUART_FWEIGHTONS_D = nothing
const CQUART_FNAMES_D = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS_D)
const CQUART_FGROUPS_D = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_D)

const CQUART_FCONTROLS_ALT_D = :Lmkat
const CQUART_FNAMES_ALT_D = (s::Symbol->Symbol(s,"_ac")).(CQUART_FNAMES_D)
const CQUART_FGROUPS_ALT_D = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_ALT_D)

###managerial vs market changes
const CQUART_FVALS_RMD = [:RDmkflev, :MDmkflev, :RDmkliab, :MDmkliab]
const CQUART_FCONTROLS_RMD = :LLret12m
const CQUART_FWEIGHTONS_RMD = nothing
const CQUART_FNAMES_RMD = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS_RMD)
const CQUART_FGROUPS_RMD = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_RMD)

const CQUART_FCONTROLS_ALT_RMD = :Lmkat
const CQUART_FNAMES_ALT_RMD = (s::Symbol->Symbol(s,"_ac")).(CQUART_FNAMES_RMD)
const CQUART_FGROUPS_ALT_RMD = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_ALT_RMD)

###BRS measures
const CQUART_FVALS_BRS_D = [:brsd, :brse, :brsx]
const CQUART_FCONTROLS_BRS_D = :LLret12m
const CQUART_FWEIGHTONS_BRS_D= nothing
const CQUART_FNAMES_BRS_D = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS_BRS_D)
const CQUART_FGROUPS_BRS_D = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_BRS_D)

#leverage + vol increases
const CQUART_FNAMES_FD = (s::Symbol->Symbol("FP_",s)).(CQUART_FVALS_D)
const CQUART_FGROUPS_FD = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_FD)

const USED_COMPA_VALUE_FIELDS = [:che, :ch, :at, :dltt, :dlc, :ceq, :txditc,
  :txdb, :pstk, :lt]
const USED_COMPQ_VALUE_FIELDS = [:rdq]
const USED_CRSP_FIELDS = [:ret, :price, :shrout]


###Retained column constants
#final columns for the decompressed file
const RETAINED_COLS_CCM = [:gvkey, :lpermno, :linkeffdate, :linkenddate, :cusip]

const RETAINED_COLS_COMP = [:gvkey; :fyear; :fyrmonth; :fdate; :Lfdate;
  :linkeffdate; :linkenddate; :Nadate; :brse; :brsd; :brsx; :netliab;
  :ddate; :adate; :lpermno; :begindate; :enddate; USED_COMPA_VALUE_FIELDS;
  :bknegcash; :cash; :netdebt; :bkequity; :bkflev; :bkliab; :mkequity; :nop;]

#not sure if we really need all of thes e fields
const RETAINED_COLS_UNIV = [RETAINED_COLS_COMP;
  :compid; :permno; :date; :ret; :mktcap; :net]

#Only keep what we need, also serves as a check to make sure we ahve what we think
const RETAINED_COLS_DATA_SERIES = unique([:gvkey; :fyear; :fyrmonth;
  :linkeffdate; :linkenddate; :Nadate; :brse; :brsd; :brsx;
  :ddate; :adate; :lpermno; :begindate; :enddate; :mkequity; :Lmkequity; :mkat;:lmkequity; :Llmkequity;
  CQUART_FCONTROLS; CQUART_FCONTROLS_ALT;
  CQUART_FVALS; (s->Symbol(:G_,s)).(CQUART_FNAMES); (s->Symbol(:G_,s)).(CQUART_FNAMES_ALT);
  CQUART_FCONTROLS_D; CQUART_FCONTROLS_ALT_D;
  CQUART_FVALS_D; (s->Symbol(:G_,s)).(CQUART_FNAMES_D); (s->Symbol(:G_,s)).(CQUART_FNAMES_ALT_D);
  CQUART_FCONTROLS_RMD; CQUART_FCONTROLS_ALT_RMD;
  CQUART_FVALS_BRS_D; (s->Symbol(:G_,s)).(CQUART_FNAMES_BRS_D); CQUART_FCONTROLS_BRS_D;
  CQUART_FVALS_RMD; (s->Symbol(:G_,s)).(CQUART_FNAMES_RMD); (s->Symbol(:G_,s)).(CQUART_FNAMES_ALT_RMD);
  CQUART_FGROUPS_FD;
  :vol252d; :Lvol252d; :LLvol252d; :LDvol252d;
  :Dlat; :LDlat; :ret12m; :Lret12m; :LLret12m; :bm; :Lbm; :nop; :Lnop;
  (s->Symbol(:LD, s)).(CQUART_FVALS)])

#Table controls
const CONTROL0 = Vector{Symbol}()
const T1_CONTROL4 = [:Lbm; :Llmkequity; :LDlat; :Lnop]
const T1_CONTROL10 = [T1_CONTROL4; :L12ret12m; :L24ret12m;
  :fyret12m; :L12vol252d; :L24vol252d; :fyvol252d]
############################IVO constants

#the following constants are used in units of millions
const DENOM_FLOOR = 0.1
const UNIT_TOL = 0.001
const MIN_AT = 10.
const MIN_BKEQUITY = 1.
const MIN_MKEQUITY = 1.
#############################Multi-threading

#available lock if needed- use lock(spinlock) unlock(spinlock) etc
const spinlock = Threads.SpinLock()

#####################REFERENCE CONSTANTS
const PARALLEL = Ref{Bool}(true)
const REPLICATION_TYPE = Ref{Symbol}(:daily)
const DIST_TYPE = Ref{Symbol}(:full)
const QRTYPE = Ref{Type}(Matrix{Float64})
const YEAR_RANGE = Ref{UnitRange{Int}}(1961:2018)
const OUT_SUFFIX = Ref{String}("all")

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
#=IN_CSV_STREAM(p::String) = LZ4SafeDecompressorStream(open(p))
IN_CSV_STREAM(F::Function, p::String) = open(F, LZ4SafeDecompressorStream, p)

IN_JLS_STREAM(p::String) = LZ4SafeDecompressorStream(open(p))
IN_JLS_STREAM(F::Function, p::String) = open(F, LZ4SafeDecompressorStream, p)

OUT_JLS_STREAM(p::String) = LZ4FastCompressorStream(open(p, "w"))
OUT_JLS_STREAM(F::Function, p::String) = open(F, LZ4FastCompressorStream, p, "w")=#

IN_CSV_STREAM(p::String) = ZstdDecompressorStream(open(p))
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
function OUT_BIN_STREAM(p::String, obj::Any; level=1)

  io = ZstdCompressorStream(open(p, "w"), level=level)
  try
    serialize(io, obj)
  finally
    close(io)
  end
end



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
include("Analysis\\FactorTests.jl")
include("Tables\\Table1.jl")
include("Tables\\Table2.jl")
include("Tables\\Table3.jl")
include("Tables\\Table4.jl")
include("Tables\\Table5.jl")
include("Tables\\Table6.jl")
include("Tables\\Table7.jl")
include("Tables\\Table8.jl")
include("Figures\\Figure1.jl")
include("Figures\\Figure2.jl")
include("Tables\\Table9.jl")
include("Figures\\Figure3.jl")
include("Tables\\Table10.jl")
include("Figures\\Figure4.jl")
include("Tables\\TableA1.jl")
include("Tables\\TableA2.jl")
include("Analysis\\Replicate.jl")

include("Event\\SEO.jl")
include("Event\\Distributions.jl")
include("Event\\FormEvent.jl")
include("Event\\FormCumex.jl")

end # module
