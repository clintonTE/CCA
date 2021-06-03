#utilities for QUtilities
#some of these are candidates for Finometrics
const CRSP_MISSINGS = [
  "A", "B", "C", "D", "E", "S", "T", "P",
  -55.0, -66.0, -77.0, -88.0, -99.0,
  -55, -66, -77, -88, -99,
  "-55.0", "-66.0", "-77.0", "-88.0", "-99.0"]

const CRSP_PARSED_MISSINGS = (-55.0, -66.0, -77.0, -88.0, -99.0)
const CRSP_PARSED_MISSINGS = (-55.0, -66.0, -77.0, -88.0, -99.0)

const USED_COMPA_VALUE_FIELDS = [:at]
const USED_COMPQ_VALUE_FIELDS = [:atq, :ltq, :ibq, :txditcq, :pstkrq,
  :seqq, :ceqq, :pstkq]
const USED_CRSP_FIELDS = [:date, :ret, :mktcap]

#final columns for the decompressed file
const RETAINED_COLS_CRSP = [:permno, :date, :ret, :lret, :mktcap, :Lmktcap]

const RETAINED_COLS_COMP = [:begindate, :enddate,
  :gvkey, :adateq, :Nadateq,
  :fyear, :fqtr, :fyearqtr, :fdate, :fdateq,
    :cyearqtr, :ddate, :ddateq, :dyearq,
  :ia, :iaq,  :roeq]

const RETAINED_COLS_CCM = [:gvkey, :lpermno, :linkeffdate, :linkenddate]

#not sure if we really need much of these
const RETAINED_COLS_UNIV = [:rid, :permno, :date, :adateq, :ddateq, :fyear,
  :fyearqtr, :ret, :lret, :mktcap, :Lmktcap, :ia, :iaq, :roeq]

const RETAINED_VWCRSP_COLS = [:date, :vwindd]
const RETAINED_RFCRSP_COLS = [:date, :t30ind]



#these are the fields from which to build the factors
const Fia = :iaq
const Fmktcap = :Lmktcap
const Froe = :roeq

const USE_LOGS = false

const WEIGHT_ON = :Lmktcap
const RETURN_FIELD = USE_LOGS ? :lret : :ret

#NOTE: Below are portfolio construction parameters
const FACTOR_FIELDS = [Fmktcap, Fia, Froe]
const FACTOR_N_GROUPS = Dict(Fmktcap=>2, Fia=>3, Froe=>3)
const FACTOR_BREAKPOINTS = Dict(Fmktcap=>(0.5,1.0),
  Fia=>(0.3,0.7,1.0), Froe=>(0.3,0.7,1.0))
const FACTOR_GROUP_WEIGHTS = Dict(Fmktcap=>(1/9,-1/9),
  Fia=>(1/6, 0.0, -1/6), Froe=>(-1/6,0.0,1/6))
const FACTOR_VALUE_TYPES = Dict(Fia=>MFloat64, Froe=>MFloat64, Fmktcap=>MInt)

#NOTE: The below adjusts for logged vs unlogged returns
#const RETAINED_FACT_COLS = (( RETURN_FIELD == :lret) ?
#  [RETURN_FIELD; FACTOR_FIELDS; :lrfr, :lmkt] :
#  [RETURN_FIELD; FACTOR_FIELDS; :rfr, :mkt]

#converts gzip to lz4
function gz2lz4(file::String; filetype::String = "csv", location::String = "",
  gzextension::String = "gz", lz4extension::String = "lz4",
  deleteold::Bool = false)

  infile::String = "$location\\$file.$filetype.$gzextension"
  outfile::String = "$location\\$file.$filetype.$lz4extension"

  @assert isfile(infile)


  #Below doesnt work
  indata = open(infile)
  outdata = open(outfile, "w")
  s = GzipDecompressorStream(LZ4CompressorStream(outdata))
  write(s, indata)
  close(s)

  #run(``)

end

FILES_TO_RECOMPRESS = ["COMP-A", "COMP-Q", "CCM", "CRSP-M", "CRSP-D"]
function recompressqfactorfiles(files = FILES_TO_RECOMPRESS; datapath=DATA_PATH)::Nothing

  for f ∈ files
    gz2lz4(f, location=datapath)
  end

  return nothing
end

recompressqfactorfile(file; datapath=DATA_PATH)::Nothing = recompressqfactorfiles(
  [file]; datapath=datapath)
