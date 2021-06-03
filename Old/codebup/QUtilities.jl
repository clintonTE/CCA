#utilities for QUtilities
#some of these are candidates for Finometrics
const CRSP_MISSINGS = [
  "A", "B", "C", "D", "E", "S", "T", "P",
  -55.0, -66.0, -77.0, -88.0, -99.0,
  -55, -66, -77, -88, -99,
  "-55.0", "-66.0", "-77.0", "-88.0", "-99.0"]

const CRSP_PARSED_MISSINGS = (-55.0, -66.0, -77.0, -88.0, -99.0)

const USED_COMPA_VALUE_FIELDS = [:at]
const USED_COMPQ_VALUE_FIELDS = [:atq, :ltq, :ibq, :txditcq, :pstkrq, :seqq, :ceqq, :pstkq]

const RETAINED_COLS_COMP_Q = [:gvkeyq, :dateq, :fyearq]

#NOTE: uncomment when we have the missings const COMP_PARSED_MISSINGS = (nothing)

#lags a variable in the time series
#=function lagwithin!(df::DataFrame, targets::Vector{Symbol};
  period::Symbol = :fisyr,
  lags::Int=1,
  laggednames::Vector{Symbol} = (lags ≠ 1 ?
    (s::Symbol->Symbol(:L, lags, s)).(targets) : (s::Symbol->Symbol(:L, s)).(targets)),
  group=:gvkey,
  sorted::Bool=false)::Nothing
end=#

#WARNING quick and dirty hack, could create issues if multiple versions of Int8
const QQuarter = Int8
const QYear = Int16

#converts gzip to lz4
function gz2lz4OLD(infile::String; outfile::String = infile[1:(length(infile)-4)] * lowercase(
  replace(uppercase(infile[(length(infile)-3):end]),".GZ"=>".LZ4")),
  deleteold::Bool = false)

  @assert isfile(infile)

  #Below doesnt work
  indata = open(infile, "r")
  outdata = open(outfile, "w")
  s = GzipDecompressorStream(LZ4CompressorStream(outdata))
  write(s, indata)
  close(s)

  #run(``)

end

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
