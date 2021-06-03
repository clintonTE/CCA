#utilities for QUtilities
#some of these are candidates for Finometrics
const CRSP_MISSINGS = [
  "A", "B", "C", "D", "E", "S", "T", "P",
  -55.0, -66.0, -77.0, -88.0, -99.0,
  -55, -66, -77, -88, -99,
  "-55.0", "-66.0", "-77.0", "-88.0", "-99.0"]

const CRSP_PARSED_MISSINGS = (-55.0, -66.0, -77.0, -88.0, -99.0)

const FRET = :ret
const CQUART_FVALS = [:mknegcash, :bknegcash, :mkflev, :bkflev, :mkliab, :bkliab]
const CQUART_FCONTROLS = :LLret12m
const CQUART_FWEIGHTONS = nothing
const CQUART_FNAMES = (s::Symbol->Symbol("P_",s)).(CQUART_FVALS)
const CQUART_FGROUPS = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES)

const CQUART_FCONTROLS_ALT = :Lmkat
const CQUART_FNAMES_ALT = (s::Symbol->Symbol(s,"_ac")).(CQUART_FNAMES)
const CQUART_FGROUPS_ALT = (s::Symbol->Symbol("G_",s)).(CQUART_FNAMES_ALT)

const USED_COMPA_VALUE_FIELDS = [:che, :ch, :at, :dltt, :dlc, :ceq, :txditc,
  :txdb, :pstk, :lt]
const USED_COMPQ_VALUE_FIELDS = [:rdq]
#const USED_COMP_VALUE_FIELDS = [USED_COMPA_VALUE_FIELDS; USED_COMPQ_VALUE_FIELDS;]
const USED_CRSP_FIELDS = [:ret, :price, :shrout]

#final columns for the decompressed file
const RETAINED_COLS_CCM = [:gvkey, :lpermno, :linkeffdate, :linkenddate]

const RETAINED_COLS_COMP = [:gvkey; :fyear; :fyrmonth;
  :linkeffdate; :linkenddate; :Nadate;
  :ddate; :adate; :lpermno; :begindate; :enddate; USED_COMPA_VALUE_FIELDS;
  :bknegcash; :cash; :netdebt; :bkequity; :bkflev; :bkliab; :mkequity; :op;]


#not sure if we really need much of these
const RETAINED_COLS_UNIV = [RETAINED_COLS_COMP;
  :compid; :permno; :date; :ret; :mktcap; :net]

#Only keep what we need, also serves as a check to make sure we ahve what we think
const RETAINED_COLS_DATA_SERIES = unique([:gvkey; :fyear; :fyrmonth;
  :linkeffdate; :linkenddate; :Nadate;
  :ddate; :adate; :lpermno; :begindate; :enddate; :mkequity; :Lmkequity; :mkat;
  CQUART_FCONTROLS; CQUART_FCONTROLS_ALT;
  CQUART_FVALS; (s->Symbol(:G_,s)).(CQUART_FNAMES); (s->Symbol(:G_,s)).(CQUART_FNAMES_ALT);
  :vol252d; :Lvol252d; :LLvol252d;
  :Dat; :LDat; :ret12m; :Lret12m; :LLret12m; :bm; :Lbm; :op; :Lop])

#Table controls
const CONTROL0 = Vector{Symbol}()
const T1_CONTROL4 = [:Lbm; :Lmkequity; :LDat; :Lop]
const T1_CONTROL10 = [T1_CONTROL4; :L12ret12m; :L24ret12m;
  :fyret12m; :L12vol252d; :L24vol252d; :fyvol252d]

#@inline Base.isfinite(::Missing)=missing

#NOTE: The below adjusts for logged vs unlogged returns
#const RETAINED_FACT_COLS = (( RETURN_FIELD == :lret) ?
#  [RETURN_FIELD; FACTOR_FIELDS; :lrfr, :lmkt] :
#  [RETURN_FIELD; FACTOR_FIELDS; :rfr, :mkt]

#consider promiting to the below to finometrics
@inline obliterate!(df::DataFrame)::DataFrame = select!(df, Not(names(df)))

#checks if a column or multiple columns are a key
@inline function checkdfkey(df::AbstractDataFrame, keys::Union{Symbol, Vector{Symbol}})
  N::Int = size(df,1)
  Nunique::Int = length(groupby(df, keys))

  (N==Nunique) || error("$keys is not a key!!")

  return nothing
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
end

function archivefile(file::String; inpath::String=DATA_PATH, inext::String = "",
  archivepath::String = ARCHIVE_PATH)

  (length(inext) > 0) && (file = "$file.$inext")
  inp::String = "$inpath\\$file"
  outp::String = "$archivepath\\$file.gz"

  open(GzipCompressorStream, outp, "w") do sout
    open(inp) do sin
      write(sout, sin)
    end
  end

  return nothing
end

function unarchivefile(file::String; outpath::String=DATA_PATH, outext::String = "",
  archivepath::String = ARCHIVE_PATH)

  (length(outext) > 0) && (file = "$file.$outext")
  outp::String = "$outpath\\$file"
  inp::String = "$archivepath\\$file.gz"

  open(GzipDecompressorStream, inp) do sin
    open(outp, "w") do sout
      write(sout, sin)
    end
  end

end

function recompressfiles(files = FILES_TO_RECOMPRESS; datapath=DATA_PATH)::Nothing

  for f ∈ files
    gz2lz4(f, location=datapath)
  end

  return nothing
end

recompressfile(file; datapath=DATA_PATH)::Nothing = recompressfiles(
  [file]; datapath=datapath)
