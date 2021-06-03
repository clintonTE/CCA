#utilities
#some of these are candidates for Finometrics


@inline function finiteormissing!(df::DataFrame, s::Symbol)::Nothing
  df[!,s] .= (f::MFloat64->coalesce(isfinite(f),false) ? f : missing).(df[!,s])

  return nothing
end

@inline function finiteormissing!(df::AbstractDataFrame, s::Symbol)::Nothing
  @mpar for r ∈ eachrow(df)
    r[s] = coalesce(isfinite(r[s]),false) ? r[s] : missing
  end

  return nothing
end

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

  isfile(infile) || error("file not found. Full path: $location\\$file.$filetype.$gzextension")


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


function uniquebyyear(df::DataFrame; Frequired::Union{Nothing, Vector{Symbol}} = nothing,
  Fidentifier::Symbol = :permno, Fgroup::Symbol = :fyear)

  for sdf::SubDataFrame ∈ groupby(df, Fgroup)
    local ssdf::SubDataFrame

    if !isnothing(Frequired)
      ssdf = view(sdf, completecases(sdf[!, Frequired]), :)
    else
      ssdf = sdf
    end

    N::Int = length(unique(ssdf[!, Fidentifier]))
    N==0 && continue
    println("$(ssdf[1, Fgroup]): $N unique records by $Fidentifier ",
      isnothing(Frequired) ? "" : "containing $Frequired")
  end
end


#=@inline function fixcsv(df::DataFrame)::DataFrame
  for c::Int ∈ 1:size(df,2)
    df[!,c] = Array(df[!,c])
  end
  return df
end=#
#=@inline function fixcsv(df::DataFrame)::DataFrame
  return DataFrame(view(df,:,:))
end=#
