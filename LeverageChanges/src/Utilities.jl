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
#=function gz2lz4(file::String; filetype::String = "csv", location::String = "",
  gzextension::String = "gz", lz4extension::String = "lz4",
  deleteold::Bool = false)

  infile::String = "$location\\$file.$filetype.$gzextension"
  outfile::String = "$location\\$file.$filetype.$lz4extension"

  isfile(infile) || error("file not found. Full path: $location\\$file.$filetype.$gzextension")


  #Below doesnt work
  indata = open(infile)
  outdata = open(outfile, "w")
  s = GzipDecompressorStream(LZ4HCCompressorStream(outdata))
  write(s, indata)
  close(s)
end=#


#converts gzip to lz4
function gz2zstd(file::String; filetype::String = "csv", location::Union{String,Nothing} = nothing,
  gzextension::String = "gz", zstdextension::String = "zstd",
  deleteold::Bool = false, level=1)

  local prefix::String

  prefix = isnothing(location) ?  "" : "$location\\"
  infile::String = "$prefix$file.$filetype.$gzextension"
  outfile::String = "$prefix$file.$filetype.$zstdextension"

  isfile(infile) || error("file not found. Full path: $prefix$file.$filetype.$gzextension")


  #Below doesnt work
  indata = open(infile)
  outdata = open(outfile, "w")
  s = GzipDecompressorStream(ZstdCompressorStream(outdata, level=level))

  try
    write(s, indata)
  finally
    close(s)
  end
end


function recompressfiles(files = FILES_TO_RECOMPRESS; datapath=DATA_PATH)::Nothing

  for f ∈ files
    gz2zstd(f, location=datapath)
  end

  return nothing
end

recompressfile(file; datapath=DATA_PATH)::Nothing = recompressfiles(
  [file]; datapath=datapath)

#check if all columns in a df or subdf are equal
#extra logic is to handle missing values
@inline function eachrowequal(df::AbstractDataFrame; cols::Vector{Symbol}=propertynames(df))::Bool
  for c ∈ eachcol(view(df, :, cols))
    completec = collect(skipmissing(c))
    (length(completec) == 0) && continue #all missing case

    #check for equality
    if  (length(completec) == length(c)) && (all(c[1] .== c))
     continue
    end

    return false
  end
  return true
end

#some useful functions for conforming numbers
@inline floororx(::Missing, ::Float64 = DENOM_FLOOR) = missing
@inline floororx(x::Real, floor::Float64 = DENOM_FLOOR) = max(x, floor)

@inline boundunit(::Missing, ::Float64=UNIT_TOL) = missing
@inline boundunit(x::Real, tol::Float64=UNIT_TOL) = max(min(x, 1.0 - tol), tol)
@inline boundunit1(x::Union{Real,Missing}) = boundunit(x, 0.0)

@inline boundabsunit(::Missing, ::Float64=UNIT_TOL) = missing
@inline boundabsunit(x::Real, tol::Float64=UNIT_TOL) = max(min(x, 1.0 -tol), -1.0 + tol)
