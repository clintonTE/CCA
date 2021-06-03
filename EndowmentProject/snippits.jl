using CodecZlib,DataFrames, uCSV
using CSV
cd("C:\\Users\\Clinton Tepper\\Dropbox\\Projects\\Endowment Project\\data\\NCCS")
p="nccs_core_2015_pc_full.csv.gz"

b = open(GzipDecompressorStream, p) do f
  b = IOBuffer("$(readline(f))\n$(readline(f))")
  (b)
end

df = CSV.read(b,# strings=:raw,
  categorical=false,  silencewarnings=true, strict=false)
nCols = size(df,2)
dfNames = (s::Symbol->Symbol(lowercase(string(s)))).(names(df))

# NOTE: Below throws an error
@time begin
  b = open(GzipDecompressorStream, p) do iStream
    b = IOBuffer(read(iStream))
    (b)
  end

  df = CSV.File(b) |> DataFrame
end

#NOTE: Below throws an error
CSV.File("nccs_core_2015_pc_full.csv") |> DataFrame


#limit=209_715

#=local df::DataFrame

try
  df = CSV.File(b) |> DataFrame
catch
  chunksize = 100_000

  df = CSV.File(b, limit=100_000)
  rowsread = 100_000

  oldsize = -1
  newsize = size(df,1)
  while oldsize ≠ newsize
    oldsize=newsize

  end
end=#




types = Dict(i=>Union{String,Missing} for i ∈ 1:nCols)

#=types = Dict(i=>Vector for i ∈ 1:nCols)

b = open(GzipDecompressorStream, p) do iStream
  b = IOBuffer(read(iStream))
  (b)
end
uCSV.read(b, coltypes=types, allowmissing=true, escape='\\', quotes='\"',
  encondings=Dict(""=>missing))=#


#=@time begin
  df = DataFrame(load(File(format"CSV", p)))
  #DataFrame(load(File(format"CSV", p), colparsers=types))
end=#

#=open("nccs_core_2013_pc_full.csv.gz") do io
    load(Stream(format"CSV", GzipDecompressorStream(io))) |> DataFrame
end=#
