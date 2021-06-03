
#=function parsestring(T::Type, s::AbstractString)
  v::Union{Nothing,T} = tryparse(T,s)

  return v
end

parsestring(Int, "3")
parsestring(Float64, "3.4")

function parsestringFloat64(s::AbstractString)
  v::Union{Nothing,Float64} = tryparse(Float64,s)

  return v
end

function parsestringInt(s::AbstractString)
  v::Union{Nothing,Int} = tryparse(Int,s)

  return v
end

parsestringInt("3")
parsestringFloat64("3.4")

function parsestring(s::AbstractString, ::Val{T}) where T
  v::Union{Nothing,T} = tryparse(T,s)

  return v
end

parsestring("3", Val(Int))
parsestring("3.4", Val(Float64))=#


#=function testreadcomp()
  df = GzipDecompressorStream(open("$(pwd())\\data\\CRSP-M.csv.gz")) |> CSV.read |> DataFrame
  println("compressed-size: $(size(df))")
end

function testread()
  df = CSV.read("$(pwd())\\data\\CRSP-M.csv") |> DataFrame
  #ostream = GzipCompressorStream(open("$(pwd())\\data\\CRSP-M.jls.gz", "w"))
  #serialize(ostream,df)
  #close(ostream)
  println("direct-size: $(size(df))")
end

function testreadbin()
  df = deserialize("$(pwd())\\data\\CRSP-M.jls")
  println("direct-size: $(size(df))")
end

function testreadbincomp()


b = open(GzipDecompressorStream, "$(pwd())\\data\\CRSP-M.jls.gz") do iStream
  b = IOBuffer(read(iStream))
  (b)
end

df = deserialize(b)

  println("direct-size: $(size(df))")
end

@time testreadcomp()
@time testread()
@time testreadbin()
@time testreadbincomp()=#
#@time gz2lz4("data\\CS-A.csv.gz")



#=
using CodecLz4, CodecZlib

function gz2lz4mwe(infile::String, outfile::String)

  @assert isfile(infile)

  indata = open(infile, "r")
  outdata = open(outfile, "w")

  s = GzipDecompressorStream(LZ4CompressorStream(outdata))

  write(s, indata)
  close(s)
end

gz2lz4mwe("data\\CS-A.csv.gz", "data\\CS-A.csv.lz4")=#

using Revise
using DataFrames, CodecLz4, Serialization

if  pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

function lz4mwe(infile::String, outfile::String)
  @assert isfile(infile)
  df = deserialize(infile)
  println("size after read: $(size(df))")
  outdata = open(outfile, "w")

  s = LZ4CompressorStream(open(outfile, "w"))

  serialize(s, df)
  close(s)

  s = LZ4DecompressorStream(open(outfile))

  df = deserialize(s)

  close(s)
  println("size after decompression: $(size(df))")

end


#lz4mwe("data\\CRSP-D.jls", "data\\CRSP-D.csv.lz4")
