using Revise
module BloombergData
using Revise
using BLPData, DataFrames, Dates, CSV, BenchmarkTools, Serialization, CodecZstd

#note- expand on the below
PARAM = Dict{Symbol, Any}(
  :inputpath=>"$(pwd())\\input",
  :outputpath=>"$(pwd())\\output",
)

function IN_BIN_STREAM(p::String)
  local obj
  open(ZstdDecompressorStream, p) do io
      obj = deserialize(io)
  end

  return obj
end

IN_SMALL_CSV_STREAM(p::String) = open(p)
IN_SMALL_CSV_STREAM(F::Function, p::String) = open(F, p)

#write an uncompressed CSV file
const SMALL_CSV_EXTENSION = "csv"
function OUT_SMALL_CSV_STREAM(p::String, obj)
  io = open(p, "w")
  try
    obj |> CSV.write(io)
  finally
    close(io)
  end
end


#write a compressed CSV file
const CSV_EXTENSION = "csv.zstd"
function OUT_CSV_STREAM(p::String, obj; level=1)
  io = ZstdCompressorStream(open(p, "w"), level=level)
  try
    obj |> CSV.write(io)
  finally
    close(io)
  end
end

const BIN_EXTENSION = "bin.zstd"
function OUT_BIN_STREAM(p::String, obj, level=1)

  io = ZstdCompressorStream(open(p, "w"), level=level)
  try
    serialize(io, obj)
  finally
    close(io)
  end
end

alock = ReentrantLock()


include("scripts\\CapacityFutures.jl")
include("scripts\\CapacityIntra.jl")
include("Session.jl")

end # module
