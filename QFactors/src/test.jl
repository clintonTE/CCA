

###############################
#NOTE: Experiments with parsing
##############################
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


###############################
#NOTE: Experiments with conversion
##############################
#=using Revise
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

end=#


###############################
#NOTE: Experiments with testing values in a df
###############################
#=
using DataStructures, DataFrames, Finometrics, BenchmarkTools

function testDS(NSet::Int = 1000, NSearched::Int = 1_000_000)
  df::DataFrame = DataFrame(cint = Vector{MInt}(rand(1:10^5, NSearched)),
    cfloat = Vector{MFloat64}(rand(1:10^5, NSearched)))

  local found::Int

  basefvec::Vector{Float64} = collect(1:NSet)
  baseivec::Vector{Int} = Vector{Int}(basefvec)

  #sftuple::NTuple{NSet, Float64} = NTuple{NSet, Float64}(basefvec)
  #situple::NTuple{NSet, Int} = NTuple{NSet, Int}(baseivec)

  sfset::Set{Float64} = Set{Float64}(basefvec)
  siset::Set{Int} = Set{Int}(baseivec)

  sfsortedset::SortedSet{Float64} = SortedSet{Float64}(basefvec)
  sisortedset::SortedSet{Int} = SortedSet{Int}(baseivec)

  print("\n")
  println("Testing Set")
  found=0
  @time found += sum((f::MFloat64->f ∈ sfset).(df.cfloat))
  @time found += sum((f::MInt->f ∈ siset).(df.cint))
  println("Set test complete found $found")

  println("Testing SortedSet")
  found=0
  @time found += sum((f::MFloat64->f ∈ sfsortedset).(df.cfloat))
  @time found += sum((f::MInt->f ∈ sisortedset).(df.cint))
  println("Set test complete found $found")

  println("Testing base")
  found=0
  @time found += sum((f::MFloat64->f ∈ basefvec).(df.cfloat))
  @time found += sum((f::MInt->f ∈ baseivec).(df.cint))
  println("Base test complete found $found")
  return nothing
end

testDS()=#


#=function mwe()
  println("first run:")
  for i ∈ 1:100:401
    basevec::Vector{Float64} = collect(1:i)
    print(i, ": ");
    @time NTuple{i, Float64}(basevec)
  end

  println("second run:")
  for i ∈ 1:100:401
    basevec::Vector{Float64} = collect(1:i)
    print(i, ": ");
    @time NTuple{i, Float64}(basevec)
  end

  println("uniontype first run")
  for i ∈ 1:18
    basevec::Vector{Union{Nothing,Float64}} = collect(1:i)
    print(i, ": ");
    @time NTuple{i, Union{Nothing,Float64}}(basevec)
  end

  println("uniontype second run")
  for i ∈ 1:18
    basevec::Vector{Union{Nothing,Float64}} = collect(1:i)
    print(i, ": ");
    @time NTuple{i, Union{Nothing,Float64}}(basevec)
  end
end
mwe()=#




#lz4mwe("data\\CRSP-D.jls", "data\\CRSP-D.csv.lz4")

##############
#NOTE: functions related to multi-threaded dataframes
################
#=using DataFrames
#WARNING: DO NOT RUN THIS

function tsmwecorrupt(N=100_000)
  df = DataFrame(rand(N,100))
  df.grpcol = (i->i%50).(1:N)

  Threads.@threads for sdf ∈ groupby(df, :grpcol)
    sdf.x3 .= -1.
  end

  println(sum(df.x3))
end

tsmwecorrupt()

function tsmwe(N=100_000)
  df = DataFrame(rand(N,100))
  df.grpcol = (i->i%50).(1:N)

  Threads.@threads for r ∈ eachrow(df)
    r.x3 = -1.
  end

  println(sum(df.x3))
end

tsmwe()


function tsmwe2(N=100_000)
  df = DataFrame(rand(N,100))
  df.grpcol = (i->i%50).(1:N)

  Threads.@threads for sdf ∈ collect(groupby(df, :grpcol))
    sdf.x3 .= -1.
  end

  println(sum(df.x3))
end


tsmwe2()=#



######################
#NOTE: Macro experiments
###################
using Random

macro tt(args...)
    println(args)

    na = length(args)
    if na != 1
        throw(ArgumentError("wrong number of arguments in @threads"))
    end
    ex = args[1]
    if !isa(ex, Expr)
        throw(ArgumentError("need an expression argument to @threads"))
    end
    if ex.head === :for
        return "ok"
    else
        throw(ArgumentError("unrecognized argument to @threads"))
    end
end


PAR = false
Random.seed!(11)
condg = true
macro mpar(cond, expr)
  quote
    if $(esc(cond))
        :($(Threads.@threads($expr)))
    else
        :($($expr))
    end
  end

end



function testmacro()

  v::Vector{Float64} = collect(1:100_000_000)
  cond::Bool = false

  @time @mpar true for i ∈ 1:length(v)
    v[i] = log10(i)
    end
  #=@mpar for i ∈ 1:length(v)
    v[i] = 1
end=#

  #println(s)
  println(sum(v))
end

testmacro()
