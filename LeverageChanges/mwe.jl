using DataFrames, CSV

#IN_CSV_STREAM(p::String) = LZ4DecompressorStream(open(p))
#IN_CSV_STREAM(F::Function, p::String) = open(F, LZ4DecompressorStream, p)
#OUT_CSV_STREAM(p::String) = LZ4CompressorStream(open(p, "w"))
#OUT_CSV_STREAM(F::Function, p::String) = open(F, LZ4CompressorStream, p, "w")

function mwe(N)
  local dfcsv::DataFrame
  #df = DataFrame(x1 = rand((string).(1:N),N), x2 = rand((string).(1:N),N))
  f = "test.csv"

  #sending the df to a CSV first breaks the code
  #open(LZ4CompressorStream, f, "w") do s
  #  CSV.write(s, df)
  #end
  #df |> CSV.write |> LZ4CompressorStream(open(f, "w"))

  #open(LZ4DecompressorStream, "$(pwd())\\data\\COMP-A.csv.lz4") do s
  #  dfcsv = CSV.read(s, threaded=true) |> DataFrame
  #end
  #df |> CSV.write(f)

  #open("$(pwd())\\data\\COMP-A.csv") do s
  #  dfcsv = CSV.read(s, threaded=false) |> DataFrame
  #end

  #dfcsv[1:75_000,8:9]  |> CSV.write(f)
  open(f) do s
    dfcsv = CSV.read(s, threaded=false) |> DataFrame
  end
  @info "unthreaded read seems to work fine"
  open(f) do s
    dfcsv = CSV.read(s, threaded=true) |> DataFrame
  end
  @info "Never gets here because the threaded read fails"

  #dfcsv = CSV.read("test.csv", copycols=false, threaded=true) |> DataFrame
  println("size dfcsv: ", size(dfcsv))

end

function mwe2()
  local dfcsv::DataFrame
  f = "test.csv"

  open(f) do s
    dfcsv = CSV.read(s, threaded=false) |> DataFrame
  end
  @info "unthreaded read seems to work fine"
  open(f) do s
    dfcsv = CSV.read(s, threaded=true) |> DataFrame
  end
  @info "Never gets here because the threaded read fails"

  println("size dfcsv: ", size(dfcsv))

end

@inline function fixcsv(df::DataFrame)::DataFrame
  #=for c::Int ∈ 1:size(df,2)
    df[!,c] = Array(df[!,c])
  end=#
  return DataFrame(view(df, :, :))
end

function testtrailing(
  ::Val{Teom}=Val(false);#allows for static-time dispath, maybe improves performance
  offset::Int=0,
  minpointsperwindow::Int = 0, iter::Int = 100_000) where Teom

  #valTeom::Int = Teom
  #println(typeof(valTeom))
  testvec::Vector{Int} = [Teom, 2]
  tot = sum(rand(testvec, iter))
  println("tot: $tot")
  return nothing
end

function testtrailingbase(
  Teom::Bool = false;#allows for static-time dispath, maybe improves performance
  offset::Int=0,
  minpointsperwindow::Int = 0, iter::Int = 100_000)

  #valTeom::Int = Teom
  #println(typeof(valTeom))
  testvec::Vector{Int} = [Teom, 2]
  tot = sum(rand(testvec, iter))
  println("tot: $tot")
  return nothing
end


@time testtrailing(iter=100_000_000)
