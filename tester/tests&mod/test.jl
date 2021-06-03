
using DataFrames, GZip, CSV


function test()::Void
  d::DataFrame = DataFrame([String, Int, Float64], [ Symbol("v",i) for i in 1:3], 0)
  push!(d, ["a" 2 3.1])
  push!(d, ["b" 33 1.1])
  println(d)

  #alternatively, append 2 dataframes
  d = DataFrame([String, Int, Float64], [ Symbol("v",i) for i in 1:3], 0)
  d = [d; DataFrame([["a","b"],[2,33],[3.1,1.1]], [:v1,:v2,:v3])]
  println(d)

  #or for best performance, use in-place broadcast loops
  d = DataFrame([String, Int, Float64], [ Symbol("v",i) for i in 1:3], 2)
  d[:, :v1] .= ["a","b"]
  d[:, :v2] .= [2,33]
  d[:, :v3] .= [3.1,1.1]
  println(d)

  #Can also just do it all in a constructor
  d = DataFrame(v1 = ["a","b"]::Vector{String},
    v2 = [2,33]::Vector{Int},
    v3 = [3.1,1.1]::Vector{Float64})
  println(d)




  writetable("outfile.csv.gz", d)

  return nothing
end

test() #No explicit delete function, but
#wrapping the dataframe in a function will delete in from memory
#if needed via the garbage collector when the function closes
#could also use other scoping mechanisms
