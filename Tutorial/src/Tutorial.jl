module Tutorial
using Revise #useful when revising code frequently

#prerequisites:
# in the directory where you want to store projects, from the REPL type ] generate "Tutorial"
# Install the following packages by typing the following into the REPL
# ] add Distributions LinearAlgebra Random Revise CUDA GLM BenchmarkTools DataFrames

#load well established
using Distributions, LinearAlgebra, DataFrames, Random, CUDA
sleep(0.1)


#NOTE: Always wrap Julia code in functions

function basics()
  #hello world
  println("Hello world")

  #The info macro is specifically designed to echo info back to the user
  @info "Hello again"

  #Simple math and assignment
  x = 3 + 4

  #string interpolation
  @info "x=$x"

  #compound assignment
  x *= 2
  @info "now x=$x" #14

  #Many ways to create vectors
  v1 = ones(5)
  v2 = [1.0 for i ∈ 1:5]
  v3 = (_->1.0).(1:5) #the dot here means we are executing the function for each value in 1:5
                      #_->1.0 is an anonymous function that takes a value as input and returns 1.
  v4 = Vector{Float64}(undef, 5); v4 .= 1.0 #The dot here is used for element wise operations

  #Use assertion to verify something is true- it will thrown an error if it is not
  @assert v1 == v2 == v3 == v4
  @info "Since we got here, v1 == v2 == v3 == v4 = $(v1 == v2 == v3 == v4)"

  #conditionals are easy
  if iseven(x)
    @info "x is $x so it is even"
  elseif isodd(x)
    @info "x is $x so it is odd"
  else
    @assert false #this pattern is useful for debugging
  end

  #The compiler works best when each variable maintains its type
  #we can enforce this by explicitly declaring a variable's type
  local M::Matrix{Float64} #declare M is a matrix of Float64 values
  local v::Vector{Float64} #declare v is a vector of Float64 values

  #assign random values to M
  M = rand(10_000,1_000) #M is a 10^6x10 matrix of random values, drawn from the std uniform
  v = rand(Normal(), 10_000) #v is a 10^6 matrix of random values,drawn from the std normal

  #Matrix algebra is easy
  Mtvsum = M'*v |> sum #|> is a pipe operator that sends output to a single argument function
  @info "M'*v |> sum: $(Mtvsum)"

  #we can also do element wise operations.
  #The .* operator expands the dimensions of v, if possible, to allow for the element wise operatoin
  Movsum = M .* v |> sum
  @info "M .* v |> sum $(Movsum)"
  @assert Mtvsum ≈ Movsum #we can use the ≈ operator to check if the answers are within rounding error

  sqrmat(X) = X'*X #functions are easy to write and will adapt to the input type
  MtMsum = M |> sqrmat |> sum
  vtvsum = v |> sqrmat |> sum #this works even though v is a vector, not a matrix
  x2 = x |> sqrmat #it even works for scalars
  @info "MtMsum: $MtMsum, vtvsum: $vtvsum, x2=$x2"

  #Suppose we wanted to define it for strings
  try # use  atry block to catch exceptions
    sqrmat("7")
  catch err
    @warn "sqrtmat(\"7\") failed with error $err"
  end

  #we can define a method which handles trhis by specifying the type in the argument
  #note that we need to explicitly set the parsing type, so we set an optional type argument
  #this is a very high performance pattern for coding
  sqrmat(X::String, ::Type{T}=Float64) where T = parse(T, X) |> sqrmat
  floatstring2 = sqrmat("7.0")
  intstring2 = sqrmat("7", Int)
  @info "sqrmat(\"7.0\"): $floatstring2, sqrmat(\"7\", Int): $intstring2"

  #suppose we wanted to square 5 matrices
  #time the results with @time (note @btime from BenchmarkTools library does a better job)
  print("Time squaring 20 matrices")
  M20 = [rand(Float32, 10_000,1000) for i ∈ 1:20]
  M20² = (Mᵢ->similar(Mᵢ,(1000,1000))).(M20) #prellocation can improve performance
  @time begin #use begin/end to create scoping blocks
    for (m20², m20ᵢ) ∈ zip(M20², M20) #loops are fast!
      m20² .= m20ᵢ |> sqrmat
    end
  end

  #for faster performnace, we can tell Julia to multi-thread
  print("Time squaring 20 matrices (parallel)")
  M20²mt = M20² |> deepcopy#use deepcopy to recursively copy contents of objects
  @time begin #use begin/end to create scoping blocks
    Threads.@threads for i ∈ 1:20 #only simple indexing allowed in loops like this
      M20²mt[i] .= M20[i] |> sqrmat
    end
  end

  print("Time squaring 20 matrices (gpu)")
  dM20² = M20² |> deepcopy .|> CuMatrix{Float32}
  dM20 = CUDA.@sync M20 .|> CuMatrix{Float32}
  @time CUDA.@sync begin #use begin/end to create scoping blocks
    for i ∈ 1:20 #only simple indexing allowed in loops like this
      dM20²[i] .= dM20[i] |> sqrmat
    end
  end

  #The DataFrames package is essential for working with data
  df = DataFrame(x = [rand() > 0.5 ? rand() : missing for i ∈ 1:1000],
                 g = rand([1,2,3,4], 1000))
  #it basically functions like a combination of other data packages
  @info "mean(df.x): $(mean(df.x)), mean(df.x |> skipmissing): $(mean(df.x |> skipmissing))"

  #sum by group
  println(combine(groupby(df, :g), [:x,:g] .=> (c)->sum(skipmissing(c))))


  #https://dataframes.juliadata.org/stable/
  #see more advanced DataFrames features at
  #
  # https://github.com/bkamins/Julia-DataFrames-Tutorial (requires Jupiter notebooks)

  #NOT covered in this tutorial but worth learning
  # structs, constructors, parametric types - see Julia docs
  # Serialization, Dates (see Julia docs- standard library)
  #Many packages- see their respective githubs:
  # CUDA (gpu)
  # Zygote (analytical derivitvies)
  # GLM (least squares)
  # BenchmarkTools (benchmarking)
  # Optim (optimizaiton)
  # CSV (text file i/o)
  # Weave (markdown)
  #  Plots, Gadfly, VegaLite (these are various plotting packages, don't need to know all)





end


basics()
end #module end
