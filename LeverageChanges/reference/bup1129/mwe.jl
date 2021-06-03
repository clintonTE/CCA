using DataFrames, CSV

function mwe()
  df = DataFrame(rand(1000,5))

  #seems to work without the csv usage

  println(Threads.nthreads())
  dfnocsv = deepcopy(df)
  filter!(r->r.x1<0.5, dfnocsv)
  println("size dfnocsv: ", size(dfnocsv))

  #sending the df to a CSV first breaks the code
  df |> CSV.write("test.csv")
  dfcsv = CSV.read("test.csv", copycols=false, threaded=false) |> DataFrame
  println(typeof(dfcsv.x1))
  filter!(r->r.x1<0.5, dfcsv)
  println("size dfcsv: ", size(dfcsv))

end

mwe()

@inline function fixcsv(df::DataFrame)::DataFrame
  #=for c::Int ∈ 1:size(df,2)
    df[!,c] = Array(df[!,c])
  end=#
  return DataFrame(view(df, :, :))
end
