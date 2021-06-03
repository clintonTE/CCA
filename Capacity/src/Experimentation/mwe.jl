using Revise

using DataFrames, Dates, Random

#=function mwe(;N)
  alldates = collect(Date(1950,1,1):Day(1):Date(2020,12,31))
  df = reduce(vcat, [DataFrame(;grp = i, date=alldates[1:(rand(1:length(alldates)))]) for i ∈ 1:N])
  df.yearmonth = (year).(df.date) .+ (month).(df.date)
  df.val = rand([1,2,3,missing],nrow(df))
  sort!(df, [:grp, :yearmonth])
  dfeom = combine(groupby(df, [:grp, :yearmonth]), propertynames(df) .=> last, nrow => :ndays)

  dfeom.tokeep = trues(nrow(dfeom))
  dfbyeom = groupby(dfeom, :yearmonth)
  Threads.@threads for i ∈ length(dfbyeom)
    sdf = dfbyeom[i]
    mindays = floor(median(sdf.ndays)) |> Int
    sdf.tokeep .= (sdf.ndays .≥ mindays) .& (mindays ≥ 2)
  end

  dfeom = dfeom[dfeom.tokeep, :]


  #=print("size(df): $(size(df))")
  for s ∈ propertynames(df)
    (!(eltype(df[!,s]) <: Union{Real, Missing})) && continue
    print(" $s $(mean(skipmissing(df[!, s]))),")
  end
  print("\n")=#

  return dfeom
end

function testmwe(;loops, N)
  Random.seed!(1111)
  base = mwe(;N)
  for i ∈ 1:loops
    Random.seed!(1111)
    iter = mwe(;N)
    for s ∈ propertynames(base)
      if !all( iter[!,s] .=== base[!,s])
        throw("Column mismatch! i=$i, col=$s")
      end
    end
  end
end

@time testmwe(loops=10, N=1000)=#
#=function mweold()
  df = DataFrame(G=[:A,:B,:C,:A,:B,:C], X = rand(6))
  local combinedlast

  println(df)

  try
    combineddf = combine(last, groupby(df, :G))
  catch err
    @warn "(:) .=> last failed. Error: $err"
  end

  #workaround
  combineddf = combine(groupby(df, :G), setdiff!(propertynames(df), [:G]) .=> last)

  println(combineddf)
end

function mwe()
  df = DataFrame(X = rand(4), Y=rand(4))

  local transformed

  try
    transformed = transform(cumsum, df)
  catch err
    @warn "(:) .=> last failed. Error: $err"
  end

  #workaround
  transformed = transform(df, propertynames(df) .=> cumsum)

  println(transformed)
end


mweold()=#


using DataFrames, BenchmarkTools, Dates, StatsBase
function mwedates()
  #build the sample
  dts = reduce(vcat, [[Date(2011,11,11) + Day(i) for j in 1:10^4] for i in 1:100])
  mdts = dts |> Vector{Union{Date, Missing}}
  id = reduce(vcat, [[j for j in 1:10^4] for i in 1:100])
  df = DataFrame(date=dts, mdate = mdts, id=id)

  #shuffle
  df = df[randperm(10^6), :]


  print("sort date and id: ")
  @btime sort($df, [:date, :id])
  print("sort date(with missings) and id: ")
  @btime sort($df, [:mdate, :id])
  print("work around performance: ")
  @btime begin
    $df.mdateconverted = $df.mdate |> Vector{Date}
    sort($df, [:mdateconverted, :id])
  end
end
