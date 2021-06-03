
#useful macro for conditionally running things in parallel
#=macro mpar(cond, expr)
  quote
    if $(esc(cond))
        :($(Threads.@threads($expr)))
    else
        :($($expr))
    end
  end
end



const spinlock = Threads.SpinLock()

using Dates, DataFrames



function turtlesort!(::Type{T}, df::AbstractDataFrame, sorts::Vector{Symbol},
  threaded::Bool=false, parentdf::DataFrame = df,
  lt::Function = (a::Tuple{T, SubArray{Int}},b::Tuple{T, SubArray{Int}})->a[1]<b[1]) where T


  target::Symbol = sorts[1]
  gdf::GroupedDataFrame = groupby(df, sorts[1])

  #map the index portions to the sub-arrays
  Ntarget::Int = length(gdf)
  sturtleidxs::Vector{SubArray{Int}} = (i->view(parentdf.turtleidx, gdf[i].turtleidx)).(1:Ntarget)

  #Sort any other sort columns
  if length(sorts) > 1
    @mpar threaded for i ∈ 1:length(gdf)
      turtlesort!(eltype(df[!, sorts[2]]), gdf[i], sorts[2:end],  false, parentdf)
    end
  end

  #now do our part of the sort
  targetv::Vector{T} = unique((i::Int->gdf[i][1,target]).(1:Ntarget))

  #map the values to the index array
  targetvindex::Vector{Tuple{T, SubArray{Int}}} = collect(zip(targetv, sturtleidxs))

  #let Julia do this part of the sort, probably can do it better than me
  sort!(targetvindex, lt = lt)

  #next part is the key
  df.turtleidx .= [(vindex::Tuple{T, SubArray{Int}}->vindex[2]).(targetvindex)...;]

  return nothing
end
function turtlesort(df::DataFrame, sorts::Vector{Symbol}; threaded::Bool=false)
  #turtleidx::Vector{Int} = collect(1:size(df,1))
  df.turtleidx = collect(1:size(df,1))
  turtlesort!(eltype(df[:,sorts[1]]), df, sorts, threaded)

  dfout::DataFrame = df[df.turtleidx,:]
  select!(dfout, Not(:turtleidx))
  return dfout
end

function turtlesort!(df::DataFrame, sorts::Vector{Symbol}; threaded::Bool=false)
  df.turtleidx = collect(1:size(df,1))
  turtlesort!(eltype(df[:,sorts[1]]), df, sorts, threaded)

  #the turtle is upside down for our purposes! Flip it back
  #(e.g. turtleidx is a correctly ordered array of row pointers,
  #but we need to associate each row pointed at with its correct order
  df.turtlebiject = df.turtleidx[df.turtleidx]
  sort!(df, :turtlebiject)
  select!(df, Not([:turtlebiject, :turtleidx]))

  return df
end



#------------------------------------------

using DataFrames, Dates, Finometrics
function testsorts(N::Int = 10_000_000, threaded::Bool = false)
  local df::DataFrame
  local dfs1::DataFrame
  local dfs2::DataFrame
  local dfs3::DataFrame
  local dfs4::DataFrame

  df = DataFrame(ones(N,100))
  df.rid = rand(1:25000,N)
  df.dt = rand(Date(1980,1,1):Day(1):Date(2018,12,31),N)
  df.rid2 = rand(1:1000,N)

  dfs1 = deepcopy(df)
  print("\ntime sort!: ")
  @time sort!(dfs1, [:rid, :dt, :rid2])

  Base.GC.gc()

  dfs2 = deepcopy(df)
  print("time sort: ")
  @time dfs2 = sort(dfs2, [:rid, :dt, :rid2])

  Base.GC.gc()

  dfs3 = deepcopy(df)
  print("time turtlesort!: ")
  @time turtlesort!(dfs3, [:rid, :dt, :rid2], threaded=threaded)

  Base.GC.gc()

  dfs4 = deepcopy(df)
  print("time turtlesort: ")
  @time dfs4 = turtlesort(dfs4, [:rid, :dt, :rid2], threaded=threaded)


  @assert sum(dfs1.rid .== dfs2.rid)==N
  @assert sum(dfs1.dt .== dfs2.dt)==N
  @assert sum(dfs1.rid2 .== dfs2.rid2)==N

  @assert sum(dfs1.rid .== dfs3.rid)==N
  @assert sum(dfs1.dt .== dfs3.dt)==N
  @assert sum(dfs1.rid2 .== dfs3.rid2)==N

  @assert sum(dfs1.rid .== dfs4.rid)==N
  @assert sum(dfs1.dt .== dfs4.dt)==N
  @assert sum(dfs1.rid2 .== dfs4.rid2)==N

  return nothing

end

testsorts(10_000)

#using StatsModels, DataFrames, GLM
#implicit_intercept(::Type{<:Any}) = true

#=function mwe(N=20)
  df = DataFrame(x = rand(N), z = (i->Symbol(:s, i)).(rand(1:3, N)))
  LHS::Nothing, RHS = nothing, :(x+z)

  #f = term(0) ~ StatsModels.terms!(StatsModels.parse!(:(x+y)))
  f = @eval(@formula(x ~ $RHS))
  #f = @eval(@formula(0 ~ $RHS))
  f = apply_schema(f, schema(f,df), StatisticalModel)
  m = modelcols(f.rhs, df)
  #mf::ModelFrame = ModelFrame(f, df)
  #m = ModelMatrix(mf).m
  display(m)
end

mwe()=#

function mwe(N=10)
  df = DataFrame(x = rand(N),
    zm = Vector{Union{Missing, Int}}(rand(1:4,N)),
    z = Vector{Int}(rand(1:4,N)))


  println("W/out missings")
  f = @eval(@formula(x ~ z))
  f = apply_schema(f, schema(f,df), StatisticalModel)
  m = modelcols(f.rhs, df)
  display(m)

  println("With missings")
  f = @eval(@formula(x ~ zm))
  f = apply_schema(f, schema(f,df), StatisticalModel)
  m = modelcols(f.rhs, df)
  display(m)
end

mwe()
=#


#makie
#=@time begin
  using AbstractPlotting, StatsMakie, DataFrames, RDatasets
  using StatsMakie: linear, smooth

function testmakie()
  N = 1000
  a = rand(1:2, N) # a discrete variable
  b = rand(1:2, N) # a discrete variable
  x = randn(N) # a continuous variable
  y = @. x * a + 0.8*randn() # a continuous variable

  scatter(x, y, markersize = 0.2)
end

testmakie()
end=#

#gadfly
#=@time begin
  using DataFrames, Gadfly, Cairo, Fontconfig, Compose

function testgadfly()
  N = 1000
  a = rand(1:2, N) # a discrete variable
  b = rand(1:2, N) # a discrete variable
  x = randn(N) # a continuous variable
  y = @. x * a + 0.8*randn() # a continuous variable

  p = plot(x=x, y=y, Geom.point)
  draw(PDF("test.pdf", 9inch, 7inch),p)
end

testgadfly()
end=#

using DataFrames, VegaLite, Revise
@time begin


  function testvega()
  N = 1000
  a = rand(1:2, N) # a discrete variable
  b = rand(1:2, N) # a discrete variable
  x = randn(N) # a continuous variable
  y = @. x * a + 0.8*randn() # a continuous variable

  df = DataFrame(x=x, y=y)

  p = df |> @vlplot(:point, x=:x, y=:y)
  save("test.pdf", p)
  end

  testvega()
end
