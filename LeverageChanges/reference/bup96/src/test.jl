
#useful macro for conditionally running things in parallel
macro mpar(cond, expr)
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



function hypersort!(df::AbstractDataFrame, sorts::Vector{Symbol}, hyperidx::Vector{Int},
  threaded::Bool=false)

  #assume the other levels are sorted
  if length(sorts) > 1
    hypersort!(df, sorts[2:end], hsidx, threaded)
  end

  target::Symbol = sorts[1]
  T::Type = eltype(df[!,target])
  gdf::Vector{SubDataFrame} = collect(groupby(df, sorts[1]))
  targetv::Vector{Float64} = (sdf::SubDataFrame->sdf[1,target]).(gdf)
  Ntarget::Int = length(targetv)

  #map the values to the index array
  cumgroups::Vector{Int} = (sdf::SubDataFrame->size(sdf,1)).(gdf)
  cumsum!(cumgroups) #this is the upperbound on the index rows
  shyperidxs::Vector{SubArray{Int}} = (i->
    view(hyperidx, (i==1 ? 1 : cumgroups[i-1]):cumgroups[i])).(1:Ntarget)
  rowmapping::Dict{T, SubArray{Int}} = Dict{T, SubArray{Int}}(zip(targetv, hsidxs))

end
function hypersort(df::AbstractDataFrame, sorts::Vector{Symbol}; threaded::Bool=false)
  hyperidx::Vector{Int} = collect(1:size(df,1))
  hypersort!(df, sorts, hyperidx, threaded)
end



#------------------------------------------
#=function sortdf3!(df::DataFrame, sorts::Vector{Symbol}; threaded::Bool = false)

  firstsort::Symbol = sorts[1]
  sort!(df, firstsort)
  (length(sorts) == 1) && (return df)

  finalsorts::Vector{Symbol} = sorts[2:end]

  dfout::DataFrame = DataFrame(view(df, falses(size(df, 1)),:))
  sdfs = groupby(df, firstsort)
  dfouts::Vector{DataFrame} = (i->deepcopy(dfout)).(1:length(sdfs))
  @mpar threaded for i ∈ 1:length(sdfs)
    sdf::SubDataFrame = sdfs[i]
    dfouts[i] = sort(sdf, finalsorts)
  end

  (d->append!(dfout,d)).(dfouts)

  return dfout
end

using DataFrames, Dates, Finometrics
function testsorts(N::Int = 10_000_000)
  local df::DataFrame
  local dfs1::DataFrame
  local dfs2::DataFrame
  local dfs3::DataFrame

  df = DataFrame(ones(N,100))
  df.rid = rand(1:25000,N)
  df.dt = rand(Date(1980,1,1):Day(1):Date(2018,12,31),N)
  df.rid2 = rand(1:1000,N)

  dfs1 = deepcopy(df)
  print("\ntime sort1: ")
  @time sort!(dfs1, [:rid, :dt, :rid2])

  Base.GC.gc()

  dfs2 = deepcopy(df)
  print("time sort2: ")
  @time begin
    sort!(dfs2, :rid)
    dfs2s = groupby(dfs2, :rid)
    Threads.@threads for i ∈ 1:length(dfs2s)
      sdfs2::SubDataFrame = dfs2s[i]
      sdfs2 .= sort(sdfs2, [:dt, :rid2])
    end
  end

  dfs3 = deepcopy(df)
  print("time sort3: ")
  @time dfs3 = sortdf3!(dfs3, [:rid, :dt, :rid2], threaded=true)

  dfs4 = deepcopy(df)
  print("time sort4: ")
  @time dfs4 .= sort(dfs4, [:rid, :dt, :rid2])

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

testsorts(100_000)=#

using StatsModels, DataFrames, GLM
implicit_intercept(::Type{<:Any}) = true

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
