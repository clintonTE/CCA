using Revise, Distributions, StatsBase, Random, BenchmarkTools, Finometrics

Random.seed!(1111)

mbsampleblockstart(N, runlength) = sample(1:(N-runlength+1))


#samplerunlength(p::Float64) = rand(Geometric(p)) + 1

#draws a non-parametric bootstrap array
function mbbootstraparray!(dat::Vector{Tdat}, sim::Vector{Tsim}; w::Int) where {Tdat, Tsim}
  simrunstart::Int = 1
  N = length(dat[1])

  #need to shift by 1 due to the way the distribution is parametrized
  while simrunstart≤N
    #draw the length of the run from the target distribution
    simrunend = min(simrunstart+w-1,N)
    runlength = simrunend - simrunstart + 1

    #draw a set of bootstrap indices
    datstart = mbsampleblockstart(N,runlength)

    #for each vector, copy the relevant indices
    for (v,b) ∈ zip(dat,sim)
      copyto!(b, simrunstart, v, datstart, runlength)
    end

    simrunstart += runlength

  end

  return nothing
end

#run the aggregation functions on the bootstrap array
function mbbootstrappass!(dat::Vector{Tv}, statfunc::Tstatfunc, sim::Vector{Tv} = deepcopy(dat);
  w::Int,
  ) where {Tv<:AbstractVector, Tstatfunc}


  @assert (length.(dat) .=== length.(sim)) |> all

  mbbootstraparray!(dat, sim; w)

  return statfunc(sim)
end

#this runs the full bootstrap
#individual passes are NOT stored
function mbbootstrap(dat::Vector{Tv}, statfunc::Tstatfunc, aggfunc::Taggfunc;
    B::Int, w::Int) where {Tv, Tstatfunc, Taggfunc}

  #pre-allocate vectors used for the bootstrap
  sims = [deepcopy(dat) for t ∈ 1:Threads.nthreads()]
  results = Dict(s => Vector{MFloat64}(undef, B) for s ∈ keys(statfunc(dat)))


  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for b ∈ 1:B
    t = Threads.threadid()
    res = mbbootstrappass!(dat, statfunc, sims[t]; w)

    for s ∈ keys(res)
      @assert results[s][b] === missing
      results[s][b] = res[s]
    end
  end

  return (; Dict(s=>aggfunc(results[s]) for s ∈ keys(results))...)
end

nonparametricp(v::AbstractVector{T}) where T<:Real = sum(v .≤ zero(T))/length(v)
nonparametricp(v::AbstractVector) = (
  sum(v .|> (x)->coalesce(x≤0.0, 0.0))/sum(v .!== missing))


function testmbbootstrap(N, B; K=2)
  dat = [rand(Normal(1,1), N) |> Vector{Union{Missing, Float64}} for i ∈ 1:K]


  #a placeholder aggregation function
  statfunc(sim) = Dict{Symbol, Float64}(
    :doublemean => (sum(sim[1]) + sum(sim[2]))/(length(sim[1]) + length(sim[2])),
    :pear=>cor(sim[1],sim[2]),
    #:kend=>corkendall(sim[1] |> Vector{Float64}, sim[2] |> Vector{Float64}),
    :spear=>corspearman(sim[1] |> Vector{Float64},sim[2] |> Vector{Float64}))

  aggfunc(v) = std(skipmissing(v))
  #aggfunc = nonparametricp


  w = 5

  sim = deepcopy(dat)
  for v ∈ sim
    v .= missing
  end

  mbbootstrappass!(dat, statfunc, sim; w)
  #@code_warntype bootstrappass!(dat, statfunc, sim; p)
  for b ∈ sim
    @assert all(b .!== missing)
  end

  results = mbbootstrap(dat, statfunc, aggfunc; B, w)
  isapprox(results.doublemean, 1/(2N)^0.5, atol=1e-2) || @warn(
    "!isapprox(results.doublemean, 1/(2N)^0.5, atol=1e-2) - could be an issue if this is unexpected")

  #create an autocorrelated series
  vauto = Vector{MFloat64}(undef, N)
  v = rand(Normal(2,2), N)
  for i ∈ 1:N
    vauto[i] = mean(v[max((i-11),1):i])
  end

  for v ∈ dat
    v .+= vauto
  end

  w=7
  nwΣ = neweywestΣfunc(5)
  X::Matrix{Float64} = hcat(dat[1] .* std(dat[2]) ./ std(dat[1]), ones(N))
  lm = FMLM(X, dat[2] |> Vector{Float64})
  pearhacz = nwΣ(lm)[1,1]^0.5



  results = mbbootstrap(dat, statfunc, aggfunc;B, w)
  @info "newey HAC: $pearhacz,
    moving se: $(results.pear),
    delta: $(abs(pearhacz-results.pear))
    cor: $(lm.β[1])
    std: $(((N-2)/(1-lm.β[1]^2))^(-0.5))"

  print("Time per bootstrap:")
  @btime mbbootstrap($dat, $statfunc, $aggfunc; B=$B, w=$w)

  return nothing
end


#@time testmbbootstrap(300,10000)

###################################
#The below is a variant that uses random blocks- shouldn't make much of a difference
#choose the moving block bootstrap since it seems more accepted in literature and is less complex
#draws a single run of an array
function sbsampleblockindices(N::Int, runlength::Int)::Vector{Int}
  datrunstart = sample(1:N)
  local out::Vector{Int}
  if datrunstart + runlength - 1 ≤ N
    out = datrunstart:(datrunstart + runlength - 1) |> collect
  else
    toplength = N - datrunstart + 1
    out = [datrunstart:N; 1:(runlength-toplength);]
  end

  #@assert length(out) == runlength "runlength=$runlength but out=$out"
  return out
end



samplerunlength(p::Float64) = rand(Geometric(p)) + 1

#draws a non-parametric bootstrap array
function sbbootstraparray!(dat::Vector{Tdat}, sim::Vector{Tsim}; p::Float64) where {Tdat, Tsim}
  simrunstart::Int = 1
  N = length(dat[1])

  #need to shift by 1 due to the way the distribution is parametrized
  while simrunstart≤N
    #draw the length of the run from the target distribution
    targetrunlength = samplerunlength(p)
    simrunend = min(simrunstart+targetrunlength-1,N)
    runlength = simrunend - simrunstart + 1

    #draw a set of bootstrap indices
    inds = sbsampleblockindices(N,runlength)

    #for each vector, copy the relevant indices
    for (v,b) ∈ zip(dat,sim)
      copyto!(b, simrunstart, @view v[inds])
    end

    simrunstart += runlength

  end

  return nothing
end


#run the aggregation functions on the bootstrap array
function sbbootstrappass!(dat::Vector{Tv}, statfunc::Tstatfunc, sim::Vector{Tv} = deepcopy(dat);
  p::Float64,
  ) where {Tv<:AbstractVector, Tstatfunc}


  @assert (length.(dat) .=== length.(sim)) |> all

  sbbootstraparray!(dat, sim; p)

  return statfunc(sim)
end

#this runs the full bootstrap
#individual passes are NOT stored
function sbbootstrap(dat::Vector{Tv}, statfunc::Tstatfunc, aggfunc::Taggfunc;
    B::Int, p::Float64) where {Tv, Tstatfunc, Taggfunc}

  #pre-allocate vectors used for the bootstrap
  sims = [deepcopy(dat) for t ∈ 1:Threads.nthreads()]
  results = Dict(s => Vector{MFloat64}(undef, B) for s ∈ keys(statfunc(dat)))


  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for b ∈ 1:B
    t = Threads.threadid()
    res = sbbootstrappass!(dat, statfunc, sims[t]; p)

    for s ∈ keys(res)
      @assert results[s][b] === missing
      results[s][b] = res[s]
    end
  end

  return (; Dict(s=>aggfunc(skipmissing(results[s])) for s ∈ keys(results))...)
end
