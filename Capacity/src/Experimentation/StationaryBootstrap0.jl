using Revise, Distributions, StatsBase, Random, BenchmarkTools, Finometrics

Random.seed!(1111)

#draws a single run of an array
function sampleblockindices(N::Int, runlength::Int)::Vector{Int}
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
function bootstraparray!(dat::Vector{Tdat}, sim::Vector{Tsim}; p::Float64) where {Tdat, Tsim}
  simrunstart::Int = 1
  N = length(dat[1])

  #need to shift by 1 due to the way the distribution is parametrized
  while simrunstart≤N
    #draw the length of the run from the target distribution
    targetrunlength = samplerunlength(p)
    simrunend = min(simrunstart+targetrunlength-1,N)
    runlength = simrunend - simrunstart + 1

    #draw a set of bootstrap indices
    inds = sampleblockindices(N,runlength)

    #for each vector, copy the relevant indices
    for (v,b) ∈ zip(dat,sim)
      copyto!(b, simrunstart, @view v[inds])
    end

    simrunstart += runlength

  end

  return nothing
end

#run the aggregation functions on the bootstrap array
function bootstrappass!(dat::Vector{Tv}, aggfunc::Taggfunc, sim::Vector{Tv} = deepcopy(dat);
  p::Float64,
  ) where {Tv<:AbstractVector, Taggfunc}


  @assert (length.(dat) .=== length.(sim)) |> all

  bootstraparray!(dat, sim; p)

  return aggfunc(sim)
end

#this runs the full bootstrap
#individual passes are NOT stored
function stationarybootstrap(dat::Vector{Tv}, aggfunc::Taggfunc, B::Int; p) where {Tv, Taggfunc}

  #pre-allocate vectors used for the bootstrap
  sims = [deepcopy(dat) for t ∈ 1:Threads.nthreads()]
  results = Dict(s => Vector{MFloat64}(undef, B) for s ∈ keys(aggfunc(dat)))


  bootstrappass!(dat, aggfunc, sims[1]; p)
  #draw a sample and compute the aggregate function, sotring the reuslts in an array
  Threads.@threads for b ∈ 1:B
    t = Threads.threadid()
    res = bootstrappass!(dat, aggfunc, sims[t]; p)

    for s ∈ keys(res)
      @assert results[s][b] === missing
      results[s][b] = res[s]
    end
  end

  return (; Dict(s=>std(skipmissing(results[s])) for s ∈ keys(results))...)
end


function teststationarybootstrap(N, B; K=2)
  dat = [rand(Normal(1,1), N) |> Vector{Union{Missing, Float64}} for i ∈ 1:K]


  #a placeholder aggregation function
  aggfunc(sim) = Dict{Symbol, Float64}(
    :doublemean => (sum(sim[1]) + sum(sim[2]))/(length(sim[1]) + length(sim[2])),
    :pear=>cor(sim[1],sim[2]),
    :kend=>corkendall(sim[1] |> Vector{Float64}, sim[2] |> Vector{Float64}),
    :spear=>corspearman(sim[1] |> Vector{Float64},sim[2] |> Vector{Float64}))

  p=0.143

  sim = deepcopy(dat)
  for v ∈ sim
    v .= missing
  end

  bootstrappass!(dat, aggfunc, sim; p)
  #@code_warntype bootstrappass!(dat, aggfunc, sim; p)
  for b ∈ sim
    @assert all(b .!== missing)
  end

  results = stationarybootstrap(dat, aggfunc, B; p)
  isapprox(results.doublemean, 1/(2N)^0.5, atol=1e-2) || @warn(
    "!isapprox(results.doublemean, 1/(2N)^0.5, atol=1e-2) - could be an issue if this is unexpected")

  #create an autocorrelated series
  vauto = Vector{MFloat64}(undef, N)
  v = rand(Normal(10,10), N)
  for i ∈ 1:N
    vauto[i] = mean(v[max((i-11),1):i])
  end

  for v ∈ dat
    v .+= vauto
  end

  nwΣ = neweywestΣfunc(4)
  X::Matrix{Float64} = hcat(dat[1] .* std(dat[2]) ./ std(dat[1]), ones(N))
  lm = FMLM(X, dat[2] |> Vector{Float64})
  pearhacz = nwΣ(lm)[1,1]^0.5



  results = stationarybootstrap(dat, aggfunc, B; p)
  @info "newey HAC: $pearhacz,
    stationarybs: $(results.pear),
    delta: $(abs(pearhacz-results.pear))
    cor: $(lm.β[1])
    std: $(((N-2)/(1-lm.β[1]^2))^(-0.5))"
  print("Time per bootstrap:")
  @btime stationarybootstrap($dat, $aggfunc, $B; p=$p)

  return nothing
end


@time teststationarybootstrap(300,10000)
