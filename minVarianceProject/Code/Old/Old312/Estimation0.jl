abstract type AbstractEstimation end

#contains the estimation problem
mutable struct Estimation <: AbstractEstimation
  Θ::Parameters
  Π::Priors
  F::Posteriors

  G::Portfolio
  TG::TrisectedPortfolio
  P::Portfolio

  nburn::Int
  nsamples::Int
  ncomplete::Int

  outparameters::Matrix{Float64}
  outweights::Matrix{Float64}
end

#creates the initial estimation object
function Estimation(U::Universe; Θ::Parameters = Parameters(),
  Π::Priors = Priors(),
  F::Posteriors = Posteriors(),
  nburn::Int = N_BURN,
  nsamples::Int = N_SAMPLES)::Estimation

  G::Portfolio = Portfolio(U)
  TG::TrisectedPortfolio = TrisectedPortfolio(G)
  P::Portfolio = initializetestportfolio(G)

  ncomplete::Int = 0


  outparameters::Matrix{Float64} = Matrix{Float64}(undef, length(copy2vector(Θ)), nsamples)
  outweights::Matrix{Float64} = Matrix{Float64}(undef, U.M, nsamples)

  return Estimation(Θ, Π, F, G, TG, P, nburn, nsamples, ncomplete, outparameters, outweights)
end

function Estimation(;datapath = DATA_PATH,
  dataname = DATA_NAME,
  datebegin::Date = DATE_BEGIN,
  dateend::Date = DATE_END,
  Θ::Parameters = Parameters(),
  Π::Priors = Priors(),
  F::Posteriors = Posteriors(),
  nburn::Int = N_BURN,
  nsamples::Int = N_SAMPLES)::Estimation

  return Estimation(Universe(datapath=datapath, dataname=dataname,
    datebegin=datebegin, dateend=dateend),
  Θ=Θ, Π=Π, F=F, nburn=nburn, nsamples=nsamples)
end

function mergeestimations!(E::Estimation...)::Estimation
  if length(E) == 1
    return E[1]
  elseif length(E) == 2
    E[1].nsamples += E[2].nsamples
    E[1].outparameters = [E[1].outparameters E[2].outparameters]
    E[1].outweights = [E[1].outweights E[2].outweights]
    return E[1]
  end

  #pivot and recurse
  pivot::Int = length(E) ÷ 2
  Eout::Estimation = mergeestimations!(
    mergeestimations!(E[1:pivot]...),
    mergeestimations!(E[(pivot+1):end]...))

  #print("length(E): $(length(E)) result: $(sum(Eout.outparameters))")
  return Eout
end


#runs a single estimation pass
function estimationpass!(E::Estimation)::Nothing

  redrawtestportfolio!(E.P)
  #E.Θ.SG = var(E.G) #mitigates error propogation
  E.Θ.SGP = cov(E.G, E.P) #where should I do this? try it here

  rotatetrisected!(E.TG)

  E.F.σ²G!(E.Θ, E.Π) #draw from the posteriors
  E.F.ζ²G!(E.Θ, E.Π)
  E.F.ζ²P!(E.Θ, E.Π)
  E.F.SG!(E.Θ, E.Π)
  E.F.SGP!(E.Θ, E.Π)

  #compute the implied weights
  ω::Vector{Float64} = getω(E.Θ, E.TG, E.P)
  rescaletrisected!(E.TG, ω)

  #record the output if we are past the burn-in
  E.ncomplete += 1
  if E.ncomplete > E.nburn
    sampleidx::Int = E.ncomplete - E.nburn

    E.outparameters[:, sampleidx] .= copy2vector(E.Θ)
    E.outweights[:, sampleidx] .= E.G.wₜ
  end

  return nothing
end

function recordestimation(E::Estimation,
  outname = OUT_NAME, outpath::String = OUT_PATH, record::Bool = true,
  timestamp::Bool = false)::Nothing

  if timestamp
    parameterspath::String = "$outpath\\$(outname)_par_$(now()).csv"
    weightspath::String = "$outpath\\$(outname)_w_$(now()).csv"
  else
    parameterspath = "$outpath\\$(outname)_par.csv"
    weightspath = "$outpath\\$(outname)_w.csv"
  end

  dfΘ::DataFrame = DataFrame(pass = collect(1:E.nsamples),
    sigma2G = E.outparameters[1,:], zeta2G = E.outparameters[2,:], zeta2P = E.outparameters[3,:],
    SG = E.outparameters[4,:], SGP = E.outparameters[5,:])
  dfw::DataFrame = DataFrame(E.outweights')

  ((oldn::Symbol, newn::Symbol)->rename!(dfw, oldn => newn)).(names(dfw), setdiff(names(E.G.U.df), [:date]))

  dfw = [DataFrame(pass = collect(1:E.nsamples)) dfw]

  CSV.write(parameterspath, dfΘ, missingstring="NA")
  CSV.write(weightspath, dfw, missingstring="NA")

  return nothing

end

#launches and performes the estimation
function estimate!(E::Estimation; verbose::Bool = true)::Estimation

  npasses::Int = E.nburn + E.nsamples

  verbose && println("Begining estimation....")
  for i::Int ∈ 1:npasses
    estimationpass!(E)
    verbose && (i % 10_000 == 0) && println("Completed $i passes out of $npasses")
  end
  verbose && println("Estimation complete and successful!")

  return E
end

#helper function which starts from the investment universe
function estimate(U::Universe; verbose::Bool = true, nburn::Int = N_BURN,
  nsamples::Int = N_SAMPLES)::Estimation

  E::Estimation = Estimation(U, nburn=nburn, nsamples=nsamples)
  estimate!(E, verbose=verbose)

  return E
end

function estimateminvariance(;datapath = DATA_PATH,
  dataname = DATA_NAME,
  datebegin::Date = DATE_BEGIN,
  dateend::Date = DATE_END,
  nburn::Int = N_BURN,
  nsamples::Int = N_SAMPLES,
  parallel::Bool = (nsamples > length(workers())*nburn), #simple heuristic for deciding if we will parllelize
  recordestimate::Bool=true)::Nothing

  U::Universe = Universe(datapath=datapath, dataname=dataname,
    datebegin=datebegin, dateend=dateend)
  local E::Estimation

  if !parallel
    E = estimate(U, nburn=nburn, nsamples=nsamples)
    recordestimate && recordestimation(E::Estimation)
    return E
  end

  #parllelized MCMC
  local pids = workers()
  local np::Int = length(pids)
  local targetassignmentsize::Int = nsamples ÷ np
  local assignments::Vector{Int} = Vector{Int}(undef, np)

  assignments[1:(end-1)] .= targetassignmentsize #divy up the assignments
  assignments[end] = nsamples - sum(assignments[1:(end-1)])

  nestimations::Int = length(assignments)
  futures::Vector{Future} = Vector{Future}()
  sizehint!(futures, np)

  #deploy the workers
  for i ∈ 1:np
    push!(futures,
      (@spawnat pids[i] estimate(U, nburn=nburn, nsamples=assignments[i])))
  end

  estimations::Vector{Estimation} = (fetch).(futures)
  estimationstest::Vector{Estimation} = deepcopy(estimations)

  E = mergeestimations!(estimations...)
  recordestimate && recordestimation(E)



  return nothing

end
