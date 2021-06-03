


##############optimsamin
mutable struct IterBBState{TΘ⁺<:VolumePartsIter, T} <: AbstractIterState{T}
  Θ⁺::TΘ⁺
  λ⁺::T

  starttime::Float64
  lasttime::Float64
  besttime::Float64
  maxelapsedsincebest::Float64
  maxelapsed::Float64

  steps::Int
  beststeps::Int
  maxstepssincebest::Int
  maxsteps::Int

  #lastreporttime::Float64
  #minreporttime::Float64
  stop::Union{Nothing,Symbol}
end

IterBBState(Θ::TΘ, ::Type{T};
  Θ⁺::TΘ = Θ |> deepcopy,
  λ⁺::T = T(Inf),

  starttime::Float64=time(),
  lasttime::Float64=time(),
  besttime::Float64=time(),
  maxelapsedsincebest::Float64=PARAM[:bbmaxelapsedsincebest],
  maxelapsed::Float64=PARAM[:bbmaxelapsed],

  #NOTE: using steps as a limit only makes sense within a particular method
  steps::Int=0,
  beststeps::Int=0,
  maxstepssincebest::Int=PARAM[:bbmaxstepsincebest],
  maxsteps::Int=PARAM[:bbmaxsteps],

  #lastreporttime::Float64
  #minreporttime::Float64
  stop::Union{Nothing,Symbol} = nothing,
  ) where {TΘ,T} = IterBBState{TΘ,T}(Θ⁺, λ⁺,
    starttime, lasttime, besttime, maxelapsedsincebest, maxelapsed,
    steps, beststeps, maxstepssincebest, maxsteps,
    #lastreporttime, minreporttime,
    stop)

#updates the problem state given the solution for an iteration
function callbackbb!(bb, Θ::AbstractVolumeParts{<:Any, <:Any, elTΘ},
  rΨ::Ref{IterBBState{TΘ,elTΨ}}) where {elTΘ, TΘ, elTΨ}

  #get the current time
  t::Float64 = time()

  #this next code is useful for comparing methods
  if (bb.start_time > rΨ[].lasttime) &&  (rΨ[].steps>0)#reset the state if the op solver itself reset
    rΨ[] = IterBBState(Θ, elTΨ, starttime=bb.start_time)
  end

  #now update the state
  Ψ = rΨ[]
  elapsed::Float64 = t - Ψ.starttime
  Ψ.lasttime = t

  Ψ.steps::Int = BlackBoxOptim.num_steps(bb)
  stepsperssecond::Float64 = round(Ψ.steps/elapsed, digits=1)
  stepssincelast::Int = bb.num_steps_since_last_report
  Δnumbetter::Int = bb.num_better_since_last_report
  improvementpercent::Int = Δnumbetter ÷ stepssincelast*100

  Ψ.steps += stepssincelast

  λ::elTΨ = bb |> BlackBoxOptim.best_fitness |> elTΨ
  if λ < Ψ.λ⁺ #keep track of the best solution
    copyto!(Ψ.Θ⁺.G, bb |> best_candidate |> Vector{elTΘ})
    Ψ.λ⁺ = λ
    Ψ.besttime=t
    Ψ.beststeps=Ψ.steps
  end

  #use the solver for stopping criteria
  #run through the stopping criteria
  elapsedsincebest::Int = round(t - Ψ.besttime)
  stepssincebest::Int = Ψ.steps - Ψ.beststeps

  println("bestloss=$(Ψ.λ⁺)," *
    "time=$(round(elapsed))/$(Ψ.maxelapsed), " *
    "timesincebest=$elapsedsincebest/$(Ψ.maxelapsedsincebest), " *
    "steps=$(Ψ.steps), " *
    "improvements=$(Δnumbetter), " *
    "stepssincebest≈$stepssincebest, " *
    "steps/s=$stepsperssecond")
  bb.num_better_since_last_report = 0
  bb.num_steps_since_last_report = 0

  if !(Ψ.stop===nothing)
    throw("Previous shutdown call failed to cut the optimizer. Manual abort!")
  elseif elapsed > Ψ.maxelapsed
    @info "STOP: max time $(round(elapsed))/$(Ψ.maxelapsed) reached"
    Ψ.stop=:maxtime
  elseif elapsedsincebest>Ψ.maxelapsedsincebest
    @info "STOP: max time with minimal improvement $elapsedsincebest/$(Ψ.maxelapsedsincebest) reached"
    Ψ.stop=:maxtimesincebest
  elseif Ψ.steps ≥ Ψ.maxsteps
    @info "STOP: max steps $(Ψ.steps)/$(Ψ.maxsteps) reached." *
      "(NOTE- the meaning of this depends on the model)"
    Ψ.stop=:maxiter
  elseif stepssincebest  ≥ Ψ.maxstepssincebest
    @info "STOP: max steps since best $stepssincebest/$(Ψ.maxstepssincebest) reached. " *
      "(NOTE- the meaning of this depends on the model)"
    Ψ.stop=:maxstepssincebest
  end

  Ψ.stop===nothing || BlackBoxOptim.shutdown_optimizer!(bb)

  return nothing
end


function BBOptions(callback, Θ::VolumePartsIter;
    bbadditionaloptions = [:PopulationSize, :TraceMode])

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))

  #construct the required parameters
  method = PARAM[:bbmethod]
  searchrange = (PARAM[:bblowerbound], PARAM[:bbupperbound])
  numdimensions=(dims.T-1)*dims.K
  callbackinterval=PARAM[:iterminreporttime]


  #make sure we only exit via the callback, setting others to be ignored
  maxfuncevals = 0
  maxsteps = 0 #some large number

  #construct the
  bboptionsdict = Dict(
    :Method=>method,
    :SearchRange=>searchrange,
    :NumDimensions=>numdimensions,
    :MaxFuncEvals=>maxfuncevals,
    :MaxSteps=>maxsteps,
    :CallbackFunction=>callback,
    :CallbackInterval=>callbackinterval) # we willl run the callback function each loop

  #load additional optional parameters
  loadedoptionalparams::Vector{Symbol} = Vector{Symbol}()
  for k ∈ bbadditionaloptions
    #parse the option name into its PARAM name
    paramname::Symbol = k |> string |> lowercase |> (s)->"bb$s" |> Symbol
    if !(PARAM[paramname] === nothing)
      bboptionsdict[k] = PARAM[paramname]
      push!(loadedoptionalparams, k)
    end
  end
  @info "Loaded bb optional params $loadedoptionalparams"

  return (; bboptionsdict...)
end

function solveA!(loss, Θ::TΘ,
    Xv::AbstractXYIter{TM, TV, T}, _Xv::XYIterAlloc{TM, TV, T},
    SolveType::Val{:bb};
    RegressionType::Val = Val(PARAM[:iterregressiontype]),
    GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TΘ, TM, TV, T}

    dims = (K=size(Θ.A,1), T=size(Θ.A,2))

    rΨ::Ref{IterBBState{TΘ,Float64}} = Ref(IterBBState(Θ,Float64))
    callback(opt) = callbackbb!(opt,Θ,rΨ)

    bboptions = BBOptions(callback, Θ)

    x::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType) #provides an initial vector for the optimizer

    #create an objective stub

    function obj(x::Vector)
      λ=updateA₁!(x, loss, Θ, Xv, _Xv, RegressionType, GrowthType) |> Float64
      return isfinite(λ) ? λ : Inf
    end

    #test the objective
    obj(x)
    local res
    try
      res = bboptimize(obj; bboptions...)
    catch err
      @warn "Aborted optimization with error $err. "
      Ψ.stop = :aborted
    finally

      t=time()
      Ψ=rΨ[]
      elapsed::Int = round(time() - Ψ.starttime)
      elapsedsincebest::Int = round(t - Ψ.besttime)
      stepssincebest::Int = Ψ.steps - Ψ.beststeps
      Θ.G .= Ψ.Θ⁺.G
      updateA₁!(loss, Θ, Xv, _Xv, RegressionType)
      updateAfromGAₜ!(Θ,Xv, _Xv,1) #update to the final values of A

      b = IOBuffer()
      write(b, "***Results info (BlackBoxOptim):\n" *
        "stop_reason: $(Ψ.stop)\n" *
        "iterations: $(Ψ.steps)/$(Ψ.maxsteps)\n" *
        "elapsed_at_completion: $(elapsed)/$(Ψ.maxelapsed)\n" *
        "elapsed since best: ~$(elapsedsincebest)/$(Ψ.maxelapsedsincebest)\n" *
        "iter since best: ~$(stepssincebest)/$(Ψ.maxstepssincebest)\n" *
        "Best resid: $(Ψ.λ⁺)\n")

      write(b, "\n\n***Other Simulation parameters:\n" *
        "growth type: $GrowthType\n" *
        "regression method: $(RegressionType)\n" *
        "report interval: $(PARAM[:iterminreporttime])\n" *
        "Types (Matrix, Vector, T): $TM, $TV, $T\n" *
        "bboptions:\n")
      for n ∈ propertynames(bboptions)
        occursin("callback", n |> string |> lowercase) && continue #callback options covered elsewhere
        write(b, "..$n=>$(bboptions[n])\n")
      end

      iterinfo = String(take!(b))
      println(iterinfo)
      open("$(PARAM[:outputpath])\\$(iterresultsname(Ψ.λ⁺))_info.txt", "w+") do f
        write(f, iterinfo)
      end

      return Ψ.λ⁺
    end

    @assert false #should never get here

  end

  #######BBOptim
