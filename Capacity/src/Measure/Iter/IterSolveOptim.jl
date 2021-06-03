



##############optim
mutable struct IterOptimState{TΘ⁺<:VolumePartsIter, Tmodel, T} <: AbstractIterState{T}
  Θ⁺::TΘ⁺
  λ⁺::T

  model::Tmodel

  starttime::Float64
  besttime::Float64
  maxelapsedsincebest::Float64
  maxelapsed::Float64

  iter::Int
  bestiter::Int
  maxitersincebest::Int
  maxiter::Int

  solvetype::Symbol
  lastreporttime::Float64
  minreporttime::Float64
  stop::Union{Nothing,Symbol}
end

IterOptimState(Θ::VolumePartsIter, ::Type{T};
  callback = (x)->false, #noop, false implies continuation of the optimization

  Θ⁺::VolumePartsIter = Θ |> deepcopy,
  λ⁺::T = T(Inf),

  model=error("model is required"),

  starttime::Float64=time(),
  besttime::Float64=time(),
  maxelapsedsincebest::Float64=PARAM[:optimmaxelapsedsincebest],
  maxelapsed::Float64=PARAM[:optimmaxelapsed],

  iter::Int=0,
  bestiter::Int=0,
  maxitersincebest::Int=PARAM[:optimmaxitersincebest],
  maxiter::Int=PARAM[:optimmaxiter],

  solvetype::Symbol=PARAM[:itersolvetype],
  lastreporttime::Float64=time(),
  minreporttime::Float64=PARAM[:iterminreporttime],
  stop::Union{Nothing,Symbol} = nothing,
  ) where T = IterOptimState(Θ⁺, λ⁺,
    model,
    starttime, besttime, maxelapsedsincebest,maxelapsed,
    iter, bestiter, maxitersincebest,maxiter,
    solvetype, lastreporttime, minreporttime, stop)


#reports the current status of the solver
function status(opt, Ψ::IterOptimState)
  (time()<Ψ.minreporttime+Ψ.lastreporttime) && return nothing
  Ψ.lastreporttime=time()

  elapsed::Int = round(time() - Ψ.starttime)
  elapsedsincebest::Int = round(time() - Ψ.besttime)

  itersincebest::Int = Ψ.iter - Ψ.bestiter

  println("iter=$(Ψ.iter), optimiter=$(opt.iteration), iter_since_best=$itersincebest," *
    "optimadjloss=$(opt.value), bestloss(unadj)=$(Ψ.λ⁺), " *
    "time=$elapsed, " * "timesincebest = $elapsedsincebest")

  return nothing
end

#reports status and initiates shutdown if the state has a stop flag
function callbackoptim!(opt, Ψ::IterOptimState{<:Any,<:Any,T}) where {T}

  status(opt, Ψ)

  if !(Ψ.stop===nothing)
    @info "STOP: Ψ.stop=$(Ψ.stop)"
    error("Working around optim issue #834 by throwing error")
  end

  return !(Ψ.stop===nothing) #doesn't stop the optimization due to issue #834
end

#main update function for optim
#NOTE: only called for evaluations of f- gradient call doesn't update state
function updateΨ!(x::AbstractVector{T}, loss, Θ::AbstractVolumePartsIter{TM,TV,T},
  Xv::AbstractXYIter{TM, TV, T},
  RegressionType::Val, Ψ::IterOptimState, GrowthType::Val) where {TM, TV, T}

  #first get the loss
  λ = updateA₁!(x, loss,Θ,Xv,RegressionType,GrowthType)

  #update the state and check additional stopping conditions
  #@info (λ, typeof(λ))
  updateΨ!(Θ, Ψ, λ)

  return λ
end

function updateΨ!(Θ::VolumePartsIter, Ψ::IterOptimState{<:Any,<:Any,T}, λ::Float64) where T
  Ψ.iter+=1
  t=time()

  if λ < Ψ.λ⁺ #keep track of the best solution
    copyto!(Ψ.Θ⁺.G, Θ.G)
    Ψ.λ⁺ = λ
    Ψ.besttime=t
    Ψ.bestiter=Ψ.iter
  end

  #use the solver for stopping criteria
  #run through the stopping criteria
  elapsed::Int = round(time() - Ψ.starttime)
  elapsedsincebest::Int = round(time() - Ψ.besttime)
  itersincebest::Int = Ψ.iter - Ψ.bestiter

  #set a stop flag if desired
  if  elapsed > Ψ.maxelapsed
    Ψ.stop=:maxtime
  elseif elapsedsincebest>Ψ.maxelapsedsincebest
    Ψ.stop=:maxtimesincebest
  elseif Ψ.iter ≥ Ψ.maxiter
    Ψ.stop=:maxiter
  elseif itersincebest  ≥ Ψ.maxitersincebest
    Ψ.stop=:maxitersincebest
  end

  return nothing
end


#handles processing of the line searhc parameter

function linesearchsetup(linesearchalg::Symbol = PARAM[:optimlinesearchalg],
    linesearchargs=PARAM[:optimlinesearchargs],
    linesearchinitial::Symbol=PARAM[:optimlinesearchinitial])
  alg=Dict(
    :hagerzhang=>LineSearches.HagerZhang,
    :morethuente=>LineSearches.MoreThuente,
    :backtracking=>LineSearches.BackTracking)

  initial=Dict(
    :initialprevious=>LineSearches.InitialPrevious,
    :initialstatic=>LineSearches.InitialStatic,
    :initialhagerzhang=>LineSearches.InitialHagerZhang,
    :initialquadratic=>LineSearches.InitialQuadratic,
    :initialconstantchange=>LineSearches.InitialConstantChange)

  return (alg=alg[linesearchalg](; linesearchargs...),
    initial=initial[linesearchinitial]())
end

#LBFGS solver
function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
    Xv::AbstractXYIter{TM, TV, T},
    SolveType::Val{:optimlbfgs};
    RegressionType::Val = Val(PARAM[:iterregressiontype]),
    GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TM, TV, T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  ls = linesearchsetup()
  lbfgs = Optim.LBFGS(
    m=PARAM[:lbfgsm],
    linesearch = ls.alg,
    alphaguess = ls.initial,
    )


  Ψ = IterOptimState(Θ,Float64, model=lbfgs)
  lowerbound::T = PARAM[:optimlowerbound]
  upperbound::T = PARAM[:optimupperbound]
  lower::Vector{T} = fill(lowerbound, (dims.T-1) * dims.K)
  upper::Vector{T} = fill(upperbound, (dims.T-1) * dims.K)

  obj(x::Vector{T}) = updateΨ!(x, loss, Θ, Xv, RegressionType, Ψ, GrowthType)
  x::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType) #provides an initial vector for the optimizer
  #println("x: $x")
  ∇obj! = ∇lossforzygote(Θ, Xv, GrowthType)
  runoptimize(options) = Optim.optimize(obj, ∇obj!, lower, upper, x, Fminbox(Ψ.model), options)
  #test
  obj(x)
  ∇obj!(deepcopy(x),x)

  return solveA!(loss, Θ, Xv, RegressionType, Ψ, runoptimize)
end


#conjugate gradient solver
#Note on line searches- HZ much better than backtracking
function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
    Xv::AbstractXYIter{TM, TV, T},
    SolveType::Val{:optimcg};
    RegressionType::Val = Val(PARAM[:iterregressiontype]),
    GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TM, TV, T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  ls = linesearchsetup()
  cg = Optim.ConjugateGradient(
    eta=PARAM[:cgeta],
    linesearch = ls.alg,
    alphaguess = ls.initial,
    )


  Ψ = IterOptimState(Θ,Float64, model=cg)
  lowerbound::T = PARAM[:optimlowerbound]
  upperbound::T = PARAM[:optimupperbound]
  lower::Vector{T} = fill(lowerbound, (dims.T-1) * dims.K)
  upper::Vector{T} = fill(upperbound, (dims.T-1) * dims.K)

  obj(x::Vector{T}) = updateΨ!(x, loss, Θ, Xv, RegressionType, Ψ, GrowthType)
  x::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType) #provides an initial vector for the optimizer
  #println("x: $x")
  ∇obj! = ∇lossforzygote(Θ, Xv, GrowthType)
  runoptimize(options) = Optim.optimize(obj, ∇obj!, lower, upper, x, Fminbox(Ψ.model), options)
  #test
  obj(x)
  ∇obj!(deepcopy(x),x)

  return solveA!(loss, Θ, Xv, RegressionType, Ψ, runoptimize)
end

####2 step solve
function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
  Xv::AbstractXYIter{TM, TV, T},
  SolveType::Val{:optimcg2step};
  RegressionType::Val = Val(PARAM[:iterregressiontype]),
  GrowthType::Val = Val(PARAM[:itergrowthtype]),
  conditionaliter=PARAM[:iterconditionaliter]) where {TM, TV, T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  lsconditional = linesearchsetup()
  cgconditional = Optim.ConjugateGradient(
    eta=PARAM[:cgeta],
    linesearch = lsconditional.alg,
    alphaguess = lsconditional.initial,
    )

  #this is just for the initial heuristic burnin
  Ψconditional = IterOptimState(Θ,Float64, model=cgconditional, maxiter=conditionaliter)
  lowerbound::T = PARAM[:optimlowerbound]
  upperbound::T = PARAM[:optimupperbound]
  lower::Vector{T} = fill(lowerbound, (dims.T-1) * dims.K)
  upper::Vector{T} = fill(upperbound, (dims.T-1) * dims.K)

  objconditional(x::Vector{T}) = updateΨ!(x, loss, Θ, Xv, Val{:none}(), Ψconditional, GrowthType)
  #provides an initial vector for the optimizer
  xconditional::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType)

  #println("x: $x")
  ∇objconditional! = ∇lossforzygote(Θ, Xv, GrowthType, RegressionType=Val{:none}())
  runoptimizeconditional(options) = Optim.optimize(objconditional, ∇objconditional!, lower, upper,
    xconditional, Fminbox(Ψconditional.model), options)
  #test
  objconditional(xconditional)
  ∇objconditional!(deepcopy(xconditional),xconditional)

  λconditional = solveA!(loss, Θ, Xv, Val{:none}(), Ψconditional, runoptimizeconditional)
  @info "Burnin complete with λconditional=$λconditional"

  #******now run the real pass
  ls = linesearchsetup()
  cg = Optim.ConjugateGradient(
    eta=PARAM[:cgeta],
    linesearch = ls.alg,
    alphaguess = ls.initial,
    )

  Ψ = IterOptimState(Θ,Float64, model=cg)
  obj(x::Vector{T}) = updateΨ!(x, loss, Θ, Xv, RegressionType, Ψ, GrowthType)
  x::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType, initializeto1=false)
  ∇obj! = ∇lossforzygote(Θ, Xv, GrowthType)
  runoptimize(options) = Optim.optimize(obj, ∇obj!, lower, upper, x, Fminbox(Ψ.model), options)
  λ = solveA!(loss, Θ, Xv, RegressionType, Ψ, runoptimize)

  @info "Second pass complete with λ=$λ"
  return λ
end

#experimental
function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
    Xv::AbstractXYIter{TM, TV, T},
    SolveType::Val{:optimadhoc};
    RegressionType::Val = Val(PARAM[:iterregressiontype]),
    GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  ls = linesearchsetup()
  lbfgs = Optim.LBFGS(
    m=PARAM[:lbfgsm],
    linesearch = ls.alg,
    alphaguess = ls.initial,
    )

  oaccel = Optim.OACCEL(nlprecon=lbfgs)

  Ψ = IterOptimState(Θ,Float64, model=oaccel)
  lowerbound::T = PARAM[:optimlowerbound]
  upperbound::T = PARAM[:optimupperbound]
  lower::Vector{T} = fill(lowerbound, (dims.T-1) * dims.K)
  upper::Vector{T} = fill(upperbound, (dims.T-1) * dims.K)

  obj(x::Vector{T}) = updateΨ!(x, loss, Θ, Xv, RegressionType, Ψ, GrowthType)
  x::Vector{T} = initializeΘgvector(Θ.g, T, GrowthType) #provides an initial vector for the optimizer
  ∇obj! = ∇lossforzygote(Θ, Xv, GrowthType)
  runoptimize(options) = Optim.optimize(obj, ∇obj!, lower, upper, x, Fminbox(Ψ.model), options)

  return solveA!(loss, Θ, Xv, RegressionType, Ψ, runoptimize)
end


function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
  Xv::AbstractXYIter{TM, TV, T},
  RegressionType::Val, Ψ::IterOptimState, runoptimize;
  callback=(opt)->callbackoptim!(opt,Ψ),
  GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TM, TV, T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))

  #can't put the below into the additoinal state due to the self-reference of the additional state
  optimoptions=Optim.Options(
    f_abstol = PARAM[:optimf_abstol],
    x_abstol = PARAM[:optimx_abstol],
    f_reltol = PARAM[:optimf_reltol],
    x_reltol = PARAM[:optimx_reltol],

    g_abstol = PARAM[:optimg_tol],
    g_reltol = 0.0,
    outer_g_reltol = 0.0,
    outer_g_abstol = PARAM[:optimg_tol], #doesn't seem to work?

    #NOTE: these limits should never trigger- iterations and time limited in callback
    iterations = 10^8,#PARAM[:optimiterations],
    outer_iterations = 10^8,
    time_limit = 10^8,#PARAM[:optimtime_limit],

    store_trace = PARAM[:optimstore_trace],
    extended_trace = PARAM[:optimextended_trace],
    allow_f_increases = PARAM[:optimallow_f_increases],
    allow_outer_f_increases=PARAM[:optimallow_f_increases],
    callback = callback)

  local res
  #res = runoptimize(optimoptions)
  #error("stop")

  local completed::Bool = false
  local res

  #short circuit for debugging purposes
  if PARAM[:runstacktraceoniter]
   res = runoptimize(optimoptions)
   error("No error to trace!:\nres: $res")
  end

  try
    res = runoptimize(optimoptions) #main optimziaiton entry point
    completed=true
  catch err
    errmessage = "$err"
    if length(errmessage) > 10_000 #sometimes we get really long error messages
      errmessage = errmessage[1:10_000]
    end

    if Ψ.stop === nothing
      Ψ.stop = :aborted
      @warn "Aborted optimization with error $errmessage"
    else #this means we intentionally threw the error to break the optimization
      @info "Completed optimization and terminated via workaround (error message $errmessage)"
    end

  finally

    elapsed::Int = round(time() - Ψ.starttime)
    elapsedsincebest::Int = round(time() - Ψ.besttime)
    itersincebest::Int = Ψ.iter - Ψ.bestiter
    #Θ.G .= Ψ.Θ⁺.G

    #make sure everything is up to date
    updateA₁!(Θ.g|>deepcopy, loss, Θ, Xv, RegressionType, Val{:identity}())
    updateA₁!(Ψ.Θ⁺.g|>deepcopy, loss, Ψ.Θ⁺, Xv, RegressionType, Val{:identity}())

    #compare the current candidate against the stored best candidate in Ψ
    #state candidate should be at least as good as curent candidate
    @info "updateAfromGAₜ!(Θ,Xv, 1)"
    updateAfromGAₜ!(Θ,Xv, 1) #update to the final values of A

    λ::T = loss(Θ, Xv, RegressionType)
    (λ > T(Ψ.λ⁺)) || @warn("Current estimated loss λ
      better than \"best\" λ⁺ estimate of loss. If the analysis was aborted part-way,
      this can be ignored, otherwise figure out why this happened")

    if T(Ψ.λ⁺) < λ
      #update the current candidate with the best candidate from Ψ
      updateA₁!(Ψ.Θ⁺.g|>deepcopy, loss, Θ, Xv, RegressionType, Val{:identity}())
      @assert (Ψ.Θ⁺.g ≈ Θ.g) && (Ψ.Θ⁺.G ≈ Θ.G)

      updateAfromGAₜ!(Θ,Xv, 1) #update to the final values of A
      λ = loss(Θ, Xv, RegressionType)

      #integrity checks- makes sure the best and current loss values are consistent
      updateA₁!(Ψ.Θ⁺.g|>deepcopy, loss, Ψ.Θ⁺, Xv, RegressionType, Val{:identity}())
      updateAfromGAₜ!(Ψ.Θ⁺,Xv, 1) #update to the final values of A
      @assert (Ψ.Θ⁺.g ≈ Θ.g) && (Ψ.Θ⁺.G ≈ Θ.G)
      (Ψ.Θ⁺.A ≈ Θ.A) || error("Ψ.Θ⁺.A ≈/ Θ.A
        *******Ψ.Θ⁺.A[:,1:5]: $(Ψ.Θ⁺.A[:,1:5])
        *******Θ.A[:,1:5]: $(Θ.A[:,1:5])")
      λ⁺ = loss(Ψ.Θ⁺, Xv, RegressionType)
      (λ ≈ Ψ.λ⁺) || error("λ ≠ Ψ.λ⁺!!! λ=$λ and Ψ.λ⁺=$(Ψ.λ⁺) and λ⁺=$λ⁺")
    end
    @assert  (λ ≤ T(Ψ.λ⁺)) || (λ ≈ T(Ψ.λ⁺))

    #why did we stop?
    #TODO- clean this up- definately a work in progress
    if Ψ.stop === nothing
      @assert Optim.converged(res)
      @info "STOP: Assuming convergence though converged=$(Optim.converged(res))"
      Ψ.stop=:converged
    end

    if completed

      #Ψ.Θ⁺.g .= conditiongrowth(res.minimizer, GrowthType)
      updateA₁!(res.minimizer, loss, Ψ.Θ⁺, Xv, RegressionType, GrowthType)

      λ⁺ = loss(Ψ.Θ⁺, Xv, RegressionType)
      λ = loss(Θ, Xv, RegressionType)

      if λ⁺ < λ
        updateA₁!(Ψ.Θ⁺.g|>deepcopy, loss, Θ, Xv, RegressionType, Val{:identity}())
      else
        updateA₁!(Θ.g|>deepcopy, loss, Ψ.Θ⁺, Xv, RegressionType, Val{:identity}())
      end

      @assert loss(Θ, Xv, RegressionType) ≈ loss(Ψ.Θ⁺, Xv, RegressionType)
      updateAfromGAₜ!(Θ,Xv, 1)
      #=((Θ.G |> Matrix |> vec ≈ resminimizer .|> T)
        || (Ψ.Θ⁺.G |> Matrix |> vec ≈ resminimizer .|> T)) || error(
        "res minimizer=\n$(resminimizer)\n ... while Θ.G=\n$(Θ.G)"
      )=#

      @info "printing res (Experimental):"
      println(res)
    end

    b = IOBuffer()
    write(b, "***Results info ($(Ψ.solvetype)):\n" *
      "stop_reason: $(Ψ.stop)\n" *
      "iterations: $(Ψ.iter)/$(Ψ.maxiter) " *
        "(optim counted=$(completed ? Optim.iterations(res) : "NA aborted"))\n" *
      "elapsed_at_completion: $(elapsed)/$(Ψ.maxelapsed)\n" *
      "elapsed since best: ~$(elapsedsincebest)\n" *
      "iter since best: ~$(itersincebest)\n" *
      "Best resid: $(Ψ.λ⁺)\n")

    write(b, "\n\n***Other Simulation parameters:\n" *
      "measuretype: $(PARAM[:measuretype])\n" *
      "regressionmethod: $(PARAM[:regressionmethod])\n" *
      "controlsbykt: $(PARAM[:controlsbykt])\n" *
      "controlsbyt: $(PARAM[:controlsbyt])\n" *
      "controlsglobal: $(PARAM[:controlsglobal])\n" *
      "growth type: $GrowthType\n" *
      "momentum weight spec: $(PARAM[:momentumweightspec])\n" *
      "sorted default thresholds: $(PARAM[:sorteddefaultthresholds])\n" *
      "regression method: $(RegressionType)\n" *
      "report interval: $(PARAM[:iterminreporttime])\n" *
      "Types (Matrix, Vector, T): $TM, $TV, $T\n" *
      "bounds: [$(PARAM[:optimlowerbound]), $(PARAM[:optimupperbound])]\n" *
      "linesearch (args): $(PARAM[:optimlinesearchalg]) ($(PARAM[:optimlinesearchargs]))\n" *
      "linesearch initial: $(PARAM[:optimlinesearchinitial])\n," *
      "optimoptions: $(optimoptions)\n" *
      "model config: $(Ψ.model)\n")
    write(b, "\n$(formatparameters())\n")

    iterinfo = String(take!(b))
    #println(iterinfo)
    open("$(PARAM[:outputpath])\\$(iterresultsname(Ψ.λ⁺))_info.txt", "w+") do f
      write(f, iterinfo)
    end
    return Ψ.λ⁺
  end

  @assert false #should never get here

end

  #######BBOptim
