

abstract type AbstractIterState{T<:Real} end
####Fixed point methods
#contains various state variables for running the solver
mutable struct IterFPState{TΘ⁺<:VolumePartsIter,T} <: AbstractIterState{T}
  Θ⁺::TΘ⁺
  λ⁺::T
  ftol::T

  griddelta::T
  gridstart::T
  gridend::T
  remainingpregriddelta::Vector{T}
  currentgriddelta::T

  starttime::Float64
  besttime::Float64
  maxelapsedsincebest::Float64
  maxelapsed::Float64

  iter::Int
  bestiter::Int
  maxitersincebest::Int
  maxiter::Int

  lastreporttime::Float64
  stop::Union{Nothing,Symbol}
end

#updates the problem state given the solution for an iteration
function updatefpstate!(Ψ::IterFPState{TΘ⁺,T}, Θ::VolumePartsIter, λ::T) where {TΘ⁺,T}
  Ψ.iter+=1
  t=time()
  if abs(Ψ.λ⁺ - λ) > Ψ.ftol
    Ψ.besttime=t
    Ψ.bestiter=Ψ.iter
  end

  if λ < Ψ.λ⁺ #keep track of the best solution, even if it doesn't meet ftol criteria
    Ψ.Θ⁺.Π .= Θ.Π
    Ψ.λ⁺ = λ
  end

  #adjusts the grid delta if a pre-conditioning glide path was provided
  if length(Ψ.remainingpregriddelta) > 0
    Ψ.currentgriddelta = popfirst!(Ψ.remainingpregriddelta)
    if length(Ψ.remainingpregriddelta) == 0
      Ψ.currentgriddelta = Ψ.griddelta #the last point should be the current grid
      @info "Grid delta pre-conditioning glide path complete"
    end
    Ψ.besttime=t
    Ψ.bestiter=Ψ.iter
    return nothing #this short-circuits the stopping criteria
  end

  #run through the stopping criteria
  elapsed::Int = round(time() - Ψ.starttime)
  elapsedsincebest::Int = round(time() - Ψ.besttime)
  itersincebest::Int = Ψ.iter - Ψ.bestiter
  if  elapsed > Ψ.maxelapsed
    @info "STOP: max time $elapsed/$(Ψ.maxelapsed) reached"
    Ψ.stop=:maxtime
  elseif elapsedsincebest>Ψ.maxelapsedsincebest
    @info "STOP: max time with minimal improvement $elapsedsincebest/$(Ψ.maxelapsedsincebest) reached"
    Ψ.stop=:maxtimesincebest
  elseif Ψ.iter ≥ Ψ.maxiter
    @info "STOP: max iter $(Ψ.iter) reached."
    Ψ.stop=:maxiter
  elseif itersincebest  ≥ Ψ.maxitersincebest
    @info "STOP: max iter since best $itersincebest reached."
    Ψ.stop=:maxitersincebest
  end

  return nothing
end

#contains some defaults for creating teh initial solver state
IterFPState(Θ::VolumePartsIter, ::Type{T};
  Θ⁺::VolumePartsIter = deepcopy(Θ),
  λ⁺::T=T(Inf),
  ftol::T=max(T(PARAM[:fpftol]), zero(T)),

  griddelta::T = PARAM[:fpgriddelta] |> T,
  gridstart::T = PARAM[:fpgridstart] |> T,
  gridend::T = PARAM[:fpgridend] |> T,
  remainingpregriddelta::Vector{T}= T.(PARAM[:fppregriddelta]),
  currentgriddelta::T= (
    length(remainingpregriddelta) > 0 ? popfirst!(remainingpregriddelta) : griddelta) |> T,

  starttime::Float64=time(),
  besttime::Float64=time(),
  maxelapsedsincebest::Float64=PARAM[:fpmaxelapsedsincebest],
  maxelapsed::Float64=PARAM[:fpmaxelapsed],

  iter::Int=0,
  bestiter::Int=0,
  maxitersincebest::Int=PARAM[:fpmaxitersincebest],
  maxiter::Int=PARAM[:fpmaxiter],

  lastreporttime::Float64=time(),
  stop::Union{Nothing,Symbol} = nothing,
  ) where T = IterFPState(Θ⁺, λ⁺, ftol,
    griddelta, gridstart, gridend, remainingpregriddelta, currentgriddelta,
    starttime, besttime, maxelapsedsincebest, maxelapsed,
    iter, bestiter, maxitersincebest, maxiter,
    lastreporttime, stop)

#reports the current status of the solver
function status(Ψ::IterFPState, λ; minreporttime::Float64 = PARAM[:iterminreporttime])
  (time()<minreporttime+Ψ.lastreporttime) && return nothing
  Ψ.lastreporttime=time()

  elapsed::Int = round(time() - Ψ.starttime)
  elapsedsincebest::Int = round(time() - Ψ.besttime)

  itersincebest::Int = Ψ.iter - Ψ.bestiter

  println("iter=$(Ψ.iter)/$(Ψ.maxiter), itersincebest=$itersincebest/$(Ψ.maxitersincebest)," *
    "loss=$λ, bestloss=$(Ψ.λ⁺), " *
    "time=$elapsed/$(Ψ.maxelapsed), " *
    "timesincebest = $elapsedsincebest/$(round(Ψ.maxelapsedsincebest)), " *
    "griddelta=$(Ψ.currentgriddelta)")

  return nothing
end



function solveA!(iterloss, Θ::VolumePartsIter,
    Xv::AbstractXYIter{TM,TV,T}, _Xv::XYIterAlloc{TM,TV,T}, RegressionType::Val, ::Val{:fixedpoint}
    ) where {TM,TV,T}


  Ψ::IterFPState = IterFPState(Θ, T)

  try
    while true #main solving loop

      #update the growths and values
      λ::T = updateA!(iterloss, Θ,Xv,_Xv,RegressionType, ΔG=Ψ.currentgriddelta)

      #update the problem state
      updatefpstate!(Ψ,Θ,λ)

      #check the stopping criteria and report status
      status(Ψ, λ)
      (!(Ψ.stop === nothing)) && break

    end
  catch err
    @warn "Aborted solver with error $err"
    Ψ.stop = :aborted
  finally #do this as a finallizer in case there is an error

    #ensure that we capture the best λ
    Θ.Π .= Ψ.Θ⁺.Π

    b = IOBuffer()
    elapsed::Int = round(time() - Ψ.starttime)
    elapsedsincebest::Int = round(time() - Ψ.besttime)
    itersincebest::Int = Ψ.iter - Ψ.bestiter
    write(b, "***Results info:\n" *
      "stop_reason: $(Ψ.stop)\n" *
      "iterations: $(Ψ.iter)/$(Ψ.maxiter)\n" *
      #"start_time: $(Ψ.starttime)\n" *
      "elapsed_at_completion: $(elapsed)/$(Ψ.maxelapsed)\n" *
      "elapsed since best: ~$(elapsedsincebest)/$(Ψ.maxelapsedsincebest)\n" *
      "iter since best: ~$(itersincebest)/$(Ψ.maxitersincebest)\n" *
      "Best resid: $(Ψ.λ⁺)\n")

    write(b, "\n\n***Other Simulation parameters:\n" *
      "regression method: $(RegressionType)\n" *
      "report interval: $(PARAM[:iterminreporttime])\n" *
      "Types (Matrix, Vector, T): $TM, $TV, $T\n" *
      "grid (excludes preconditioning) start:delta:end: $(Ψ.gridstart):$(Ψ.griddelta):$(Ψ.gridend)\n" *
      "currentgriddelta: $(Ψ.currentgriddelta)\n")

    if length(PARAM[:iterpregriddelta]) > 0
      write(b,
      "precondition max: $(maximum(PARAM[:iterpregriddelta]))\n" *
      "precondition length: $(length(PARAM[:iterpregriddelta]))\n"
      )
    end
    #compute the full jacobian with AD

    #∇Θver = gradient(()->loss(Θver,X,y),Flux.params(Θver))
    #rmsver = sqrt(mean(vec∇(Θver, ∇Θver).^2))

    #write(b, "\n\nlst squares gradient rms: $rmsver\n")

    iterinfo = String(take!(b))
    println(iterinfo)
    open("$(PARAM[:outputpath])\\$(iterresultsname(Ψ.λ⁺))_info.txt", "w+") do f
      write(f, iterinfo)
    end

    return Ψ.λ⁺
  end

  @assert false #should never hit this point due to the finalizer
end
