

#these map model types to model functions
fluxmodels = Dict{Symbol, Function}(
  :levelmle => fluxlevelmlemodel,
  :levellsq => fluxlevellsqmodel,
)



fluxresultsname(λ⁺) = "results_$(PARAM[:crspdatalabel])" *
  "_$(PARAM[:fluxopt])" *
  "_$(PARAM[:fluxmodel])" *
  "_$(PARAM[:fluxpredictiontype])" *
  "_$(PARAM[:fluxresultssuffix])" *
  "_$λ⁺"

#primary function for estimating the measure using flux
function formfluxregression(panel::DataFrame, ms::MeasureSpec;
  fluxmodel::Function = fluxmodels[PARAM[:fluxmodel]])

  issorted(panel, [:permno, :date]) || error("panel must be sorted")

  zs::ZSpec = condition!(panel, ms)

  #testing
  if PARAM[:testmeasure]
    testvolumepartsindiceslsq()
    testvolumepartsindicesmle()
  end

  if PARAM[:refreshmeasure]
    λ⁺, results::DataFrame = fluxmodel(panel, zs)
    destandardize!(results, zs)

    #this will be the name of the results file
    resultname::String = fluxresultsname(λ⁺)

    results |> CSV.write("$(PARAM[:outputpath])\\$(resultname).csv")

  end

end

#pre-run on optimimization to find optimal learning rate
function findηmax(optgen, ::Type{T}=Float32;
  loss::Function = error("loss function is required"),
  Θ::AbstractVolumeParts = error("original VolumeParts model is required"),
  X::NamedTuple= error("X is required"),
  y = error("y is required"),
  ηsearchmax::T = T(PARAM[:fluxsearchmax]),
  ηstart::T = T(PARAM[:fluxsearchstart]),
  Δη::T = T(PARAM[:fluxsearchdelta]),
  searchiter::Int=PARAM[:fluxsearchiter],
  reportinterval::NInt = nothing) where T<:Real

  Θ₀ = copytrainable(Θ)  #the best model
  Π = Flux.params(Θ) #formatted parameters

  #record η here
  ηmax::T=ηstart
  λmax::T = -T(Inf)
  η::T=ηstart
  λmin::T = T(Inf)

  #set up an iterator so we can use the built-in training loop
  repdat = Iterators.repeated((X,y), searchiter)

  #callback for reporting
  reportcb() = println("η: $η\tηmax: $ηmax\t ηsupermax: $ηsearchmax λmin: $λmin")
  throttledreportcb = reportinterval === nothing ? reportcb : throttle(reportcb, reportinterval)

  @info "Beginning search for ηmax"
  while η ≤ ηsearchmax
    opt = optgen(η)
    Flux.train!(loss, Π, repdat, opt)
      #cb=throttle(()->println("Searching...loss=", loss(ws, RLws,V)), reportinterval))
    λ::T = loss(X,y)

    if λ < λmin
      λmin = λ
      ηmax = η
    end

    if λ > λmax
      λmax = λ
    end
    throttledreportcb()

    if λ === T(Inf)
      @warn "infinite loss. Assuming ηmax ($ηmax) < η($η) and breaking."
      break
    end

    #reset and repeat
    updatexfromy!(Θ, Θ₀)
    η *= exp(Δη)
  end

  updatexfromy!(Θ, Θ₀)
  ((λmax - λmin)/(abs(λmin) + abs(λmax)) ≤ 0.01) && error("
    λmin ($λmin) nearly equal to λmax ($λmax), consider expanding search bounds")
  (λmin === T(Inf)) && error("no ηmax found!!")

  @info "ηmax search complete."
  reportcb()
  return ηmax
end

function solvevol!(optgen, ::Type{T}=Float32;
    loss::Function = error("loss function is required"),
    Θ::AbstractVolumeParts = error("original vol parts model is required"),
    X::NamedTuple = error("X is required"),
    y = error("y is required"),
    debugcb::Function = ()->nothing,
    maxiter = PARAM[:fluxmaxiter],
    maxiternoimprovement = PARAM[:fluxmaxiternoimprovement],
    reportinterval = PARAM[:fluxreportinterval],
    ηsupermindenom = T(PARAM[:fluxlearningsupermindenom]),
    Δη⁺ = T(PARAM[:fluxlearningdeltaup]),
    Δη⁻ = T(PARAM[:fluxlearningdeltadown]),
    cyclemaxdenom = T(PARAM[:fluxcyclemaxdenom]),
    cyclemindenom = T(PARAM[:fluxcyclemindenom]),
    cycledeltadenom = T(PARAM[:fluxcycledeltadenom]),
    cycles = PARAM[:fluxcycles],
    ∇rmstol = PARAM[:fluxgradrmstol],
    stopnoimprovement = PARAM[:fluxstopnoimprovement],
    stopsmallgradient = PARAM[:fluxstopsmallgradient]
  ) where T<:Real

  #get the data
  Π = Flux.params(Θ) #formatted parameters
  Θ⁺ = copytrainable(Θ)  #the best model
  Θtemp = copytrainable(Θ)  #a temporary storage location

  PARAM[:fluxopt] == :adadelta && error("udpate not supported with :adadelta")

  #initialize the current state
  local ∞::T = T(Inf)
  local λ ::T= ∞
  local λ⁺::T = ∞ #the best loss

  local itersincebest::Int = 0
  local iter::Int = 0

  #cycle parameters
  local ηmax::T = findηmax(optgen, T, loss=loss, Θ=Θ, X=X, y=y, reportinterval=reportinterval)
  ηmax = ηmax / cyclemaxdenom
  local ηmin::T = ηmax / cyclemindenom

  local ηsupermin::T = ηmin / ηsupermindenom
  (ηmax ≤ ηmin) && error("ηmax $ηmax ≤ ηmin $ηmin")
  local Δcycle::T = ηmin / cycledeltadenom

  #learning parameters
  η::T = ηmin
  local η⁺::T = η #learnign rate at best
  local ηupmult::T = exp(Δη⁺) #multiply η by this when we decrease the learning rate
  local ηdownmult::T = exp(-Δη⁻) #multiply η by this when we reduce the learning rate
  local isbestλ::Bool

  #this will yield an updated optimizer given a learning rate
  #updateopt(η) = (η, optgen(η))
  @inline function updateopt!(η, opt)
    opt.eta = η
  end
  updateopt(η) = (η, optgen(η))

  #initialize the optimizer
  opt = optgen(η)
  opt⁺ = deepcopy(opt)


  #main gradient function
  get∇Θ() = gradient(Π) do
      λ = loss(X,y) #store the loss so that we don't have to calculate it again

      if λ < λ⁺ #check if this is the best we have
        λ⁺ = λ
        η⁺ = η
        itersincebest = 0
        isbestλ = true
      else
        itersincebest += 1
        isbestλ=false
      end

      return λ
    end

  #Report progress to the user
  function reportcb(∇Θ)

    #get the gradients for monitoring
    rmsnow = ∇rms(Θ, ∇Θ)
    updatexfromy!(Θtemp, Θ) #temp=Θ
    updatexfromy!(Θ, Θ⁺) #Θ=Θ⁺
    ∇Θ = get∇Θ()
    rmsbest = ∇rms(Θ, ∇Θ)
    updatexfromy!(Θ, Θtemp) #Θ=temp

    println("loss: $λ\tbestloss: $(λ⁺)\tη:$η\tcurrent∇rms: $(rmsnow)\tbest∇rms: $(rmsbest)\t",
      iter==0 ? "" : "iter: $iter/$maxiter\titersincebest: $itersincebest/$maxiternoimprovement")

  end
  throttledreportcb = reportinterval === nothing ? reportcb : throttle(reportcb, reportinterval)

  #three exit conditions: 1) maxiter 2) no improvement for some time 3) an error
  @info "beginning cycle conditioning. ηmax=$ηmax, ηmin=$ηmin, η=$η"
  @progress for cycle ∈ 0:cycles
    @info "beginning cycle $cycle of $(cycles). Est. $((ηmax-ηmin)/Δcycle*2.0) iter/cycle"
    (cycle==cycles) && @info "YAY! ITS THE FINAL CYCLE"
    while η ≤ ηmax
      ∇Θ = get∇Θ()
      (λ==∞) && debugcb() #infinity means something went wrong
      Flux.update!(opt, Π, ∇Θ)
      if isbestλ
        updatexfromy!(Θ⁺, Θ)
        opt₊ = deepcopy(opt)
      end
      η, opt = updateopt(η + Δcycle)
      #η = updateopt!(η + Δcycle, opt)
      throttledreportcb(∇Θ) #keep user appraised
    end
    while η ≥ ηmin
      ∇Θ = get∇Θ()
      (λ==∞) && debugcb() #infinity means something went wrong
      Flux.update!(opt, Π, ∇Θ)
      if isbestλ
        updatexfromy!(Θ⁺, Θ)
        opt₊ = deepcopy(opt)
      end
      η, opt = updateopt(η - Δcycle)
      #η = updateopt!(η - Δcycle, opt)
      throttledreportcb(∇Θ) #keep user appraised
    end
  end

  @info "Final cycle complete. Beginning final convex optimization from best point"
  updatexfromy!(Θ, Θ⁺)
  itersincebest = 0
  isbestλ = false

  opt = opt⁺
  @progress for i ∈ 1:maxiter
    iter = i #for reporting purposes
    ∇Θ = get∇Θ()

    isbestλ && updatexfromy!(Θ⁺, Θ)
    if λ==∞ #reset and lower learning rate
      debugcb()  #this is an error scenario where learning rate can't decrease
    elseif isbestλ && (η ≤ ηmax) #if its the best answer, note model and accelerate learning
      #η, opt = updateopt(η * ηupmult)
      η = updateopt!(η * ηupmult, opt)
    elseif (!isbestλ) && η ≥ ηsupermin #if its not, decellerate learning
      #η, opt = updateopt(η * ηdownmult)
      η = updateopt!(η * ηdownmult, opt)
    end

    #run through the stopping criteria
    if stopnoimprovement && (itersincebest ≥ maxiternoimprovement) #this is the main stopping criteria
      @info "No improvement in $maxiternoimprovement iterations."
      reportcb(∇Θ) #run the callback once more
      break
    elseif stopsmallgradient && ∇rms(Θ, ∇Θ) < ∇rmstol
      @info "GOOD NEWS: ∇rms ($(∇rms(Θ, ∇Θ))) < tol ($∇rmstol)."
      reportcb(∇Θ) #run the callback once more
      break
    elseif iter == maxiter
      @warn "max iter reached, returning best result"
      reportcb(∇Θ)
    else
      throttledreportcb(∇Θ) #keep user appraised
      Flux.update!(opt, Π, ∇Θ)
    end
  end

  ∇Θ = get∇Θ()
  printτ = :τ ∈ propertynames(Θ) #check if τ (risk or precision) is used in the optimization process
  open("$(PARAM[:outputpath])\\$(fluxresultsname(λ⁺))_info.txt", "w+") do f
    println(f, "******Loss: loss: $λ\tbestloss: $(λ⁺)\tη:$η\tgradrms: $(∇rms(Θ, ∇Θ))",
    "\n******Current A:\n", ∇Θ[Θ.A]',
    "\n******Current V₀:\n", ∇Θ[Θ.V₀],
    (printτ ? "\n******Current τ:\n$(display(∇Θ[Θ.τ]'))" : ""),
    "\n++++++Switch to best: ", isnothing(updatexfromy!(Θ,Θ⁺)),
    "\n++++++Best A:\n", ∇Θ[Θ.A]',
    "\n++++++Best V₀:\n", ∇Θ[Θ.V₀],
    (printτ ? "\n++++++Best τ:\n$(display(∇Θ[Θ.τ]'))" : ""),
    )
  end
  return λ⁺
end
