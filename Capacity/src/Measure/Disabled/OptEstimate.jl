

#=#these map model types to model functions
optmodels = Dict{Symbol, Function}(
  :levellsq => optlevellsqmodel,
  :levellsqbb => optlevellsqbbmodel,
)=#



optresultsname(λ⁺, ::Val{:levellsq}) = "results_$(PARAM[:crspdatalabel])" *
  boundrangestr() *
  "_$(PARAM[:predictiontype])" *
  "_$(PARAM[:optresultssuffix])" *
  "_$λ⁺"

optresultsname(λ⁺, ::Val{:levellsqbb}) = "results_$(PARAM[:crspdatalabel])" *
  boundrangestr() *
  "_$(PARAM[:predictiontype])" *
  "_$(PARAM[:optresultssuffix])" *
  "_$λ⁺" *
  "_$(Dates.format(now(),"yyyymmdd_HHMMSS"))"

optresultsname(::Val{:levellsqbb}) = "results_$(PARAM[:crspdatalabel])" *
  boundrangestr() *
  "_$(PARAM[:predictiontype])" *
  "_$(PARAM[:optresultssuffix])" *
  "_$(Dates.format(now(),"yyyymmdd_HHMMSS"))"
optresultsname(λ⁺) = optresultsname(λ⁺, Val(PARAM[:optmodel]))
optresultsname() = optresultsname(Val(PARAM[:optmodel]))
#primary function for estimating the measure using optimizations
#note that this still uses the same volume parts construct as flux
function formoptregression(panel::DataFrame, ms::MeasureSpec,
  OptModel::Val = Val(PARAM[:optmodel]))

  issorted(panel, [:permno, :date]) || error("panel must be sorted")


  #in contrast to other methods, drop any column where we don't have
  #the lag weight (improves performance I think)

  #conditionforturing!(panel::SubDataFrame, ms::MeasureSpec)

  zs::ZSpec = condition!(panel, ms)

  #testing
  if PARAM[:testmeasure]
    testdlsq(panel, zs.Zms)
  end

  if PARAM[:refreshmeasure]
    λ⁺, results::DataFrame = optmodel(panel, zs, OptModel)
    destandardize!(results, zs)

    #this will be the name of the results file
    resultname::String = optresultsname(λ⁺)

    results |> CSV.write("$(PARAM[:outputpath])\\$(resultname).csv")

  end

end

#mostly taken from /BlackBoxOptim.jl/src/opt_controller.jl


mutable struct BBAdditionalState
  lastresid::Float64
  lasttime::Float64
  besttime::Float64
  timesincebest::Float64
  stepssincebest::Int
  stopnoimprovement::Float64
end

function BBAdditionalState()
  lastresid=Inf
  lasttime=0
  stepssincebest=0
  besttime = 0.0
  timesincebest=0.0
  stopnoimprovement=PARAM[:optbbstopnoimprovement]

  return   BBAdditionalState(lastresid, lasttime, besttime, timesincebest, stepssincebest, stopnoimprovement)
end

function bbcallback(op, rΨ::Ref{BBAdditionalState})
  t = BlackBoxOptim.elapsed_time(op)
  if t < rΨ[].lasttime #special case when comparing multiple models
    rΨ[] = BBAdditionalState() #reinitialize
  end
  Ψ::BBAdditionalState = rΨ[]


  steps::Int = BlackBoxOptim.num_steps(op)
  stepsps = steps/t
  resid = BlackBoxOptim.best_fitness(op)
  stepssincelast = op.num_steps_since_last_report
  Δnumbetter = op.num_better_since_last_report

  if resid < Ψ.lastresid #this means we found a new best solution
    Ψ.lastresid = resid
    Ψ.besttime = t
    Ψ.timesincebest=0.0
    Ψ.stepssincebest = 0
  else
    Ψ.stepssincebest += stepssincelast
    Ψ.timesincebest = t - Ψ.besttime
  end
  Ψ.lasttime=t
  #reset the built-in state
  op.num_better_since_last_report = 0
  op.num_steps_since_last_report = 0

  form(x) = round(x, digits=1)
  @info ("t=$(form(t))s, steps=$steps, steps/s=$(form(stepsps)) ," *
    " resid=$(form(resid)), Δimprovsteps=$Δnumbetter, stepssincebest=$(Ψ.stepssincebest)" *
    " timesincebest=$(form(Ψ.timesincebest))")

  #BlackBoxOptim.shutdown_optimizer!(op)
  if Ψ.timesincebest > Ψ.stopnoimprovement
    BlackBoxOptim.shutdown_optimizer!(op)
  end

  return nothing
end

function optbbestimate(::Type{T}=PARAM[:opttype],
    ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{T} : Matrix{T};
  bounds=error("bounds (tuple or Vector of tuples) is required"),
  loss::Function = error("loss function is required"),
  Θ::AbstractVolumeParts = error("original vol parts model is required"),
  X::NamedTuple = error("X is required"),
  y::TV = error("y is required"),
  bbalg::Symbol = PARAM[:optbbalg],
  maxtime::Float64 = PARAM[:optbbmaxtime],
  reportinterval = PARAM[:optbbreportinterval],
  ) where {
    T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}


  rΨ::Ref{BBAdditionalState} = Ref(BBAdditionalState())


  opt=bboptimize(loss,
    SearchRange=bounds,
    Method=bbalg,
    MaxTime=maxtime,
    TraceMode=:silent,
    CallbackFunction=(op)->bbcallback(op, rΨ),
    CallbackInterval = reportinterval)

  return opt, rΨ[]
end

#this function writes information about the simulation run
function optbbresults(::Type{T}=PARAM[:opttype];
    loss=error("loss function is required"),
    opt=error("opt is required"),
    Ψ::BBAdditionalState=error("Ψ is required"),
    Θ::AbstractVolumeParts=error("Θ is required"),
    Θver::AbstractVolumeParts=error("Θver (AD-valid Θ) is required"),
    X::NamedTuple=error("X is required"),
    y::AbstractVector=error("y is required")) where T<:Real

  dims = (T=size(Θ.A,2), K=size(Θ.A,1))
  λ⁺ = loss(Θ, X,y)

  b = IOBuffer()
  write(b, "***Results info:\n" *
    "iterations: $(opt.iterations)\n" *
    "start_time: $(opt.start_time)\n" *
    "time_to_completion: $(opt.start_time)\n" *
    "maxtime: $(PARAM[:optbbmaxtime])\n" *
    "time since best: ~$(Ψ.timesincebest)\n" *
    "steps since best: ~$(Ψ.stepssincebest)\n" *
    "Best resid: $λ⁺\n")

  write(b, "\n\n***Other Simulation parameters:\n" *
    "method: $(opt.method)\n" *
    "stop_reason (if not timesincebest): $(opt.stop_reason)\n" *
    "report interval: $(PARAM[:optbbreportinterval])\n" *
    "eval many models: $(PARAM[:optbbcomparemodels])")
  #compute the full jacobian with AD

  ∇Θver = gradient(()->loss(Θver,X,y),Flux.params(Θver))
  rmsver = sqrt(mean(vec∇(Θver, ∇Θver).^2))

  write(b, "\n\nlst squares gradient rms: $rmsver\n")

  optinfo = String(take!(b))
  println(optinfo)
  open("$(PARAM[:outputpath])\\$(optresultsname(λ⁺))_info.txt", "w+") do f
    write(f, optinfo)
  end

  return nothing
end

function optbbcomparemodels(::Type{T}=PARAM[:opttype],
    ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{T} : Matrix{T};
  bounds=error("bounds (tuple or Vector of tuples) is required"),
  loss::Function = error("loss function is required"),
  Θ::AbstractVolumeParts = error("original vol parts model is required"),
  X::NamedTuple = error("X is required"),
  y::TV = error("y is required"),
  maxtime::Float64 = PARAM[:optbbmaxtime],
  reportinterval = PARAM[:optbbreportinterval],
  comparealgs::Vector{Symbol} = PARAM[:comparealgs])   where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}


  rΨ::Ref{BBAdditionalState} = Ref(BBAdditionalState())

  #make sure all of the methods we are going to run are valid
  for alg ∈ comparealgs
    @assert haskey(BlackBoxOptim.SingleObjectiveMethods, alg)
  end


  maxtime>1800 && error("$maxtime too high given that evalmodels=true")


  results = compare_optimizers(loss,
    SearchRange=bounds,
    MaxTime=maxtime,
    TraceMode=:silent,
    CallbackFunction=(op)->bbcallback(op, rΨ),
    CallbackInterval = reportinterval,
    Methods=comparealgs)

  #this section is based on code at  BlackBoxOptim.jl/src/compare_optimizers.jl
  sorted = sort(results, by = (t) -> t[4])
  comparedf = DataFrame(N=1:length(sorted),
    method=(v->v[1]).(sorted),
    status=(v->v[2]).(sorted),
    fitness=(v->v[4]).(sorted),
    time = (v->v[5]).(sorted))

  comparedf|>CSV.write("$(PARAM[:outputpath])\\$(optresultsname())_compare.csv")

end
