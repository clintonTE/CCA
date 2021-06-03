module  ChernovA1P2

#this is useful if we have mdoules
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, ForwardDiff

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "m/d/y"
const SP500_NAME = "GSPC"

const DEF_NUM_MOMENTS = 6


function preProcessSP500()::Void
  dateFormat::String = DATE_FORMAT_STR
  sp500DF::DataFrame = readtable("$DATA_PATH\\$SP500_NAME.csv")

  #fix the dates
  rename!(sp500DF, :Date, :DATE)
  sp500DF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(sp500DF[:,:DATE])::DataVector{Date}
  sort!(sp500DF, cols = [:DATE])
  sp500DF[:,:returns] = similar(sp500DF[:,:Adj_Close])
  sp500DF[2:end,:returns] = log.(sp500DF[2:end,:Adj_Close] ./  sp500DF[1:(end-1),:Adj_Close])
  #showcols(sp500DF)

  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls", "w")
  serialize(stream, sp500DF)
  close(stream)

  return nothing
end

#this function gets the sample moments of a given vector of data
function getSampleMoments!(x::Vector{Float64}, N::Int, moments::Vector{Float64} = Vector{Float64}(N))

  #for each moment we are getting
  for i::Int ∈ 1:N
    moments[i] = sum(x.^Float64(i))/Float64(length(x))
  end

  println(moments)
  println("Var: ")

  return moments
end

function getCumulantsOfJumpNorm(α::Real, σ²::Real, θ::Real, δ²::Real, ω::Real)::Vector{Real}
  κ::Vector{Real} = Vector{Real}(DEF_NUM_MOMENTS)

  #δ²::Real = δ*δ #helpers
  θ²::Real = θ*θ
  θ³::Real = θ*θ²

  κ[1] = α + ω*θ
  κ[2] = σ²+ω*(δ²+θ²) #σ^2 +ω*(δ^2+θ^2)
  κ[3] = 3.*ω*δ²*θ+ω*θ³ #3*ω*δ^2*θ+ω*θ^3
  κ[4] = ω*(6.*δ²*θ²+θ²*θ²+3.*δ²*δ²) #6*ω*δ^2*θ^2+ω*θ^4 +3*ω*δ^4
  κ[5] = ω*(15.*θ*δ²*δ²+10.*δ²*θ³+θ³*θ²) #θ*ω*(15*δ^4+10*δ^2*θ^2+θ^4)
  κ[6] = ω*(15.*δ²*δ²*δ²+45.*δ²*δ²*θ²+15.*δ²*θ²*θ²+θ³*θ³) #ω*(15*δ^6+45*δ^4*θ^2+15*δ^2*θ^4+θ^6)

  return κ
end

function get6CumulantMoments(κ...)::Vector{Real}
  μ::Vector{Real} = Vector{Real}(DEF_NUM_MOMENTS)

  κ1²::Real = κ[1] * κ[1]
  κ1³::Real = κ[1] * κ1²
  κ2²::Real = κ[2] * κ[2]

  μ[1] = κ[1]
  μ[2] = κ[2] + κ1² #κ[2]+κ[1]^2
  μ[3] = κ[3] + 3.*κ[2]*κ[1] + κ1³
  μ[4] = κ[4] + 4.*κ[3]*κ[1] + 3.*κ2² + 6.*κ[2]*κ1²+κ1²*κ1²
  μ[5] = κ[5] + 5.*κ[4]*κ[1] + 10.*κ[3]*κ[2] + 10.*κ[3]*κ1² + 15.*κ2² *κ[1] + 10.*κ[2]*κ1³+κ1³*κ1²
  μ[6] = κ[6] + 6.*κ[5]*κ[1] + 15.*κ[4]*κ[2] + 15.*κ[4]*κ1² + 10.*κ[3]*κ[3]+60*κ[3]*κ[2]*κ[1] +
    20.*κ[3]*κ1³+15.*κ2²*κ[2] + 45.*κ2²*κ1² + 15.*κ[2]*κ1²*κ1²+κ1³*κ1³

  return μ
end


function getParametersOfJumpNorm(x::Vector{Float64})::Void

  numMoments::Int = DEF_NUM_MOMENTS

  #get the sample moments
  M::Vector{Float64} = getSampleMoments!(x, numMoments)

  #this is our objective function
  function obj(α::Real, σ²::Real, θ::Real, δ²::Real, ω::Real)::Real

    κ::Vector{Real} = getCumulantsOfJumpNorm(α, σ², θ, δ², ω)
    μ::Vector{Real} = get6CumulantMoments(κ...)
    d::Real = sum((M .- μ) .* (M .- μ))

    return d
  end


  αBest::Float64=0.0; σ²Best::Float64=0.0; θBest::Float64=0.0
  δ²Best::Float64=0.0; ωBest::Float64=0.0; dBest::Float64=100.0

  solveFails::Int = 0
  solveAttempts::Int = 0
  failed::Bool = false

  #encapsulate the routing for searching a grid
  #Algorithms: :LD_LBFGS, :LN_COBYLA, :LN_BOBYQA
  #  :LD_SLSQP, :LD_TNEWTON, :LD_VAR2/:LD_VAR1,
  #  :LD_MMA, :GN_ESCH, :GN_DIRECT, GN_DIRECT_L, :GD_STOGO
  # GN_ISRES, :LN_NELDERMEAD, :LN_SBPLX,
  #GN_CRS2_LM
  function solveProblem(αMin::Float64, σ²Min::Float64, θMin::Float64, δ²Min::Float64, ωMin::Float64, Δ::Float64)
    solveAttempts += 1
    mod = Model(solver=NLoptSolver(algorithm=:LD_SLSQP, maxtime=1))

    #JuMP.build(mod)
    #maxtime!(internalmodel(mod), 1.0)

    @variable(mod, αMin <= α <= αMin + Δ)
    @variable(mod, σ²Min <= σ² <= σ²Min + Δ)
    @variable(mod, θMin <= θ <= θMin + Δ)
    @variable(mod, δ²Min <= δ² <= δ²Min + Δ)
    @variable(mod, ωMin <= ω <= ωMin + Δ)

    #register the outside functions
    JuMP.register(mod, :getCumulantsOfJumpNorm, 5, getCumulantsOfJumpNorm, autodiff=true)
    JuMP.register(mod, :get6CumulantMoments, 6, get6CumulantMoments, autodiff=true)
    JuMP.register(mod, :obj, 5, obj, autodiff=true)

    @NLobjective(mod, Min, obj(α, σ², θ, δ², ω))

    #set the starting point at the midpoint of the grid
    setvalue(α , (αMin + Δ/2.0))
    setvalue(σ², (σ²Min + Δ/2.0))
    setvalue(θ , (θMin + Δ/2.0))
    setvalue(δ², (δ²Min + Δ/2.0))
    setvalue(ω , (ωMin + Δ/2.0))

    #indicator
    failed = false
    try #solve the model
      output = solve(mod)
    catch
      solveFails+=1
      failed = true
    end

    objValue::Float64 = obj(getvalue(α), getvalue(σ²), getvalue(θ), getvalue(δ²), getvalue(ω))
    if !failed && objValue<dBest #if both statements are true, update the best parameters
      dBest = objValue
      αBest, σ²Best, θBest, δ²Best, ωBest = getvalue(α), getvalue(σ²), getvalue(θ), getvalue(δ²), getvalue(ω)
    end
  end

  #iterate the parameter space
  Δ::Float64=0.1
  @fastmath for αTry::Float64 ∈ 0.0:Δ:(0.1-Δ),
     σ²Try::Float64 ∈ 0.0:Δ:(0.2-Δ),
      θTry::Float64 ∈ -1.0:Δ:(0.0-Δ),
      δ²Try::Float64 ∈ 0.0:Δ:(0.2-Δ),
      ωTry::Float64 ∈ 0.0:Δ:(1.0-Δ)
        solveProblem(αTry, σ²Try, θTry, δ²Try, ωTry, Δ)
  end

  #= Output code here
  println("Trys: $solveAttempts")
  println("Fails: $solveFails")
  println("α=$(αBest), σ²=$(σ²Best), θ=$(θBest), δ²=$(δ²Best), ω=$(ωBest)")
  println("Also: σ=$(σ²Best^0.5) annualized: $((σ²Best)^0.5*255.0^0.5)")
  println("Objective: $dBest")
  println("For Reference: x̄=$(mean(x)), std(x)= $(std(x))")
  =#

  return nothing
end

#the main method for running the optimziation routines
function problemScript(;refreshData::Bool = true)::Void
  if !isfile("$DATA_PATH\\$SP500_NAME.jls") || refreshData
    preProcessSP500()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$SP500_NAME.jls")
  sp500DF::DataFrame = deserialize(stream)
  close(stream)

  x::Vector{Float64} = sp500DF[2:end, :returns]
  getParametersOfJumpNorm(x)

  return nothing
end

@time problemScript()
end
