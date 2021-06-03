module TimeSeriesA3P1

#Use this to print
#=
weave(Pkg.dir("$(pwd())\\FamaMacbeth.jl"),
  informat="script",
  out_path = "$(pwd())\\FamaMacbethOUT.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  ParallelAccelerator, Measures, Formatting

importall CTIO, CTReg, CTStat

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const YAHOO_DATE_FORMAT = "m/d/yyyy"
const FRED_DATE_FORMAT = "m/d/yyyy"
const SP500_NAME = "GSPC"
const DBV_NAME = "DBV"
const FRED_NAME = "CPI"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
const PlotContainer = Union{Plot,Gadfly.Compose.Context} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]



#preprocess inflation data
function preProcessFredDF()::Void

  fredDF::DataFrame = readtable("$DATA_PATH\\$FRED_NAME.csv")


  #read and pre-process the data
  dateFormat::String = FRED_DATE_FORMAT
  if :Date ∈ names(fredDF)
    rename!(fredDF, :Date, :DATE)
  elseif :_Date ∈ names(fredDF)
    rename!(fredDF, :_Date, :DATE)
  end

  fredDF[:,:DATE] = ((s::String)->Date(s, dateFormat)::Date).(fredDF[:,:DATE])::DataVector{Date}
  sort!(fredDF, cols = [:DATE])

  fredDF[:,:lCPI] = log.(fredDF[:,:CPI])

  stream::IOStream = open("$DATA_PATH\\$FRED_NAME.jls", "w")
  serialize(stream, fredDF)
  close(stream)

  return nothing
end



#simple container function for the autocorrelation work
function A3BuildAutoρDF(fredDF::DataFrame; lagRange::Range = 1:36)::DataFrame

  #make the dataframe
  ρDF::DataFrame = DataFrame(lag=lagRange, actualρ=autoρ(fredDF[:,:inflation], lagRange))
  ρ1::Float64 = autoρ(fredDF[:,:inflation], 1)
  ρDF[:,:AR1ρ] = (ρ1).^(1:36)

  return ρDF

end

function estimateA3P5Once(sampleMoments::NTuple{4,Float64},
    lowerBound::NTuple{4,Float64}, upperBound::NTuple{4,Float64},
    initialValues::NTuple{4,Float64}, alg::Symbol, timeLimit::Float64)


  function obj(πₐ::A, σ²ᵢ::B, ρ::C, σ²ₓ::D) where {A<:Number, B<:Number, C<:Number, D<:Number}
    oneMinusρSq = (1. - ρ^2)
    SxSi = (σ²ᵢ*σ²ₓ)^0.5
    σ² = (σ²ₓ+σ²ᵢ*oneMinusρSq) / oneMinusρSq
    ρ1 = ((ρ*σ²ₓ+SxSi*oneMinusρSq)/oneMinusρSq)/σ²
    ρ2 = ρ*ρ1

    return (πₐ - sampleMoments[1])^2 + (σ² - sampleMoments[2])^2 +
      (ρ1 - sampleMoments[3])^2 + (ρ2 - sampleMoments[4])^2
  end

  #build the optimizer
  #Algorithms: :LD_LBFGS (fails), :LN_COBYLA (bad), :LN_BOBYQA (bad)
  #  :LD_SLSQP (fails), :LD_TNEWTON (fails), :LD_VAR2/:LD_VAR1 (fails),
  #  :LD_MMA (bad), :GN_ESCH (hang), :GN_DIRECT (questionable), GN_DIRECT_L (questionable),
  # GN_ISRES (bad), :LN_NELDERMEAD(bad), :LN_SBPLX(bad), :LD_LBFGS (fails), :GD_STOGO (fails)
  #GN_CRS2_LM (OK)

  #beginning to solve model
  mod = Model(solver=NLoptSolver(algorithm=alg, maxtime=timeLimit))

  #large boundaries
  @variable(mod, lowerBound[1] ≤ πₐ ≤ upperBound[1])
  @variable(mod, lowerBound[2] ≤ σ²ᵢ ≤ upperBound[2])
  @variable(mod, lowerBound[3] ≤ ρ ≤ upperBound[3])
  @variable(mod, lowerBound[4] ≤ σ²ₓ ≤ upperBound[4])

  setvalue(πₐ, initialValues[1])
  setvalue(σ²ᵢ, initialValues[2])
  setvalue(ρ, initialValues[3])
  setvalue(σ²ₓ, initialValues[4])

  JuMP.register(mod, :obj, 4, obj, autodiff=true)
  @NLobjective(mod, Min, obj(πₐ, σ²ᵢ, ρ, σ²ₓ))

  #println("Running Solve: ")
  output = solve(mod)
  πₐAns::Float64, σ²ᵢAns::Float64, ρAns::Float64, σ²ₓAns::Float64 =
  getvalue(πₐ), getvalue(σ²ᵢ), getvalue(ρ), getvalue(σ²ₓ)

  println("Dennis Obj: $(obj(.0287445, 1.271375E-03^2,.7237045,  2.919860E-03^2))")
  println(" Dennis Inf: $(((1.271375E-03^2+2.919860E-03^2*(1-.7237045^2)) / (1-.7237045^2))^0.5)")


  return (obj(πₐAns, σ²ᵢAns, ρAns, σ²ₓAns), πₐAns, σ²ᵢAns, ρAns, σ²ₓAns)
end

function estimateA3P5(inflation::Vector{Float64}; Δx::Float64 = 0.5,
    lowerBound::NTuple{4,Float64} = (0.0, 0.0, 0.0, 0.0),
    upperBound::NTuple{4,Float64} = (0.01, 0.01, 0.9999999, 0.01),
    initialValue::NTuple{4,Float64} = (upperBound .+ lowerBound)./2,
    hybridGridSearch::Bool = true, alg::Symbol = :LN_SBPLX, timeLimit::Float64 = 2.)

  #first calculate the sample moments
  μ::Float64 = mean(inflation)
  σ²::Float64 = var(inflation)
  ρ1::Float64 = autoρ(inflation, 1)
  ρ2::Float64 = autoρ(inflation, 2)
  sampleMoments::NTuple{4,Float64} = (μ, σ², ρ1, ρ2)
  bestObj::Float64 = 10.^ 10.
  runCtr::Int = 0
  failureCtr::Int = 0

  ΔxAbs::NTuple{4,Float64} = Δx .* (upperBound .- lowerBound)

  if hybridGridSearch
    #this does the grid serach
    for πₐInit::Float64 ∈ (lowerBound[1]+ΔxAbs[1]/2):ΔxAbs[1]:upperBound[1]-ΔxAbs[1]/2,
      σ²ᵢInit::Float64 ∈ (lowerBound[2]+ΔxAbs[2]/2):ΔxAbs[2]:upperBound[2]-ΔxAbs[2]/2,
      ρInit::Float64 ∈ (lowerBound[3]+ΔxAbs[3]/2):ΔxAbs[3]:upperBound[3]-ΔxAbs[3]/2,
      σ²ₓInit::Float64 ∈ (lowerBound[4]+ΔxAbs[4]/2):ΔxAbs[4]:upperBound[4]-ΔxAbs[4]/2

      runCtr += 1

      initialValues::NTuple{4, Float64} = (πₐInit, σ²ᵢInit, ρInit, σ²ₓInit)
      candidateObj::Float64 = 10^10.

      try
        candidateObj,  πₐAns::Float64, σ²ᵢAns::Float64, ρAns::Float64, σ²ₓAns::Float64 =
          estimateA3P5Once(sampleMoments, initialValues .- ΔxAbs./2, initialValues .+ ΔxAbs./2,
            initialValues, alg, timeLimit)
      catch
        failureCtr+=1
      end

      if candidateObj < bestObj
        bestObj=candidateObj
        bestValues::NTuple{4, Float64} = (πₐAns, σ²ᵢAns, ρAns, σ²ₓAns)
      end
    end

  else
    bestObj,  πₐAns, σ²ᵢAns, ρAns, σ²ₓAns =
      estimateA3P5Once(sampleMoments, lowerBound, upperBound, (lowerBound .+ upperBound) ./ 2,
        alg, timeLimit)
    bestValues = (πₐAns, σ²ᵢAns, ρAns, σ²ₓAns)
  end

  println("Simulation Results:")
  println("πₐ = $(bestValues[1]), σᵢ = $(bestValues[2]^0.5), ρ = $(bestValues[3]), σ²ₓ=$(bestValues[4]^0.5)")
  println("Objective: $(bestObj)")
  println("Implied Sigma: $(((bestValues[4]+bestValues[2]*(1-bestValues[3]^2))/(1-bestValues[3]^2))^0.5)")
  println("Sample moments: μ=$μ, σ=$(σ²^0.5), ρ1=$ρ1, ρ2=$ρ2)")
  println("$failureCtr failures out of $runCtr")

  return nothing
end


#problem specific graphs
function visualizeA3P1(fredDF::DataFrame, ρDF::DataFrame)::Void
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()
  lagRange::Range = 1:size(ρDF,1)

  Gadfly.push_theme(:default)
  outDF::DataFrame = deepcopy(melt(fredDF[:,[:DATE,:inflationx12,:inflationYoY]], [:DATE]))
  completecases!(outDF)

  push!(plotNames, "Inflation Over Time")
  push!(plots, plot(outDF, x=:DATE, y=:value, color=:variable,  Geom.line,
    Guide.xlabel("Date"),Guide.ylabel("Inflation"),
    Guide.title("Inflation over Time")))

  push!(plotNames, "Inflation auto ρ")
  push!(plots, plot(
    layer(x=collect(lagRange), y=ρDF[:,:actualρ],  Geom.line,
      Theme(default_color=COLORS[1])),
    layer(x=collect(lagRange), y=ρDF[:,:AR1ρ],  Geom.line,
      Theme(default_color=COLORS[2])),
    Guide.manual_color_key("Legend", ["Actual", "AR1"], COLORS[1:2]),
    Guide.xlabel("Lag"),Guide.ylabel("Autocorrelation"),
    Guide.title("AutoCorrelation")))

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing
end


#problem specific code
function A3P1Script(; refreshData::Bool = true)
  fredDF::DataFrame = getDF("$DATA_PATH\\$FRED_NAME.jls", refreshData, preProcessFredDF)

  #some inflation math
  fredDF[:,:inflation] = similar(fredDF[:,:CPI])
  fredDF[2:end,:inflation] = fredDF[2:end,:lCPI] .-  fredDF[1:(end-1),:lCPI]

  fredDF[:,:inflationx12] = fredDF[:,:inflation] .* 12.
  fredDF[:,:inflationYoY] =  similar(fredDF[:,:inflation])
  fredDF[13:end,:inflationYoY] .= fredDF[13:end, :lCPI] .- fredDF[1:(end-12), :lCPI]

  #build the autocorrelation dataframe
  ρDF::DataFrame = A3BuildAutoρDF(fredDF)

  inflation::Vector{Float64} = fredDF[2:end,:inflation]

  #perform the LB test
  TLBTest::Int = round(log(length(inflation)))
  LBTest::UnivariateTest = UnivariateTest(ljungBox, Chisq(TLBTest))

  Q::Float64, p::Float64 = test(LBTest, inflation)

  println("Ljung Box Test Completed with statistic of $Q and p value of $p")

  #problem 5 solve
  #Algorithms: :LD_LBFGS (fails), :LN_COBYLA (bad), :LN_BOBYQA (bad)
  #  :LD_SLSQP (fails), :LD_TNEWTON (fails), :LD_VAR2/:LD_VAR1 (fails),
  #  :LD_MMA (bad), :GN_ESCH (hang), :GN_DIRECT (questionable), GN_DIRECT_L (questionable),
  # GN_ISRES (bad), :LN_NELDERMEAD(bad), :LN_SBPLX(bad), :LD_LBFGS (fails), :GD_STOGO (fails)
  #GN_CRS2_LM (OK)
  #GN_DIRECT .060, LN_SBPLX .03, MMA .278 (Slow, 0.5), CRS2 0.3 (Slow, 0.5), SLSQP 0.08 (Rocky)

  estimateA3P5(inflation, hybridGridSearch=false, Δx = 0.3, alg=:LD_SLSQP, timeLimit = 2.)

  visualizeA3P1(fredDF, ρDF)

end

@time begin
  A3P1Script()
end

end
