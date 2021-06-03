
module  A3

#Use this to print
#=
weave(Pkg.dir("$(pwd())\\A3.jl"),
  informat="script",
  out_path = "$(pwd())\\A3.html",
  doctype = "md2html")
=#
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, ForwardDiff

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "m/d/y"
const SP500_NAME = "GSPC"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"

const DEFAULT_NUM_SIMS = 10000
const PlotContainer = Union{Plot,Gadfly.Compose.Context}


#this holds the parameters of an AR1 Process
struct AR1Process
  α::Float64
  β::Float64
  params::Vector{Float64}
  ε::UnivariateDistribution
end

#main constructor for an AR1 with sensible defaults
AR1Process(α::Float64, β::Float64; params::Vector{Float64}=Vector{Float64}(),
  dist::Type = Normal)::AR1Process =
  AR1Process(α, β, params, Normal(params...))

#constructor to get a normal AR1 Process
AR1Process(α::Float64, β::Float64, μ::Float64, σ::Float64)::AR1Process =
          AR1Process(α, β, params=[μ,σ], dist=Normal)


#this function gets a series given an AR1 process
function getAR1Series!(p::AR1Process, T::Int,
  p₀::Float64=(p.β==1.0)?rand(p.ε)+p.α:rand(Normal(p.α,(p.params[2]^2/(1-p.β^2)))),
  series::Vector{Float64} = Vector{Float64}(T+1))::Vector{Float64}

  series[1] = p₀
  series[2:(T+1)] = rand(p.ε, T)

  for i::Int ∈ 1:T
    series[1+i] += p.α + p.β*series[i]
  end

  return series
end

function getAR1ForOLS!(p::AR1Process, T::Int;
                      p₀::Float64=rand(p.ε)+p.α,
                      AR1DF::DataFrame = DataFrame([Float64, Float64], [:Y, :X], T),
                      series::Vector{Float64} = Vector{Float64}(T+1))::DataFrame

  getAR1Series!(p, T, p₀, series)

  AR1DF[:,:X] .= series[1:(end-1)]
  AR1DF[:,:Y] .= series[2:end]

  return AR1DF
end

#this holds information from the simulation
mutable struct AR1OLSSim
  AR1::AR1Process
  α::Vector{Float64}
  β::Vector{Float64}
  σ::Vector{Float64}
  #Σ::Vector{Matrix{Float64}}
end

#this function constructs and runs the OLS sim
function AR1OLSSim(p::AR1Process, T::Int, N::Int)::AR1OLSSim

  #pre-allocate information from the regression
  α::Vector{Float64} = Vector{Float64}(N)
  β::Vector{Float64} = Vector{Float64}(N)
  σ::Vector{Float64} = Vector{Float64}(N)
  #Σ::Vector{Matrix{Float64}} = [Matrix{Float64}(2,2) for i::Int ∈ 1:N]

  #use this to store the data for a given OLS (ONLY FOR SINGLE THREADED RUNS)
  #series::Vector{Float64} = Vector{Float64}(T+1)
  #AR1DF::DataFrame = DataFrame([Float64, Float64], [:Y, :X], T)

  #Use some cores
  @fastmath Threads.@threads for i::Int ∈ 1:N
    AR1DF::DataFrame = getAR1ForOLS!(p,T)

    reg::CTLM = CTLM(AR1DF,  :X, :Y, XNames = [:intercept, :X], YName=:Y,
        eliminateNulls=false)

      #get the coefficients
      α[i], β[i] = reg.β
      σ[i] = std(reg.ε)

    #getModWhiteΣ!(reg, Σ[i])
  end

  return AR1OLSSim(deepcopy(p),α, β, σ)
end

function makeSummaryTables(sims::Vector{AR1OLSSim})
  ########################make the summary data table
  summaryDecimals::Int = 4
  numCols = length(sims)
  summaryCols::Vector{String} = ["(\$\\beta=1.0, T=50\$)",
    "(\$\\beta=1.0, T=600\$)", "(\$\\beta=0.95, T=50\$)"]
  summaryRowNames::Vector{String} =
    ["\$\\beta_{true}\$", "\$Mean(\\beta_{est})\$", "\$\\beta\$ Standard Deviation",
    "\$\\alpha_{true}\$", "\$Mean(\\alpha_{est})\$", "\$\\alpha\$ Standard Deviation",
    "\$\\sigma_{true}\$", "\$Mean(\\sigma_{est})\$", "\$\\sigma\$ Standard Deviation",
      "T-Stat (\$\\beta_{est}-\\beta_{true}\$)", "\$\\beta_{true}\$ - 5\\% T-Tail",
      "\$\\beta_{true}\$ - 1\\% T-Tail", "N"]
  summaryNumRows::Int = length(summaryRowNames)

  βTrues::Vector{Float64} = ((i::Int)->(sims[i].AR1.β)).(1:numCols)
  βEstimates::Vector{Float64} = ((sim::AR1OLSSim)->mean(sim.β)).(sims)
  βStdDevs::Vector{Float64} = ((sim::AR1OLSSim)->std(sim.β)).(sims)
  βTStats::Vector{Float64} = (βEstimates .- βTrues) ./ βStdDevs
  Ns::Vector{Int} = ((sim::AR1OLSSim)->length(sim.β)).(sims)
  βT5Tails = βTrues .+ βStdDevs .* ((n::Int)->quantile(TDist(n-2), 0.05)).(Ns)
  βT1Tails = βTrues .+ βStdDevs .* ((n::Int)->quantile(TDist(n-2), 0.01)).(Ns)

  αTrues::Vector{Float64} = ((i::Int)->(sims[i].AR1.α)).(1:numCols)
  αEstimates::Vector{Float64} = ((sim::AR1OLSSim)->mean(sim.α)).(sims)
  αStdDevs::Vector{Float64} = ((sim::AR1OLSSim)->std(sim.α)).(sims)

  σTrues::Vector{Float64} = ((i::Int)->(sims[i].AR1.params[2])).(1:numCols)
  σEstimates::Vector{Float64} = ((sim::AR1OLSSim)->mean(sim.σ)).(sims)
  σStdDevs::Vector{Float64} = ((sim::AR1OLSSim)->std(sim.σ)).(sims)

  #make the content matrix
  summaryContent::Vector{Vector{String}} = Vector{Vector{String}}()


  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βTrues)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βEstimates)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βStdDevs)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(αTrues)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(αEstimates)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(αStdDevs)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(σTrues)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(σEstimates)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(σStdDevs)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βTStats)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βT5Tails)))
  push!(summaryContent, string.(((x::Float64)->round(x,summaryDecimals)).(βT1Tails)))
  push!(summaryContent, string.(Ns))

  tableSummary::String = texTable( "Summary of Results",
    """See Tex File""", #caption
    [summaryCols], #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    summaryRowNames, #descRowNames
    summaryContent, #descContent
    Vector{String}(),
    columnSepPt = -5)
end

#problem specific script
function simulationScript(; N::Int=10000,
  T::Int = 50,
  α::Float64 = 0.0,
  β::Float64 = 1.0,
  μ::Float64 = 0.0,
  σ::Float64 = 0.2,
  staticsPoints::Int = 50,
  discountForStatics::Float64 = 0.1)::AR1OLSSim

  #main specification
  mainProblem::AR1Process = AR1Process(α, β, μ, σ)
  mainSim::AR1OLSSim = AR1OLSSim(mainProblem, T, N)
  println("""α mean: $(mean(mainSim.α))
            β mean: $(mean(mainSim.β)))""")

  #comparitive statics specifications
  σStatics::Vector{AR1Process} =  AR1Process.(α, β, μ, [i / staticsPoints
    for i::Int ∈ 1:staticsPoints]::Vector{Float64})
  αStatics::Vector{AR1Process} =  AR1Process.([i / staticsPoints - 0.5
    for i::Int ∈ 1:staticsPoints]::Vector{Float64}, β, μ, σ)
  simσStatics::Vector{AR1OLSSim} = AR1OLSSim.(σStatics, T, Int(floor(N*discountForStatics)))
  simαStatics::Vector{AR1OLSSim} = AR1OLSSim.(αStatics, T, Int(floor(N*discountForStatics)))

  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  push!(plotNames, "Statics_Beta$(β)_N$(N)_T$(T)")
  push!(plots, vstack(
    plot(x=((ar1::AR1Process)->ar1.α).(αStatics),
      y=((sim::AR1OLSSim)->mean(sim.β)).(simαStatics), Geom.point,
      Guide.ylabel("β"), Guide.xlabel("α"),
      Guide.title("Estimates of β vs. α (True β=1.0, N=$(Int(floor(discountForStatics*N))))")),
    plot(x=((ar1::AR1Process)->ar1.params[2]).(σStatics),
      y=((sim::AR1OLSSim)->mean(sim.β)).(simσStatics), Geom.point,
      Guide.ylabel("β"), Guide.xlabel("σ"),
      Guide.title("Estimates of β vs. σ (True β=1.0, N=$(Int(floor(discountForStatics*N))))"))))

      push!(plotNames, "Statics_AlphaSigma$(β)_N$(N)_T$(T)")
      push!(plots, vstack(
        plot(x=((ar1::AR1Process)->ar1.α).(αStatics),
          y=((sim::AR1OLSSim)->mean(sim.σ)).(simαStatics), Geom.point,
          Guide.ylabel("σ"), Guide.xlabel("α"),
          Guide.title("Estimates of σ vs. α (True σ=0.2, N=$(Int(floor(discountForStatics*N))))")),
        plot(x=((ar1::AR1Process)->ar1.params[2]).(σStatics),
          y=((sim::AR1OLSSim)->mean(sim.α)).(simσStatics), Geom.point,
          Guide.ylabel("α"), Guide.xlabel("σ"),
          Guide.title("Estimates of α vs. σ (True α=0.0, N=$(Int(floor(discountForStatics*N))))"))))



    for i ∈ 1:length(plots)
      draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 7.5inch, 9inch),plots[i])
    end

  return mainSim
end

function A3Script()::Void

  sims::Vector{AR1OLSSim} = Vector{AR1OLSSim}()

  push!(sims, simulationScript())
  push!(sims, simulationScript(T=600))
  push!(sims, simulationScript(β=0.95))

  writeTables2File([makeSummaryTables(sims)],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "SummaryTable.tex")

  return nothing
end

#enable to run
@time A3Script()


end
