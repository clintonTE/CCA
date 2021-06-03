using Revise
using DataFrames, CSV, Statistics, Distributions, Gadfly, Cairo, Fontconfig

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end
using CTCopy

const PLOT_NAME = "plot"
const TAX_FIELD = :interesttax
const CORP_FIELD = :corplt
const MUNI_FIELD = :muni20
const TREASURY_FIELD = :tnote20
const N_SIMS = 10_000

const LHS = [:creditspread, :taxtnote20, TAX_FIELD]
const RHS = (s->Symbol(s,"L1")).(LHS)
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}}


#gets the file
readPlot(; plotName=PLOT_NAME)::DataFrame = CSV.read("$plotName.csv")

#make the intermediate variables
function getvar!(df; taxfield = TAX_FIELD, corpfield=CORP_FIELD, munifield = MUNI_FIELD,
  treasuryfield = TREASURY_FIELD)::VAR

  df[:creditspread] = df[treasuryfield] .- df[corpfield]
  df[:munispread] = df[treasuryfield] .- df[munifield]
  df[:taxtnote20] = df[treasuryfield] .* df[taxfield]
  df[:resid1t] = df[:munispread] .- (df[:creditspread] .* (1.0 .- df[taxfield]) .+ df[:taxtnote20])
  df[:resid1] = df[:munispread] .- (df[:creditspread] .+ df[:taxtnote20])

  #labels for the VAR


  (s->lagDF!(df,s)).(LHS)

  var::VAR = (VAREstimate(df, LHS, RHS)).var
  display(var.ϕ₀)
  display(var.ϕ₁)

  return var
end

#simulate the VAR once
function simulateVAR!(var::VAR, σε::Float64, n::Int; initial=zeros(length(var.ϕ₀)), k::Int = length(var.ϕ₀),
    outmatrix = Matrix{Float64}(undef, k, n),
    simdf = DataFrame(
      creditspread = Vector{MFloat64}(missing, n),
      taxtnote20 = Vector{MFloat64}(missing, n),
      taxfield = Vector{MFloat64}(missing, n),
      residuals=Vector{MFloat64}(missing, n),
      munispread=Vector{MFloat64}(missing, n)))

  k::Int = length(var.ϕ₀)

  outmatrix::Matrix{Float64} = Matrix{Float64}(undef, k, n)
  outmatrix[:,1] .= initial

  shocks::Vector{Float64} = Vector{Float64}(undef, k)

  for t ∈ 1:(n-1)
    shocks .= rand(Normal(), k)
    outmatrix[:,t+1] = var.ϕ₀ .+ var.ϕ₁ * outmatrix[:,t] + var.L*shocks
  end

  simdf[:creditspread] .= outmatrix[1,:]
  simdf[:taxtnote20] .= outmatrix[2,:]
  simdf[:taxfield] .= outmatrix[3,:]
  simdf[:residuals] .= rand(Normal(0.0, σε), n)
  simdf[:munispread1t] .= simdf[:creditspread] .* (1.0 .- simdf[:taxfield]) .+
    simdf[:taxtnote20] .+ simdf[:residuals]
  simdf[:munispread] .= simdf[:creditspread] .+ simdf[:taxtnote20] .+ simdf[:residuals]

  return simdf
end

#runs many regressions
function runRegressions(df::DataFrame, var::VAR; residualfield=:resid1, nsims::Int = N_SIMS)
  σε::Float64 = std(skipmissing(df[residualfield]))
  n::Int = sum(completecases(df[RHS]))
  k::Int = length(var.ϕ₀)

  #do some pre-allocation
  outmatrix::Matrix{Float64} = Matrix{Float64}(undef, k, n)
  simdf = DataFrame(
    creditspread = Vector{MFloat64}(missing, n),
    taxtnote20 = Vector{MFloat64}(missing, n),
    taxfield = Vector{MFloat64}(missing, n),
    residuals=Vector{MFloat64}(missing, n),
    munispread1t=Vector{MFloat64}(missing, n),
    munispread=Vector{MFloat64}(missing, n))
  resultsdf = DataFrame(b0 = Vector{MFloat64}(missing, nsims),
    bcreditspread = Vector{MFloat64}(missing, nsims),
    btaxtnote20 = Vector{MFloat64}(missing, nsims))

  #set the regression spec
  xspec::CTExpr = Meta.parse("creditspread + taxtnote20")
  xnames::Vector{Symbol} = [:intercept, :creditspread, :taxtnote20]
  yspec::Symbol = :munispread

  #results
  for i ∈ 1:nsims
     simdf = simulateVAR!(var, σε, n, outmatrix=outmatrix, simdf=simdf)
     reg::CTLM = CTLM(simdf, xspec, yspec)
     resultsdf[i,:b0] = reg.β[1]
     resultsdf[i,:bcreditspread] = reg.β[2]
     resultsdf[i,:btaxtnote20] = reg.β[3]
   end

   return resultsdf
end

function makehistograms(resultsdf::DataFrame)

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()


  push!(plotNames, "intercept")
  push!(plots, plot(resultsdf,
    x=:b0,
    Guide.title("$(plotNames[end])"),
    Guide.xlabel("value"), Guide.ylabel("count"),
    Geom.histogram(bincount=22, density=false)))

  push!(plotNames, "betacreditspread")
  push!(plots, plot(resultsdf,
    x=:bcreditspread,
    Guide.title("$(plotNames[end])"),
    Guide.xlabel("value"), Guide.ylabel("count"),
    Geom.histogram(bincount=22, density=false)))

  push!(plotNames, "betatautnote20")
  push!(plots, plot(resultsdf,
    x=:btaxtnote20,
    Guide.title("$(plotNames[end])"),
    Guide.xlabel("value"), Guide.ylabel("count"),
    Geom.histogram(bincount=22, density=false)))

    for i ∈ 1:length(plots) #write the graphs
      draw(PDF("output\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
      #println("Graph $(plotNames[i]) written.")
    end

  show(describe(resultsdf, stats=[:mean, :std, :min, :q25, :median, :q75, :max]))

  println("\nGraphs written.")


  return nothing
end

function runPlacebo()::Nothing
  df::DataFrame = readPlot()
  var::VAR = getvar!(df)
  resultsdf = runRegressions(df, var)
  makehistograms(resultsdf)
  return nothing
end

@time runPlacebo()
