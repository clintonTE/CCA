module TimeSeriesA3P6

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
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

const DEFAULT_BURN = 1_000
const NUM_SIMS = 10_000
const DEFAULT_LAG = 5


#generates the VAR data
function getVARData!(T::Int,
    y::Vector{Float64} = Vector{Float64}(T),
    x::Vector{Float64} = Vector{Float64}(T))::NTuple{2, Vector{Float64}}

  #pregenerate RVs
  ε::Vector{Float64} = rand(Normal(), T)
  η::Vector{Float64} = rand(Normal(), T)


  y[1] = 0.0
  x[1] = 0.0

  @fastmath for i::Int ∈ 2:T
    y[i] = x[i-1] + ε[i]
    x[i] = x[i-1] + η[i]
  end

  return (y, x)
end

#hold the regression data
mutable struct VARReg
  N::Int
  T::Int
  lag::Int
  βs::Matrix{Float64}
  SEs::Matrix{Float64}
  Ts::Matrix{Float64}
end

function VARReg(N::Int, T::Int; lag::Int=DEFAULT_LAG, verbose::Bool = false,
    βs::Matrix{Float64} = Matrix{Float64}(N,2),
    SEs::Matrix{Float64} = Matrix{Float64}(N,2),
    Ts::Matrix{Float64} = Matrix{Float64}(N,2))::VARReg

  #preallocate regression components
  Z::Matrix{Float64} = ones(T,2)
  Y::Vector{Float64} = Vector{Float64}(T)
  z2::Vector{Float64} = Vector{Float64}(T) #holds VAR results for Z

  ZNames::Vector{Symbol} = [:intercept, :Z]
  YName::Symbol = :Y

  #container for holding the X data, even though we don't need it
  X::Vector{Float64} = Vector{Float64}(T) #holds VAR results for X

  for i::Int ∈ 1:N
    getVARData!(T, z2, X)
    Z[:, 2] .= z2
    getVARData!(T, Y, X)
    reg::CTLM = CTLM(Z, Y, XNames=ZNames, YName=YName)
    βs[i,:] .= reg.β
    SEs[i,:] .= (diag(getNeweyWest!(reg, lag))).^0.5
    Ts[i,:] .= βs[i,:] ./ SEs[i,:]

    if verbose
      println("NW: $(SEs[i,:]), NWSlow: $((diag(getNeweyWestSlow(reg, lag))).^0.5)")
    end
  end

  return VARReg(N, T, lag, βs, SEs, Ts)
end


#convenience method for the above
VARReg(;N::Int=NUM_SIMS, T::Int=50, lag=DEFAULT_LAG,
    βs::Matrix{Float64} = Matrix{Float64}(N,2),
    SEs::Matrix{Float64} = Matrix{Float64}(N,2),
    Ts::Matrix{Float64} = Matrix{Float64}(N,2))::VARReg =
  VARReg(N, T, lag=lag, βs=βs, SEs=SEs, Ts=Ts)

function visualizeA3P6(VRs::Vector{VARReg})

  #convenience containers
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)
  #convenience method for number formatting
  n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

  push!(plotNames, "Betas")
  push!(plots, Vector{Plot}())
  for VR::VARReg ∈ VRs
    push!(plots[end],
      plot(x=VR.βs[:,2], Geom.histogram(bincount=25, density=true, position=:dodge),
        Guide.xlabel("β"),Guide.ylabel("Density"),
        Guide.title("T=$(n2s(VR.T)), " *
          "μ(β)=$(n2s(mean(VR.βs[:,2]))), σ(β)=$(n2s(std(VR.βs[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "SEs")
  push!(plots,Vector{Plot}())
  for VR::VARReg ∈ VRs
    push!(plots[end],
      plot(x=VR.SEs[:,2], Geom.histogram(bincount=25, density=true, position=:dodge),
        Guide.xlabel("SE"),Guide.ylabel("Density"),Coord.Cartesian(xmin=0.0, xmax=1.0),
        Guide.title("T=$(n2s(VR.T)), " *
          "μ(SE)=$(n2s(mean(VR.SEs[:,2]))), σ(SE)=$(n2s(std(VR.SEs[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "Ts")
  push!(plots,Vector{Plot}())
  for VR::VARReg ∈ VRs
    push!(plots[end],
      plot(x=VR.Ts[:,2], Geom.histogram(bincount=25, density=true, position=:dodge),
        Guide.xlabel("t-stat (Newey-West Lag=5)"),Guide.ylabel("Density"),
        Guide.title("T=$(n2s(VR.T)), " *
          "μ(t)=$(n2s(mean(VR.Ts[:,2]))), σ(t)=$(n2s(std(VR.Ts[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
  end

end

#script for the problem
function A3P6Script()::Void
  y::Vector{Float64}, x::Vector{Float64} = getVARData!(10_000)
  println("Autocorrelation: $(autoρ(y,1))")
  println("Autocorrelation 1D: $(autoρ(y[2:end]-y[1:(end-1)],1))")
  #VARReg(100, 120, verbose=true)
  @time VRs::Vector{VARReg} = [ VARReg(N=10_000, T=120), VARReg(N=10_000, T=600)]
  visualizeA3P6(VRs)

  return nothing
end

@time begin
  A3P6Script()
end

end
