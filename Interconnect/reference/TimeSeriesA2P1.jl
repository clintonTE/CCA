module TimeSeriesA2P1

#Use this to print
#=
using Weave
codeName = "TimeSeriesA2P1"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTIO, CTReg, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  ParallelAccelerator, Measures, Formatting

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const YAHOO_DATE_FORMAT = "m/d/yyyy"
const SP500_NAME = "GSPC"
const DBV_NAME = "DBV"

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]

const DEFAULT_BURN = 1_000
const NUM_SIMS = 10_000
const DEFAULT_LAG = 5

BLAS.set_num_threads(1)


mutable struct MCSpec
  T::Int
  ρ::Float64
  burn::Int
end

#generates the VAR data
function getVARData!(MC::MCSpec,
    y::Vector{Float64} = Vector{Float64}(MC.T),
    z::Vector{Float64} = Vector{Float64}(MC.T))::NTuple{2, Vector{Float64}}

  #pregenerate RVs
  εy::Vector{Float64} = rand(Normal(), MC.T+MC.burn)
  εz::Vector{Float64} = rand(Normal(), MC.T+MC.burn)

  #first do the burn
  yburn::Float64 = 0.0
  zburn::Float64 = 0.0

  @fastmath for i::Int ∈ 1:(MC.burn+1)
    yburn = MC.ρ*yburn + εy[i]
    zburn = MC.ρ*zburn + εz[i]
  end

  y[1] = yburn
  z[1] = zburn

  @fastmath for i::Int ∈ 2:MC.T
    y[i] = MC.ρ*y[i-1] + εy[i+MC.burn]
    z[i] = MC.ρ*z[i-1] + εz[i+MC.burn]
  end

  return (y,z)
end

mutable struct VARReg
  N::Int
  MC::MCSpec
  lag::Int
  βs::Matrix{Float64}
  SEs::Matrix{Float64}
  Ts::Matrix{Float64}
end

function VARReg(N::Int, MC::MCSpec; lag::Int=DEFAULT_LAG, verbose::Bool = false,
    βs::Matrix{Float64} = Matrix{Float64}(N,2),
    SEs::Matrix{Float64} = Matrix{Float64}(N,2),
    Ts::Matrix{Float64} = Matrix{Float64}(N,2))::VARReg

  #preallocate regression components
  X::Matrix{Float64} = ones(MC.T,2)
  Y::Vector{Float64} = Vector{Float64}(MC.T)
  XNames::Vector{Symbol} = [:intercept, :X]
  YName::Symbol = :Y

  #container for holding the Z data
  Z::Vector{Float64} = Vector{Float64}(MC.T)

  for i::Int ∈ 1:N
    getVARData!(MC, Y, Z)
    X[:, 2] .= Z
    reg::CTLM = CTLM(X, Y, XNames=XNames, YName=YName)
    βs[i,:] .= reg.β
    SEs[i,:] .= (diag(getNeweyWest!(reg, lag))).^0.5
    Ts[i,:] .= βs[i,:] ./ SEs[i,:]

    if verbose
      println("NW: $(SEs[i,:]), NWSlow: $((diag(getNeweyWestSlow(reg, lag))).^0.5)")
    end

  end

  return VARReg(N, MC, lag, βs, SEs, Ts)
end

#convenience method for the above
VARReg(;N::Int=NUM_SIMS, T::Int=50, ρ::Float64=1.0, burn::Int = DEFAULT_BURN, lag=DEFAULT_LAG,
    βs::Matrix{Float64} = Matrix{Float64}(N,2),
    SEs::Matrix{Float64} = Matrix{Float64}(N,2),
    Ts::Matrix{Float64} = Matrix{Float64}(N,2))::VARReg =
  VARReg(N, MCSpec(T, ρ, burn), lag=lag, βs=βs, SEs=SEs, Ts=Ts)

function visualizeA2P1(VRs::Vector{VARReg})

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
      plot(x=VR.βs[:,2], Geom.histogram(bincount=20, density=true, position=:dodge),
        Guide.xlabel("β"),Guide.ylabel("Density"),
        Guide.title("ρ=$(n2s(VR.MC.ρ)), T=$(n2s(VR.MC.T)), " *
          "μ(β)=$(n2s(mean(VR.βs[:,2]))), σ(β)=$(n2s(std(VR.βs[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "SEs")
  push!(plots,Vector{Plot}())
  for VR::VARReg ∈ VRs
    push!(plots[end],
      plot(x=VR.SEs[:,2], Geom.histogram(bincount=20, density=true, position=:dodge),
        Guide.xlabel("SE"),Guide.ylabel("Density"),
        Guide.title("ρ=$(n2s(VR.MC.ρ)), T=$(n2s(VR.MC.T)), " *
          "μ(SE)=$(n2s(mean(VR.SEs[:,2]))), σ(SE)=$(n2s(std(VR.SEs[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  push!(plotNames, "Ts")
  push!(plots,Vector{Plot}())
  for VR::VARReg ∈ VRs
    push!(plots[end],
      plot(x=VR.Ts[:,2], Geom.histogram(bincount=20, density=true, position=:dodge),
        Guide.xlabel("t-stat (Newey-West Lag=5)"),Guide.ylabel("Density"),
        Guide.title("ρ=$(n2s(VR.MC.ρ)), T=$(n2s(VR.MC.T)), " *
          "μ(t)=$(n2s(mean(VR.Ts[:,2]))), σ(t)=$(n2s(std(VR.Ts[:,2])))")))
  end
  plots[end] = gridstack(reshape(plots[end], :, 2))

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
  end

end

#########################getNeweyWest ##########OLS Only
function getNeweyWest!(X::Matrix{Float64}, xqr::CTQR, ε::Vector{Float64}, lag::Int,
  Σ::Matrix{Float64} = Matrix{Float64}(xqr.K, xqr.K))::Matrix{Float64}

  #pre-allocate for the spectral matrix

  Rv::Matrix{Float64} = Matrix{Float64}(xqr.K, xqr.K) #pre-allocate working matrix
  RRInv::Matrix{Float64} = BLAS.gemm('N', 'T', xqr.RInv, xqr.RInv) #this is equivelent to [X'X]^-1

  #need to multiply through by the error
  Xe::Matrix{Float64} = X .* ε
  ST::Matrix{Float64} = BLAS.gemm('T','N',1.0/xqr.N, Xe, Xe)
  for v::Int ∈ 1:lag
    #overwrites Rv with (1/N)R'R
    BLAS.gemm!('T', 'N', 1.0/xqr.N, view(Xe, (v+1):(xqr.N), :),view(Xe, 1:(xqr.N-v), :), 0.0, Rv)
    #Rv .= view(Xe, (v+1):(xqr.N), :)' * view(Xe, 1:(xqr.N-v), :) .* (1.0/xqr.N)
    ST .+= (lag+1.-v)/(lag+1.) .* (Rv .+ Rv')
  end

  #this is [X'X]^-1S=[R'R]^-1S
  RRInvS::Matrix{Float64} = BLAS.gemm('N', 'N', RRInv, ST)

  #finally we have T[X'X]^-1S[X'X]^-1
  BLAS.gemm!('N','N',Float64(xqr.N), RRInvS, RRInv, 0.0, Σ)
  return Σ
end


function getNeweyWest!(lin::CTLM, lag::Int,
  Σ::Matrix{Float64} = Matrix{Float64}(lin.K, lin.K))::Matrix{Float64}

  return getNeweyWest!(lin.X, lin.xqr, lin.ε, lag, Σ)
end

#helper function in case the input requires a single argument function
getNeweyWest!(lag::Int) = (lin::CTLM)-> getNeweyWest!(lin, lag)

#test function to work
function testSalil()::Void
  Y::Vector{Float64} = SALIL_Y
  X::Matrix{Float64} = ones(length(SALIL_Z),2)
  X[:,2] .= SALIL_Z
  reg::CTLM = CTLM(X, Y)

  println("βs: $(reg.β)")
  println("σs: $((diag(getNeweyWest!(reg, 5))).^0.5)")

  return nothing
end

function A2P1Script(refreshData::Bool = false)::Void
  @time VRs::Vector{VARReg} = [VARReg(N=10_000, T=50, ρ=1.0),
      VARReg(N=10_000, T=50, ρ=0.9),
      VARReg(N=10_000, T=50, ρ=0.5),
      VARReg(N=10_000, T=1000, ρ=1.0),
      VARReg(N=10_000, T=1000, ρ=0.9),
      VARReg(N=10_000, T=1000, ρ=0.5)]

  visualizeA2P1(VRs)

  return nothing
end

#uncomment to run
@time begin
  A2P1Script()
end


end
