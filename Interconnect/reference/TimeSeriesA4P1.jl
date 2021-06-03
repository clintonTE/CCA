module TimeSeriesA4P1

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

using DataFrames, Distributions, StatsBase, GZip, JLD, ForwardDiff,
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

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)

struct AR2
  ϕ₀::Float64
  ϕ₁::Float64
  ϕ₂::Float64
  dist::T where T<:UnivariateDistribution
end

AR2(ϕ₁::Float64, ϕ₂::Float64; ϕ₀::Float64 = 0.0,
    dist::T where T<:UnivariateDistribution = Normal())::AR2 = AR2(ϕ₀, ϕ₁, ϕ₂, dist)

#recursive formulation of AR2 autocorrelation
#the shared memory structure makes this O(n) as opposed to O(2ⁿ)
function AR2ρ!(Θ::AR2, k::Int,
    ρₖ::Vector{Float64}=Vector{Float64}(k), solved::Vector{Bool}=fill(false,k))::Float64

  if k == 0
    ans::Float64 = 1.
  elseif solved[k]
    ans = ρₖ[k]
  elseif k == 1
    ans = Θ.ϕ₁/(1.-Θ.ϕ₂)
    ρₖ[1] = ans
    solved[k] = true
  else
    ans = Θ.ϕ₁*AR2ρ!(Θ, k-1, ρₖ, solved) + Θ.ϕ₂*AR2ρ!(Θ, k-2, ρₖ, solved)
    ρₖ[k] = ans
    solved[k] = true
  end

  return ans
end

#helper which creates the shocks
AR2Returns(Θ::AR2, k::Int; rM1::Float64 = 0.0, rM2::Float64 = 0.0)::Vector{Float64} =
  AR2Returns!(Θ, rand(Θ.dist, k), [Vector{Float64}(k-2)])

#fills the AR2
function AR2Returns!(Θ::AR2, ε::Vector{T}, rM1::Float64=0.0, rM2::Float64 = 0.0,
    r::Vector{T} = zeros(T, length(ε)))::Vector{T} where T<:Real
  k::Int = length(r)

  #fill the special cases
  r[1] = Θ.ϕ₀ + Θ.ϕ₁*rM1 + Θ.ϕ₂*rM2 + ε[1]
  if k > 1
    r[2] = Θ.ϕ₀ + Θ.ϕ₁*r[1] + Θ.ϕ₂*rM1 + ε[2]
  end

  #fill the rest
  for i::Int ∈ 3:k
    r[i] = Θ.ϕ₀ + Θ.ϕ₁*r[i-1] + Θ.ϕ₂*r[i-2] + ε[i]
  end

  return r
end

#gets the dynamic multiplier as a function
function dynamicAR2(Θ::AR2, k::Int; debug::Bool = false)::Float64

  #AR2 as function of shocks
  fAR2(ε::Vector{T} where T<:Real) = (AR2Returns!(Θ, ε))[end]
  ∇AR2::Function = ε::Vector -> ForwardDiff.gradient(fAR2, ε)

  if debug
    println("Case 1: $(∇AR2(rand(Θ.dist,k)))")
    println("Case 2: $(∇AR2(rand(Θ.dist,k)))")
    println("Case 3: $(∇AR2(rand(Θ.dist,k)))")

    ε::Vector{Float64} = rand(Θ.dist,k)
    Δh::Float64 = 10.^-5.
    εPΔh::Vector{Float64} = copy(ε)
    εPΔh[1] += Δh
    println("Case fd: $((fAR2(εPΔh)[1] - fAR2(ε)[1])/Δh)")
  end

  return ∇AR2(rand(Θ.dist,k))[1]
end




#problem specific graphs
function visualizeA4P1(Θ::Vector{AR2}, ks::Vector{Int}, ρₖ::Vector{Vector{Float64}})::Void
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)

  for i::Int ∈ 1:length(Θ)
    push!(plotNames, "ρₖ for AR2, (ϕ₁=$(n2s(Θ[i].ϕ₁)), ϕ₂=$(n2s(Θ[i].ϕ₂)))")
    push!(plots, plot(x=ks, y=ρₖ[i], Geom.line, Geom.point,
      Guide.xlabel("k"),Guide.ylabel("ρₖ"),
      Guide.title("ρₖ for ϕ₁=$(n2s(Θ[i].ϕ₁)), ϕ₂=$(n2s(Θ[i].ϕ₂))")))
  end

  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #display(plots[i])
  end

  return nothing
end

#script for the problem
function A4P1Script()::Void
  ϕ₁::Float64 = 1.1
  ϕ₂::Float64 = -0.25

  Θ::Vector{AR2} = [AR2(1.1,-0.25), AR2(0.9,0.8)]
  maxk::Int = 50
  ks::Vector{Int} = Vector{Int}(1:maxk)
  ρₖ::Vector{Vector{Float64}} = (i::Int->Vector{Float64}(maxk)).(1:length(Θ))
  println("Calculating auto correlations and dynamic multipliers")

  #calculating dynamic multiplier

  for i::Int ∈ 1:length(Θ)
    AR2ρ!(Θ[i],maxk,ρₖ[i])
    λ::Float64 = dynamicAR2(Θ[i],7,debug=true)
    println("Dynamic multiplier for ϕ₁=$(n2s(Θ[i].ϕ₁)), ϕ₂=$(n2s(Θ[i].ϕ₂)): $λ")
  end

  visualizeA4P1(Θ, ks, ρₖ)



  return nothing
end

@time begin
  A4P1Script()
end

end
