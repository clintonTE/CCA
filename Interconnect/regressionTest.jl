module scratchpad

#Use this to print
#=
using Weave
codeName = "TimeSeriesA6"
weave(Pkg.dir("$(pwd())\\$(codeName).jl"),
  informat="script",
  out_path = "$(pwd())\\output\\$(codeName)_Appendix.html",
  doctype = "md2html")
=#

#NOTE: To prepare the data, conform missing data to "". That means eliminating
#missing codes "N/A", "missing", and "--"

if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, CSV, StatsBase, GZip, ForwardDiff, Query,
  Gadfly, NLopt, HCubature, StaticArrays, Measures, Formatting


importall CTIO, CTReg, CTStat, CTVAR, CTNumerical

#############parameters
const VOL_PERIODS = 1

#path constants
const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT = "m/d/yyyy"
const CHINA_DS_NAME = "ChinaDS1.1"
const CHINA_DS_NAME_L = "$(CHINA_DS_NAME)_L"
const CHINA_DS_NAME_W = "$(CHINA_DS_NAME)_W"
const R_TEST_PATH = "C:\\Users\\Clinton\\Dropbox\\Projects\\InterconnectRTester\\trainTester\\trainTester"
const ACCEPTABLE_RANGE = #use this to check if values are in acceptable ranges
    Dict(:value=>(-10., 10.), :marketCap=>(.1*10.^9.,Inf))
const DROP_OUTSIDE_RANGE = true
const MIN_SZSECOM_MARKETCAP = 6.*10.^9.
const KEEP_MISSING_MARKETCAP = false

const MIN_SAMPLES_IN_INDUSTRY = 5
#const MIN_SAMPLES_INTERACTED = 50 #obsolete
const MIN_COMPANY_QUARTERS = 1 #threshold for dropping companies
const SET_MISSING_IF_BELOW = true
const KEEP_LARGEST_WEIGHT_IF_DUP = true # keep the largest weight if multiple shareclasses

const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const DECIMALS = 2
const NUM_FORMAT = "{:.$(DECIMALS)f}"
#const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs

const DATE_SYMBOLS =  Symbol("12/31/2010"), Symbol("3/31/2011"), Symbol("6/30/2011"),
  Symbol("9/30/2011"), Symbol("12/31/2011"), Symbol("3/31/2012"), Symbol("6/30/2012"),
  Symbol("9/30/2012"), Symbol("12/31/2012"), Symbol("3/31/2013"), Symbol("6/30/2013"),
  Symbol("9/30/2013"), Symbol("12/31/2013"), Symbol("3/31/2014"), Symbol("6/30/2014"),
  Symbol("9/30/2014"), Symbol("12/31/2014"), Symbol("3/30/2015"), Symbol("6/30/2015"),
  Symbol("9/30/2015"), Symbol("12/31/2015"), Symbol("3/30/2016"), Symbol("6/30/2016"),
  Symbol("9/30/2016"), Symbol("12/31/2016"), Symbol("3/30/2017"), Symbol("6/30/2017"),
  Symbol("9/30/2017"), Symbol("12/31/2017"), Symbol("3/31/2018")

const FOCAL_DATES = Date("12/31/2010",  DATE_FORMAT):Dates.Month(3):Date("3/31/2018",  DATE_FORMAT)
const OUTCOMES = ["NORMALIZED_ACCRUALS_CF_METHOD"]

const MFloat64 = Union{Missing, Float64}
const MInt = Union{Missing, Int}
const MBool = Union{Missing, Bool}
const MSymbol = Union{Missing, Symbol}

n2s(x::T where T<:Real) = num2Str(x, 3, Ints=true)


function script(N::Int = 10^4)
X::Matrix{Float64}= [ones(N) rand(N) rand(N).*10000.]
Y::Vector{Float64} = X[:,2].*10 + rand(N).*5

println("First regression results:")
lm1::CTLM = CTLM(X[:,1:2],Y)
println("Coef: $(lm1.β)")
println("SE: $(diag(getWhiteΣ!(lm1)).^0.5)")

println("Noisy regression results:")
lm2::CTLM = CTLM(X,Y)
println("Coef: $(lm2.β)")
println("SE: $(diag(getWhiteΣ!(lm2)).^0.5)")

end

script()

end
