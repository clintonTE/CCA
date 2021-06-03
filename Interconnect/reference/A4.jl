module A4


#Use this to print
#=
weave(Pkg.dir("$(pwd())\\A4.jl"),
  informat="script",
  out_path = "$(pwd())\\A4.html",
  doctype = "md2html")
=#
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD,
  Gadfly, CTModCopy, JuMP, NLopt, ForwardDiff

const DATA_PATH = pwd() * "\\data" #path where we store the data
const OUTPUT_PATH = pwd() * "\\output"
const DATE_FORMAT_STR = "yyyy-mm-dd"
const SP500_NAME = "GSPC"
const TBILL_NAME = "t-bill"
const TERM_NAME = "termStructure"
const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
const HEADER_NAME = "header.tex"
const MAX_YEAR = 2016
const IN_SAMPLE_CUTOFF_YEAR = 1975

const DEFAULT_NUM_SIMS = 10000
const DEFAULT_LOGμ_TOL = 10.0^-10.0
const PlotContainer = Union{Plot,Gadfly.Compose.Context}
const CSRegSpec = Tuple{Symbol, Symbol, Int}
const CPRegSpec = Tuple{Symbol, Int}


#preprocess terms structure
#the default year cutoffes ensure the largest data set with a complete year
function preProcessTermStructure(;maxYear::Int = 2016, minYear::Int = 1962, decimalYields = true)::Vector{Symbol}
  dateFormat::String = DATE_FORMAT_STR
  termDF::DataFrame = readtable("$DATA_PATH\\$TERM_NAME.csv")

  #showcols(termDF)
  termDF[:,:date] = ((s::String)->Date(s, dateFormat)::Date).(termDF[:,:date])::DataVector{Date}
  sort!(termDF, cols = [:date])
  termDF[:,:year] = Dates.year.(termDF[:,:date])
  termDF = termDF[(termDF[:,:year] .≤ maxYear) && (termDF[:,:year] .≥ maxYear), :]

  #we will need this for the predictive regressions
  termDF[:,:dayOfQuarter] = zeros(Int,size(termDF,1))
  for i::Int ∈ 2:size(termDF,1)
    if Dates.firstdayofquarter(termDF[i,:date]) ≠ Dates.firstdayofquarter(termDF[i-1,:date])
      termDF[i,:dayOfQuarter] = 1
    else
      termDF[i,:dayOfQuarter] = termDF[i-1,:dayOfQuarter] + 1
    end
  end

  fieldSymbols = [((i::Int)->Symbol("SVENY", i<10?"0":"",i)).(collect(1:30));
    ((i::Int)->Symbol("SVENF", i<10?"0":"",i)).(collect(1:30))]


  if decimalYields
    for s::Symbol ∈ fieldSymbols
      termDF[:,s] *= 0.01
          #  termDF[:,s] ./= 250.0
    end

    for s::Symbol ∈ fieldSymbols[31:60]
      termDF[:,s] = exp.(termDF[:,s]) - 1.0
    end

  end

  stream::IOStream = open("$DATA_PATH\\$TERM_NAME.jls", "w")
  serialize(stream, termDF)
  close(stream)

  return fieldSymbols
end


mutable struct BYParam
    β::Float64
    g::Float64
    γ1::Float64
    ϕg::Float64
    v::Float64
    ϕv::Float64
    ν0::Float64

    α::Float64
    ρ::Float64

    logμ::Float64

end
#=
From paper: (δ = 0.998, μ = μd = 0.0015, ρ = 0.979, σ = 0.0078, φ =3,
 ϕe = 0.044, and ϕd = 4.5), the parameters of the stochastic volatility process are ν1 = 0.987
and σw = 0.23 × 10−5

Notation Mapping:
β = δ
g = μ
ϕv = ν1
γ1 = ϕe
v = σ^2
ϕg = ρ
ν0 = σw

θ = (1-γ)/(1-1/ψ)
ρ= -θ/ψ+1
α = ρ -1 - θ



=#

#Use BY2004 calibration as the default
function BYParam(;
    β::Float64 = 0.998,
    g::Float64 = 0.0015,
    γ1::Float64 = 0.044,
    ϕg::Float64 = 0.979,
    v::Float64 = 0.0078 * 0.0078,
    ϕv::Float64 = 0.987,
    ν0::Float64 = 0.23 * 10.0^-5.0,
    γ::Float64 = 10.0 ,
    ψ::Float64 = 1.5)::BYParam

    θ::Float64 = (1.0-γ)/(1.0-1.0/ψ)
    ρ::Float64 = 1.0-θ/((θ-1)*ψ)
    α::Float64 = ρ*(θ - 1.0) +ρ

    logμ::Float64 = 0.05

    return BYParam(β, g, γ1, ϕg, v, ϕv, ν0, α, ρ, logμ)
end


b1(p::BYParam) = p.β*exp(p.ρ*p.logμ)/((1-p.β)+p.β*exp(p.ρ*p.logμ))
#b1Check(p::BYParam) = p.β * exp(p.ρ*p.logμ)/((1-p.β)+p.β*exp(p.ρ*p.logμ))
b0(p::BYParam) = 1/p.ρ*log(1-p.β+p.β*exp(p.ρ*p.logμ))-b1(p)*p.logμ
#b0Check(p::BYParam) = 1/p.ρ*log((1-p.β)+p.β*exp(p.ρ*p.logμ)) -b1(p)*p.logμ

px(p::BYParam)::Float64 = b1(p)/(1-p.ϕg*b1(p))
#pxCheck(p::BYParam)::Float64 = b1(p)/(1-p.ϕg*b1(p))

 pv(p::BYParam)::Float64 = b1(p)*p.α*(1+px(p)^2*p.γ1^2)/(2*(1-p.ϕv*b1(p)))
#pvCheck(p::BYParam)::Float64 = b1(p)*p.α*(1+px(p)^2*p.γ1^2)/(2*(1-p.ϕv*b1(p)))

u(p::BYParam)::Float64 = (b0(p)+b1(p)*(p.g+pv(p)*(1-p.ϕv)*p.v+p.α/2*pv(p)^2*p.ν0^2))/(1-b1(p))
#uCheck(p::BYParam)::Float64 = (b0(p)+b1(p)*(p.g+pv(p)*(1-p.ϕv)*p.v+p.α/2*pv(p)^2*p.ν0^2))/(1-b1(p))

ERP(p::BYParam)::Float64 = 0.5*(p.ρ-1)^2.0*p.v+0.5*(p.α-p.ρ)^2*(p.v*(1+px(p)^2*p.γ1^2)+pv(p)^2*p.ν0^2)

Em(p::BYParam)::Float64 = log(p.β)+(p.ρ-1)*(p.g+p.v)-p.α/2*(p.v*(1+px(p)^2*p.γ1^2)+pv(p)^2*p.ν0^2)*(p.α-p.ρ)
#EmCheck(p::BYParam)::Float64 = log(p.β)+(p.ρ-1)*(p.g+p.v)-p.α/2*(p.v*(1+px(p)^2*p.γ1^2)+pv(p)^2*p.ν0^2)*(p.α-p.ρ)

Erfr(p::BYParam)::Float64 = -Em(p)-ERP(p)

ELogμt(p::BYParam)::Float64 =
  u(p)+p.g+p.v*(p.α/2*(1+px(p)^2*p.γ1^2)+pv(p)*p.ϕv)+pv(p)*(1-p.ϕv)*p.v+p.α/2*pv(p)^2*p.ν0^2


function solveForLogμ!(p::BYParam, tol = DEFAULT_LOGμ_TOL, limit=10^7)
  ctr::Int = 0
  last::Float64 = 9999.0


  while abs(p.logμ-last) > tol && ctr< limit
    last = p.logμ
    p.logμ = ELogμt(p)
    ctr+=1
    #print(" ",p.logμ)
  end


  if ctr == limit
    println("WARNING: Iteration Limit Reached")
  end

end


function A4P1Script()
  p::BYParam = BYParam()
  solveForLogμ!(p)

  println("\n",p)
  println("ERP: ", ERP(p))
  println("rfr: ", Erfr(p))

  println("\nAlternate Configuration ϕv=0.999")
  p = BYParam(ϕv=0.999)
  solveForLogμ!(p)
  println("ERP: ", ERP(p))
  println("rfr: ", Erfr(p))


end


#assume yns is in order
function setupCS!(termDF::DataFrame, yns::Vector{Symbol}, ns::Vector{Int}=collect(1:length(yns)))

  T::Int =  size(termDF, 1)
  seriesNames::Vector{CSRegSpec} = Vector{CSRegSpec}(length(yns)-1)

  daysPerPeriod = size(termDF,1) ÷ (maximum(termDF[:,:year]) - minimum(termDF[:,:year] + 1))
  println("Days per period: $daysPerPeriod")

  for i::Int ∈ 2:length(yns)
    sDif = Symbol(yns[i],"_$(ns[i])D")
    sSlope = Symbol(yns[i],"_$(ns[i])S")

    termDF[:, sDif] = 0.0
    termDF[:, sDif] = NA

    termDF[:, sSlope] = 0.0
    termDF[:, sSlope] = NA

    termDF[1:(end-daysPerPeriod), sDif] =
      termDF[(daysPerPeriod+1):end, yns[i-1]] .- termDF[1:(end-daysPerPeriod), yns[i]]

    termDF[1:(end-daysPerPeriod), sSlope] = (termDF[1:(end-daysPerPeriod), yns[i]] .-
      termDF[1:(end-daysPerPeriod), yns[1]]) ./ (ns[i] - 1.0)

    seriesNames[i-1] = (sDif, sSlope, ns[i])
  end

  return seriesNames
end

function runCS(termDF::DataFrame, seriesNames::Vector{CSRegSpec})

  #run the regression

  regs::Vector{CTLM} = Vector{CTLM}(length(seriesNames))

  for i::Int ∈ 1:length(seriesNames)
    regs[i] = CTLM(termDF, seriesNames[i][2], seriesNames[i][1],
      XNames = [:intercept, :phi], YName=:Y)
  end

  return regs
end

function getCSTable(termDF::DataFrame, seriesNames::Vector{CSRegSpec};
  titleCaption="Campbell Shiller Regressions")

  regs::Vector{CTLM} = runCS(termDF, seriesNames)


  descRowNames::Vector{String} = ["N"]
  descContent::Vector{Vector{String}} =
      [Vector{String}(length(regs)) for i::Int ∈ 1:length(descRowNames)]

  for i ∈ 1:length(regs)
    descContent[1][i] = "$(regs[i].N)"
  end

  tableText::String = texTable(regs, getModWhiteΣ!, [:intercept; :phi],
      titleCaption = titleCaption,
      colNames = [((t::CSRegSpec)->"$(t[3])").(seriesNames)],
      contentRowNames = ["intercept", "\$\\phi\$"],
      descRowNames = descRowNames,
      descContent = descContent,
      decimalDigits = 2,
      columnSepPt = 0,
      scaling = [1.0, 1.0],
      caption = "")

  return tableText
end

function setupCP!(termDF::DataFrame, yns::Vector{Symbol}, fns::Vector{Symbol})

  T::Int =  size(termDF, 1)
  seriesNames::Vector{CPRegSpec} = Vector{CPRegSpec}(length(yns)-1)
  ns::Vector{Int} = collect(1:length(yns))
  numMaturies::Float64 = length(seriesNames)

  for i::Int ∈ 2:length(yns)
    syx = Symbol(yns[i],"_$(ns[i])rxp1")

    termDF[:, syx] = 0.0
    termDF[:, syx] = NA
    termDF[1:(end-1), syx] = termDF[2:end, yns[i]] .- termDF[2:end, yns[1]]

    seriesNames[i-1] = (syx, ns[i])
  end

  termDF[:, :rxBarp1] = 0.0
  termDF[:, :rxBarp1] = NA

  # now need to get the average excess return
  @fastmath for i::Int = 1:T
    tot::Float64 = 0.0
    foundNA::Bool = false #will not record an entry if there is a missing value

    for s ∈ seriesNames
      if isna(termDF[i,s[1]])
        foundNA = true
      else
        tot += termDF[i,s[1]]
      end
    end

    termDF[i, :rxBarp1] = tot/numMaturies
  end

  return seriesNames
end

function getγ!(termDF::DataFrame, fns::Vector{Symbol},
  rowCutoff::Int = size(termDF,1))::Vector{Float64}

  T::Int = size(termDF,1)
  γNames = ((i::Int)->Symbol("gamma",i)).(0:length(fns))

  γf::DataVector{Float64} = DataVector{Float64}(Vector{Float64}(T))
  γf .= zeros(T)

  XExpr::CTExpr = parse(vec2String(fns, "+"))

  #println(XExpr)
  #run the regression
  reg::CTLM= CTLM(termDF[1:rowCutoff,:], XExpr, :rxBarp1, XNames = γNames, YName=:rxBarp1)

  γf .+= reg.β[1]

  for i ∈ 1:length(fns)
    γf .+= reg.β[i+1] .* termDF[:,fns[i]]
  end

  termDF[:, :γf] = γf

  return reg.β
end

#runs the CP regression
function runCP(termDF::DataFrame, seriesNames::Vector{CPRegSpec}, rowCutoff::Int = size(termDF,1))::Vector{CTLM}

  #run the regression

  regs::Vector{CTLM} = Vector{CTLM}(length(seriesNames))

  XExpr::CTExpr = parse("γf+0")
  for i::Int ∈ 1:length(seriesNames)
    regs[i] = CTLM(termDF[1:rowCutoff,:], XExpr, seriesNames[i][1],
      XNames = [:bn], YName=:Y)
  end

  return regs
end

function getCPbns(termDF::DataFrame, seriesNames::Vector{CPRegSpec}, rowCutoff::Int = size(termDF,1))::Vector{Float64}
  regs::Vector{CTLM} = runCP(termDF, seriesNames, rowCutoff)

  return ((i::Int)->(regs[i].β[1])).(1:length(seriesNames))
end


function getCPTable(termDF::DataFrame, seriesNames::Vector{CPRegSpec}, fns::Vector{Symbol};
  titleCaption="Cochrane Piazzesi Regressions (longest maturity=$(length(seriesNames)+1))")

  γ::Vector{Float64} = getγ!(termDF, fns)
  caption::String = vec2String(((i::Int)->"\$\\gamma_{$(i-1)}=$(round(γ[i],1))\$ ").(1:length(γ)),"")
  println(caption)
  regs::Vector{CTLM} = runCP(termDF, seriesNames)

  descRowNames::Vector{String} = ["\$R^{2}\$", "N"]
  descContent::Vector{Vector{String}} =
      [Vector{String}(length(regs)) for i::Int ∈ 1:length(descRowNames)]

  for i ∈ 1:length(regs)
    descContent[1][i] = "$(round(getR(regs[i])^2.0,2))"
        descContent[2][i] = "$(regs[i].N)"
  end

  Gadfly.push_theme(:default)
  p = plot(x=collect(0:(length(γ)-1)),
    y=γ, Geom.line,
    Guide.ylabel("γ (restricted)"), Guide.xlabel("gamma value"))

  draw(SVG("$OUTPUT_PATH\\gammaPlot_$(length(γ)-1).svg", 7.5inch, 4.5inch),p)

  tableText::String = texTable(regs, getModWhiteΣ!, [:bn],
      titleCaption = titleCaption,
      colNames = [((t::CPRegSpec)->"$(t[2])").(seriesNames)],
      contentRowNames = ["\$b_n\$"],
      descRowNames = descRowNames,
      descContent = descContent,
      decimalDigits = 2,
      columnSepPt = 20,
      scaling = [1.0],
      caption = caption)

  return tableText
end


#run predictive regressions
function CPPredictive!(termDF::DataFrame, inSampleCutoffYear::Int, seriesNames::Vector{CPRegSpec},
  yns::Vector{Symbol}, fns::Vector{Symbol};
  predictedSyms::Vector{Symbol} = ((i::Int)->Symbol(:predicted, i)).(1:length(seriesNames)))

  T::Int = size(termDF,1)
  N::Int = length(seriesNames)

  #pre-allocate
  γ::Vector{Float64} = Vector{Float64}(length(fns)+1)
  bns::Vector{Float64} = Vector{Float64}(N)

  #set the training sample
  NTrain::Int = size(termDF[termDF[:,:year] .≤ inSampleCutoffYear,:],1)
  startRow::Int = NTrain+1

  #allocate the following column
  for s::Symbol ∈ predictedSyms
    termDF[:,s] = 0.0
    termDF[:,s] = NA
  end

  #run the simulations
  @fastmath for i::Int ∈ startRow:size(termDF,1)

    if termDF[i,:dayOfQuarter] == 1 #train the data if its the beginning of the quarter
      γ .= getγ!(termDF, fns, i-1) #Only include data through the end of the prior quarter
      bns .= getCPbns(termDF, seriesNames, i-1)
    end

    #now get the predictions
    for n::Int ∈ 1:N
      termDF[i,predictedSyms[n]] = bns[n] * termDF[i,:γf]
    end

  end

  #graph the predictions
  palette = ["Black", "DodgerBlue"]
  plottedDF = termDF[startRow:end,[:date; predictedSyms; [seriesNames[i][1] for i ∈ 1:N]]]

  for n::Int ∈ 1:N

    RMSN::Float64 = length(dropna(
      (plottedDF[:,predictedSyms[n]] .- plottedDF[:,seriesNames[n][1]])))

    RMSDeviation::Float64 = (sum(dropna(
      (plottedDF[:,predictedSyms[n]] .- plottedDF[:,seriesNames[n][1]]).^2))/(RMSN-1))^0.5
    p = plot(
      layer(plottedDF[completecases(plottedDF[:,[:date,predictedSyms[n]]]),:], x="date", y=predictedSyms[n],
        Geom.line,Theme(line_width=0.5pt, default_color=palette[1])),
      layer(plottedDF[completecases(plottedDF[:,[:date,seriesNames[n][1]]]),:], x="date", y=seriesNames[n][1],
        Geom.line,Theme(line_width=0.5pt, default_color=palette[2])),
      Guide.ylabel("Excess Return"), Guide.xlabel(nothing),
      style(key_position = :bottom, key_max_columns = 6, key_label_font_size=8pt, key_label_font="courier"),
      Guide.title("Predicted sv Realized n=$(n+1), longest maturity=$(N+1) (Tracking σ = $(round(RMSDeviation*100.0,4))% points)"),
      Guide.manual_color_key("Legend",  ["actual", "predicted"], palette))

    draw(SVG("$OUTPUT_PATH\\predicted_$(N+1)_$(n+1).svg", 7.5inch, 3.5inch),p)
  end



end




function A4P2Script(;refreshData=true, runCP=true)::Void

  CSRegRange::UnitRange = 1:7
  const CSRegSymbols = ((i::Int)->Symbol("SVENY", i<10?"0":"",i)).(CSRegRange)


  if !isfile("$DATA_PATH\\$TERM_NAME.jls") || refreshData
    preProcessTermStructure()
  end

  #read the binary file
  stream::IOStream = open("$DATA_PATH\\$TERM_NAME.jls")
  termDF::DataFrame = deserialize(stream)
  close(stream)

  fieldSymbols = [((i::Int)->Symbol("SVENY", i<10?"0":"",i)).(collect(1:30));
    ((i::Int)->Symbol("SVENF", i<10?"0":"",i)).(collect(1:30))]

  #run the full sample CS regressions
  #termDFCopy::DataFrame = deepcopy(termDF)
  seriesNamesCS::Vector{CSRegSpec} = setupCS!(termDF, CSRegSymbols#=, ((i::Int)->(12*i)).(CSRegRange)=#)
  CSTable::String = getCSTable(termDF, seriesNamesCS)

  if runCP
    #####get first five forward rates
    fns::Vector{Symbol} = [fieldSymbols[1]; fieldSymbols[31:34]]
    CPRegSymbols::Vector{Symbol} = CSRegSymbols[1:5]
    #run the full sample CP regressions
    seriesNamesCP::Vector{CPRegSpec} = setupCP!(termDF, CPRegSymbols, fns)
    CPTable::String = getCPTable(termDF, seriesNamesCP, fns)

    #Now do the predictive regressions
    CPPredictive!(termDF, IN_SAMPLE_CUTOFF_YEAR, seriesNamesCP,
      CPRegSymbols, fns)

    #####Now do it with 7
    #termDF = termDFCopy
    fns = [fieldSymbols[1]; fieldSymbols[31:36]]
    CPRegSymbols = CSRegSymbols[1:7]
    #run the full sample CP regressions
    seriesNamesCP = setupCP!(termDF, CPRegSymbols, fns)
    CPTableAll::String = getCPTable(termDF, seriesNamesCP, fns)

    #Now do the predictive regressions
    CPPredictive!(termDF, IN_SAMPLE_CUTOFF_YEAR, seriesNamesCP,
    CPRegSymbols, fns)
  else
    CPTable = ""
    CPTableAll = ""
  end

  writeTables2File([CSTable, CPTable, CPTableAll],
      HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
      outName = "RegTables_$A4.tex")

  #Run the full-sample Campbell-Shiller Regressions

  return nothing
end

@time A4P2Script(refreshData=false, runCP=true)
@time A4P1Script()

end
