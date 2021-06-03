module TimeSeriesA2P2

#Use this to print
#=
using Weave
codeName = "TimeSeriesA2P2"
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
const DECIMALS = 3
const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const OVERRIDE_STAR_STRINGS = ["\\ostar","\\dstar","\\tstar"]

#asignment specific constants
const FACTOR_NAME = "FMFactors"
const ASSET_NAME = "PORTSBM55"
const CUTOFF_DATE = 1963_07
const FACTOR_SYMS = [:Mkt_RF, :SMB, :HML]
const RFR_SYM = :RF
const MKT_SYM = :Mkt


BLAS.set_num_threads(1)

  n2s(x::T where T<:Real) = num2Str(x, DECIMALS, Ints=true, scaleHurdle=99.)

#pre-process the asset data
function preProcessMonthlyFA()::Void

  #read in the data (note the non-standard NA values)
  assetDF::DataFrame = readtable("$DATA_PATH\\$ASSET_NAME.csv", nastrings=["-999", "-99.99"])
  factorDF::DataFrame = readtable("$DATA_PATH\\$FACTOR_NAME.csv", nastrings=["-999", "-99.99"])

  factorDF[:,:Mkt] = factorDF[:,:Mkt_RF] .+ factorDF[:,:RF]

  #merge the dataframes
  DF::DataFrame = join(factorDF, assetDF, on = :DATE, kind=:inner)

  #write out a JLS binary file
  stream::IOStream = open("$DATA_PATH\\$(FACTOR_NAME)_$(ASSET_NAME).jls", "w")
  serialize(stream, DF)
  close(stream)
  return nothing
end

#gets the dataframe as a serial file
function getDF(pathStr::String, refreshData::Bool, preProcessFunc::Function)::DataFrame

  #re-process the data if needed or desired
  if !isfile(pathStr) || refreshData
    preProcessFunc()
  end

  #read the binary file
  stream::IOStream = open(pathStr)
  DF::DataFrame = deserialize(stream)
  close(stream)

  return DF
end

mutable struct PCA
  λ::Vector{Float64}
  ν::Vector{Vector{Float64}}
  svd::Base.LinAlg.SVD
  labels::Vector{Symbol}
  labelsInd::Dict
  normConstant::Float64
end

function PCA(Σ::Matrix{Float64},
    labels::Vector{Symbol} = ((i::Int)->Symbol("v$i")).(1:size(Σ,1)))::PCA
  #setup and factorize
  M::Int, N::Int = size(Σ)
  labelsInd::Dict = Dict(labels[i]=>i for i::Int ∈ 1:M)

  Σsvd::Base.LinAlg.SVD = svdfact(Σ)

  #recover the eigenvalues
  λ::Vector{Float64} = Σsvd[:S]
  normConstant::Float64 = sum(λ)
  λ ./= normConstant

  #recover the eigenvectors
  ν::Vector{Vector{Float64}} =
    ((i::Int)->(Σsvd[:U][:,i]./sum((abs2).(Σsvd[:U][:,i])))).(1:N)

  return PCA(λ, ν, Σsvd,labels,labelsInd,normConstant)
end

function PCA(DF::DataFrame, labels::Vector{Symbol})
  D::Matrix{Float64} = Matrix{Float64}(size(DF[:,labels]))
  D .= Matrix{Float64}(DF[:,labels])
  D .-= mean(D,1)
  return PCA(BLAS.gemm('T','N',1.0, D,D), labels) #(This is simply D'D)
end

#adds principal componenets to an existing dataframe
function PC(pca::PCA, DF::DataFrame;
    labels::Vector{Symbol}=pca.labels,
    dataLabels::Vector{Symbol}=labels,
    PCLabels::Vector{Symbol} = ((i::Int)->(Symbol("PC$(i)"))).(1:length(labels)),
    )::Tuple{DataFrame, Vector{Symbol}}

    D::Matrix{Float64} = DF[:,dataLabels]
    for i::Int ∈ 1:length(PCLabels)
      DF[:,PCLabels[i]] = BLAS.gemv('N', 1.0, D, pca.ν[pca.labelsInd[labels[i]]])
    end

    return (DF, PCLabels)
end

#Calculates the summary statistics
function summaryA2P2(FADF::DataFrame, PCLabels::Vector{Symbol}, cols::Vector{Symbol})
  #label the table
  colNames::Vector{Vector{String}} =
    [((s::Symbol)->replace(String(s),"_","")).(cols)]
  numCols = length(colNames[end])

  rowNames::Vector{String} = ((s::Symbol)->replace(String(s),"_","")).(PCLabels)

  numRows::Int = length(rowNames)
  content::Vector{Vector{String}} = [Vector{String}(numCols) for i::Int ∈ 1:numRows]

  #run the tests and record the results
  for r::Int ∈ 1:numRows, c::Int ∈ 1:numCols
    content[r][c] = n2s(cor(FADF[:,PCLabels[r]], FADF[cols[c]]))
  end

  #write the table
  summaryTable::String = texTable( "Correlations of Principal Components",
    """See Tex File""", #caption
    colNames, #colNames
    Vector{String}(),#contentRowNames
    Vector{Matrix{String}}(), #content
    rowNames, #descRowNames
    content, #descContent
    Vector{String}(), #notes
  )

  return summaryTable
end

#problem specific visualizations
function visualizeA2P2(pca::PCA, FADF::DataFrame, PCLabels::Vector{Symbol})::Void

  #convenience containers
  K::Int = length(pca.λ)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  Gadfly.push_theme(:default)

  #some basic data processing for the eigenvalues
  λLabels::Vector{Symbol} = ((i::Int)->(Symbol("λ$i"))).(1:K)
  λDataLabels::Vector{String} = (n2s).(pca.λ)
  λCum::Vector{Float64} = ((i::Int)->sum(pca.λ[1:i])).(1:K)
  λCumDataLabels::Vector{String} = (n2s).(λCum)

  push!(plotNames, "eigenvalues")
  push!(plots, plot(y=pca.λ, x=λLabels, label=λDataLabels, Geom.bar, Geom.label(position=:above),
    Guide.xlabel("λ"),Guide.ylabel("λᵢ / sum(λ)"),Coord.Cartesian(ymin=0., ymax=1.05),
    Guide.title("Eigenvalues λ")))

  push!(plotNames, "cumulative eig")
  push!(plots, plot(y=λCum, x=λLabels, label=λCumDataLabels,
    Geom.point, Geom.line, Geom.label(position=:above),
    Guide.xlabel("λ"),Guide.ylabel("Cumulative λᵢ / sum(λ)"),Coord.Cartesian(ymin=0., ymax=1.05),
    Guide.title("Cumulative Eigenvalues λ")))

  PC3DF::DataFrame = melt(DataFrame(labels=pca.labels,
    ν1=pca.ν[1], ν2=pca.ν[2], ν3=pca.ν[3]),:labels)
  push!(plotNames, "eigenvectors")
  push!(plots, plot(PC3DF, y=:value, x=:labels, color=:variable,
    Geom.bar(position=:dodge),# Geom.line, Geom.point,
    Guide.xlabel("λ"), Guide.ylabel("νᵢ/sum(ν.²)"),
    Guide.title("Eigenvectors ν")))

  push!(plotNames, "eigenvectors ν₁ vs ν₂")
  push!(plots,
    plot(x=pca.ν[1], y=pca.ν[2], label=(string).(pca.labels),
      Geom.point,Geom.label,
      Guide.xlabel("ν₁"), Guide.ylabel("ν₂"),
      Coord.Cartesian(ymin=-0.5, ymax=0.5, xmin=-0.5, xmax=0.5),
      Guide.title("Eigenvectors ν₁ vs ν₂")))
  push!(plotNames, "eigenvectors ν₂ vs ν₃")
  push!(plots,
    plot( x=pca.ν[2], y=pca.ν[3], label=(string).(pca.labels),
      Geom.point,Geom.label,# Geom.line, Geom.point,
      Guide.xlabel("ν₂"), Guide.ylabel("ν₃"),
      Coord.Cartesian(ymin=-0.5, ymax=0.5, xmin=-0.5, xmax=0.5),
      Guide.title("Eigenvectors ν₂ vs ν₃")))

  push!(plotNames, "eigenvectors ν₃ vs ν₁")
  push!(plots,
    plot(x=pca.ν[3], y=pca.ν[1], label=(string).(pca.labels),
      Geom.point,Geom.label,# Geom.line, Geom.point,
      Guide.xlabel("ν₃"), Guide.ylabel("ν₁"),
      Coord.Cartesian(ymin=-0.5, ymax=0.5, xmin=-0.5, xmax=0.5),
      Guide.title("Eigenvectors ν₃ vs ν₁")))



  #write the graphs
  for i::Int ∈ 1:length(plots)
    draw(SVG("$OUTPUT_PATH\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
  end

  return nothing

end

#problem specific code
function A2P2Script(refreshData::Bool = false)::Void
  FADF::DataFrame = getDF("$DATA_PATH\\$(FACTOR_NAME)_$(ASSET_NAME).jls",
    refreshData, preProcessMonthlyFA)
  FADF = FADF[FADF[:,:DATE].≥196307,:]

  #rename the columns to a common pattern
  rename!(FADF, [:SMALL_LoBM, :SMALL_HiBM, :BIG_LoBM, :BIG_HiBM],
    [:ME1_BM1, :ME1_BM5, :ME5_BM1, :ME5_BM5])

  #get portfolio handles
  FSyms::Vector{Symbol} = FACTOR_SYMS
  RFRSym::Symbol = RFR_SYM
  RSyms::Vector{Symbol} =
    names(FADF[:,setdiff(names(FADF), [:DATE; RFR_SYM; MKT_SYM; FACTOR_SYMS])])

  #do the eigenvalue decomposition
  Σpca::PCA = PCA(FADF, RSyms)

  #get the principal components
  (FADF, PCLabels::Vector{Symbol}) = PC(Σpca, FADF)
  visualizeA2P2(Σpca, FADF, PCLabels)

  summaryTable::String = summaryA2P2(FADF,PCLabels[1:5],[RFR_SYM; MKT_SYM; FACTOR_SYMS])

  writeTables2File([summaryTable],
    HEADER_NAME, FOOTER_NAME, path=OUTPUT_PATH,
    outName = "TimeSeries A2P2 corTable.tex")

  return nothing
end

#uncomment to run
@time begin
  A2P2Script()
end


end
