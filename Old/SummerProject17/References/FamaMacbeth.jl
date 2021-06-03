module CTFM

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
  Gadfly, CTModCopy, JuMP, NLopt, HCubature, StaticArrays, CTNumerical,
  ParallelAccelerator, Measures

  #=const WORKERS = 3 #seems like logical cores * (3/4) is the magic value
  const CLEAN_WORKERS = false

  if nworkers() < WORKERS || (nworkers()==1 && WORKERS == 1)
    addprocs(nworkers()==1?WORKERS-nworkers()+1:WORKERS-nworkers())
  end

  @everywhere using DataFrames, Distributions, StatsBase, GZip, JLD,
    Gadfly, CTModCopy, JuMP, NLopt, HCubature, StaticArrays, CTNumerical=#

  const constDef = true
  const DATA_PATH = pwd() * "\\data" #path where we store the data
  const OUTPUT_PATH = pwd() * "\\output"
  const DATE_FORMAT_STR = "yyyymmdd"
  const SP500_NAME = "GSPC"
  const TBILL_NAME = "t-bill"
  const FACTOR_NAME = "FMFactors"
  const ASSET_NAME = "PORTSBM1010"
  #const ASSET_NAME = "PORTSBM23"
  #const ASSET_NAME = "PORTPI55"
  const FOOTER_NAME = "footer.tex" #these are just headers for the tex tables
  const HEADER_NAME = "header.tex"

  const PlotContainer = Union{Plot,Gadfly.Compose.Context}
  const CSRegSpec = Tuple{Symbol, Symbol, Int}
  const CPRegSpec = Tuple{Symbol, Int}

  const FACTOR_SYMS = [:Mkt_RF, :SMB, :HML]
  #const FACTOR_SYMS = [:Mkt_RF]
  const RFR_SYM = :RF

  const T_LIMIT = 99999
  #const T_LIMIT = 100


  BLAS.set_num_threads(1) #need this since we are already multi-threading at a higher level

  #pre-process the asset data
  function preProcessFA()::Void

    #read in the data (note the non-standard NA values)
    dateFormat::String = DATE_FORMAT_STR
    assetDF::DataFrame = readtable("$DATA_PATH\\$ASSET_NAME.csv", nastrings=["-999", "-99.99"])
    assetDF[:,:DATE] = ((i::Int)->Date("$i", dateFormat)::Date).(assetDF[:,:DATE])::DataVector{Date}

    factorDF::DataFrame = readtable("$DATA_PATH\\$FACTOR_NAME.csv", nastrings=["-999", "-99.99"])
    factorDF[:,:DATE] = ((i::Int)->Date("$i", dateFormat)::Date).(factorDF[:,:DATE])::DataVector{Date}
    factorDF = factorDF[:,[:DATE; FACTOR_SYMS; RFR_SYM]]

    #merge the dataframes
    DF::DataFrame = join(factorDF, assetDF, on = :DATE, kind=:inner)

    #write out a JLS binary file
    stream::IOStream = open("$DATA_PATH\\$(FACTOR_NAME)_$(ASSET_NAME).jls", "w")
    serialize(stream, DF)
    close(stream)
    return nothing
  end



struct CTFMWindow
  FADF::DataFrame
  FSyms::Vector{Symbol}
  RSyms::Vector{Symbol}

  β::Matrix{Float64} #(NxK)
  α::DataMatrix{Float64} #(NxT)
  λ::Matrix{Float64} #(KxT)

  T::Int #T time slices
  N::Int #N assets
  K::Int #k factors
end

#basic constructor to get coefficients
function CTFMWindow(FADF::DataFrame, FSyms::Vector{Symbol}, RSyms::Vector{Symbol})

  BLAS.set_num_threads(1) #we will do our own multi-threading here
  N::Int = length(RSyms)
  K::Int = length(FSyms)
  T::Int = min(size(FADF,1), T_LIMIT)
  FADF = FADF[(end - T + 1):(end),:]

  #println("size: ",size(FADF,1))

  #prepare the first stage
  #form the RHS
  FExpr::CTExpr = parse(join((String).(FSyms), " + "))
  RHSNames::Vector{Symbol} = [:intercept; FSyms]

  #run the N 1st stage regressions
  β = Matrix{Float64}(N,K)
  @fastmath Threads.@threads for i::Int ∈ 1:N
    β[i,:] .= ((CTLM(FADF, FExpr, RSyms[i], XNames=RHSNames)).β)[2:end]
  end

  #prepare the 2nd stage
  βDF = DataFrame(β)

  names!(βDF, FSyms)
  TSyms::Vector{Symbol} = ((i::Int)->Symbol("T$(i)")).(1:T)

  #transpose the R Matrix and append to the right of the factor results
  FADFRows = eachrow(FADF[:,RSyms])
  vToTranspose = DataVector{Float64}(Vector{Float64}(N))
  for t::Int ∈ 1:T
    βDF[:,TSyms[t]] = DataVector{Float64}(Vector{Float64}(N)) #allocate space
    for i::Int ∈ 1:N # copy the entries 1 by 1
      vToTranspose[i] = FADFRows[t][i]
    end
    βDF[:,TSyms[t]] .= vToTranspose
  end

    #showall(βDF)

  α::DataArray{Float64,2} = DataArray{Float64,2}(Matrix{Float64}(N,T))
  α .= NA

  #run the 2nd stage regressions
  cλ::Matrix{Float64} = Matrix{Float64}(K+1,T)
  @fastmath Threads.@threads for t::Int ∈ 1:T
    cλ[:,t] .= CTLM(βDF,FExpr,TSyms[t],XNames=RHSNames).β
  end

  @fastmath Threads.@threads for t::Int ∈ 1:T
    α[:,t] = βDF[:,TSyms[t]] .- cλ[1,t]
  end

  α -= β * cλ[2:end,:]

#=println("cλ1: ", cλ[1,1])
println("βDF[:,TSyms[1]]: ", βDF[:,TSyms[1]])
println("β * cλ[2:end,t]:", β * cλ[2:end,1])
println("α1: ", α[:,1])=#

  #store the resutls in an object and return it
  return CTFMWindow(FADF, FSyms, RSyms, β, α, cλ[2:end,:], T, N, K)
end

#provides the mean risk premia
meanλ(fm::CTFMWindow)::Vector{Float64} = vec(mean(fm.λ,2))

#get the standard errors of the estimates
function σλ(fm::CTFMWindow)::Vector{Float64}
  λ₀::Vector{Float64} = meanλ(fm)
  return ((1.0/fm.T^2) .* ((k::Int)->sum((fm.λ[k,:]-λ₀[k]).^2)).(1:fm.K)).^0.5
end

#provides the mean risk premia across time
αi(fm::CTFMWindow)::Vector{Float64} = ((i::Int)->mean(dropna(fm.α[i,:]))).(1:fm.N)

#performs the cochrane pricing error test and returns (χ²,p)
function χ²α(fm::CTFMWindow)::Tuple{Float64,Float64}
  α₀::Vector{Float64} = αi(fm)
  #BLAS.set_num_threads(8)

  #first build the covariance matrix


  #pre-allcoate the space
  αtMα₀::DataVector{Float64} = DataVector{Float64}(Vector{Float64}(fm.N))
  covα::Matrix{Float64} = zeros(fm.N, fm.N)
  covαt::DataArray{Float64,2} = DataArray{Float64,2}(Matrix{Float64}(fm.N,fm.N)) #temporary for holding NAs
  TMat::Matrix{Int} = zeros(Int,fm.N, fm.N) #counts the number of non-missing values
  NAMat::BitArray{2} = trues(fm.N, fm.N)

  for t::Int ∈ 1:fm.T
    αtMα₀ .= fm.α[:,t] .- α₀
    covαt .= αtMα₀ * (αtMα₀)'
    NAMat .= (isna).(covαt)
    TMat .+= 1 .- NAMat
    covαt[NAMat] .= 0.0
    covα .+= covαt
  end

  covα .= covα ./ (TMat.^2)

  #calculate the inverse via the cholesky process
  covαInv::Matrix{Float64} = covα \ eye(fm.N)
  χ²::Float64 = α₀' * covαInv * α₀
  p::Float64 = cdf(Chisq(fm.N-1), χ²)
  return χ², p
end

#get the standard errors of the pricing errors
function σα(fm::CTFMWindow)::Vector{Float64}
  αi₀::Vector{Float64} = αi(fm)
  σ::Vector{Float64} = Vector{Float64}(fm.N)
  for i::Int ∈ 1:fm.N
    T::Float64 = length(fm.α[i,:])
    σ[i] = (sum((dropna(fm.α[i,:]).-αi₀[i]).^2) / T^2)^0.5
  end
  return σ #((1.0/fm.T^2) .* ((i::Int)->sum((dropna(fm.α[i,:]).-αi₀[i]).^2)).(1:fm.N)).^0.5
end

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

#an encapsulating function for making the graphs
function makeFMScriptGraphs(fm::CTFMWindow, label::Symbol; sorts::Tuple{Symbol, Symbol}=(:ME, :BM),
    sortDims::Tuple{Int,Int} = (10,10))::Void
  pathStr::String = "$OUTPUT_PATH\\$(FACTOR_NAME)_$(ASSET_NAME)_$label"
  FSyms::Vector{Symbol} = FACTOR_SYMS
  RFRSym::Symbol = RFR_SYM

  #this structure allows us to get a vector of symbols corresponding to a single
  #dimension of the sort. So LinkedRSyms[1][4][3] corresponds to ME3
  formToGrid(v::Vector)::Matrix = reshape(v, sortDims[2], sortDims[1])

  RSymMat::Matrix{Symbol} = formToGrid(fm.RSyms)
  RSymDict::Dict = Dict(fm.RSyms[i] => i for i::Int ∈ 1:length(fm.RSyms))

  #turn the betas into a matrix
  βMat::NTuple{fm.K, Matrix{Float64}} = (((k::Int)->Matrix{Float64}(sortDims[2], sortDims[1])).(1:(fm.K))...)
  for k = 1:fm.K
    for j = 1:sortDims[1]
      for i = 1:sortDims[2]
        βMat[k][i,j] = fm.β[RSymDict[RSymMat[i,j]],k]
      end
    end
  end

  Gadfly.push_theme(:default)
  plots::Vector{PlotContainer} = Vector{PlotContainer}()
  plotNames::Vector{String} = Vector{String}()

  push!(plotNames, "SMB1st")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[2], 1)),
    Geom.line,Geom.point,
    Guide.ylabel("β (SMB)"), Guide.xlabel("SMB Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage SMB Beta by SMB Decile")))

  push!(plotNames, "HML1st")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[3], 2)),
    Geom.line,Geom.point,
    Guide.ylabel("β (HML)"), Guide.xlabel("HML Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage HML Beta by HML Decile")))

  push!(plotNames, "SMB2nd")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[2], 2)),
    Geom.line,Geom.point,
    Guide.ylabel("β (SMB)"), Guide.xlabel("HML Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage SMB Beta by HML Decile")))

  push!(plotNames, "HML2nd")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[3], 1)),
    Geom.line,Geom.point,
    Guide.ylabel("β (HML)"), Guide.xlabel("SMB Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage HML Beta by SMB Decile")))


  push!(plotNames, "MRF1")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[1], 1)),
    Geom.line,Geom.point,
    Guide.ylabel("β (HML)"), Guide.xlabel("SMB Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage M-Rf Beta by SMB Decile")))

  push!(plotNames, "MRF2")
  push!(plots, plot(x=collect(1:10), y=vec(mean(βMat[1], 1)),
    Geom.line,Geom.point,
    Guide.ylabel("β (HML)"), Guide.xlabel("HML Portfolio Decile"),
    Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("First Stage M-Rf Beta by HML Decile")))

  push!(plotNames, "ErrorHist")
  push!(plots, plot(x=σα(fm)./100*255^.5,
    Geom.histogram(bincount=20,density=true),
    Guide.ylabel("Freq"), Guide.xlabel("σ (ann.)"),
    #Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
    Guide.title("Frequency of Time Series Avg Pricing Errors")))


    push!(plotNames, "ErrorHistT")
    push!(plots, plot(x=((αi(fm)/100+1.0).^255-1.0) ./ (σα(fm)./100*255^.5),
      Geom.histogram(bincount=20,density=true),
      Guide.ylabel("Freq"), Guide.xlabel("T"),
      #Coord.Cartesian(ymin=-0.5,ymax=1.5, xmin=0, xmax=10),
      Guide.title("Frequency of Avg T Statistics")))

  for i ∈ 1:length(plots)
    draw(SVG("$(pathStr)_$(plotNames[i]).svg", 4inch, 4inch), plots[i])
  end




  #βMat::NTuple{fm.K, Matrix{Float64}} = (((k::Int)->formToGrid(fm.β[:,k])).(1:(fm.K))...)

  return nothing
end

function A7Script(; refreshData::Bool = true, rerunReg::Bool = true,
    label::Symbol=Symbol(""), calcχ²::Bool=true)::Void
  FADF::DataFrame = getDF("$DATA_PATH\\$(FACTOR_NAME)_$(ASSET_NAME).jls", refreshData, preProcessFA)

  #rename non-standard column names
  rename!(FADF, [:SMALL_LoBM, :SMALL_HiBM, :BIG_LoBM, :BIG_HiBM],
    [:ME1_BM1, :ME1_BM10, :ME10_BM1, :ME10_BM10])
    #=rename!(FADF, [:LoOP_LoINV, :LoOP_HiINV, :HiOP_LoINV, :HiOP_HiINV],
      [:OP1_INV1, :OP1_INV10, :OP10_INV1, :OP10_INV10])=#
  #build the vector of return symbols
  #showcols(FADF)
  #setup the factor symbols
  FSyms::Vector{Symbol} = FACTOR_SYMS
  RFRSym::Symbol = RFR_SYM

  RSyms::Vector{Symbol} = names(FADF[:,setdiff(names(FADF), [:DATE; RFR_SYM; FACTOR_SYMS])])

  #net out the risk-free rate
  for s::Symbol ∈ RSyms
    FADF[:,s] = FADF[:,s] .- FADF[:,RFRSym]
  end

  #run the FM regressions
  if rerunReg
    gstream::GZipStream = gzopen("$DATA_PATH\\REG_OUT$label.jls.gz", "w")
    serialize(gstream, CTFMWindow(FADF, FSyms, RSyms))
    close(gstream)
  end

  gstream = gzopen("$DATA_PATH\\REG_OUT$label.jls.gz")
  fm::CTFMWindow = deserialize(gstream)
  close(gstream)

  makeFMScriptGraphs(fm, label)

  λ₀::Vector{Float64} = meanλ(fm)
  σλ₀::Vector{Float64} = σλ(fm)

  println(fm.FSyms, "\n", "λ(Ann.): $((1.0.+λ₀/100).^255-1)", "\n", "λ: $λ₀",
    "\n", "σ(λ): $σλ₀", "\n", "T(λ): $(λ₀./σλ₀)")

  if calcχ²
    println("\nCalculating χ² statistics")
    χ²::Float64, p::Float64 = χ²α(fm)
    println("χ²($(fm.N-1))): $χ², p=$p")
  end

  return nothing
end

#uncomment to run code
#=@time begin
  A7Script(refreshData=false, rerunReg=false, calcχ²=true, label=Symbol("$(Dates.today())"))
  #cleanworkers()
end

gc()=#

end
