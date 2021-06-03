module  CumEx

#this is useful if we have mdoules
if pwd() ∉ LOAD_PATH
  push!(LOAD_PATH,pwd())
end

using DataFrames, Distributions, StatsBase, GZip, JLD, Gadfly

#=#NOTE##############
Some companies have multiple PERMNOs. The returns are different, so I am leaving them as is.
=####################
const DATA_PATH = pwd() * "\\data" #path where we store the data
const DATA_FILE_NAME = "ReturnsDivOnly1926-2017" #data file name
const DATE_FORMAT_STR = "yyyymmdd"
const WINSOR_MAX = 0.2
const NET_OF_MARKET = true
const MIN_DIVIDEND = 10.0^-4.0 #1/100 of a penny
const MIN_PRICE = 1.0
const EXC_SIC_CODES = nothing

#useful type
const CTSym = Union{Symbol,Void}

#This function converts a string datavector to a real number
#IN: A datavector string
#OUT: A data vector with all string numbers converted to numbers,
# and everything else converted to an NA
function parseDataVector(strings::DataVector{String})::DataVector{Float64}
  N::Int = length(strings)

  vals::DataVector{Float64} = DataVector{Float64}(Vector{Float64}(N),trues(N))

  #threading makes this much faster
  @fastmath @inbounds Threads.@threads for i::Int ∈ 1:N
    parsed::Union{Symbol, Float64} = parse(strings[i])
    vals[i] = typeof(parsed)<:Float64?parsed:NA
  end

  return vals
end

#in this case, do nothing
parseDataVector(vals::DataVector{Float64})::DataVector{Float64}  = vals


#This function winsorizes a data column
#IN: A datavector string
#OUT: A data vector with all string numbers converted to numbers,
# and everything else converted to an NA
winsorize(col::DataVector{Float64}, minVal::Float64, maxVal::Float64)::DataVector{Float64} =
  min.(max.(col, minVal),maxVal)

#same as above but works at the element level
winsorize(val::T, minVal::Float64, maxVal::Float64) where T<:Union{NAtype,Float64} =
  min.(max.(val, minVal),maxVal)

#this function pre-processes the data files
#IN: A switch for reading or re-reading the csv
#OUT: nothing. A clean and compressed julia serial file has all of the data
function processStockData(;newBinary::Bool = true)::Void
  dtFormat::DateFormat = DateFormat(DATE_FORMAT_STR) ## a date format object

  if !isfile("$DATA_PATH\\$DATA_FILE_NAME.jls.gz") || newBinary


  #extract from the zip file, package into dataframe and serialize into zip
  #Saves a lot of time not to have to do this each time
    igStream::GZipStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME.csv.gz")
    ogStream::GZipStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME.jls.gz","w")
    serialize(ogStream, readtable(igStream))

    close(igStream)
    close(ogStream)
  end

  igStream = gzopen("$DATA_PATH\\$DATA_FILE_NAME.jls.gz")
  stockDF::DataFrame  = deserialize(igStream)
  close(igStream)

  #unfortunately, some of these values have error codes
  stockDF = stockDF[completecases(stockDF[:,[:RET,:RETX]]),:]
  stockDF[:,:RET] = parseDataVector(stockDF[:,:RET])
  stockDF[:,:RETX] = parseDataVector(stockDF[:,:RETX])
  stockDF = stockDF[completecases(stockDF[:,[:RET,:RETX,:PRC]]),:] #can get more NAs from before

  stockDF = stockDF[stockDF[:,:DIVAMT] .≥ MIN_DIVIDEND,:] #Filter out extremely small dividends
  stockDF = stockDF[stockDF[:,:PRC] .≥ MIN_PRICE,:] #Filter out penny stocks

  if EXC_SIC_CODES ≠ nothing
    stockDF = stockDF[(stockDF[:,:SICCD] .> maximum(EXC_SIC_CODES)) | (
      stockDF[:,:SICCD] .< minimum(EXC_SIC_CODES)),:] #Filter out financials
  end

  stockDF = deepcopy(stockDF[:,[:PERMNO,:date,:PERMCO, :DIVAMT,:PRC,:SHROUT, :RET,:RETX,:vwretd]])
  stockDF[:,:date] = ((x::Int)->Date("$x",dtFormat)).(stockDF[:,:date])
  N::Int = size(stockDF,1)

  #generate the derived columns
  stockDF[:,:divYield] = stockDF[:,:RET] .- stockDF[:,:RETX]
  stockDF[:,:RETXW] = winsorize(stockDF[:,:RETX],-1.0*WINSOR_MAX,WINSOR_MAX)
  stockDF[:,:divYieldW] = winsorize(stockDF[:,:divYield],-1.0*WINSOR_MAX,WINSOR_MAX)


  stockDF[:,:marketCap] = stockDF[:,:SHROUT] .* stockDF[:,:PRC]

  stockDF[:,:year] = (Dates.year).(stockDF[:,:date])
  stockDF[:,:month] = (Dates.month).(stockDF[:,:date])
  stockDF[:,:quarter] = (Dates.quarterofyear).(stockDF[:,:date])
  stockDF[:,:quarterOfYear] = (Dates.firstdayofquarter).(stockDF[:,:date])
  stockDF[:,:monthOfYear] = (Dates.firstdayofmonth).(stockDF[:,:date])


  ogStream = gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_clean.jls.gz","w")
  serialize(ogStream, stockDF)
  close(ogStream)

  #showcols(stockDF)
  return nothing
end

#this function is meant to be used internally in aggDivCumEx to build the time series
#IN: a dataframe or subdataframe of input data and a weighting vector
#OUT: the weighted sum of the input vector. This is a weighted mean if Σw=1
function getNegROverDNet(df::T, weightVec::Vector{Float64})::Float64 where {T<:AbstractDataFrame}

  N::Int = size(df, 1)

  ans::Float64 = sum(-(df[:,:RETXW].-df[:,:vwretd])./df[:,:divYieldW] .* weightVec)
  ans /= sum(weightVec)
  return ans
end


function getDivLevels(df::T, weightVec::Vector{Float64})::Float64 where {T<:AbstractDataFrame}

  ans::Float64 = sum(df[:,:DIVAMT] .* weightVec)
  #ans /= sum(weightVec)

  return ans
end


#aggregates the Elton-Gruber 1970 statistic based on the supplied groupBy column and weight
#IN: The dataframe, the column to split by, weights, and an aggregation function
#OUT: A dataframe aggregated according to the input parameters
#NOTE: Default options supplied for testing purposes
function aggDivCumEx(stockDF::DataFrame, weightCol::Symbol; byCol::Symbol=:year,
  statFunc::Function = getNegROverDNet)::DataFrame

  #apply the supplied function to each subsection
  outDF::DataFrame = by(stockDF, byCol) do subDF::SubDataFrame
    DataFrame(EGStat = statFunc(subDF,Vector{Float64}(subDF[:,weightCol])))
  end

  return outDF
end


#this version allows for a vector of weighting columns.
#The results are merged into a single dataframe
#IN: the source data column, the weighting columns, the by column for aggregation, and
# the aggregation function
#OUT: The merged data frame.
#NOTE: Each weighting column can be void, and if so this creates a default equal-weighted column
function aggDivCumEx(stockDF::DataFrame, weightCols::Vector{T}=[nothing]; byCol::Symbol=:year,
  statFunc::Function = getNegROverDNet)::DataFrame where T<:CTSym

  #set unit weights if no other weights
  if nothing ∈ weightCols
    weightCol = :equalWeights
    stockDF[:,:equalWeights] = ones(size(stockDF,1))
  end

  #replace the nothing column
  weightCols =
    [weightCols[i]==nothing?(:equalWeights):(weightCols[i]) for i::Int ∈ 1:length(weightCols)]

  #Create the dataframe from the first aggregation
  outDF = aggDivCumEx(stockDF, :equalWeights, byCol=byCol, statFunc = getDivLevels)
  rename!(outDF, :EGStat, :Levels) #need to set the name for each column

  #if more than one weight is provided, execute again and append to the dataframe

  for weightCol::Symbol ∈ weightCols[1:end]
    outDF =
      [outDF aggDivCumEx(stockDF, weightCol, byCol=byCol, statFunc = getNegROverDNet)[:,:EGStat]]
    #println(outDF)
    rename!(outDF, :x1, Symbol(:EGStat,"_", weightCol))
  end


  return outDF
end


#This is a script for running the aggregation. Uses the output of processStockData
#IN: Nothing, but reads in a dataframe as a serial file
#OUT: Nothing, but writes the aggregated tables to the disk
function getCumExTablesScript()::Void
  igStream = gzopen("$DATA_PATH\\$(DATA_FILE_NAME)_clean.jls.gz")
  stockDF::DataFrame  = deserialize(igStream)
  close(igStream)

  #set the dimensions of the cum-ex aggregation
  weightCols::Vector{CTSym} = [nothing, :divYieldW, :marketCap]
  byCols::Vector{Symbol} = [:year, :quarterOfYear, :monthOfYear]
  aggTables::Vector{DataFrame} = Vector{DataFrame}()

  #get the dataframes
  for byCol::Symbol ∈ byCols
    push!(aggTables, aggDivCumEx(stockDF, weightCols, byCol=byCol, statFunc = getNegROverDNet))
  end

  #write out the dataframes
  for i::Int ∈ 1:length(aggTables)
    pth::String = "$DATA_PATH\\$(DATA_FILE_NAME)_$(byCols[i]).csv"
    println("Writing: $pth")
    writetable(pth, aggTables[i])
  end

  println("All tables writtens successfully.")

  return nothing
end



#run the appropriate functions
@time begin
processStockData(newBinary = false)
getCumExTablesScript()
end

end
