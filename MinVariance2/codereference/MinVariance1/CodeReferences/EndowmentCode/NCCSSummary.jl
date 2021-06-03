#This file contains a wrapper object and some helper functions
#relating to endowment files
#creates a table of summary values
function summarizeNCCSFile(data::NCCSFile)::DataFrame


  dataDF::DataFrame = data.df #for convenience since we are just reading, not writing
  summaryDF::DataFrame = metaTable(data.meta) #get the metadata in tabular format
  fields::Vector{Symbol} = summaryDF[:fieldNames]
  #the below is an UGLY hack for CSV compat. Fix it later!
  summaryDF[:longDescriptions] .= (s::String->strip("""γ,γ$s""")).(summaryDF[:longDescriptions])

  NData::Int = size(dataDF,1)
  NFields::Int = length(fields)
  #allocate the new columns
  fileRows::Vector{MInt} = fill(NData, NFields)
  N::Vector{MInt} = missings(NFields)
  pMissing::Vector{MFloat64} = missings(NFields)
  Unique::Vector{MInt} = missings(NFields)
  Mean::Vector{MFloat64} = missings(NFields)
  MeanSq::Vector{MFloat64} = missings(NFields)
  Std::Vector{MFloat64} = missings(NFields)
  Skew::Vector{MFloat64} = missings(NFields)
  Kurtosis::Vector{MFloat64} = missings(NFields)
  Min::Vector{MFloat64} = missings(NFields)
  p05::Vector{MFloat64} = missings(NFields)
  p25::Vector{MFloat64} = missings(NFields)
  p50::Vector{MFloat64} = missings(NFields)
  p75::Vector{MFloat64} = missings(NFields)
  p95::Vector{MFloat64} = missings(NFields)
  Max::Vector{MFloat64} = missings(NFields)

  for i::Int ∈ 1:NFields
    #first get some type info
    T::Type = data.meta.dictionary[fields[i]]
    TSub::Type = typeintersect(T, Union{Float64, Bool, Int, String, CategoricalValue}) #get the nonmissable type
    TSub = CategoricalValue <: TSub || TSub <: CategoricalValue ? Symbol : TSub

    #compute if the type is numeric, then get the column type
    numeric::Bool = (T <: MFloat64) || (T <: MInt)
    VecT = (T <: CategoricalValue) || (CategoricalValue <: T) ? Vector{Union{Symbol,Missing}} : Vector{T}

    #form the rows and do the applicable computations
    col::VecT = VecT(undef, NData)
    col .= VecT(dataDF[fields[i]])
    sansMissing::Vector{TSub} = collect(skipmissing(col))
    if length(sansMissing) == 0
      N[i] = 0.0
      pMissing[i] = 1.0
    elseif numeric
      N[i] = length(sansMissing)
      pMissing[i] = 1.0 - length(sansMissing) / length(col)
      Unique[i] = length(unique(col))
      Mean[i] = mean(sansMissing)
      MeanSq[i] = mean(sansMissing .* sansMissing)
      Std[i] = std(sansMissing)
      Skew[i] = skewness(sansMissing)
      Kurtosis[i] = kurtosis(sansMissing)
      Min[i] = minimum(sansMissing)
      p05[i] = percentile(sansMissing,0.05)
      p25[i] = percentile(sansMissing,0.25)
      p50[i] = percentile(sansMissing,0.50)
      p75[i] = percentile(sansMissing,0.75)
      p95[i] = percentile(sansMissing,0.95)
      Max[i] = maximum(sansMissing)
    else
      N[i] = length(sansMissing)
      pMissing[i] = 1.0 - length(sansMissing) / length(col)
      Unique[i] = length(unique(col))
    end
  end

  #finish building the dataframe
  summaryDF = hcat(summaryDF, DataFrame(fileRows=fileRows, N=N, pMissing=pMissing, Unique=Unique, Mean=Mean,
    MeanSq=MeanSq, Std=Std, Skew=Skew, Kurtosis=Kurtosis, Min=Min, p05=p05, p25=p25,
    p50=p50, p75=p75, p95=p95, Max=Max))

  return summaryDF
end

#convenience method since loading many files is memory intensive
summarizeNCCSFile(meta::MetaDictionary) =
  summarizeNCCSFile(loadNCCSFile(meta, reserializeData=false, returnNCCSFile=true))::DataFrame

#vectorized method for creating summary tables
function summarizeNCCSFiles(metadata::Vector{MetaDictionary}; mode=:serial)::DataFrame
  #get the summary tables

  if mode == :serial #terminal branch, summarize the data
    NData::Int = length(metadata)
    #println("typeof: $(summarizeNCCSFile(metadata[1]))")
    summaryTable::DataFrame = vcat((summarizeNCCSFile).(metadata)...)
    return summaryTable

  elseif mode == :parallel #relaunch the dedicated parallel function
    numFiles::Int = length(metadata)
    pids = workers()
    np::Int = length(pids)

    #divide up the work
    assignments::Vector{Vector{Int}} = ((i::Int)->Vector{Int}()).(1:np)
    for i::Int ∈ 1:numFiles
      push!(assignments[i%np+1], i) #assigns a file to a pid
    end

    #set up a vector for pushing and receiving assignments
    futures::Vector{Future} = Vector{Future}()
    sizehint!(futures, np)

    for i::Int ∈ 1:np
      if length(assignments[i]) > 0
        push!(futures,
          (@spawnat pids[i] summarizeNCCSFiles(metadata[assignments[i]], mode=:serial)))
      end
    end
    outDFs::Vector{DataFrame} = (fetch).(futures)
    outDF::DataFrame = vcat(outDFs...) #concatenate the results
    return outDF
  else
    error("Incorrect mode for function summarizeNCCSFiles")
  end

end


#NOTE this ugly hack works around an issue in CSV.write
#if deleting this, be sure to get the flag where it is inserted elsehwere in the code
function csvWorkaroundWrite(filePath::String, df::DataFrame)

    b = IOBuffer()
    CSV.write(b, df)
    #b = Stream(format"CSV", IOBuffer())
    #save(b, df)
    outString::String = String(take!(copy(b)))
    outString = replace(outString, "\r\n"=>"InnerNewLineρη")
    outString = replace(outString, "\n"=>"\r\n")
    outString = replace(outString, "InnerNewLineρη"=>"\n")
    outString = replace(outString, "γ,γ"=>"")
    write(filePath, outString)
    #end ugly hack

    return nothing
end

#this function constructs and write a summary data file
function writeNCCSFileSummary(metadata::Vector{MetaDictionary};
  documentationPath::String =DOCUMENTATION_PATH,  summaryName = SUMMARY_NAME)::DataFrame

  summaryDF::DataFrame = summarizeNCCSFiles(metadata, mode=:serial)
  csvWorkaroundWrite("$documentationPath\\$summaryName.csv", summaryDF) #workaround on CSV

  #also make a serialized version
  oStream::IOStream = open("$documentationPath\\$summaryName.jls", "w")
  serialize(oStream, summaryDF)
  close(oStream)

  println("Summary table written.")

  return summaryDF
end

#consolidates the summary table
function consolidateSummary(metadata::Vector{MetaDictionary};
  refreshSummary::Bool = true, documentationPath::String = DOCUMENTATION_PATH,
  summaryName::String = SUMMARY_NAME, consolidatedName::String = CONSOLIDATED_NAME)

  #load and/or generate the table
  if !isfile("$documentationPath\\$summaryName.jls") || refreshSummary
    summaryDF::DataFrame = writeNCCSFileSummary(metadata)
  else
    iStream::IOStream = open("$documentationPath\\$summaryName.jls")
    summaryDF = deserialize(iStream)
    close(iStream)
  end


  summaryDF = summaryDF[(!).(summaryDF[:fullFiles]),:]
  NSummary::Int = size(summaryDF,1)

  #copy and do some initial filtering
  uniqueDF::DataFrame = deepcopy(summaryDF)
  deletecols!(uniqueDF, [:Unique, :Skew, :Kurtosis, :p05, :p25, :p50, :p75, :p95])
  uniqueDF[[:fileRows, :N, :pMissing, :Mean, :MeanSq, :Std, :Min, :Max]] =  missing
  NUnique::Int = size(uniqueDF,1)

  #create a smaller version with only unique columns
  uniqueDF[:keep] = trues(NUnique)
  uniqueDFBy::GroupedDataFrame = groupby(uniqueDF, [:fieldNames, :nonprofitTypes])
  for df::AbstractDataFrame in uniqueDFBy
    maxYear::Int = maximum(df[:fileYears])

    df[:keep] = (maxYear .== df[:fileYears])
    if sum(df[:keep]) > 1
      @warn "Multiple rows kept for field $(df[1, :fieldNames]), type $(df[1, :nonprofitTypes]),
        fullFile $(df[1, :fullFiles]), num kept $(sum(df[:keep])) in file $(df[1, :fileNames])"
    end
  end

  #make a dictionary to get the address of unique values

  uniqueDF = uniqueDF[uniqueDF[:keep], :]
  NUnique = size(uniqueDF, 1)
  uniqueIndex::Dict = Dict((uniqueDF[i, :fieldNames], uniqueDF[i, :nonprofitTypes])=>i for i ∈ 1:NUnique)

  uniqueDF[:fileYearsCovered] = Vector{MInt}(missing,NUnique)
  uniqueDF[:fileYearsMin] = Vector{MInt}(missing,NUnique)
  uniqueDF[:fileYearsMax] = Vector{MInt}(missing,NUnique)
  uniqueDF[:fileYearsGap] = Vector{MInt}(missing,NUnique)
  uniqueDF[:UniqueMin] = Vector{MInt}(missing,NUnique)
  uniqueDF[:UniqueMax] = Vector{MInt}(missing,NUnique)
  uniqueDF[:inconsistentType] = Vector{MBool}(missing,NUnique)
  uniqueDF[:inconsistentNumeric] = Vector{MBool}(missing,NUnique)

  #now get the aggregate summary stats
  #first make sure we only select columns for which we have information
  #validSummaryCols::Vector{Bool} =  ((f::Symbol, t::Symbol) -> haskey(uniqueIndex, (f,t))).(
  #  summaryDF[:fieldNames], summaryDF[:nonprofitTypes])
  summaryDFBy = groupby(summaryDF, [:fieldNames, :nonprofitTypes])
  for df::SubDataFrame in summaryDFBy
    #by(summaryDF, [:fieldNames, :nonprofitTypes]) do df::AbstractDataFrame

    row::Int = uniqueIndex[(df[1, :fieldNames], df[1, :nonprofitTypes])]

    #now aggregate the summary info
    uniqueDF[row,:fileYearsCovered] = size(df,1)
    uniqueDF[row,:fileYearsMin] = minimum(skipmissing(df[:fileYears]))
    uniqueDF[row,:fileYearsMax] = maximum(skipmissing(df[:fileYears]))
    uniqueDF[row, :fileYearsGap] =  #count the missing years
      (uniqueDF[row,:fileYearsMax]-uniqueDF[row,:fileYearsMin] + 1) - uniqueDF[row,:fileYearsCovered]

    #check if the metadata is at least somewhat consistent

    #make sure we don't reduce over 0D
    if length(df[(!ismissing).(df[:Unique]),:Unique]) > 0
      uniqueDF[row, :UniqueMin] = maximum(skipmissing(df[:Unique]))
      uniqueDF[row, :UniqueMax] = minimum(skipmissing(df[:Unique]))
    end
    uniqueDF[row,:inconsistentType] = length(unique(df[:nonprofitTypes])) > 1
    uniqueDF[row,:inconsistentNumeric] = length(unique(df[:numeric])) > 1

    #Get information on the quantity of missing values
    uniqueDF[row, :fileRows] = sum(skipmissing(df[:fileRows]))
    uniqueDF[row, :N] = sum(skipmissing(df[:N]))
    uniqueDF[row, :pMissing] =
      (uniqueDF[row, :N] > 0) ? (uniqueDF[row, :fileRows] - uniqueDF[row, :N]
      ) / uniqueDF[row, :fileRows] : 0.0

    #do the numeric processing
    if uniqueDF[row, :numeric]
      uniqueDF[row, :Mean] = sum(skipmissing(df[:Mean] .* df[:N])) / uniqueDF[row, :N]
      uniqueDF[row, :MeanSq] = sum(skipmissing(df[:MeanSq] .* df[:N])) / uniqueDF[row, :N]
      uniqueDF[row, :Std] = (uniqueDF[row, :MeanSq] - uniqueDF[row, :Mean]^2)^0.5
      uniqueDF[row, :Min] = minimum(skipmissing(df[:Min]))
      uniqueDF[row, :Max] = maximum(skipmissing(df[:Max]))
    end
  end

  csvWorkaroundWrite("$documentationPath\\$consolidatedName.csv", uniqueDF)


end
