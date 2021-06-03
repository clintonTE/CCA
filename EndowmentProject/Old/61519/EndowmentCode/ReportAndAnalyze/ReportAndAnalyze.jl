

#oneoff tasks
function oneOffCategories(data::NCCSData)
  metadf::DataFrame = deepcopy(data.df[data.df[:fisyr] .== 2013, [:category, :categoryname, :naics2name, :naics2]])
  dropmissing!(metadf)
  outdf::DataFrame = by(metadf, [:category, :categoryname, :naics2name, :naics2],
    [:category, :categoryname, :naics2name, :naics2] =>
      x -> (n=length(x.category)))

  CSV.write("$WORKING_PATH\\categoryMapping.csv",outdf)

  return nothing
end


function oneOffTasks(data::NCCSData)
  println(data.df[1:10, :datefiscal])

  return nothing
end

#conducts main data processing and analysis
#results are cached in the serial file. Useful to avoid redoing computation
function constructPrimaryOutput(; refreshFilter::Bool = true,
    primaryOutputName::String = PRIMARY_OUTPUT_NAME,
    writePrimaryCSV::Bool = false,
    dataPath::String = DATA_PATH,
    refreshoutsidedata::Bool = false,
    filternoncentral::Bool = true,
    refreshopen990data=true)::NCCSData

  local data::NCCSData
  local outPath = "$dataPath\\$(primaryOutputName)_long.csv"
  if refreshFilter
    data = deserializeNCCSData()
    data = processAll!(data,
      refreshoutsidedata=refreshoutsidedata,
      filternoncentral=filternoncentral,
      refreshopen990data=refreshopen990data)

    refreshopen990data && prepsensitivityregressions!(data)

    serializeNCCSData(data, dataName = primaryOutputName)
  else
    data = deserializeNCCSData(dataName = primaryOutputName)
  end

  #scenario where we only want to rerun the open990 data
  if (!refreshFilter) && refreshopen990data
    prepsensitivityregressions!(data)
    serializeNCCSData(data, dataName = primaryOutputName)
  end

  #for debugging
  #uCSV.write(outPath, data.df[400_000:450_000,:])

  writePrimaryCSV && uCSV.write(outPath, data.df) #write to a csv if desired

  return data

end


function writedata(df::DataFrame; rows::AbstractRange = 1:size(df,1),
  cols::Vector{Symbol}=names(df),
  testname::String = "dataout", workingpath::String = WORKING_PATH)

  path::String = "$workingpath\\dataout.csv"
  CSV.write(path, df[rows,cols])
end


#launching point for the descriptive output
function reportandanalyze(;refreshFilter::Bool = true,
    refreshDescriptiveFigures::Bool = true,
    refreshDescriptiveTables::Bool = true, writePrimaryCSV::Bool = false,
    doOneOffTasks::Bool= false,
    refreshβ::Bool = true,
    refreshoutsidedata::Bool = false,
    refreshpersistence::Bool = false,
    refreshtopinstitutions::Bool = false,
    refreshopen990data::Bool = true,
    refreshsensitivity::Bool = true,
    refreshtails=true
    )


    @time data::NCCSData = constructPrimaryOutput(refreshFilter=refreshFilter,
      writePrimaryCSV=writePrimaryCSV,
      refreshoutsidedata=refreshoutsidedata,
      refreshopen990data=refreshopen990data)

    #NOTE: Begin debug code
    sort!(data.df, [:ein, :fisyr])

    #=writedata(data.df[data.df[:fisyr] .≥ 2009,[:ein, :fisyr, :lreturn,
      :Llcontributions, :lcontributions, :contributions,
      :Lpublicsupport, :publicsupport, :publicsupport_avg,
      :Llprogramexpense, :lprogramexpense, :programexpense,
      :programexpense_avg, :lprogramexpense_avg,
      :LpublicsupportXLlcont,  :LlprogramexpenseXLlcont, :publicsupport_avgXLlcont,
      :lprogramexpense_avgXLlcont
      ]])=#

    #end debug code


    refreshDescriptiveFigures && makeDescriptiveFigures(data)
    refreshDescriptiveTables && makeDescriptiveTables(data)
    refreshtopinstitutions && getTopInstitutions(data, sorted=true) #sorted=true ok after descriptive tables
    refreshβ && processβ(data)
    refreshpersistence && analyzepersistence(data)

    #no need to re-run this if it was run in the filter
    refreshopen990data && (!refreshFilter) && prepsensitivityregressions!(data)
    refreshsensitivity && olssensitivitytables(data)

    doOneOffTasks && oneOffTasks(data)
    refreshtails && summarizetails(data)


    return nothing
end
