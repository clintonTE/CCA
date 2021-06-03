
#run the time series regressions
function analyzeqfactors(;Ftest::Vector{Symbol} = FF_PORT_VALS,
  Ffactors::Vector{Symbol} = FACTOR_FIELDS,
  factname::String = FACT_NAME, outpath::String = OUT_PATH)

  #read in the data nd make the requisite structure
  df = CSV.read("$outpath\\$factname.csv") |> DataFrame
  fs::FactorSpec = FactorSpec(df, :date, Ffactors, Ftest)

  #not strictly necessary, but helps stop doubling up on the multi-threading
  PARALLEL[] && (BLAS.set_num_threads(1))

  #run the regressions
  results::Vector{FMLM} = timeseriesregressions(fs)

  #make the output table

end
