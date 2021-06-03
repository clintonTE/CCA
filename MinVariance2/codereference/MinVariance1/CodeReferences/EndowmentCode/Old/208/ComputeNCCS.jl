


const ALL_COLS = [:name, :charitytype, :ein, :fisyr, :fileyear, :invrealized,
  :invrealizeddirect, :totrev, :totexp,
  :invrealizedindirect, :name, :netassets,
  :bookassets, :bookliabilities,
  :pfexpbook, :pfexpcharity, :totinv, :category,
  :returnamt, :returncheckamt,:netincome, :naics]

const CATEGORY_CONSOLIDATED_THRESHOLD = 100_000

#=const :Arts, :Education, :Environmental, :AnimalWelfare, :Health, :MentalHealth,
:MedicalDiseases, :MedicalResearch, :CrimeAndLegal, :Employment, :FoodAndAgriculture,
:Housing, :PublicSafety, :Recreation, :YouthDevelopment, :HumanServices, :ForeignAffairs,
:SocialAction, :CommunityImprovement, :GrantmakingFoundations, :ScienceandTechResearch,
:SocialScienceResearch, :SocietyBenefitAndMultipurpose, :ReligionRelated,
:MutualBenefitOrg, :Unknown, :Other=#

function processReturns!(data::NCCSData)

  local N::Int = size(data.df,1) #length of dataframe
  local revenueCutoff::Float64 = REVENUE_CUTOFF
  local largeDifferenceCutoff::Float64 = LARGE_DIFFERENCE_CUTOFF #max discrepenacy between methods


  println("data.df N1: $N")
  println("data.df pf N1: $(sum(data.df[:charitytype] .== :pf))")
  #filter out organizations which are too small
  data.df[:tokeep] = trues(size(data.df,1))
  dfByein::GroupedDataFrame = groupby(data.df, :ein)
  for subdf::SubDataFrame ∈ dfByein
    subN::Int = size(subdf,1)
    maxRevenue::Float64 = maximum(collect(skipmissing(subdf[:totrev])))
    if (maxRevenue < revenueCutoff)
      subdf[:tokeep] = falses(subN)
    end

    if subN ≠ length(unique(subdf[:fisyr])) # check if we have duplicate filings
      deleteOrg::Bool = false
      subdfByfisyr = groupby(subdf, :fisyr)
      for ssubdf::SubDataFrame ∈ subdfByfisyr
        nrows::Int = size(ssubdf, 1)
        if nrows > 1 # if we are at the org-year with the dup
          fileYears::Vector{Int} = unique(ssubdf[:fileyear])
          if length(fileYears) == nrows #if we can reconcile by fileYear
            maxFileYear::Int = maximum(fileYears) # keep the newest version
            ssubdf[ssubdf[:fileyear] .≠ maxFileYear, :tokeep] = false
          else
            deleteOrg = true #If we cannot correct, throw out the org
          end
        end
      end
      if deleteOrg #delete the org if we need to
        subdf[:tokeep] = false
      end
    end
  end

  #delete the invalid rows
  data.df = data.df[data.df[:tokeep] .== true, :]
  N = size(data.df,1)

  println("data.df N2: $N")

  #println("N2-pf: $(sum(data.df[:charitytype] .== :pf)) \nN2-pf-Return Amt: $(
  #  sum((data.df[:charitytype] .== :pf) .& ((!ismissing).(data.df[:returnamt]))))")

  #data.df[:returnamt] = Vector{Union{Float64,Missing}}(missings(N)) #don't need this since its done during mapping
  data.df[:return] = Vector{Union{Float64,Missing}}(missings(N))

  #data.df[:returncheckamt] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:returncheck] = Vector{Union{Float64,Missing}}(missings(N))

  data.df[:largeDifference] = Vector{Union{Bool,Missing}}(missings(N))

  sort!(data.df, [:ein, :fisyr])
  #ctr::Int = 0
  dfByein = groupby(data.df, :ein)
  for subdf::SubDataFrame ∈ dfByein
    subN = size(subdf,1)

    if (subdf[1,:charitytype] == :co) || (subdf[1,:charitytype] == :pc)

      for i ∈ 2:subN
        if subdf[i-1, :fisyr] == subdf[i, :fisyr] - 1

          #NOTE: returns = unrealized + inv income s.t. unrealized = Δassets - netincome
          returnCandidateAmt::MFloat64 = (subdf[i, :netassets] .- subdf[i-1, :netassets]
            ) .+ subdf[i, :invrealizeddirect] .- subdf[i, :netincome]
          returnCandidate::MFloat64 = returnCandidateAmt / subdf[i-1, :netassets]

          returnCheckCandidateAmt::MFloat64 = (subdf[i, :netassets] .- subdf[i-1, :netassets]
            ) .+ subdf[i, :invrealizedindirect] .- subdf[i, :netincome]
          returnCheckCandidate::MFloat64 = returnCheckCandidateAmt / subdf[i-1, :netassets]

          if ((!ismissing(returnCandidate)) && returnCandidate < 2. && returnCandidate > -0.9 &&
              abs(returnCandidate - returnCheckCandidate) ≤ largeDifferenceCutoff)# && abs(subdf[i, :totinv]) ≥ 1. #make sure we have valid numbers here
            subdf[i, :returnamt] = returnCandidateAmt
            subdf[i, :return] = returnCandidate

            subdf[i, :returncheckamt] = returnCheckCandidateAmt
            subdf[i, :returncheck] = returnCheckCandidate
          end
        end
      end
    elseif subdf[1,:charitytype] == :pf
      for i ∈ 2:subN
        if subdf[i-1, :fisyr] == subdf[i, :fisyr] - 1

          #NOTE: returns = unrealized + inv income s.t. unrealized = Δassets - netincome
          returnCandidateAmt::MFloat64 = subdf[i, :returnamt]
          returnCandidate::MFloat64 = returnCandidateAmt / subdf[i-1, :netassets]

          returnCheckCandidateAmt::MFloat64 = subdf[i, :returncheckamt]
          returnCheckCandidate::MFloat64 = returnCheckCandidateAmt / subdf[i-1, :netassets]

          if ( (!ismissing(returnCandidate)) && (returnCandidate < 2.) && (returnCandidate > -0.9) &&
              (!ismissing(returnCheckCandidate)) &&
              abs(returnCandidate - returnCheckCandidate) ≤ largeDifferenceCutoff)# && abs(subdf[i, :totinv]) ≥ 1. #make sure we have valid numbers here
            subdf[i, :returnamt] = returnCandidateAmt
            subdf[i, :return] = returnCandidate

            subdf[i, :returncheckamt] = returnCheckCandidateAmt
            subdf[i, :returncheck] = returnCheckCandidate
          #=elseif ctr < 50
            ctr+=1
            println("###\nreturn: $(returnCandidate)")
            println("returnchk: $(returnCheckCandidate)")
            println("returnamt: $(returnCandidateAmt)")
            println("returnchk: $((returnCheckCandidateAmt))")=#
          end
        end
      end
    end
  end

  sort!(data.df, [:ein, :fisyr])

  return data
end

function processWealth!(data::NCCSData)

  data.df[:lnetassets] = (log).(data.df[:netassets])
  data.df[:ladjnetassets] = (log).(data.df[:adjnetassets])

  return data
end

#data processing related to the categorization
function processCategories!(data::NCCSData;
  minPointsPerCategory::Int = MIN_POINTS_PER_CATEGORY,
  categoryConsolidatedThreshold::Int = CATEGORY_CONSOLIDATED_THRESHOLD)::Nothing

  data.df[:categoryname] = (s::MSymbol->((ismissing(s) || (!haskey(category, s))) ?
    missing : category[s])).(Vector{MSymbol}(data.df[:category]))

  data.df[:tokeep] .= true
  data.df[:categoryfiltered] = deepcopy(data.df[:category])
  data.df[:categorynamefiltered] = deepcopy(data.df[:categoryname])

  dfBycategoryfisyr::GroupedDataFrame = groupby(data.df, [:category, :fisyr])
  for subdf::SubDataFrame ∈ dfBycategoryfisyr
    if size(subdf,1) < minPointsPerCategory
      subdf[:categoryfiltered] = :ZZ
      subdf[:categorynamefiltered] = category[:ZZ]
    end
  end

  data.df[:categoryconsolidated] = deepcopy(data.df[:categoryfiltered])
  data.df[:categorynameconsolidated] = deepcopy(data.df[:categorynamefiltered])

  dfBycategorynamefiltered::GroupedDataFrame = groupby(data.df, :categorynamefiltered)
  for subdf::SubDataFrame ∈ dfBycategorynamefiltered
    if size(subdf,1) < categoryConsolidatedThreshold
      subdf[:categoryconsolidated] = :ZZ
      subdf[:categorynameconsolidated] = category[:ZZ]
    end
  end

  data.df[:naics2] = (s::MString->(!ismissing(s) ? Symbol(s[1:2]) : missing)).(data.df[:naics])
  data.df[:naics2name] = (s::MSymbol->(naicsCodes[s])).(data.df[:naics2])
  for subdf::SubDataFrame ∈ groupby(data.df, :categoryname)
    cname::MSymbol = subdf[1,:categoryname]
    subdf[:naics2name] = (s::MSymbol->
      (ismissing(s) || s == :Unknown) ? ntee2naics2names[cname] : s).(subdf[:naics2name])
  end

  return nothing
end

#adjusts for inflation
#currently just does wealth, eventually can use this for return fields
function adjustInflation!(data::NCCSData; fredPath::String = FRED_PATH,
  fredInflationFile::String = FRED_INFLATION_FILE, inflationSeries::Symbol = :CPI,
  baseYear::Int = BASE_YEAR, fredDateFormat::String = FRED_DATE_FORMAT,
  adjustCutoff::Bool = true)::Nothing

  #read in the inflation info
  freddf::DataFrame = CSV.read("$fredPath\\$fredInflationFile.csv")

  #make a dictionary
  cpiIndex::Dict =
    Dict(Dates.year(Date(freddf[i,:DATE], fredDateFormat))=>freddf[i, inflationSeries]
      for i ∈ 1:size(freddf,1))

  #get vectors for performance purpuses
  nrows::Int = size(data.df, 1)
  baseValue::Float64 = cpiIndex[baseYear]
  netassets::Vector{Float64} = data.df[:netassets]
  totrev::Vector{Float64} = data.df[:totrev]
  adjnetassets = Vector{MFloat64}(undef, nrows)
  adjtotrev = Vector{MFloat64}(undef, nrows)
  fisyr::Vector{Int} = data.df[:fisyr]
  @inbounds @simd for i ∈ 1:nrows
    if haskey(cpiIndex, fisyr[i])
      adjnetassets[i] = netassets[i] / (cpiIndex[fisyr[i]] / baseValue)
      adjtotrev[i] = totrev[i] / (cpiIndex[fisyr[i]] / baseValue)
    end
  end

  data.df[:adjnetassets] = adjnetassets
  data.df[:adjtotrev] = adjtotrev
  if adjustCutoff #this is to account for the cutoff inflating forward
    inflatedRevenueCutoff::Float64 = REVENUE_LOWER_BOUND * baseValue/cpiIndex[minimum(data.df[:fisyr])]
    println("Inflated revenue cutoff of $inflatedRevenueCutoff")
    data.df = data.df[
      (f::MFloat64-> (!ismissing(f)) && (f ≥ inflatedRevenueCutoff)).(data.df[:adjtotrev]), :]

    inflatedWealthCutoff::Float64 = WEALTH_LOWER_BOUND * baseValue/cpiIndex[minimum(data.df[:fisyr])]
    data.df = data.df[
      (f::MFloat64-> (!ismissing(f)) && (f ≥ inflatedWealthCutoff)).(data.df[:adjnetassets]), :]
  end


  return nothing
end

#does an initial filter to get down the file size
function processAll!(data::NCCSData; allcols::Vector{Symbol} = ALL_COLS,
  revenueCutoff::Float64 = REVENUE_LOWER_BOUND,
  wealthCutoff::Float64 = WEALTH_LOWER_BOUND,
  )

  data = removeFields!(data, setdiff(names(data.df), allcols))
  data.df[:rows] = collect(1:size(data.df,1))
  println("N0: $(size(data.df,1))")
  deleterows!(data.df, data.df[(!).(completecases(data.df[[:fisyr, :netassets, :totrev, :totexp]])),:rows])

  adjustInflation!(data)
  processReturns!(data)
  processWealth!(data)
  processCategories!(data)

  println("N3: $(size(data.df,1)) N3-Return: $(sum((!ismissing).(data.df[:return])))")
  #=println("N3-pf: $(sum(data.df[:charitytype] .== :pf)) \nN3-pf-Return: $(
    sum((data.df[:charitytype] .== :pf) .& ((!ismissing).(data.df[:return]))))\nN3 check return: $(
    sum((data.df[:charitytype] .== :pf) .& ((!ismissing).(data.df[:returncheckamt]))))\nN3 net income: $(
    sum((data.df[:charitytype] .== :pf) .& ((!ismissing).(data.df[:netincome]))))")=#

  data.df = data.df[((!ismissing).(data.df[:adjtotrev])) .&
    (data.df[:adjtotrev] .≥ revenueCutoff), :]
  data.df = data.df[((!ismissing).(data.df[:adjnetassets])) .&
    (data.df[:adjnetassets] .≥ wealthCutoff), :]

  println("N4: $(size(data.df,1)) N4-Return: $(sum((!ismissing).(data.df[:return])))")

  return data
end


function constructPrimaryOutput(; refreshFilter::Bool = true,
    primaryOutputName::String = PRIMARY_OUTPUT_NAME,
    writePrimaryCSV::Bool = false,
    dataPath::String = DATA_PATH)::NCCSData

  local data::NCCSData
  local outPath = "$dataPath\\$(primaryOutputName)_long.csv"
  if refreshFilter
    data = deserializeNCCSData()
    data = processAll!(data)
    serializeNCCSData(data, dataName = primaryOutputName)
  else
    data = deserializeNCCSData(dataName = primaryOutputName)
  end

  #for debugging
  #uCSV.write(outPath, data.df[400_000:450_000,:])

  writePrimaryCSV && uCSV.write(outPath, data.df) #write to a csv if desired

  return data

end
