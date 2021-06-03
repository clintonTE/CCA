


const ALL_COLS = [:name, :charitytype, :ein, :fisyr, :fileyear, :invrealized,
  :invrealizeddirect, :totrev, :totexp,
  :invrealizedindirect, :name, :netassets,
  :bookassets, :bookliabilities, :contributions,
  :pfexpbook, :pfexpcharity, :totinv, :category,:datefiscal,
  :returnamt, :returncheckamt,:netincome, :naics]

const CATEGORY_CONSOLIDATED_THRESHOLD = 100_000

#=const :Arts, :Education, :Environmental, :AnimalWelfare, :Health, :MentalHealth,
:MedicalDiseases, :MedicalResearch, :CrimeAndLegal, :Employment, :FoodAndAgriculture,
:Housing, :PublicSafety, :Recreation, :YouthDevelopment, :HumanServices, :ForeignAffairs,
:SocialAction, :CommunityImprovement, :GrantmakingFoundations, :ScienceandTechResearch,
:SocialScienceResearch, :SocietyBenefitAndMultipurpose, :ReligionRelated,
:MutualBenefitOrg, :Unknown, :Other=#

#processes the returns in a reasonable manner
function processReturns!(data::NCCSData; maxReturnCutoff::Float64 = MAX_RETURN_CUTOFF,
  minReturnCutoff::Float64 = MIN_RETURN_CUTOFF)

  local N::Int = size(data.df,1) #length of dataframe
  local revenueCutoff::Float64 = REVENUE_CUTOFF
  local largeDifferenceCutoff::Float64 = LARGE_DIFFERENCE_CUTOFF #max discrepenacy between methods


  println("data.df N1: $N")
  println("data.df pf N1: $(sum(data.df[:charitytype] .== :pf))")
  #filter out organizations which are too small
  data.df[:tokeep] = trues(size(data.df,1))

  #dfByein::GroupedDataFrame = groupby(data.df, :ein)
  for subdf::SubDataFrame ∈ groupby(data.df, :ein)
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

  data.df[:return] = Vector{Union{Float64,Missing}}(missings(N))

  #data.df[:returncheckamt] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:returncheck] = Vector{Union{Float64,Missing}}(missings(N))

  data.df[:largeDifference] = Vector{Union{Bool,Missing}}(missings(N))

  data.df[:lagnetassets] = Vector{MFloat64}(missing, N)
  data.df[:lag2netassets] = Vector{MFloat64}(missing, N)
  data.df[:lagcontributions] = Vector{MFloat64}(missing, N)

  sort!(data.df, [:ein, :fisyr])
  #ctr::Int = 0
  #dfByein = groupby(data.df, :ein)
  for subdf::SubDataFrame ∈ groupby(data.df, :ein)
    subN = size(subdf,1)

    if (subdf[1,:charitytype] == :co) || (subdf[1,:charitytype] == :pc)

      for i ∈ 2:subN

        if subdf[i-1, :fisyr] == subdf[i, :fisyr] - 1
          subdf[i, :lagnetassets] = subdf[i-1, :netassets]
          subdf[i, :lagcontributions] = subdf[i-1, :contributions]
          (i≥3) && (subdf[i, :lag2netassets] = subdf[i-2, :netassets])

          #NOTE: returns = unrealized + inv income s.t. unrealized = Δassets - netincome
          returnCandidateAmt::MFloat64 = (subdf[i, :netassets] .- subdf[i, :lagnetassets]
            ) .+ subdf[i, :invrealizeddirect] .- subdf[i, :netincome]
          returnCandidate::MFloat64 = returnCandidateAmt / subdf[i, :lagnetassets]

          returnCheckCandidateAmt::MFloat64 = (subdf[i, :netassets] .- subdf[i, :lagnetassets]
            ) .+ subdf[i, :invrealizedindirect] .- subdf[i, :netincome]
          returnCheckCandidate::MFloat64 = returnCheckCandidateAmt / subdf[i, :lagnetassets]

          if ((!ismissing(returnCandidate)) && (returnCandidate < maxReturnCutoff) &&
              (returnCandidate > minReturnCutoff) &&
              (abs(returnCandidate - returnCheckCandidate) ≤ largeDifferenceCutoff)
               && abs(subdf[i, :totinv]) ≥ 1.) #make sure we have valid numbers here
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
          subdf[i, :lagnetassets] = subdf[i-1, :netassets]
          subdf[i, :lagcontributions] = subdf[i-1, :contributions]
          (i≥3) && (subdf[i, :lag2netassets] = subdf[i-2, :netassets])

          #NOTE: returns = unrealized + inv income s.t. unrealized = Δassets - netincome
          returnCandidateAmt::MFloat64 = subdf[i, :returnamt]
          returnCandidate::MFloat64 = returnCandidateAmt / subdf[i, :lagnetassets]

          returnCheckCandidateAmt::MFloat64 = subdf[i, :returncheckamt]
          returnCheckCandidate::MFloat64 = returnCheckCandidateAmt / subdf[i, :lagnetassets]

          if ( (!ismissing(returnCandidate)) && (returnCandidate < maxReturnCutoff) &&
              (returnCandidate > minReturnCutoff) &&
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

#cleans the wealth if possible
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


function mergeSP500!(data::NCCSData;
    fiscaldatestr::String = FISCAL_DATE_FORMAT, refreshSP500::Bool = false,
    minReturnCutoff = MIN_RETURN_CUTOFF, maxReturnCutoff = MAX_RETURN_CUTOFF)

  #get the S&P500 data and index it
  local spdf = getWRDSSP500(refreshSP500 = refreshSP500)
  local idxbegin::Union{Int,Missing}
  local idxend::Union{Int,Missing}
  local spdfindex::Dict{Union{Date,Missing},MInt} =
    Dict{Union{Date,Missing},MInt}(spdf[i,:date]=>i for i ∈ 1:(size(spdf,1)))
  spdfindex[missing] = missing
  #println(spdf[1:5,:])

  #make date fields
  rename!(data.df, :datefiscal=>:datefiscalstr)
  data.df[:datefiscal] =
    (s->ismissing(s) ? missing : Date(s,fiscaldatestr)).(data.df[:datefiscalstr])

  #preallocate and make a lookup table for efficiency
  data.df[:sp500ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:lagsp500ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:t30ret] = Vector{MFloat64}(missing, size(data.df,1))
  data.df[:monthscovered] = Vector{MFloat64}(missing, size(data.df,1))


  printein=false
  printctr = 0

  #correct bad fiscal dates
  data.df[:datefiscal] = ((fisyr,datefiscal)->
    (!ismissing(datefiscal)) && (Year(datefiscal) != fisyr) ?
      Date(fisyr, Dates.value(Month(datefiscal)), 1) : datefiscal
      ).(data.df[:fisyr], data.df[:datefiscal])

  sort!(data.df, [:ein, :datefiscal])
  for subdf in groupby(data.df,:ein)
    N::Int = size(subdf,1)

    for i::Int ∈ 1:N
      idxend = spdfindex[subdf[i,:datefiscal]]
      if !ismissing(idxend)
        #assume eom dates (fiscal date given as yyyy-mm)
        #uses 12 month lag if no prior date is available
        idxbegin = (((i==1) || (ismissing(subdf[i-1,:datefiscal]))) ||
          subdf[i-1,:fisyr] != subdf[i,:fisyr] - 1 ?
          (idxend - 11) : (spdfindex[subdf[i-1,:datefiscal]] + 1))
        #(idxbegin > idxend) &&  @warn "ComputeNCCS df not sorted msg2341. I-1 dt
          #$(subdf[i-1,:datefiscal]) I dt $(subdf[i,:datefiscal])"
        printein = (idxbegin > idxend) || printein
        subdf[i,:sp500ret] = exp(sum(spdf[idxbegin:idxend, :lexcess]))-1.0
        (i≥2) && (subdf[i,:lagsp500ret] = subdf[i-1,:sp500ret])
        subdf[i,:t30ret] = exp(sum(spdf[idxbegin:idxend, :lt30ret]))-1.0
        subdf[i, :monthscovered] = idxend - idxbegin + 1
      end
    end
  end

  data.df[:excess] = data.df[:return] .- data.df[:t30ret]
  data.df[:pcontributions] = (data.df[:contributions]) ./ data.df[:lagnetassets]
  data.df[:pexpenses] = (data.df[:totexp]) ./ data.df[:lagnetassets]

  #trim rediculous values
  for v ∈ [data.df[:excess], data.df[:pcontributions], data.df[:pexpenses]]
    v .= (f->((ismissing(f)) || (f>maxReturnCutoff) || (f<minReturnCutoff)) ? missing : f).(v)
  end
  #data.df[:pcontributions] = (data.df[:totrev] .- data.df[:returncheckamt]) ./ data.df[:lagnetassets]

  #println(data.df[1000:1050,[:ein, :fisyr, :datefiscal, :return, :sp500ret]])

  #NOTE: Probably move this to another method
  #mod::FMLM = FMLM(data.df, Meta.parse("lagsp500ret"), :return,
  #  XNames = [:interecept, :lagsp500ret], YName=:pcontributions)
  #println(mod.β)

  return data
end

#ranks performance as a quantile
#NOTE: Returns can still be mis-aligned by year
function performanceByYear(data::NCCSData; sorted::Bool = false)
  local F::ECDF

  (!sorted) && sort!(data.df, [:ein, :fisyr])

  data.df[:preturnbyyr] = Vector{MFloat64}(undef, length(data.df[:return]))
  data.df[:preturnbyyrtype] = Vector{MFloat64}(undef, length(data.df[:return]))

  #split the dataframe and rank all of the returns
  for subdf::SubDataFrame ∈ groupby(data.df, :fisyr)
    F = ecdf(collect(skipmissing(subdf[:return])))
    subdf[:preturnbyyr] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(subdf[:return])

    for ssubdf::SubDataFrame ∈ groupby(subdf, :charitytype)
      F = ecdf(collect(skipmissing(ssubdf[:return])))
      ssubdf[:preturnbyyrtype] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(ssubdf[:return])
    end
  end

  #now do the lags
  data.df[:lpreturnbyyr] = Vector{MFloat64}(undef, length(data.df[:return]))
  data.df[:lpreturnbyyrtype] = Vector{MFloat64}(undef, length(data.df[:return]))

  for subdf::SubDataFrame ∈ groupby(data.df, :ein)
    for i ∈ 2:size(subdf,1)
      if subdf[i-1,:fisyr] == subdf[i, :fisyr] - 1 #only bring the return forward if one year gap
        subdf[i, :lpreturnbyyr] = subdf[i-1, :preturnbyyr]
        subdf[i, :lpreturnbyyrtype] =subdf[i-1, :preturnbyyrtype]
      end
    end
  end
  display(describe(data.df[[:fisyr, :preturnbyyr, :preturnbyyrtype, :lpreturnbyyr, :lpreturnbyyrtype]]))
  display(data.df[1000:1200,[:fisyr, :preturnbyyr, :preturnbyyrtype, :lpreturnbyyr, :lpreturnbyyrtype]])
end

#does an initial filter to get down the file size
function processAll!(data::NCCSData; allcols::Vector{Symbol} = ALL_COLS,
  revenueCutoff::Float64 = REVENUE_LOWER_BOUND,
  wealthCutoff::Float64 = WEALTH_LOWER_BOUND,
  refreshSP500::Bool = false
  )

  data = removeFields!(data, setdiff(names(data.df), allcols))
  data.df[:rows] = collect(1:size(data.df,1))
  println("N0: $(size(data.df,1))")
  deleterows!(data.df, data.df[(!).(completecases(data.df[[:fisyr, :netassets, :totrev, :totexp]])),:rows])

  adjustInflation!(data)
  processReturns!(data)
  performanceByYear(data, sorted=true) #NOTE: Assumes already sorted
  processWealth!(data)
  processCategories!(data)
  mergeSP500!(data, refreshSP500=refreshSP500)


  println("N3: $(size(data.df,1)) N3-Return: $(sum((!ismissing).(data.df[:return])))")

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
    dataPath::String = DATA_PATH,
    refreshSP500::Bool = false)::NCCSData

  local data::NCCSData
  local outPath = "$dataPath\\$(primaryOutputName)_long.csv"
  if refreshFilter
    data = deserializeNCCSData()
    data = processAll!(data, refreshSP500=refreshSP500)
    serializeNCCSData(data, dataName = primaryOutputName)
  else
    data = deserializeNCCSData(dataName = primaryOutputName)
  end

  #for debugging
  #uCSV.write(outPath, data.df[400_000:450_000,:])

  writePrimaryCSV && uCSV.write(outPath, data.df) #write to a csv if desired

  return data

end
