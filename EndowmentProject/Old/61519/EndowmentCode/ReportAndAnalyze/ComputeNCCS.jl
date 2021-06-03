#NOTE: We may want to convert this to monthly at some point

const ALL_COLS = [:name, :charitytype, :ein, :fisyr, :fileyear, :invrealized,
  :invrealizeddirect, :totrev, :totexp,
  :invrealizedindirect, :name, :netassets,
  :bookassets, :bookliabilities, :contributions,
  :pfexpbook, :pfexpcharity, :totinv, :category,:datefiscal,
  :returnamt, :returncheckamt,:netincome, :naics, :affiliationcode]

const CATEGORY_CONSOLIDATED_THRESHOLD = 50_000

#=const :Arts, :Education, :Environmental, :AnimalWelfare, :Health, :MentalHealth,
:MedicalDiseases, :MedicalResearch, :CrimeAndLegal, :Employment, :FoodAndAgriculture,
:Housing, :PublicSafety, :Recreation, :YouthDevelopment, :HumanServices, :ForeignAffairs,
:SocialAction, :CommunityImprovement, :GrantmakingFoundations, :ScienceandTechResearch,
:SocialScienceResearch, :SocietyBenefitAndMultipurpose, :ReligionRelated,
:MutualBenefitOrg, :Unknown, :Other=#

#processes the returns in a reasonable manner
function processReturns!(data::NCCSData; maxReturnCutoff::Float64 = MAX_RETURN_CUTOFF,
  minReturnCutoff::Float64 = MIN_RETURN_CUTOFF,
  winsorizeprop::Float64 = WINSORIZE_PROP)

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
  filter!(r->(r[:tokeep]),data.df)
  #data.df = data.df[data.df[:tokeep] .== true, :]
  N = size(data.df,1)

  println("data.df N2: $N")

  data.df[:return] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:lreturn] = Vector{Union{Float64,Missing}}(missings(N))

  #data.df[:returncheckamt] = Vector{Union{Float64,Missing}}(missings(N))
  data.df[:returncheck] = Vector{Union{Float64,Missing}}(missings(N))

  data.df[:largeDifference] = Vector{Union{Bool,Missing}}(missings(N))

  data.df[:lagnetassets] = Vector{MFloat64}(missing, N)
  data.df[:lag2netassets] = Vector{MFloat64}(missing, N)
  data.df[:lagcontributions] = Vector{MFloat64}(missing, N)

  sort!(data.df, [:ein, :fisyr])

  #double check we have no duplicate entries
  #@assert size(data.df,1) == size(unique(data.df[[:ein,:fisyr]]), 1)
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
          #write down 2 different methods
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
            #subdf[i, :lreturn] = log(1.0+returnCandidate)

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

  #winsorize!(data.df[:return], prop=winsorizeprop)
  data.df[:lreturn] = (log).(1.0 .+ data.df[:return])
  #NOTE: Already sorted, so don't need sort!(data.df, [:ein, :fisyr])

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

function filternoncentral!(data::NCCSData)::Nothing

  #conform the common affiliation codes
  data.df[:affiliationcode2] = (s::MSymbol->affiliationcodes[s]).(
    Vector{MSymbol}(data.df[:affiliationcode]))

  #filter out non-central orgs
  filter!(r->conformedaffiliationcodes[r[:affiliationcode2]]==:central, data.df)


  return nothing
end

#adjusts for inflation
#currently just does wealth, eventually can use this for return fields
function adjustInflation!(data::NCCSData;
  fredPath::String = FRED_PATH,
  fredinflationfile::String = FRED_INFLATION_FILE,
  inflationseries::Symbol = :CPI,
  baseyear::Int = BASE_YEAR,
  freddategormat::String = FRED_DATE_FORMAT,
  adjustCutoff::Bool = true,
  revenuelowerbound::Float64 = REVENUE_LOWER_BOUND,
  wealthlowerbound::Float64 = WEALTH_LOWER_BOUND,
  lowerezrevenuelimit::Float64 = LOWER_EZ_REVENUE_LIMIT,
  prefilterlowerbound::Float64 = PREFILTER_REVENUE_LOWER_BOUND)::Nothing

  #read in the inflation info
  freddf::DataFrame = CSV.read("$fredPath\\$fredinflationfile.csv")

  #make a dictionary
  cpiindex::Dict =
    Dict(Dates.year(Date(freddf[i,:DATE], freddategormat))=>freddf[i, inflationseries]
      for i ∈ 1:size(freddf,1))

  #get vectors for performance purpuses
  nrows::Int = size(data.df, 1)
  basevalue::Float64 = cpiindex[baseyear]
  #=netassets::Vector{Float64} = data.df[:netassets]
  totrev::Vector{Float64} = data.df[:totrev]
  adjnetassets =
  adjtotrev =
  fisyr::Vector{Int} = data.df[:fisyr]
  @inbounds @simd for i ∈ 1:nrows
    if haskey(cpiindex, fisyr[i])
      adjnetassets[i] = netassets[i] * (basevalue / cpiindex[fisyr[i]])
      adjtotrev[i] = totrev[i] * (basevalue / cpiindex[fisyr[i]])
    end
  end

  data.df[:adjnetassets] = adjnetassets
  data.df[:adjtotrev] = adjtotrev=#

  data.df[:adjnetassets] = Vector{MFloat64}(undef, nrows)
  data.df[:adjtotrev] = Vector{MFloat64}(undef, nrows)
  data.df[:tokeep] = trues(nrows)

  for subdf ∈ groupby(data.df, :fisyr)
    currentfisyr::Int = subdf[1,:fisyr]

    #make sure its a valid fisyr and set the inflation factors
    current2basefactor::Float64 = (haskey(cpiindex, currentfisyr) ?
      basevalue / cpiindex[currentfisyr] : 0.0)
    base2currentfactor::Float64 = (haskey(cpiindex, currentfisyr) ?
      cpiindex[currentfisyr] / basevalue : 0.0)

    wealthcutoff::Float64 = wealthlowerbound * base2currentfactor
    revenuecutoff::Float64 = revenuelowerbound * base2currentfactor
    (lowerezrevenuelimit ≥ revenuecutoff) &&  @warn (
      "$currentfisyr revenue limit ($revenuecutoff) below IRS minimum" *
      "$(lowerezrevenuelimit). Expect selection bias.")

    (prefilterlowerbound ≥ revenuecutoff) && @warn (
      "$currentfisyr revenue limit ($revenuecutoff) below inflation adjusted revenue limit" *
      "$(prefilterlowerbound). Expect moving limit in real terms.")

    #adjust assets and revenue to the base year
    subdf[:adjnetassets] .= (f::MFloat64->f*current2basefactor).(subdf[:netassets])
    subdf[:adjtotrev] .= (f::MFloat64->f*current2basefactor).(subdf[:totrev])

    #delete rows below the limit
    subdf[:tokeep] .= ((fna::MFloat64, ftr::MFloat64) ->
      (!ismissing(fna)) && (!ismissing(ftr)) &&
      (wealthcutoff ≤ fna) && (revenuecutoff ≤ ftr)).(
      subdf[:netassets], subdf[:totrev])

#    println("$currentfisyr Minimum assets: $(minimum(subdf[subdf[:tokeep],:totrev])) " *
#      "Adj: $(minimum(subdf[subdf[:tokeep],:adjtotrev]))")
#    println("$currentfisyr Minimum assets: $(minimum(subdf[subdf[:tokeep],:netassets])) " *
#      "Adj: $(minimum(subdf[subdf[:tokeep],:adjnetassets]))")
  end


  data.df = data.df[data.df[:tokeep],:]
  #=if adjustCutoff #this is to account for the cutoff inflating forward
    inflatedrevenuecutoff::Float64 = REVENUE_LOWER_BOUND * basevalue/cpiindex[minimum(data.df[:fisyr])]
    println("Inflated revenue cutoff of $inflatedRevenueCutoff")
    data.df = data.df[
      (f::MFloat64-> (!ismissing(f)) && (f ≥ inflatedrevenuecutoff)).(data.df[:adjtotrev]), :]

    inflatedwealthcutoff::Float64 = WEALTH_LOWER_BOUND * basevalue/cpiindex[minimum(data.df[:fisyr])]
    data.df = data.df[
      (f::MFloat64-> (!ismissing(f)) && (f ≥ inflatedwealthcutoff)).(data.df[:adjnetassets]), :]
  end=#


  return nothing
end

#creates quantiles given a benchmark
function relativeperformance!(data::NCCSData, period::Int, lbenchmark::NSymbol;
  sorted::Bool = false, lretsym1y::NSymbol = :lreturn)::Nothing
  (!sorted) && sort!(data.df, [:ein, :fisyr])

  local N::Int = size(data.df,1)

  #first establish the field names
  lbenchstring::String = isnothing(lbenchmark) ? "" : string("net", lbenchmark)

  retstring1y::String = string(something(lretsym1y, "zero"))

  local retsym = Symbol("$retstring1y", lbenchstring, period, "yr")
  local Lretsym = Symbol("L$(period)$retstring1y", lbenchstring, period, "yr")
  local percentilesym::Symbol = Symbol("p$retstring1y", lbenchstring, period, "yr")
  local percentilebytypesym::Symbol = Symbol("p$retstring1y", lbenchstring, period, "yrtype")
  local Lpercentilesym::Symbol = Symbol("L$(period)p$retstring1y", lbenchstring, period, "yr")
  local Lpercentilebytypesym::Symbol = Symbol("L$(period)p$retstring1y", lbenchstring, period, "yrtype")


  #this will be useful
  local uniqueyrs::Vector{Int} = collect(
    minimum(data.df[:fisyr]):maximum(data.df[:fisyr]))#unique(data.df[:fisyr])

  #pre-allocate
  data.df[retsym] = Vector{MFloat64}(undef, N)
  data.df[Lretsym] = Vector{MFloat64}(undef, N)
  data.df[percentilesym] = Vector{MFloat64}(undef, N)
  data.df[percentilebytypesym] = Vector{MFloat64}(undef, N)
  data.df[Lpercentilesym] = Vector{MFloat64}(undef, N)
  data.df[Lpercentilebytypesym] = Vector{MFloat64}(undef, N)
  data.df[:ltempreturn] = Vector{MFloat64}(undef, N)

  #compute the excess returns, net of a benchmark if there is one
  data.df[:ltempreturn] .= isnothing(lretsym1y) ? zeros(N) : data.df[lretsym1y]
  data.df[:ltempreturn] .-= isnothing(lbenchmark) ? 0.0 : data.df[lbenchmark]

  #split the dataframe and compute the cumulative reutrns
  for subdf::SubDataFrame ∈ groupby(data.df, :ein)
    for i ∈ (1+period):size(subdf,1)
      #make sure we are covering the correct amount of time
      if subdf[i,:fisyr] == subdf[i-period,:fisyr] + period
        subdf[i, retsym] = sum(subdf[collect((i-(period-1)):i), :ltempreturn])
      end

    end
  end

  #split the dataframe and rank all of the returns
  for subdf::SubDataFrame ∈ groupby(data.df, :fisyr)

    F = ecdf(collect(skipmissing(subdf[retsym])))

    if isfinite(F(0.0)) #we can shortcut this loop if there are no valid ecdf values
      subdf[percentilesym] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(subdf[retsym])

      #now rank within types
      for ssubdf::SubDataFrame ∈ groupby(subdf, :charitytype)
        F = ecdf(collect(skipmissing(ssubdf[retsym])))
        ssubdf[percentilebytypesym] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(ssubdf[retsym])
      end
    end
  end

  #finally record the lag quantile of one period where available
  for subdf::SubDataFrame ∈ groupby(data.df, :ein)
    for i ∈ (1+2 * period):size(subdf,1)
      #make sure we are covering the correct amount of time
      if subdf[i,:fisyr] == subdf[i-period,:fisyr] + period

        #copy in the lag quantile
        subdf[i, Lretsym] = subdf[i - period,retsym]
        subdf[i, Lpercentilesym] = subdf[i - period,percentilesym]
        subdf[i, Lpercentilebytypesym] = subdf[i - period,percentilebytypesym]
      end
    end
  end

  return nothing
end


#ranks performance as a quantile
#NOTE: Returns can still be mis-aligned by year
function relativeperformance!(data::NCCSData; sorted::Bool = false,
    persistenceperiods::Vector{Int} = PERSISTENCE_PERIODS,
    lbenchmarks::Vector{NSymbol} = LBENCHMARKS,
    lbmxfields = LBM_XFIELDS)::Nothing

  for s ∈ lbenchmarks
    #need to build up column names so we don't overwrite anything

    #different procedure if no benchmark or we are not beta adjusting
    if isnothing(s)
      (i::Int->relativeperformance!(data, i, s)).(persistenceperiods)
    else
      regressbynonprofit!(data, bmfield=s, YField=:lreturn, XFields=lbmxfields[s])

      (i::Int->relativeperformance!(data, i, s)).(persistenceperiods)
      #(i::Int->relativeperformance(data, i, predlbenchmark,
      #  lretsym1y = nothing)).(persistenceperiods)
    end
  end

  return nothing
end

#lags a variable in the time series
function lagwithin!(data::NCCSData, target::Symbol;
  period::Symbol = :fisyr,
  lags::Int=1,
  laggedname = lags ≠ 1 ? Symbol(:L, lags, target) : Symbol(:L, target),
  group=:ein,
  sorted::Bool=false)::Nothing

  if !sorted
    sort!(data.df, [group, period])
  end

  #pre-allocate the space
  data.df[laggedname] = Vector{Union{eltype(data.df[target]), Missing}}(undef, size(data.df, 1))

  for subdf ∈ groupby(data.df, group)
    Nsub::Int = size(subdf, 1)

    #only run these routines if we need to
    if Nsub > lags
      for i ∈ (lags+1):Nsub #iterate for all values
        fisyr::Int = subdf[i, :fisyr]

        #first see if we can lag the easy way
        if subdf[i-lags, :fisyr] == fisyr - lags
          subdf[i,laggedname] = subdf[i-lags, target]
        else #now check the hard way (finding the lagged year)
          location::NInt = findfirst(isequal(fisyr-lags), subdf[1:(i-1),:fisyr])

          #record the lagged value if it is available
          if !isnothing(location)
            subdf[i,laggedname] = subdf[location, target]
          end
        end
      end
    end
  end

  return nothing
end

#creates a differenced column
function differencewithin!(data::NCCSData, target::Symbol;
  period::Symbol = :fisyr,
  differencedname = Symbol(:D, target),
  group=:ein,
  sorted::Bool=false,
  createlag::Bool=true,
  deletelag::Bool=true,
  laggedname::Symbol = deletelag ? Symbol(:L, target, :_temp) : Symbol(:L, target)
  )::Nothing

  createlag && lagwithin!(data, target,
    period=period,
    laggedname=laggedname,
    group=group,
    sorted=sorted)

  data.df[differencedname] = data.df[target] .- data.df[laggedname]

  deletelag && deletecols!(data.df, laggedname)

  return nothing

end

#does an initial filter to get down the file size
function processAll!(data::NCCSData; allcols::Vector{Symbol} = ALL_COLS,
  revenueCutoff::Float64 = REVENUE_LOWER_BOUND,
  wealthCutoff::Float64 = WEALTH_LOWER_BOUND,
  refreshoutsidedata::Bool = false,
  filternoncentral::Bool = true,
  refreshopen990data = true
  )

  data = removeFields!(data, setdiff(names(data.df), allcols))
  data.df[:rows] = collect(1:size(data.df,1))
  println("N0: $(size(data.df,1))")
  deleterows!(data.df, data.df[(!).(completecases(data.df[[:fisyr, :netassets, :totrev, :totexp]])),:rows])

  adjustInflation!(data)
  processReturns!(data)
  processWealth!(data)
  processCategories!(data)
  mergeoutsidedatasources!(data, refreshoutsidedata=refreshoutsidedata)
  relativeperformance!(data, sorted=false) #NOTE: Assumes already sorted
  mergeopen990data!(data, refreshopen990data=refreshopen990data)

  println("N3: $(size(data.df,1)) N3-Return: $(sum((!ismissing).(data.df[:return])))")
  filternoncentral && filternoncentral!(data)
  filternoncentral && println("Filtered out intermediate organziations.")
  filter!(r->(!ismissing(r[:adjtotrev])) && (!ismissing(r[:adjnetassets])) &&
    (r[:adjtotrev] ≥ revenueCutoff) && (r[:adjnetassets] ≥ wealthCutoff) , data.df)
  #data.df = data.df[((!ismissing).(data.df[:adjtotrev])) .&
  #  (data.df[:adjtotrev] .≥ revenueCutoff), :]
  #data.df = data.df[((!ismissing).(data.df[:adjnetassets])) .&
  #  (data.df[:adjnetassets] .≥ wealthCutoff), :]

  println("N4: $(size(data.df,1)) N4-Return: $(sum((!ismissing).(data.df[:return])))")

  return data
end
