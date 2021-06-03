const WRDS_FLDS = [:gvkey, :datadate, :fyear, :tic, :cusip, :conm, :indfmt,
  :curcd, :curncd, :costat, :loc, :mkvalt, :prcc_f,
  :stalt, :stko, :idbflag, :revt, :adrr, :seq, :naics, :at]
const WRDS_OK_MISSING_FLDS = [:adrr, :stalt]


function wrdsAggDF(; wrdsPath::String = WRDS_PATH, compustatName::String = COMPUSTAT_NAME,
    outcols::Vector{Symbol} = [:naics2name, :conm, :revt, :at, :seq],
    datasuffix::String=DATA_SUFFIX)::DataFrame

  local outdf::DataFrame
  local compdf::DataFrame
  #now work on the compustat data
  compdf = CSV.read("$wrdsPath\\$compustatName.csv")
  compdf = compdf[WRDS_FLDS]
  compdf = compdf[completecases(compdf[setdiff(WRDS_FLDS, WRDS_OK_MISSING_FLDS)]), :]

  filterfunc(r::DataFrameRow) = (
    r[:indfmt] == "INDL" &&
    r[:curcd] == "USD" && #native currency and reported currency in USD
    r[:curncd] == "USD" && #ibid
    #r[:costat] == "A" && #active in db
    r[:loc] == "USA" && #incorporated or registered in US
    r[:mkvalt] > 0.1 && #has a market value > 100k
    (ismissing(r[:stalt]) ||  r[:stalt] == "AW") &&#not in reogranization/bankruptcy
    r[:stko] < 2 && #ignore very thinly traded securities (no stale quotes)
    #r[:idbflag] ≠ "I" && #incorporated or registered in US
    r[:revt] > 0.1 && #revenue is reported and > 100k
    (ismissing(r[:adrr])) && #not an adr
    r[:seq] > 0.1 #book equity reported and > 100k
  )

  compdf = compdf[((r::DataFrameRow)->filterfunc(r)).(eachrow(compdf)),:]

  #make sure we don't have dupe companies
  compdf[:tokeep] = trues(size(compdf,1))
  compdfBygvkey::GroupedDataFrame = groupby(compdf, :gvkey)
  for subdf::SubDataFrame ∈ compdfBygvkey
    if size(subdf,1) > 1
      subdf[:tokeep] = subdf[:datadate] .== maximum(subdf[:datadate])
      (sum(subdf[:tokeep]) > 1) && error("Duplicate firm remains in compustat data, add new filtering rule")
    end
  end
  compdf = compdf[compdf[:tokeep],:]
  compdf[:naics2] = (s::Int->Symbol(("$s")[1:2])).(compdf[:naics])
  compdf[:naics2name] = (s::Symbol->(naicsCodes[s])).(compdf[:naics2])

  #now write out the largest firms
  sort!(compdf, (:naics2name, order(:at,rev=true)))
  outdf = DataFrame(deepcopy(compdf[[],outcols]))

  for subdf::SubDataFrame ∈ groupby(compdf, :naics2name)
    numentries::Int = min(3,size(subdf,1))
    append!(outdf, DataFrame(deepcopy(subdf[1:numentries,outcols])))
  end
  outdf[:firmtype] = :nonprofit

  #housekeeping to match the output dataframe
  rename!(outdf, [:revt=>:revenue, :conm=>:name, :at=>:assets, :seq=>:netrentassets])
  CSV.write("$WORKING_PATH\\largestPublicForProfits_$datasuffix.csv", outdf)

  #this aggregates the industry data using new data frames 15.2+ functions
  compaggdf::DataFrame = by(
    compdf, :naics2name,
      [:naics2name,:mkvalt, :seq, :revt, :at] =>
        x -> (n = length(x.naics2name),
        revenue=sum(skipmissing(x.revt)) * 1_000_000.,
        assets=sum(skipmissing(x.at)) * 1_000_000.,
        netassets=sum(skipmissing(x.seq)) * 1_000_000.))

  #create normalized versions
  compaggdf[:drevenue] = compaggdf[:revenue] ./ sum(compaggdf[:revenue])
  compaggdf[:dassets] = compaggdf[:assets] ./ sum(compaggdf[:assets])
  compaggdf[:dnetassets] = compaggdf[:netassets] ./ sum(compaggdf[:netassets])

  compaggdf[:firmtype] = :publicforprofit

  return compaggdf
end

const SOI_FLDS = [:naics2, :n,	:assets, :netassets,	:revenue]

#this grabs the pre-aggregated SOI data (already aggregated by IRS)
function soiAggDF(irsSOIPath::String = IRS_SOI_PATH, soiName::String = SOI_NAME)::DataFrame

  compaggdf::DataFrame = CSV.read("$irsSOIPath\\$(soiName).csv")
  compaggdf = compaggdf[SOI_FLDS]

  compaggdf[:naics2name] = (i::Int->(naicsCodes[Symbol("$i")])).(compaggdf[:naics2])
  deletecols!(compaggdf, :naics2)

  #create normalized versions
  compaggdf[:drevenue] = compaggdf[:revenue] ./ sum(compaggdf[:revenue])
  compaggdf[:dassets] = compaggdf[:assets] ./ sum(compaggdf[:assets])
  compaggdf[:dnetassets] = compaggdf[:netassets] ./ sum(compaggdf[:netassets])

  compaggdf[:firmtype] = :soiforprofit

  return compaggdf
end

#make the industry distribution plots
function industryCrossDF(data::NCCSData; assetYear::Int = INDUSTRY_CROSS_YEAR,
    industryComparisonSource::Symbol=:wrds,
    outcols::Vector{Symbol} = [:naics2name, :name, :ein, :totrev, :bookassets, :bookliabilities],
    datasuffix::String=DATA_SUFFIX
    )::DataFrame
  local presubdf::SubDataFrame = view(data.df, data.df[:fisyr].==assetYear,:)
  local comparisondf::DataFrame
  local outdf::DataFrame

  #write out the largest firms by assets
  #need to copy it to use the sort command
  local preoutdf::DataFrame = DataFrame(deepcopy(presubdf[:,outcols]))
  sort!(preoutdf, (:naics2name, order(:bookassets, rev=true)))

  outdf = deepcopy(preoutdf[[],:]) #create a stub df for collecting
  for subdf::SubDataFrame ∈ groupby(preoutdf, :naics2name)
    numentries::Int = min(3,size(subdf,1))
    append!(outdf, DataFrame(deepcopy(subdf[1:numentries,outcols])))
  end
  outdf[:netassets] = outdf[:bookassets] .- outdf[:bookliabilities]
  outdf[:firmtype] = :nonprofit

  #housekeeping to match the output dataframe
  deletecols!(outdf, :bookliabilities)
  rename!(outdf, [:totrev=>:revenue, :bookassets=>:assets])
  CSV.write("$WORKING_PATH\\largestNonProfits_$datasuffix.csv", outdf)

  #this aggregates the industry data using new data frames 15.2+ functions
  local aggdf::DataFrame = by(
    presubdf[:,[:naics2, :naics2name,:totrev, :totexp, :netassets, :bookassets, :bookliabilities]],
      :naics2name,
      [:naics2name,:totrev,:totexp,:netassets, :bookassets, :bookliabilities] =>
        x -> (n = length(x.naics2name),
        revenue=sum(skipmissing(x.totrev)),
        #profits = sum(skipmissing(x.totrev .- x.totexp)),
        assets=sum(skipmissing(x.bookassets)),
        netassets=sum(skipmissing(x.bookassets .- x.bookliabilities))))
  #create normalized versions
  (assetYear ≥ 2000) && (aggdf[:drevenue] = aggdf[:revenue] ./ sum(aggdf[:revenue]))
  aggdf[:dassets] = aggdf[:assets] ./ sum(aggdf[:assets])
  aggdf[:dnetassets] = aggdf[:netassets] ./ sum(aggdf[:netassets])

  aggdf[:firmtype] = :nonprofit

  #load the comparison df based on the requested source
  if industryComparisonSource == :wrds
    comparisondf = wrdsAggDF()
  elseif industryComparisonSource == :soi
    comparisondf = soiAggDF()
  end

  assetYear < 2000 && deletecols!(comparisondf, :revenue)

  aggdf = vcat(aggdf, comparisondf)

  #missing values really correspond to 0
  for s::Symbol ∈ [:revenue, :assets, :netassets, :drevenue, :dassets, :dnetassets]
    aggdf[:,s] .= (f::MFloat64->ismissing(f) ? 0 : f).(aggdf[:,s])
  end

  #this snippit creates a 0-value row if an industry is missing
  #the purpose is to provide the placeholder value to the bar charts
  for s::Symbol ∈ unique(aggdf[:naics2name])
    for t::Symbol ∈ unique(aggdf[:firmtype])
      if t ∉ aggdf[aggdf[:naics2name] .== s, :firmtype]
        aggdf = vcat(aggdf, DataFrame(naics2name=[s], firmtype=[t],n=[0.],
          revenue=[0.], assets=[0.], netassets=[0.],
          drevenue=[0.], dassets=[0.], dnetassets=[0.]))
      end
    end
  end

  #create proportional logged results
  #each value is the corresponding percentage of the total unlogged value in each industry
  # times the logged value
  aggdf[:prevenue] = similar(aggdf[:revenue])
  aggdf[:passets] = similar(aggdf[:assets])
  aggdf[:pnetassets] = similar(aggdf[:netassets])
  aggdfBynaics2name::GroupedDataFrame = groupby(aggdf, :naics2name)
  for subdf::SubDataFrame ∈ aggdfBynaics2name

    subdf[:prevenue] = log(sum(subdf[:revenue])) .* (subdf[:revenue] ./ sum(subdf[:revenue]))
    subdf[:passets] = log(sum(subdf[:assets])) .* (subdf[:assets] ./ sum(subdf[:assets]))
    subdf[:pnetassets] = log(sum(subdf[:netassets])) .* (subdf[:netassets] ./ sum(subdf[:netassets]))
  end


  return aggdf
end

function industryCross(data::NCCSData;
        outputPath::String = GRAPH_PATH,
        industryComparisonSource::Symbol = INDUSTRY_COMPARISON_SOURCE)::Nothing

  #runs the routine for multiple comparison sources if desired
  if industryComparisonSource == :soiwrds
    industryCross(data, outputPath=outputPath, industryComparisonSource=:soi)
    industryComparisonSource = :wrds
  end

  local aggdf::DataFrame = industryCrossDF(data, industryComparisonSource=industryComparisonSource)
  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  sort!(aggdf, [:naics2name, :firmtype])

  #first make the concentration plots
  push!(plotNames, "Concentration_of_Revenue_by_Industry")
  push!(plots, plot(aggdf,
    x=:naics2name, y=:drevenue, color=:firmtype,
    Guide.title("Distribution of Revenue by Industry"),
    Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Percentage of Revenue"),
     Geom.bar(position=:dodge)))

   push!(plotNames, "Concentration_of_Assets_by_Industry")
   push!(plots, plot(aggdf,
     x=:naics2name, y=:dassets, color=:firmtype,
     Guide.title("Distribution of Assets by Industry"),
     Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Percentage of Total Assets"),
      Geom.bar(position=:dodge),
      Coord.cartesian(ymin=0.0, ymax=0.4)))

   push!(plotNames, "Concentration_of_Book_Equity_by_Industry")
   push!(plots, plot(aggdf,
     x=:naics2name, y=:dnetassets, color=:firmtype,
     Guide.title("Distribution of Equity by Industry"),
     Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Percentage of Total Book Equity"),
      Geom.bar(position=:dodge)))

  #lets make some stacked plots
  push!(plotNames, "Revenue_by_Industry")
  push!(plots, plot(aggdf,
    x=:naics2name, y=:prevenue, color=:firmtype,
    Guide.title("Revenue by Industry"),
    Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Log Revenue (Color areas are % of total bar)"),
     Geom.bar(position=:stack)))

   push!(plotNames, "Assets_by_Industry")
   push!(plots, plot(aggdf,
     x=:naics2name, y=:passets, color=:firmtype,
     Guide.title("Assets by Industry"),
     Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Log Assets (Color areas are % of total bar)"),
      Geom.bar(position=:stack)))

   push!(plotNames, "Book_Equity_by_Industry")
   push!(plots, plot(aggdf,
     x=:naics2name, y=:pnetassets, color=:firmtype,
     Guide.title("Equity by Industry"),
     Guide.xlabel("Sector (NAICS-2)"), Guide.ylabel("Log Book Equity (Color areas are % of total bar)"),
      Geom.bar(position=:stack)))

  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i])_$(industryComparisonSource)_$INDUSTRY_CROSS_YEAR.pdf", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-secitonal industry graphs drawn.")

  return nothing
end
