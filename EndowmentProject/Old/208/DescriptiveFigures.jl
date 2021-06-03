#TODO: Inflation-index the wealth distributions in figuresByCategory


##Graph types and constants
const COLORS = ["DodgerBlue", "OrangeRed", "Green",
  "Brown", "Black", "BlueViolet"]

const PlotContainer = Union{Plot,Gadfly.Compose.Context,Vector{Plot}} #holds graphs
const ASSET_YEARS = [2015, 2010, 2005, 2000, 1997]
const ASSET_YEARS_FULL = collect(1997:2015)
const WEALTH_FOCAL_FIELD = :adjnetassets

const RETURN_FOCAL_FIELD = :return
const RETURN_YEARS = ASSET_YEARS
const RETURN_YEARS_FULL = collect(2000:2015)

const CATEGORY_FIGURE_FIELD = :categorynameconsolidated
const BIN_ON_FIELDS = [:fisyr, CATEGORY_FIGURE_FIELD]


#NOTE: Create a new more consolidated category field

getNonmissingCategories(data::NCCSData; categoryField::Symbol=CATEGORY_FIGURE_FIELD) =
  unique(data.df[(!ismissing).(data.df[categoryField]), categoryField])

#const FIGURES_MIN_ASSETS = 1_500_000.

#contains the functions for writing the grpah by year
function wealthFiguresByYear(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS,
    outputPath::String = GRAPH_PATH, focalField::Symbol=WEALTH_FOCAL_FIELD,
    categoryFigureField::Symbol = CATEGORY_FIGURE_FIELD,
    categoriesUsed::Vector = getNonmissingCategories(data, categoryField = categoryFigureField))

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  #local presubdf::SubDataFrame =  view(data.df, data.df[focalField] .≥ figuresMinAssets,:)

  for y::Int ∈ assetYears #main wealth distribution
    subdf = view(data.df, data.df[:fisyr] .== y,:) #make a view for performance reasons

    push!(plotNames, "Wealth_Distribution_of_Non-Profits_($y)-dodge")
    push!(plots, plot(subdf, x=focalField, color=:charitytype,
      Guide.title("Wealth Distribution of Non-Profits ($y)-dodge"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=false, position=:dodge),
      Coord.cartesian(ymin=0.0, ymax = 2500.0, xmin=5.0)))

    push!(plotNames, "Wealth_Distribution_of_Non-Profits_($y)")
    push!(plots, plot(subdf,
      x=focalField, color=categoryFigureField,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),  Scale.x_log10,
      Geom.line, Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 3000.0, xmin=5.0)
      ))


  #=  push!(plotNames, "Wealth_Distribution_of_Private_Charities_($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
      x=focalField, color=categoryFigureField,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
      Scale.x_log10, Geom.line,
      Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 2500.0, xmin=5.0)
      ))


    push!(plotNames, "Wealth_Distribution_of_Private_Foundations_($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
      x=focalField, color=categoryFigureField,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
      Scale.x_log10, Geom.line,
      Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=5.0)
      ))=#

    #=push!(plotNames, "Log_Wealth_Distribution_of_Non-Profits_($y)")
    push!(plots, plot(subdf,
    x=Symbol('l',focalfield), color=:charitytype,
      Guide.title("Log Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=true)))=#
  end

  for i ∈ 1:length(plots) #write the graphs
    try
      draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    catch err
      @error "ERROR writing $(plotNames[i]). Error message $err"
      println("Size: $(size(subdf))")
    end
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Wealth Graphs by Year Drawn")
  return nothing
end

#contains the functions for writing the grpah by year
function returnFiguresByYear(data::NCCSData; returnYears::Vector{Int} = RETURN_YEARS,
    outputPath::String = GRAPH_PATH, focalField::Symbol=RETURN_FOCAL_FIELD,
    categoryFigureField::Symbol = CATEGORY_FIGURE_FIELD,
    categoriesUsed::Vector = getNonmissingCategories(data, categoryField = categoryFigureField))

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame = view(data.df, (!ismissing).(data.df[focalField]),:)
  presubdf = view(presubdf, (!ismissing).(presubdf[categoryFigureField]),:)

  for y::Int ∈ returnYears #main wealth distribution
    subdf = view(presubdf, presubdf[:fisyr] .== y,:) #make a view for performance reasons

    #println("\n####\nYear: $y; rows: $(size(subdf,1)),
    #  rowspf: rows: $(sum(subdf[:charitytype] .== :pf)),
    #  rowspc: rows: $(sum(subdf[:charitytype] .== :pc)))")

    push!(plotNames, "Return_Distribution_of_Non-Profits_($y)_by_Type")
    push!(plots, plot(subdf, x=focalField, color=:charitytype,
      Guide.title("Return Distribution of Non-Profits ($y) by Type"),
      Guide.xlabel("Return"),Guide.ylabel("Density"),
      Geom.histogram(bincount=50, density=true#=, position=:dodge=#),
      Coord.cartesian(xmin=-0.5, xmax=0.5)))

    push!(plotNames, "Return_Distribution_of_Non-Profits_($y)_by_Category")
    push!(plots, plot(subdf,
      x=focalField, color=categoryFigureField,
      Guide.title("Return_Distribution_of_Non-Profits_($y)_by_Category"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line, Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)
      ))


    push!(plotNames, "Return_Distribution_of_Private_Charities_($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
      x=focalField, color=categoryFigureField,
      Guide.title("Return Distribution of Non-Profits ($y)"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line,
      Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)
      ))


    push!(plotNames, "Return_Distribution_of_Private_Foundations_($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
      x=focalField, color=categoryFigureField,
      Guide.title("Return Distribution of Non-Profits ($y)"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line,
      Stat.histogram(bincount=50, density=true),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)
      ))

    #=push!(plotNames, "Log_Wealth_Distribution_of_Non-Profits_($y)")
    push!(plots, plot(subdf,
    x=Symbol('l',focalfield), color=:charitytype,
      Guide.title("Log Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=true)))=#
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Return Graphs by Category Drawn")
  return nothing

end

#contains the functions for writing the graphs by category
function wealthFiguresByCategory(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS,
    outputPath::String = GRAPH_PATH,
    minPointsPerCategory::Int = MIN_POINTS_PER_CATEGORY,
    focalField::Symbol=WEALTH_FOCAL_FIELD,
    categoryFigureField::Symbol = CATEGORY_FIGURE_FIELD,
    categoriesUsed::Vector = getNonmissingCategories(data, categoryField = categoryFigureField))

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame

  #get only the categories with valid values
  presubdf = view(data.df, ((y::Int)-> y ∈ assetYears).(data.df[:fisyr]),:)
  #presubdf = view(presubdf, presubdf[focalField] .≥ figuresMinAssets)

  push!(plotNames, "Wealth_Distribution_of_Non-Profits_by_Year")
  push!(plots, plot(presubdf,
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Non-Profits by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=50, density=true),
    Coord.cartesian(ymin=0.0, ymax = 0.8, xmin=5.)),)

  push!(plotNames, "Wealth_Distribution_of_Public_Charities_by_Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pc, :],
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Public Charities by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=50, density=true),
    Coord.cartesian(ymin=0.0, ymax = 1., xmin=5.)))

  push!(plotNames, "Wealth_Distribution_of_Private_Foundations_by_Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pf,:],
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Private Foundations by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=50, density=true),
    Coord.cartesian(ymin=0.0, ymax = 1., xmin=5.)))

  #filter out funds with no category
  presubdf = view(presubdf, (!ismissing).(presubdf[categoryFigureField]),:)


  #=for c::MSymbol ∈ categoriesUsed

    #filter to the category and exclude non-studied asset years
    subdf = view(presubdf, Vector{Bool}(presubdf[categoryFigureField] .== c), :)

    #we need to make sure we have a viable sample
    if size(subdf,1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth_Distribution_of Non-Profit_($c)")
      push!(plots, plot(subdf,
        x=focalField, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(ymin=0.0, ymax = 1000.0)
        ))
    end

    if size(subdf[subdf[:charitytype] .== :pc,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth_Distribution_of_Public Charities_($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
        x=focalField, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(ymin=0.0, ymax = 1000.0)))
    end

    if size(subdf[subdf[:charitytype] .== :pf,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth_Distribution_of_Private_Foundations_($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
        x=:adjnetassets, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=5.0)))
    end
  end=#

  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Wealth Graphs by Category Drawn.")

  return nothing
end

#contains the functions for writing the graphs by category
function returnFiguresByCategory(data::NCCSData;
    returnYears::Vector{Int} = RETURN_YEARS,
    outputPath::String = GRAPH_PATH,
    minPointsPerCategory::Int = MIN_POINTS_PER_CATEGORY,
    focalField::Symbol=RETURN_FOCAL_FIELD,
    categoryFigureField::Symbol = CATEGORY_FIGURE_FIELD,
    categoriesUsed::Vector = getNonmissingCategories(data, categoryField = categoryFigureField))



  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame

  #get only the categories with valid values
  presubdf = view(data.df, ((y::Int)-> y ∈ returnYears).(data.df[:fisyr]),:)
  presubdf = view(presubdf, (!ismissing).(presubdf[:return]), :)

  push!(plotNames, "Return_Distribution_of_Non-Profits_by_Year")
  push!(plots, plot(presubdf,
    x=focalField, color=:fisyr,
    Guide.title("Return of Non-Profits by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=50, density=false),
    Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))

  #=push!(plotNames, "Return_Distribution_of_Public Charities_by_Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pc, :],
    x=focalField, color=:fisyr,
    Guide.title("Return of Public Charities by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=50, density=false),
    Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))

  push!(plotNames, "Return_Distribution_of_Private Foundations_by_Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pf,:],
    x=focalField, color=:fisyr,
    Guide.title("Return of Private Foundations by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=50, density=false),
    Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))=#

  #filter out funds with no category
  presubdf = view(presubdf, (!ismissing).(presubdf[categoryFigureField]),:)

  #=for c::MSymbol ∈ categoriesUsed

    #filter to the category and exclude non-studied asset years
    subdf = view(presubdf, Vector{Bool}(presubdf[categoryFigureField] .== c), :)

    #we need to make sure we have a viable sample
    if size(subdf,1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return_Distribution_of_Non-Profit_($c)")
      push!(plots, plot(subdf,
        x=focalField, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))
    end

    if size(subdf[subdf[:charitytype] .== :pc,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return_Distribution_of_Public_Charities_($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
        x=focalField, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))
    end

    if size(subdf[subdf[:charitytype] .== :pf,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return_Distribution_of_Private Foundations_($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
        x=:adjnetassets, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=50, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5, ymin=0.0, ymax=1000.0)))
    end
  end=#

  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Return Graphs by Category Drawn.")

  return nothing
end

function binForFigures(data::NCCSData;
  returnYears::Vector{Int} = RETURN_YEARS_FULL,
  binOnFields::Vector{Symbol}=BIN_ON_FIELDS,
  focalField::Symbol=RETURN_FOCAL_FIELD,
  assetsField::Symbol = :adjnetassets)::DataFrame

  #pre-filter the dataframe
  local subdf::SubDataFrame
  local presubdf::SubDataFrame = view(data.df, (y::Int -> y ∈ returnYears).(data.df[:fisyr]), :)
  presubdf = view(presubdf, (!ismissing).(presubdf[focalField]), :)

  for f ∈ binOnFields #get rid of missing values
    presubdf = view(presubdf, (!ismissing).(presubdf[f]), :)
  end

  outdf::DataFrame = unique(presubdf[binOnFields])
  nrows::Int = size(outdf, 1)

  local N::Vector{MInt} = Vector{MInt}(undef, nrows)
  local equalmean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local valuemean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  #local meanm1::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  #local meanp1::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)

  #make an index for the rows- do it this way for performance
  outIndex::Dict = Dict(Tuple(outdf[i, binOnFields])=>i for i::Int ∈ 1:nrows)
  presubdfBybinOnFields::GroupedDataFrame = groupby(presubdf, binOnFields)
  for subdf::SubDataFrame ∈  presubdfBybinOnFields
    r::Int = outIndex[Tuple(subdf[1, binOnFields])]
    focal::Vector = subdf[focalField]

    #compute the appropriate statistics
    N[r] = length(focal)

    if N[r] .≥ 100 #NOTE: this is a hardcoded cutoff
      equalmean[r] = StatsBase.mean(focal)
      valuemean[r] =
        sum(processmissing(subdf[focalField] .* subdf[assetsField])) /
        sum(processmissing(subdf[(!ismissing).(subdf[focalField]), assetsField]))
      stddev[r] = StatsBase.std(focal)
      p5[r] = StatsBase.quantile(focal, .05)
      p25[r] = StatsBase.quantile(focal, .25)
      p50[r] = StatsBase.quantile(focal, 0.5)
      p75[r] = StatsBase.quantile(focal, 0.75)
      p95[r] = StatsBase.quantile(focal, 0.95)
    end

  end



  #make the dataframe
  outdf[:N] = N
  outdf[:equalmean] = equalmean
  outdf[:valuemean] = valuemean
  outdf[:stddev] = stddev #NOTE: Should this be weighted standard deviation???
  outdf[:valuemeanm1] = outdf[:valuemean] .- outdf[:stddev]
  outdf[:valuemeanp1] = outdf[:valuemean] .+ outdf[:stddev]
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95

  for f ∈ binOnFields # we want to get some sub-aggregates
    outdf[Symbol(:p5, f)] = Vector{MFloat64}(undef, nrows)
    outdf[Symbol(:p25, f)] = Vector{MFloat64}(undef, nrows)
    outdf[Symbol(:p75, f)] = Vector{MFloat64}(undef, nrows)
    outdf[Symbol(:p95, f)] = Vector{MFloat64}(undef, nrows)

    for subdf::SubDataFrame ∈ groupby(presubdf, f)
      focal = subdf[focalField]
      rows::Vector{Bool} = (outdf[f] .== subdf[1:1, f])
      outdf[rows, Symbol(:p5, f)] = StatsBase.quantile(focal, .05)
      outdf[rows, Symbol(:p25, f)] = StatsBase.quantile(focal, .25)
      outdf[rows, Symbol(:p75, f)] = StatsBase.quantile(focal, .75)
      outdf[rows, Symbol(:p95, f)] = StatsBase.quantile(focal, .95)
    end
  end
  return outdf

end

#contains the functions for writing the graphs by category
function returnFiguresBinned(data::NCCSData;
    returnYears::Vector{Int} = RETURN_YEARS_FULL,
    outputPath::String = GRAPH_PATH,
    focalField::Symbol=RETURN_FOCAL_FIELD,
    categoryFigureField::Symbol = CATEGORY_FIGURE_FIELD)


  local aggdf::DataFrame = binForFigures(data, returnYears=returnYears,
    focalField=focalField, binOnFields = [:fisyr, categoryFigureField])

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  aggdf = aggdf[((!ismissing).(aggdf[categoryFigureField])) .& ((!ismissing).(aggdf[:valuemean])), :]

  sort!(aggdf, [:fisyr, categoryFigureField])
  push!(plotNames, "Mean_Returns_by_Category")
  push!(plots, plot(aggdf,
    x=categoryFigureField, y=:valuemean, ymin = Symbol(:p5,categoryFigureField),
      ymax = Symbol(:p95,categoryFigureField),
    Stat.x_jitter(range=0.5), color=:fisyr,
    Guide.title("Mean Returns by Category"),
    Guide.xlabel("Category"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar
    #Stat.histogram(bincount=50, density=true),
    #Coord.cartesian(ymin=0.0, ymax = 3.)
  ))

  push!(plotNames, "Mean_Returns_by_Year")
  push!(plots, plot(aggdf,
    x=:fisyr, y=:valuemean,
    ymin = :p5fisyr, ymax = :p95fisyr,
    #Stat.x_jitter(range=0.5),
    color=categoryFigureField,
    Guide.title("Mean Returns by Year"),
    Guide.xlabel("Year"), Guide.ylabel("Return"),
    Geom.errorbar, Geom.line,
    #Stat.histogram(bincount=50, density=true),
    Coord.cartesian(xmin=1999, xmax=2016)))


  sort!(aggdf, [:fisyr, categoryFigureField])
  push!(plotNames, "Median_Returns_by_Category")
  push!(plots, plot(aggdf,
    x=categoryFigureField, y=:p50, ymin = :p25, ymax = :p75,
    Stat.x_jitter(range=0.5), color=:fisyr,
    Guide.title("Median Returns by Category"),
    Guide.xlabel("Category"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar
    #Stat.histogram(bincount=50, density=true),
    #Coord.cartesian(ymin=0.0, ymax = 3.)
  ))

  push!(plotNames, "Median_Returns_by_Year")
  push!(plots, plot(aggdf,
    x=:fisyr, y=:p50, ymin = :p25, ymax = :p75,
    #Stat.x_jitter(range=0.5),
    color=categoryFigureField,
    Guide.title("Median Returns by Year"),
    Guide.xlabel("Year"), Guide.ylabel("Return"),
    Geom.errorbar, Geom.line,
    #Stat.histogram(bincount=50, density=true),
    Coord.cartesian(xmin=1999, xmax = 2016)
  ))



  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Binned return graphs drawn.")

  return nothing
end

function makeDescriptiveFigures(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS)::Nothing

  Gadfly.push_theme(:default)

  #futures::Vector{Future} = Vector{Future}()

  #push!(futures, @spawn wealthFiguresByYear(data))
  #push!(futures, @spawn wealthFiguresByCategory(data))
  #push!(futures, @spawn returnFiguresByYear(data))
  #push!(futures, @spawn returnFiguresByCategory(data))
  #push!(futures, @spawn returnFiguresBinned(data))

  #(fetch).(futures)

  wealthFiguresByYear(data)
  wealthFiguresByCategory(data)
  #returnFiguresByYear(data)
  #returnFiguresByCategory(data)
  #returnFiguresBinned(data)
  industryCross(data)

  return nothing
end

#oneoff tasks
function oneOffTasks(data::NCCSData)
  metadf::DataFrame = deepcopy(data.df[data.df[:fisyr] .== 2013, [:category, :categoryname, :naics2name, :naics2]])
  dropmissing!(metadf)
  outdf::DataFrame = by(metadf, [:category, :categoryname, :naics2name, :naics2],
    [:category, :categoryname, :naics2name, :naics2] =>
      x -> (n=length(x.category)))

  CSV.write("$WORKING_PATH\\categoryMapping.csv",outdf)

  return nothing
end


#launching point for the descriptive output
function makeDescriptive(;refreshFilter::Bool = true,
    refreshDescriptiveFigures::Bool = true,
    refreshDescriptiveTables::Bool = true, writePrimaryCSV::Bool = false,
    doOneOffTasks::Bool= false)


    data::NCCSData = constructPrimaryOutput(refreshFilter=refreshFilter,
      writePrimaryCSV=writePrimaryCSV)

    refreshDescriptiveFigures && makeDescriptiveFigures(data)
    refreshDescriptiveTables && makeDescriptiveTables(data)

    doOneOffTasks && oneOffTasks(data)


    return nothing
end
