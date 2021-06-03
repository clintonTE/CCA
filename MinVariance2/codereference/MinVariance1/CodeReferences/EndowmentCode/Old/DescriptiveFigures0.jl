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

const BIN_ON_FIELDS = [:fisyr, :categorynamefiltered]

#const FIGURES_MIN_ASSETS = 1_500_000.

#contains the functions for writing the grpah by year
function wealthFiguresByYear(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS,
    outputPath::String = OUTPUT_PATH, focalField::Symbol=WEALTH_FOCAL_FIELD#=,
    figuresMinAssets::Float64 = FIGURES_MIN_ASSETS=#)

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  #local presubdf::SubDataFrame =  view(data.df, data.df[focalField] .≥ figuresMinAssets,:)

  for y::Int ∈ assetYears #main wealth distribution
    subdf = view(data.df, data.df[:fisyr] .== y,:) #make a view for performance reasons

    push!(plotNames, "Wealth Distribution of Non-Profits ($y)")
    push!(plots, plot(subdf, x=:adjnetassets, color=:charitytype,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=true, position=:dodge)))

    push!(plotNames, "Wealth Distribution of Non-Profits ($y)")
    push!(plots, plot(subdf,
      x=focalField, color=:categorynamefiltered,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),  Scale.x_log10,
      Geom.line, Stat.histogram(bincount=50, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0)))


    push!(plotNames, "Wealth Distribution of Private Charities ($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
      x=focalField, color=:categorynamefiltered,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
      Scale.x_log10, Geom.line,
      Stat.histogram(bincount=25, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0)))


    push!(plotNames, "Wealth Distribution of Private Foundations ($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
      x=focalField, color=:categorynamefiltered,
      Guide.title("Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
      Scale.x_log10, Geom.line,
      Stat.histogram(bincount=25, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0)))

    #=push!(plotNames, "Log Wealth Distribution of Non-Profits ($y)")
    push!(plots, plot(subdf,
    x=Symbol('l',focalfield), color=:charitytype,
      Guide.title("Log Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=true)))=#
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(SVG("$outputPath\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Wealth Graphs by Year Drawn")
  return nothing
end

#contains the functions for writing the grpah by year
function returnFiguresByYear(data::NCCSData; returnYears::Vector{Int} = RETURN_YEARS,
    outputPath::String = OUTPUT_PATH, focalField::Symbol=RETURN_FOCAL_FIELD)

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame = view(data.df, (!ismissing).(data.df[focalField]))
  presubdf = view(presubdf, (!ismissing).(presubdf[:categorynamefiltered]))

  for y::Int ∈ returnYears #main wealth distribution
    subdf = view(presubdf, presubdf[:fisyr] .== y,:) #make a view for performance reasons

    #println("\n####\nYear: $y; rows: $(size(subdf,1)),
    #  rowspf: rows: $(sum(subdf[:charitytype] .== :pf)),
    #  rowspc: rows: $(sum(subdf[:charitytype] .== :pc)))")

    push!(plotNames, "Return Distribution of Non-Profits ($y) by Type")
    push!(plots, plot(subdf, x=focalField, color=:charitytype,
      Guide.title("Return Distribution of Non-Profits ($y) by Type"),
      Guide.xlabel("Return"),Guide.ylabel("Density"),
      Geom.histogram(bincount=25, density=true, position=:dodge),
      Coord.cartesian(xmin=-0.5, xmax=0.5)))

    push!(plotNames, "Return Distribution of Non-Profits ($y) by Category")
    push!(plots, plot(subdf,
      x=focalField, color=:categorynamefiltered,
      Guide.title("Return Distribution of Non-Profits ($y) by Category"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line, Stat.histogram(bincount=25, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)))


    push!(plotNames, "Return Distribution of Private Charities ($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
      x=focalField, color=:categorynamefiltered,
      Guide.title("Return Distribution of Non-Profits ($y)"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line,
      Stat.histogram(bincount=25, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)))


    push!(plotNames, "Return Distribution of Private Foundations ($y)")
    push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
      x=focalField, color=:categorynamefiltered,
      Guide.title("Return Distribution of Non-Profits ($y)"),
      Guide.xlabel("Return"), Guide.ylabel("Amount"),
      Geom.line,
      Stat.histogram(bincount=25, density=false),
      Coord.cartesian(ymin=0.0, ymax = 1000.0, xmin=-0.5, xmax=0.5)))

    #=push!(plotNames, "Log Wealth Distribution of Non-Profits ($y)")
    push!(plots, plot(subdf,
    x=Symbol('l',focalfield), color=:charitytype,
      Guide.title("Log Wealth Distribution of Non-Profits ($y)"),
      Guide.xlabel("Wealth"),Guide.ylabel("Density"),
      Scale.x_log10,
      Geom.histogram(bincount=50, density=true)))=#
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(SVG("$outputPath\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Return Graphs by Category Drawn")
  return nothing

end

#contains the functions for writing the graphs by category
function wealthFiguresByCategory(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS,
    outputPath::String = OUTPUT_PATH,
    minPointsPerCategory::Int = MIN_POINTS_PER_CATEGORY,
    focalField::Symbol=WEALTH_FOCAL_FIELD#=,
    figuresMinAssets::Float64 = FIGURES_MIN_ASSETS=#)



  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame
  local categoriesUsed::Vector{Symbol}

  #get only the categories with valid values
  presubdf = view(data.df, ((y::Int)-> y ∈ assetYears).(data.df[:fisyr]),:)
  #presubdf = view(presubdf, presubdf[focalField] .≥ figuresMinAssets)

  push!(plotNames, "Wealth Distribution of Non-Profits by Year")
  push!(plots, plot(presubdf,
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Non-Profits by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=6.)),)

  push!(plotNames, "Wealth Distribution of Public Charities by Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pc, :],
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Public Charities by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=6.)))

  push!(plotNames, "Wealth Distribution of Private Foundations by Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pf,:],
    x=focalField, color=:fisyr,
    Guide.title("Wealth of Private Foundations by Year"),
    Guide.xlabel("Wealth"), Guide.ylabel("Density"),
    Scale.x_log10, Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=6.)))

  #filter out funds with no category
  presubdf = view(presubdf, (!ismissing).(presubdf[:categorynamefiltered]),:)
  categoriesUsed =
    unique(data.df[(!ismissing).(data.df[:categorynamefiltered]), :categorynamefiltered])

  for c::MSymbol ∈ categoriesUsed

    #filter to the category and exclude non-studied asset years
    subdf = view(presubdf, Vector{Bool}(presubdf[:categorynamefiltered] .== c), :)

    #we need to make sure we have a viable sample
    if size(subdf,1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth Distribution of Non-Profit ($c)")
      push!(plots, plot(subdf,
        x=focalField, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=25, density=false),
        Coord.cartesian(ymin=0.0, ymax = 100.0)))
    end

    if size(subdf[subdf[:charitytype] .== :pc,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth Distribution of Public Charities ($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
        x=focalField, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=25, density=false)#=,
        Coord.cartesian(ymin=0.0, ymax = 100.0))=#))
    end

    if size(subdf[subdf[:charitytype] .== :pf,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Wealth Distribution of Private Foundations ($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
        x=:adjnetassets, color=:fisyr,
        Guide.title("Wealth Distribution of Non-Profits ($c)"),
        Guide.xlabel("Wealth"), Guide.ylabel("Amount"),
        Scale.x_log10, Geom.line,
        Stat.histogram(bincount=25, density=false)#=,
        Coord.cartesian(ymin=0.0, ymax = 100.0))=#))
    end
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(SVG("$outputPath\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Wealth Graphs by Category Drawn.")

  return nothing
end

#contains the functions for writing the graphs by category
function returnFiguresByCategory(data::NCCSData;
    returnYears::Vector{Int} = RETURN_YEARS,
    outputPath::String = OUTPUT_PATH,
    minPointsPerCategory::Int = MIN_POINTS_PER_CATEGORY,
    focalField::Symbol=RETURN_FOCAL_FIELD)



  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame
  local presubdf::SubDataFrame
  local categoriesUsed::Vector{Symbol}

  #get only the categories with valid values
  presubdf = view(data.df, ((y::Int)-> y ∈ returnYears).(data.df[:fisyr]),:)
  presubdf = view(presubdf, (!ismissing).(presubdf[:return]))

  push!(plotNames, "Return Distribution of Non-Profits by Year")
  push!(plots, plot(presubdf,
    x=focalField, color=:fisyr,
    Guide.title("Return of Non-Profits by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=-0.5, xmax=0.5)))

  push!(plotNames, "Return Distribution of Public Charities by Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pc, :],
    x=focalField, color=:fisyr,
    Guide.title("Return of Public Charities by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=-0.5, xmax=0.5)))

  push!(plotNames, "Return Distribution of Private Foundations by Year")
  push!(plots, plot(presubdf[presubdf[:charitytype] .== :pf,:],
    x=focalField, color=:fisyr,
    Guide.title("Return of Private Foundations by Year"),
    Guide.xlabel("Return"), Guide.ylabel("Amount"),
     Geom.line,
    Stat.histogram(bincount=25, density=true),
    Coord.cartesian(ymin=0.0, ymax = 2., xmin=-0.5, xmax=0.5)))

  #filter out funds with no category
  presubdf = view(presubdf, (!ismissing).(presubdf[:categorynamefiltered]),:)
  categoriesUsed =
    unique(data.df[(!ismissing).(data.df[:categorynamefiltered]), :categorynamefiltered])

  for c::MSymbol ∈ categoriesUsed

    #filter to the category and exclude non-studied asset years
    subdf = view(presubdf, Vector{Bool}(presubdf[:categorynamefiltered] .== c), :)

    #we need to make sure we have a viable sample
    if size(subdf,1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return Distribution of Non-Profit ($c)")
      push!(plots, plot(subdf,
        x=focalField, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=25, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5)))
    end

    if size(subdf[subdf[:charitytype] .== :pc,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return Distribution of Public Charities ($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pc,:],
        x=focalField, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=25, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5)))
    end

    if size(subdf[subdf[:charitytype] .== :pf,:],1) ≥ minPointsPerCategory
      push!(plotNames, "categories\\Return Distribution of Private Foundations ($c)")
      push!(plots, plot(subdf[subdf[:charitytype] .== :pf,:],
        x=:adjnetassets, color=:fisyr,
        Guide.title("Return Distribution of Non-Profits ($c)"),
        Guide.xlabel("Return"), Guide.ylabel("Amount"),
         Geom.line,
        Stat.histogram(bincount=25, density=false),
        Coord.cartesian(xmin=-0.5, xmax=0.5)))
    end
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(SVG("$outputPath\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Cross-sectional Graphs by Category Drawn.")

  return nothing
end

function binForFigures(data::NCCSData;
  returnYears::Vector{Int} = RETURN_YEARS_FULL,
  binOnFields::Vector{Symbol}=BIN_ON_FIELDS,
  focalField::Symbol=RETURN_FOCAL_FIELD)::DataFrame

  #pre-filter the dataframe
  local sdf::SubDataFrame
  local ssdf::SubDataFrame
  local presubdf::SubDataFrame = view(data.df, (y::Int -> y ∈ returnYears).(data.df[:fisyr]), :)
  presubdf = view(presubdf, (!ismissing).(presubdf[focalField]), :)

  for f ∈ binOnFields #get rid of missing values
    presubdf = view(presubdf, (!ismissing).(presubdf[f]), :)
  end

  outdf::DataFrame = unique(presubdf[binOnFields])
  nrows::Int = size(outdf, 1)

  local N::Vector{MInt} = Vector{MInt}(undef, nrows)
  local mean::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local stddev::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  #local meanm1::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  #local meanp1::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p5::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p25::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p50::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p75::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)
  local p95::Vector{MFloat64} = Vector{MFloat64}(undef, nrows)



  sort!(outdf, binOnFields)

  local firstBinValue::Any = outdf[1, binOnFields[1]]
  local sdf = view(presubdf, presubdf[binOnFields[1]] .== firstBinValue, :)
  for i::Int ∈ 1:nrows

    #The below snippit provides an initial level of filtering for performance purposes
    if firstBinValue ≠ outdf[i,binOnFields[1]]
      firstBinValue = outdf[i,binOnFields[1]]
      sdf = view(presubdf, presubdf[binOnFields[1]] .== firstBinValue, :)
    end

    if length(binOnFields) > 1
      for f::Symbol ∈ binOnFields[2:end] #form the bin
        ssdf = view(presubdf, sdf[f] .== outdf[i,f], :)
      end
    else
      ssdf=sdf
    end

    focal::Vector = ssdf[focalField]

    #compute the appropriate statistics
    N[i] = length(focal)

    if N[i] .≥ 10 #NOTE: this is a hardcoded cutoff
      mean[i] = StatsBase.mean(focal)
      stddev[i] = StatsBase.std(focal)
      p5[i] = StatsBase.quantile(focal, .05)
      p25[i] = StatsBase.quantile(focal, .25)
      p50[i] = StatsBase.quantile(focal, 0.5)
      p75[i] = StatsBase.quantile(focal, 0.75)
      p95[i] = StatsBase.quantile(focal, 0.95)
    end

  end

  #make the dataframe
  outdf[:N] = N
  outdf[:mean] = mean
  outdf[:stddev] = stddev
  outdf[:meanm1] = outdf[:mean] .- outdf[:stddev]
  outdf[:meanp1] = outdf[:mean] .+ outdf[:stddev]
  outdf[:p5] = p5
  outdf[:p25] = p25
  outdf[:p50] = p50
  outdf[:p75] = p75
  outdf[:p95] = p95

  return outdf

end

#contains the functions for writing the graphs by category
function returnFiguresBinned(data::NCCSData;
    returnYears::Vector{Int} = RETURN_YEARS_FULL,
    outputPath::String = OUTPUT_PATH,
    focalField::Symbol=RETURN_FOCAL_FIELD)


  local aggdf::DataFrame = binForFigures(data, returnYears=returnYears,
    focalField=focalField, binOnFields = [:fisyr, :categorynamefiltered])

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  aggdf = aggdf[((!ismissing).(aggdf[:categorynamefiltered])) .& ((!ismissing).(aggdf[:mean])), :]

  sort!(aggdf, [:fisyr, :categorynamefiltered])
  push!(plotNames, "Mean Returns by Category")
  push!(plots, plot(aggdf,
    x=:categorynamefiltered, y=:mean, ymin = :meanm1, ymax = :meanp1,
    Stat.x_jitter(range=0.5), color=:fisyr,
    Guide.title("Mean Returns by Category"),
    Guide.xlabel("Category"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar
    #Stat.histogram(bincount=50, density=true),
    #Coord.cartesian(ymin=0.0, ymax = 3.)
  ))

  push!(plotNames, "Mean Returns by Year")
  push!(plots, plot(aggdf,
    x=:fisyr, y=:mean, ymin =:meanm1, ymax = :meanp1,
    Stat.x_jitter(range=0.5), color=:categorynamefiltered,
    Guide.title("Mean Returns by Year"),
    Guide.xlabel("Year"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar,
    #Stat.histogram(bincount=50, density=true),
    Coord.cartesian(xmin=1999, xmax=2016)))


  sort!(aggdf, [:fisyr, :categorynamefiltered])
  push!(plotNames, "Median Returns by Category")
  push!(plots, plot(aggdf,
    x=:categorynamefiltered, y=:p50, ymin = :p25, ymax = :p75,
    Stat.x_jitter(range=0.5), color=:fisyr,
    Guide.title("Median Returns by Category"),
    Guide.xlabel("Category"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar
    #Stat.histogram(bincount=50, density=true),
    #Coord.cartesian(ymin=0.0, ymax = 3.)
  ))

  push!(plotNames, "Median Returns by Year")
  push!(plots, plot(aggdf,
    x=:fisyr, y=:p50, ymin = :p25, ymax = :p75,
    Stat.x_jitter(range=0.5), color=:categorynamefiltered,
    Guide.title("Median Returns by Year"),
    Guide.xlabel("Year"), Guide.ylabel("Return"),
    Geom.point, Geom.errorbar,
    #Stat.histogram(bincount=50, density=true),
    Coord.cartesian(xmin=1999, xmax = 2016)
  ))



  for i ∈ 1:length(plots) #write the graphs
    draw(SVG("$outputPath\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    #println("Graph $(plotNames[i]) written.")
  end

  println("Binned return graphs drawn.")

  return nothing
end



function makeDescriptiveFigures(data::NCCSData; assetYears::Vector{Int} = ASSET_YEARS)::Nothing

  Gadfly.push_theme(:default)

  futures::Vector{Future} = Vector{Future}()

  #push!(futures, @spawn wealthFiguresByYear(data))
  #push!(futures, @spawn wealthFiguresByCategory(data))
  #push!(futures, @spawn returnFiguresByYear(data))
  #push!(futures, @spawn returnFiguresByCategory(data))
  push!(futures, @spawn returnFiguresBinned(data))

  (fetch).(futures)

  return nothing
end

#launching point for the descriptive output
function makeDescriptive(;refreshFilter::Bool = true,
    refreshDescriptiveFigures::Bool = true,
    refreshDescriptiveTables::Bool = true, writePrimaryCSV::Bool = false)


    data::NCCSData = constructPrimaryOutput(refreshFilter=refreshFilter,
      writePrimaryCSV=writePrimaryCSV)

    refreshDescriptiveFigures ? makeDescriptiveFigures(data) : nothing
    refreshDescriptiveTables ? makeDescriptiveTables(data) : nothing


    return nothing
end
