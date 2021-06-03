
#NOTE: May also want to check multiple years
const YEAR_PAIRS = [(2009,2010),(2013,2014)]

function singleYearPlots(data::NCCSData; yearPairs::Vector{Tuple{Int,Int}} = YEAR_PAIRS,
  outputPath::String = GRAPH_PATH)
  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame



  for p ∈ YEAR_PAIRS
    @assert p[2] == p[1] + 1 #only support one year difference for now
    subdf = view(data.df, data.df[:fisyr] .== p[2], :)
    subdf = view(subdf, (!ismissing).(subdf[:lpreturnbyyr]), :)
    subdf = view(subdf, (!ismissing).(subdf[:preturnbyyr]), :)

    push!(plotNames, "Performance Percentile in $(p[2]) v. $(p[1])")
    push!(plots, plot(subdf, x=:lpreturnbyyr, y=:preturnbyyr,
      #size=:lnetassets,
      Guide.xlabel("Return in $(p[1])"),Guide.ylabel("Return in $(p[2])"),
      Geom.point,
      Coord.cartesian(xmin = 0.0, xmax=1.0, ymin=0.0, ymax = 1.0)))
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(PNG("$outputPath\\$(plotNames[i]).png", 9inch, 7inch), plots[i])
    #savefig(plots[i], "$outputPath\\$(plotNames[i])$(graphsuffix).$ext")
    #println("Graph $(plotNames[i]) written.")
  end

  return nothing
end


function persistancePlots(data::NCCSData)

  singleYearPlots(data)

  println("completed persistance plots")
end
