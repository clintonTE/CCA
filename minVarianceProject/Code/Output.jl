
function makeconvergencegraphs(E::Estimation;
  convergencewindows::Vector{Int} = CONVERGENCE_WINDOWS,
  prefix::String = "",  outgraphpath::String=OUT_GRAPH_PATH, display::Bool = false)

  local plotnames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subplots::Vector{PlotContainer} = Vector{PlotContainer}()

  df::DataFrame = paramdf(E::Estimation)
  df[:zetaGP_abs] = (abs).(df[:zetaGP])
  df[:SGP_abs] = (abs).(df[:SGP])
  paramnames::Vector{Symbol} = [(Symbol).(names(E.Θ)); :zetaGP_abs; :SGP_abs]

  sort!(df, :pass)
  for f::Symbol ∈ [:sigma2G, :SG, :zeta2G, :zeta2P, :zetaGP_abs]
    for convergencewindow ∈ convergencewindows
      dflong::DataFrame = melt(df[df[:pass] .≤ convergencewindow,:], :pass, paramnames)
      df2Gplot::SubDataFrame = view(dflong, (s->s∈[:pass, f]).(dflong[:variable]),:)
      push!(subplots,
        plot(df2Gplot,
        x=:pass, color=:variable,y=:value,
        Scale.y_log10,
        #Guide.title(plotnames[end]),
        Guide.ylabel("Value"),Guide.xlabel(nothing),
        Guide.Theme(key_position=:none, line_width=0.2mm),
        Geom.line, Coord.cartesian(ymin=-9, ymax = -3, xmin=1, xmax=convergencewindow)))
      end

      push!(plotnames, "$(f)_convergence")
      push!(plots, vstack([subplots[1], subplots[2], subplots[3], subplots[4]]))

      for i ∈ 1:length(plots) #write the graphs
        draw(PNG("$outgraphpath\\$(prefix)_$(plotnames[i]).png", 9inch, 7inch), plots[i])
      end

      plots = Vector{PlotContainer}()
      subplots = Vector{PlotContainer}()
      plotnames = Vector{String}()
    end

# NOTE TODO: Uncomment THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for convergencewindow ∈ convergencewindows
    subdf::SubDataFrame  = view(df, df[:pass] .≤ convergencewindow, :)
    push!(subplots,
      plot(subdf,
      x=:sigma2G, y=:SG,
      Scale.y_log10, Scale.x_log10,
      Guide.title("Minvar sample and pop var i=$convergencewindow"),
      Guide.xlabel("Pop var"), Guide.ylabel("Minvar sample var"),
      Gadfly.Theme(line_width=0.2mm),
      Coord.cartesian(ymin=-9, ymax = -3, xmin=-9, xmax = -3),
      Geom.path))
    push!(subplots,
      plot(subdf,
      x=:zeta2G, y=:zeta2P,
      Scale.y_log10, Scale.x_log10,
      Guide.title("Var of Var i=$convergencewindow"),
      Guide.xlabel("Minvar Var of Sample Var"), Guide.ylabel("Test Port. Var of Sample Var"),
      Gadfly.Theme(line_width=0.2mm),
      Coord.cartesian(ymin=-8, ymax = -2, xmin=-8, xmax = -2),
      Geom.path))
  end

  push!(plotnames, "sample_pop_var_convergence")
  push!(plots,
    gridstack([subplots[1] subplots[3]; subplots[5] subplots[7]]))

  push!(plotnames, "var_of_var_convergence")
  push!(plots,
    gridstack([subplots[2] subplots[4]; subplots[6] subplots[8]]))


  for i ∈ 1:length(plots) #write the graphs
    draw(PNG("$outgraphpath\\$(prefix)_$(plotnames[i]).png", 9inch, 7inch), plots[i])
  end



  println("Convergence graphs drawn")
  return nothing
end

function makeconvergencegraphs(;outpath::String = OUT_PATH,
  convergencewindows::Vector{Int} = CONVERGENCE_WINDOWS,
  outname::String = OUT_ESTIMATE_NAME,
  prefix="")::Nothing

  makeconvergencegraphs(deserialize("$outpath\\$(prefix)_$(outname).jls"),
    prefix=prefix, convergencewindows=convergencewindows)


  return
end
