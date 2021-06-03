function makefigure1(dist::DataFrame;
  Fret::Symbol = DIST_RET,
  figurepath::String = FIGURE_PATH)

  local sdist::SubDataFrame
  local ssdist::SubDataFrame
  local aggstd::DataFrame

  generatepath(panelstring::String) = "$figurepath\\figure1$panelstring-$(REPLICATION_TYPE[]).pdf"
  generateagg(df::AbstractDataFrame) = aggregatebyday(x->std(x)*100, df, selectcols=[:day, Fret])
  #sdist::SubDataFrame = view(dist, completecases(dist[!, [Fret, :day]]), :)


  #display(allstdXday)

  #=panel22 = plot(allstdXday, x=:day, y=Fret, Guide.xlabel("Day"), Guide.ylabel("Std"),
    Geom.line, Geom.point, Guide.title("All Dividend Events"))
  draw(PDF(figurepath, 9inch, 7inch),panel22)
  allstdXday |> CSV.write("output\\allstdXday-$(REPLICATION_TYPE[]).csv")=#

  #=panel22 = plot(allstdXday[!,Fret], allstdXday[!,:day],
    xlabel="Event Day",
    marker=:dot,
    ylabel = "SD/day in \\%")
  savefig(panel22, figurepath)=#

  #=panel22 = allstdXday|> @vlplot(
    title="All Dividend Events",
    mark={:line,
      point={filled=false}},
    x={:day,
      title="Event Day",
      type=:quantitative},
    y={Fret,
      axis={title="SD in %",
        tickCount=5,
        grid=false},
      type=:quantitative,
      scale={zero=false},
      },
    )=#
  sdist = view(dist, (r->r.primary && (!ismissing(r[Fret]))).(eachrow(dist)), :)
  aggstd = generateagg(sdist)
  panel22 = makedividendeventplot(aggstd,
    Fret = Fret,
    title = "All Dividend Events",
    yaxisname = "SD in %")
  save(generatepath("panel22"), panel22)

  ssdist = view(sdist,sdist.yield .< 0.0075, :)
  aggstd = generateagg(ssdist)
  panel11 = makedividendeventplot(aggstd,
    Fret = Fret,
    title = "Low-Yield: <0.75% Dividend",
    yaxisname = "SD in %")
  save(generatepath("panel11"), panel11)

  ssdist = view(sdist,(f::Float64-> 0.0075≤f && f<.015).(sdist.yield), :)
  aggstd = generateagg(ssdist)
  panel12 = makedividendeventplot(aggstd,
    Fret = Fret,
    title = "Mid-Yield: 0.75 to 1.5% Dividend",
    yaxisname = "SD in %")
  save(generatepath("panel12"), panel12)

  ssdist = view(sdist,0.015 .≤ sdist.yield, :)
  aggstd = generateagg(ssdist)
  panel21 = makedividendeventplot(aggstd,
    Fret = Fret,
    title = "Mid-Yield: >1.5% Dividend",
    yaxisname = "SD in %")
  save(generatepath("panel21"), panel21)
end

function makedividendeventplot(df::AbstractDataFrame;title::String = error("title is required"),
  xaxisname::String = "Event Day",
  yaxisname::String = error("yaxisname is required"),
  Fret::Symbol = error("Fret is required"))

  panel = df|> @vlplot(
    title=title,
    mark={:line,
      point={filled=false}},
    x={:day,
      axis={title=xaxisname,
        grid=true,
        gridDash=[2,2]},
      type=:quantitative},
    y={Fret,
      axis={title=yaxisname,
        grid=true,
        gridDash=[2,2]},
      type=:quantitative,
      scale={zero=false},
      },
    )

  return panel
end

function aggregatebyday(F::Function, dist::AbstractDataFrame; Fgroup::Symbol = :day,
  selectcols::Vector{Symbol} = names(dist), suffix::String = "")

  sdist::SubDataFrame = view(dist, :, selectcols)
  aggdf::DataFrame = aggregate(groupby(sdist, Fgroup), F)
  names!(aggdf, (n::Symbol->Symbol(replace(string(n), "_function"=>suffix))).(names(aggdf)))

  return aggdf
end
