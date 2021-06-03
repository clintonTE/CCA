

#true dodger blue = colorant"#005A9C"
graphparams = Dict(
  :color => [colorant"black", #=colorant"#005C5C"=#colorant"dodgerblue3",],
  :height => "18cm",
  :width => "24cm")



    #=comombenchmarkcols = [
      :lG_WXplmomentum_LP366d,
      :lGP_WXplmomentum_LP366d,
      :lGSP_WXplmomentum_LP366d,]=#
 comomgraphbenchmark = "comom"
 comomgraphmeasure = "lGP_WXplmomentum_LP366d"

function comomentumgraphs(; figurepath = PARAM[:figurepath],
    measurefilename = PARAM[:tabcomommeasurefilename],
    rid = PARAM[:tabcomomrid],
    resultsname = PARAM[:tabcomomresultsname],
    analysispath = PARAM[:analysispath],
    graphname = PARAM[:figcomomgraphname]
    #=benchmarkcols = comombenchmarkcols=#)


  comom = CSV.File("$analysispath\\$measurefilename\\$rid\\$(rid)_pl_series.csv") |> DataFrame


  #dt = cg.date
  #benchmark = cg.comom


  #form a complete cases ar
  cg = comom[completecases(comom, ["date"; comomgraphmeasure; comomgraphbenchmark]),:]


  #cg.zmeasure = cg[!, measure] ./ std(cg[!, measure])
  #cg.zbenchmark = cg[!, benchmark] ./ std(cg[!, benchmark])
  ztransform(v) = (v .- mean(v)) ./ std(v)
	#ztransform(v) = (v .- mean(v)) ./ (maximum(v)-minimum(v))
  cg.zmeasure = cg[!, comomgraphmeasure] |> ztransform
  cg.zbenchmark = cg[!, comomgraphbenchmark] |> ztransform




  comomgrowth = @pgf Axis(
  {
    height=graphparams[:height],
    width=graphparams[:width],
    date_coordinates_in = "x",
    xlabel="date",
    xticklabel={"\\year"},
    x_tick_label_style = "{rotate=90}",
    ylabel="growth",
	#ymin=-4.0,
	#ymax=4.0,
  },
  PlotInc({
    no_marks,
    color=graphparams[:color][1],
    style={thick},
    },
  Table(
    {
      x="date",
      y="zmeasure",
    },
    cg,
  )),
	LegendEntry("Measure"),
  PlotInc({
    no_marks,
    color=graphparams[:color][2],
    dashed,
    style={thick},
    },
  Table(
    {
      x="date",
      y="zbenchmark",
    },
    cg,
  )),
	LegendEntry("Comomentum"),
  )


  ######### Scatter- doesn't work well
  #cgstack = stack(cg, [:zmeasure; :zbenchmark], "date")
  #cg.variable = string.(cg.variable)
  #=comomscatter = @pgf Axis(
  {
    height=graphparams[:height],
    width=graphparams[:height],
    xlabel="measure",
    ylabel="comom",
		axis_equal,
		xmin=-3.0,
		xmax=3.0,
		ymin=-3.0,
		ymax=3.0,
    "axis on top"=false,
		"axis x line"="middle",
		"axis y line"="middle",
  },
  Plot(
  {
    #scatter,
    "only marks",
		mark = "o",
    #"scatter src"="explicit symbolic",
    #"scatter/classes"=
    #color=graphparams[:color][1],
    #==#
  },
  Table(
  {
    x="zmeasure",
    y="zbenchmark",
  },
  cg)))=#


	############
	icg = cg[:, [:date, :zmeasure, :zbenchmark]]

	@assert issorted(icg, :date)
	icg = vcat(DataFrame(date=[icg.date[1] - Month(1)], zmeasure=[0.0], zbenchmark=[0.0]), icg)
	icg.izmeasure = cumsum(icg.zmeasure)
	icg.izbenchmark = cumsum(icg.zbenchmark)
	comomint = @pgf Axis(
  {
    height=graphparams[:height],
    width=graphparams[:width],
    date_coordinates_in = "x",
    xlabel="date",
    xticklabel={"\\year"},
    x_tick_label_style = "{rotate=90}",
    ylabel="Integrated Growth Z-score",
		#ymode="log",
  },
  PlotInc({
    no_marks,
    color=graphparams[:color][1],
    style={thick},
    },
  Table(
    {
      x="date",
      y="izmeasure",
    },
    icg,
  )),
	LegendEntry("Measure"),
  PlotInc({
    no_marks,
    color=graphparams[:color][2],
    dashed,
    style={thick},
    },
  Table(
    {
      x="date",
      y="izbenchmark",
    },
    icg,
  )),
	LegendEntry("Comomentum"),
  )



  #display(p)
  pgfsave("$figurepath\\$(graphname)_growth.svg", comomgrowth)
  #pgfsave("$figurepath\\$(graphname)_scatter.svg", comomscatter)
  pgfsave("$figurepath\\$(graphname)_int.svg", comomint)
  #display("image/png", p)


end



fundgraphbenchmarks = Dict(
  :umdadj=>:lG_A_aumequityXabsw_F_ff3m_umd_LP366d,
  :measureadj => :lG_A_aumequityXabsw_WXmomentum121_LP366d,
   )


function fundgraphs(;)

  for fundtype ∈ [:hf,:mf], fundadjtype = [:umdadj, :measureadj]
	  fundgraph(;fundtype, fundadjtype)
  end
end

function fundgraph(;
	fundtype::Symbol,
	fundadjtype::Symbol,
	fundgraphbenchmark::String=fundgraphbenchmarks[fundadjtype]|>string,
	fundgraphmeasure::String=fundmeasurecol|>string,
	graphname=PARAM[:figfundgraphname],
	figurepath = PARAM[:figurepath],
  measurefilename = PARAM[:tabfundmeasurefilename],
  rid = PARAM[:tabfundrid],
  resultsname = PARAM[:tabfundresultsname],
  analysispath = PARAM[:analysispath],)

	fund = CSV.File("$analysispath\\$measurefilename\\$rid\\$(rid)_$(fundtype)_series.csv") |> DataFrame


	fg = fund[completecases(fund, ["date"; fundgraphmeasure; fundgraphbenchmark]),:]
  benchmarklabel = verificationlabels[fundadjtype][Symbol(fundgraphbenchmark)]



	ztransform(v, vol=std(v)) = (v .- mean(v)) ./ vol

	fg.zmeasure = fg[!, fundgraphmeasure] |> ztransform
	fg.zbenchmark = fg[!, fundgraphbenchmark] |> ztransform




	fundgrowth = @pgf Axis(
	{
	  height=graphparams[:height],
	  width=graphparams[:width],
	  date_coordinates_in = "x",
	  xlabel="date",
	  xticklabel={"\\year"},
	  x_tick_label_style = "{rotate=90}",
	  ylabel="growth",
	  #ymin=-4.0,
	  #ymax=4.0,
	},
	PlotInc({
	  no_marks,
	  color=graphparams[:color][1],
	  style={thick},
	  },
	Table(
	  {
		x="date",
		y="zmeasure",
	  },
	  fg,
	)),
	  LegendEntry("Measure"),
	PlotInc({
	  no_marks,
	  color=graphparams[:color][2],
	  dashed,
	  style={thick},
	  },
	Table(
	  {
		x="date",
		y="zbenchmark",
	  },
	  fg,
	)),
	  LegendEntry(latexstring(benchmarklabel)),
	)



	  ############
  ifg = fg[:, [:date, :zmeasure, :zbenchmark]]

  @assert issorted(ifg, :date)
  ifg = vcat(DataFrame(date=[ifg.date[1] - Month(1)], zmeasure=[0.0], zbenchmark=[0.0]), ifg)
  ifg.izmeasure = cumsum(ifg.zmeasure)
  ifg.izbenchmark = cumsum(ifg.zbenchmark)
	fundint = @pgf Axis(
	{
	  height=graphparams[:height],
	  width=graphparams[:width],
	  date_coordinates_in = "x",
	  xlabel="date",
	  xticklabel={"\\year"},
	  x_tick_label_style = "{rotate=90}",
	  ylabel="Integrated Growth Z-score",
		  #ymode="log",
	},
	PlotInc({
	  no_marks,
	  color=graphparams[:color][1],
	  style={thick},
	  },
	Table(
	  {
		x="date",
		y="izmeasure",
	  },
	  ifg,
	)),
	  LegendEntry("Measure"),
	PlotInc({
	  no_marks,
	  color=graphparams[:color][2],
	  dashed,
	  style={thick},
	  },
	Table(
	  {
		x="date",
		y="izbenchmark",
	  },
	  ifg,
	)),
	  LegendEntry(latexstring(benchmarklabel)),
	)



	#display(p)
	pgfsave("$figurepath\\$(graphname)-$fundadjtype-$(fundtype)_growth.svg", fundgrowth)
	#pgfsave("$figurepath\\$(graphname)_scatter.svg", comomscatter)
  pgfsave("$figurepath\\$(graphname)-$fundadjtype-$(fundtype)_int.svg", fundint)
	#display("image/png", p)

end
