function makefigure4(dist::DataFrame;)

  ###Common Data and Parameters

  local sdist::SubDataFrame = view(dist, dist.primary, :)

  F4_STYLE_IDX::Dict = Dict(:textsize => 10.0,
      :linesize => 0.1,
      :shapesize => 0.25,
      :xdomain => [-12, 12],
      :width => 400,
      :height => 300,
      :xaxisname=>"Payment Evt Day",
      :legendtextsize => 9.0,
      :yaxisname => "#",
      :ydomain => [0.,15.],
      :seriescolor => :black)


  ###Panel C
  F4_FX = :day
  F4_FOCAL_FIELDS = [:lnumtrd, :ldolvol]

  figure4(sdist,
    Fx = F4_FX,
    focalfields = F4_FOCAL_FIELDS,
    styleidx = F4_STYLE_IDX)

end



function figure4(dist::AbstractDataFrame;
  focalfields::Vector{Symbol} = error("focalfields is required"),
  Fx::Symbol = error("Fx is required"),
  styleidx::Dict = error("style idx is required"))

  figure4fullpath::String = "$FIGURE_PATH\\figure4-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[]).pdf"

  agg::DataFrame = aggregatebyday(v->(mean(skipmissing(v))), dist, selectcols=[Fx; focalfields])
  long::DataFrame = stack(agg, focalfields, Fx)

  #maybe add both the top dividend rate and the marginal rate?
  figure4 = makefigure4graph(long,
    Fy = :value,
    Fx = Fx,
    Fgroup = :variable,
    styleidx=styleidx)
  save(figure4fullpath, figure4)
end

function makefigure4graph(df::DataFrame;
  Fy::Symbol = error("Fy is required"),
  Fx::Symbol = error("Fx is required"),
  Fgroup = error("Fgroup is required"),
  styleidx::Dict = error("styleidx is required"))

  ###NOTE: Begin the JSON for the plot
  #some plot parameters are below

  series = @vlplot() +
    @vlplot( #first part contains the line and circle dots
      data = df,
      width=styleidx[:width],
      height=styleidx[:height],
      title=nothing,
      mark={:line,
        clip=true,
        size=styleidx[:linesize]},
      x={Fx,
        axis={title=nothing,#styleidx[:xaxisname],
          formatType = :number,
          format= ".4",
          grid=false,
          titleFontSize=styleidx[:textsize],
          tickCount=12,
          titleFontWeight=:normal,
          labels=true,
          labelOpacity = {
            condition = {test = "datum.$(Fx) % 10 == 0",
              value = 0.0},
            value = 1.0}},
        type=:quantitative,
        scale={domain=styleidx[:xdomain]}},
      y={Fy,
        axis={title=styleidx[:yaxisname],
            grid=false,
            tickCount=4,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
          scale={zero=false,
            domain=styleidx[:ydomain]}
        },
      color={value=styleidx[:seriescolor]},
      detail={Fgroup,
        type=:nominal}
    ) +
    @vlplot(
      data = df,
      mark={:point,
        clip=true,
        size=styleidx[:shapesize]},
      y={Fy,
        type=:quantitative},
      x={Fx,
        type=:quantitative,
        scale={domain=styleidx[:xdomain]}},
      shape={Fgroup,
        type=:nominal,
        legend=nothing},
      color={value=styleidx[:seriescolor]},
    ) + @vlplot( #handles horizontal text
      data={values=[
        {xloc=-11, tx="Log Dollar Volume of Trades", yloc=3},
        {xloc=2, tx="Log Number of Trades", yloc=12}]},
      mark={:text,
        size=styleidx[:textsize],
        align=:left},
      x={:xloc,
        type=:quantitative},
      y={:yloc,
        type=:quantitative},
      color = {value=:black},
      text = {:tx, type=:nominal}
    )

  return series
end
