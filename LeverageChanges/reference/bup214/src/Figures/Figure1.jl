function makefigure1(dist::DataFrame;
  Fret::Symbol = DIST_RET,
  days::Int = DIST_DAYS)

  local sdist::SubDataFrame

  generateagg(df::AbstractDataFrame) = aggregatebyday(x->std(x)*100, df, selectcols=[:day, Fret])

  ydomain::Vector{Float64} = [2.0,2.4]
  yaxistext::String = "SD (daily cross- in %)"
  figurename="figure1"

  sdist = view(dist, (r->r.primary && (!ismissing(r[Fret]))).(eachrow(dist)), :)

  makedividendeventplots(generateagg, sdist,
    Fret=Fret,
    ydomain=ydomain,
    yaxistext=yaxistext,
    figurename=figurename
    )


end

function makedividendeventplots(generateagg::Function, dist::AbstractDataFrame;
  Fret::Symbol = error("Fret is required"),
  ydomain::Vector{Float64} = [minimum(df[!,Fret]), maximum(df[!,Fret])],
  yaxistext::String = error("yaxisname is required"),
  figurename::String=error("figurename is required"),
  figurepath::String = FIGURE_PATH
  )

  local aggstd::DataFrame
  local sdist::SubDataFrame

  generatepath(panelstring::String) = (
    "$figurepath\\$(figurename)$panelstring-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[]).pdf")

  #NOTE: "All Dividend Events"
  aggstd = generateagg(dist)
  panel22 = makedividendeventplot(aggstd,
    Fret = Fret,
    yaxisname = yaxistext,
    ydomain=ydomain)
  save(generatepath("panel22"), panel22)

  #NOTE: "Low-Yield: <0.75% Dividend"
  sdist = view(dist,dist.yield .< 0.0075, :)
  aggstd = generateagg(sdist)
  panel11 = makedividendeventplot(aggstd,
    Fret = Fret,
    yaxisname = yaxistext,
    ydomain=ydomain)
  save(generatepath("panel11"), panel11)

  #NOTE: "Low-Yield: <0.75% Dividend"
  sdist = view(dist,(f::Float64-> 0.0075≤f && f<.015).(dist.yield), :)
  aggstd = generateagg(sdist)
  panel12 = makedividendeventplot(aggstd,
    Fret = Fret,
    yaxisname = yaxistext,
    ydomain=ydomain)
  save(generatepath("panel12"), panel12)

  #NOTE: "High-Yield: >1.5% Dividend"
  sdist = view(dist,0.015 .≤ dist.yield, :)
  aggstd = generateagg(sdist)
  panel21 = makedividendeventplot(aggstd,
    Fret = Fret,
    yaxisname = yaxistext,
    ydomain=ydomain)
  save(generatepath("panel21"), panel21)
end


function makedividendeventplot(dist::AbstractDataFrame;
  xaxisname::String = "Event Day",
  yaxisname::String = error("yaxisname is required"),
  Fret::Symbol = error("Fret is required"),
  ydomain::Vector{Float64} = [minimum(df[!,Fret]), maximum(df[!,Fret])])

  ###NOTE: Begin the JSON for the plot
  #some plot parameters are below
  textsize::Float64 = 12.0
  dotsize::Float64 = 30.0
  linesize::Float64 = 0.1
  rulethickness::Float64 =4.0
  Δ::Float64 = ydomain[2] - ydomain[1]
  bottomtextlocation::Float64 = ydomain[1] + Δ * .05 #heuristic
  gridwidth::Float64 = 0.75

  panel = @vlplot() +
    @vlplot( #first part contains the line and circle dots
      data = dist,
      width=400,
      height=300,
      title=nothing,
      mark={:line,
        point={filled=false,
          strokeWidth={value=linesize},
          size=dotsize},
        clip=true,
        size=linesize},
      x={:day,
        axis={title=xaxisname,
          grid=true,
          gridDash=[1,3],
          gridColor=:gray,
          gridWidth=gridwidth,
          titleFontSize=textsize,
          tickCount=5,
          titleFontWeight=:normal},
        type=:quantitative},
      y={Fret,
        axis={title=yaxisname,
            grid=true,
            gridDash=[1,3],
            gridColor=:gray,
            gridWidth=gridwidth,
            titleFontSize=textsize,
            titleFontWeight=:normal},
          type=:quantitative,
          scale={zero=false,
            domain=ydomain},
        },
      color={value=:black}
    ) + @vlplot( #handles the grey shading
      data={values=[
        {start=-12, stop=-9},
        {start=-2, stop=2}]},
      mark=:rect,
      x={:start,
        type=:quantitative},
      x2={:stop,
        type=:quantitative},
      color = {value=:black},
      opacity = {value=0.3}
    ) + @vlplot( #handles the dark center lines
      data={values=[
        {thresh=-2}, {thresh=0}, {thresh=2}]},
      mark={:rule,
        size=rulethickness},
      x={:thresh,
        type=:quantitative},
      color = {value=:black},
    )+ @vlplot( #handles the verticle text
      data={values=[
        {xloc=-11, tx="(Around and Pre-) Declarations", yloc=bottomtextlocation},
        {xloc=-8, tx="Around Declarations", yloc=bottomtextlocation}]},
      mark={:text,
        size=textsize,
        angle=270,
        align=:left},
      x={:xloc,
        type=:quantitative},
      y={:yloc,
        type=:quantitative},
      color = {value=:black},
      text = {:tx, type=:nominal}
    ) + @vlplot( #handles horizontal text
      data={values=[
        {xloc=-7, tx="Clean Cum Div", yloc=bottomtextlocation},
        {xloc=2.5, tx="Same Win Ex div", yloc=bottomtextlocation}]},
      mark={:text,
        size=textsize,
        align=:left},
      x={:xloc,
        type=:quantitative},
      y={:yloc,
        type=:quantitative},
      color = {value=:black},
      text = {:tx, type=:nominal}
    )



  return panel
end

function aggregatebyday(F::Function, dist::AbstractDataFrame; Fgroup::Symbol = :day,
  selectcols::Vector{Symbol} = names(dist), suffix::String = "")

  sdist::SubDataFrame = view(dist, :, selectcols)
  aggdf::DataFrame = aggregate(groupby(sdist, Fgroup), F)
  rename!(aggdf, (n::Symbol->Symbol(replace(string(n), "_function"=>suffix))).(names(aggdf)))

  return aggdf
end
