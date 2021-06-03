function makefigure3(cumex::DataFrame;)

  ###Common Data and Parameters

  local scumex::SubDataFrame

  F3_STYLE_IDX::Dict = Dict(:textsize => 12.0,
      :linesize => 0.1,
      :shapesize => 0.15,
      :rulethickness => 0.25,
      :gridwidth => 0.25,
      :xdomain => [1960, 2020],
      :width => 300,
      :height => 100,
      :xaxisname=>"year",
      :legendtextsize => 9.0)


  #WARNING & NOTE: In the below notation, X and Y refer to the x and y values for the plots
  #while focal and dep refer to the regression focal and dependent variables respectivelly.

  ###Panel A
  F3_PANEL_AB_FFOCAL = :yield
  F3_PANEL_AB_FY = :coef
  F3_PANEL_AB_FX = :exyear

  F3_PANEL_A_YAXIS = "Volatility Change"
  F3_PANEL_A_YDOMAIN = [-0.5,0.8]
  F3_PANEL_A_FDEP = :Dabs
  F3_PANEL_A_NAME = "panela"
  F3_PANEL_A_COLOR = :blue

  scumex = view(cumex, cumex.primary, :)
  scumex = view(scumex, completecases(view(scumex, :, [F3_PANEL_AB_FFOCAL,F3_PANEL_A_FDEP])), :)

  figure3panelab(scumex,
    Fy = F3_PANEL_AB_FY,
    Fx = F3_PANEL_AB_FX,
    Fdep = F3_PANEL_A_FDEP,
    Ffocal = F3_PANEL_AB_FFOCAL,
    yaxisname = F3_PANEL_A_YAXIS,
    ydomain=F3_PANEL_A_YDOMAIN,
    panelname = F3_PANEL_A_NAME,
    seriescolor=F3_PANEL_A_COLOR,
    styleidx = F3_STYLE_IDX
    )

  ###Panel B
  F3_PANEL_B_YAXIS = "Average Return Change"
  F3_PANEL_B_YDOMAIN = [-1.2,0.1]
  F3_PANEL_B_FDEP = :Dret
  F3_PANEL_B_NAME= "panelb"
  F3_PANEL_B_COLOR = :red


  figure3panelab(scumex,
    Fy = F3_PANEL_AB_FY,
    Fx = F3_PANEL_AB_FX,
    Fdep = F3_PANEL_B_FDEP,
    Ffocal = F3_PANEL_AB_FFOCAL,
    yaxisname = F3_PANEL_B_YAXIS,
    ydomain=F3_PANEL_B_YDOMAIN,
    panelname = F3_PANEL_B_NAME,
    seriescolor=F3_PANEL_B_COLOR,
    styleidx = F3_STYLE_IDX)

  ###Panel C
  F3_PANEL_C_YAXIS = "dividend tax rate"
  F3_PANEL_C_YDOMAIN = [0.,100.]
  F3_PANEL_C_FY = :avg_marginal
  F3_PANEL_C_FX = :year
  F3_PANEL_C_NAME = "panelc"
  F3_PANEL_C_COLOR = :black

  figure3panelc(
    Fx = F3_PANEL_C_FX,
    yaxisname = F3_PANEL_C_YAXIS,
    ydomain=F3_PANEL_C_YDOMAIN,
    panelname = F3_PANEL_C_NAME,
    seriescolor = F3_PANEL_C_COLOR,
    styleidx = F3_STYLE_IDX)

end

generatefigure3path(panelstring::String) = (
  "$FIGURE_PATH\\figure3$panelstring-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[]).pdf")

function figure3panelab(cumex::AbstractDataFrame;
  taxespath::String = TAXES_PATH,
  taxesname::String = TAXES_NAME,
  Ffocal = error("Ffocal is required"),
  Fdep = error("Fdep is required"),
  Fy = error("Fy is required"),
  Fx = error("Fx is required"),
  yaxisname = error("yaxisname is required") ,
  ydomain = error("ydomain is required"),
  panelname = error("panelname is required"),
  seriescolor::Symbol = error("seriescolor is required"),
  styleidx::Dict = error("style idx is required"))


  regs::DataFrame = cumexregressionsbyyear(cumex,
    Fdep = Fdep,
    Ffocal = Ffocal)

  #regs.coef .*=100.

  regs |> CSV.write("output\\regs$(panelname)-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[]).csv")
  #println("Running $fullfigurepath")
  #println("Regression coef: ", mean(skipmissing(regs.coef)), "\n")
  panelab = makefigure3panelabplot(regs,
    Fy = Fy,
    Fx = Fx,
    yaxisname = yaxisname,
    ydomain=ydomain,
    seriescolor=seriescolor,
    styleidx=styleidx)
  save(generatefigure3path(panelname), panelab)
end

#runs the regressions within each year to create the plots
function cumexregressionsbyyear(cumex::AbstractDataFrame;
  Fgroup = :exyear,
  Fdep::Symbol = error("Fdep is required"),
  Ffocal::Symbol = error("Ffocal is required"),
  minpointsforreg = 3)

  local Xnames::Vector{String}

  @assert sum((ismissing).(cumex[!,Ffocal])) + sum((ismissing).(cumex[!,Fdep])) == 0

  Xnames = (string).([:intercept, Ffocal])

  #this will hold the results of the regressions
  regs::DataFrame = DataFrame(exyear = unique(cumex.exyear))
  regs.coef = Vector{MFloat64}(undef, size(regs,1))
  regidx::Dict = Dict(r.exyear=>r for r ∈ eachrow(regs))

  #run the regressions and extract the coefficients
  for scumex ∈ groupby(cumex, :exyear)
    (size(scumex,1) < minpointsforreg) && continue
    r::DataFrameRow = regidx[scumex.exyear[1]]

    reg::FMLM = FMLM(scumex, Ffocal, Fdep, Xnames=Xnames, Yname="$Fdep", containsmissings=false)
    @assert Xnames[2] == "$Ffocal"
    r.coef = reg.β[2]
  end

  return regs[completecases(regs),:]
end




function makefigure3panelabplot(df::DataFrame;
  yaxisname::String = error("yaxisname is required"),
  Fy::Symbol = error("Fy is required"),
  Fx::Symbol = error("Fx is required"),
  ydomain::Vector{Float64} = error("ydomain is required"),
  seriescolor::Symbol = error("series color is required"),
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
        axis={title=nothing,#xaxisname,
          formatType = :number,
          format= ".4",
          grid=true,
          gridDash=[1,3],
          gridColor=:gray,
          gridWidth=styleidx[:gridwidth],
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
        axis={title=yaxisname,
            grid=false,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
          scale={zero=false,
            domain=ydomain}
        },
      color={value=seriescolor}
    ) + @vlplot( #handles the dark center lines
        data={values=[{thresh=0}]},
        mark={:rule,
          size=styleidx[:rulethickness],
          strokeDash = [1,3]},
        y={:thresh,
          type=:quantitative},
        color = {value=:black})

  return series
end


function figure3panelc(;
  Fx::Symbol = error("Fx is required"),
  taxespath::String = TAXES_PATH,
  taxesname::String = TAXES_NAME,
  yaxisname = error("yaxisname is required") ,
  ydomain = error("ydomain is required"),
  panelname = error("panelname is required"),
  seriescolor::Symbol = error("seriescolor is required"),
  styleidx::Dict = error("style idx is required"),
  yearrange::UnitRange{Int} = YEAR_RANGE[])


  taxes::DataFrame = CSV.File("$taxespath\\$taxesname.csv") |> DataFrame
  rename!(taxes, [:us_dividend_marg_rate=>:avg_marginal, :top_div_rate=>:top_marginal])
  select!(taxes, [:year, :avg_marginal, :top_marginal])
  taxes = stack(taxes, [:avg_marginal, :top_marginal], [:year])
  taxes.value .*= 100.
  taxes = taxes[completecases(taxes),:]
  filter!(r::DataFrameRow->r.year ∈ yearrange, taxes) #restrict to desired years


  #maybe add both the top dividend rate and the marginal rate?
  panelc = makefigure3panelcplot(taxes,
    Fy = :value,
    Fx = Fx,
    Fgroup = :variable,
    yaxisname = yaxisname,
    ydomain=ydomain,
    seriescolor=seriescolor,
    styleidx=styleidx)
  save(generatefigure3path(panelname), panelc)
end

function makefigure3panelcplot(df::DataFrame;
  yaxisname::String = error("yaxisname is required"),
  Fy::Symbol = error("Fy is required"),
  Fx::Symbol = error("Fx is required"),
  Fgroup = error("Fgroup is required"),
  ydomain::Vector{Float64} = error("ydomain is required"),
  seriescolor::Symbol = error("series color is required"),
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
          grid=true,
          gridDash=[1,3],
          gridColor=:gray,
          gridWidth=styleidx[:gridwidth],
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
        axis={title=yaxisname,
            grid=false,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
          scale={zero=false,
            domain=ydomain}
        },
      color={value=seriescolor},
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
        legend={title=nothing,
          orient=:bottom,
          labelFontSize=styleidx[:legendtextsize]}
        },
      color={value=seriescolor},
    )

  return series
end
