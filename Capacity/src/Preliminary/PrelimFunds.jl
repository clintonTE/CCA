#graphs of assets

function formpreliminaryfunds()

  local funds::DataFrame

  #form the data to grph
  funds = readprocessfundsdata()
  series = aggregatefunds(funds)

  PRELIM_STYLE_IDX::Dict = Dict(:textsize => 12.0,
      :linesize => 0.5,
      :shapesize => 0.25,
      #:xdomain => [1926, 2018],
      :width => 400,
      :height => 300,
      :widthscatter => 1600,
      :heightscatter => 1200,
      :xaxisname=>"Year",
      :legendtextsize => 9.0,
      :xaxisscattername=>"10ˣ_assets",
      :pointsize=>8.0,
      #:ydomain => [-1.0,2.0]#,
      #=:seriescolor => :black=#)

  #println(describe(series))
  graphfundassets(series,
    styleidx=PRELIM_STYLE_IDX,
    graphname="fundassets",
    yaxisname = "Assets")

  graphscatterfunds(funds,
    styleidx=PRELIM_STYLE_IDX,
    graphname="scattercash",
    yaxisname = "% cash",
    Fy = :wpercash)

    graphscatterfunds(funds,
      styleidx=PRELIM_STYLE_IDX,
      graphname="scattereq",
      yaxisname = "% eq",
      Fy = :wpercom)

end


function graphfundassets(
  series::AbstractDataFrame;
  outputfolder::String = PARAM[:figurepath],
  styleidx::Dict = error("styleidx is required"),
  Fx=:date,
  Fy=:value,
  Fgroup=:variable,
  graphname::String = error("graph name is required"),
  yaxisname::String = error("y-axis name is required"))

 #  println(series)

  panel = @vlplot() +
    @vlplot( #first part contains the line and circle dots
      data = series,
      width=styleidx[:width],
      height=styleidx[:height],
      title=nothing,
      mark={:line,
        clip=true,
        size=styleidx[:linesize]},
      x={Fx,
        axis={title=styleidx[:xaxisname],
          grid=false,
          titleFontSize=styleidx[:textsize],
          tickCount=11,
          titleFontWeight=:normal,
          labels=true},
        type=:temporal,
        },
      y={Fy,
        axis={title=yaxisname,
            grid=false,
            tickCount=4,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
          scale={type=:log},
        },
      color={Fgroup,
        type=:nominal,
        legend = {
          title=nothing,
          orient="bottom-right"
          }
        }
      )

  save("$outputfolder\\$graphname.pdf", panel)
end

#read in the appropriate file
function readprocessfundsdata(;prelimfundspath::String = "$(PARAM[:preliminarypath])\\$(PARAM[:prelimfunddata])",
  prelimfundretspath::String = "$(PARAM[:preliminarypath])\\$(PARAM[:prelimfundretsdata])")

  local altwords::Vector{String} = PARAM[:prelimaltwords]
  local funds::DataFrame
  local fundrets::DataFrame

  #scrubs the names and dates
  function initialclean(df::DataFrame)
    rename!(df, (cleanname).(names(df)))
    df.date = (i::Int->Date("$i", PARAM[:wrdsdateformat])).(df.caldt)
    df.year = (d->year(d)).(df.date)

    return df
  end


  if PARAM[:refreshpreliminaryfunds]
    #first get the summery data
    funds = CSV.read("$(prelimfundspath).csv") |> DataFrame
    funds = initialclean(funds)
    sort!(funds, [:year, :crspfundno, :date])
    funds.tokeep = trues(size(funds,1))
    #dedup funds

    for sfunds ∈ groupby(funds, [:year, :crspfundno])
      if size(sfunds,1) > 1
        sfunds.tokeep .= false
        sfunds.tokeep[end] = true
      end
    end
    funds = funds[funds.tokeep,:]
    #select!(funds, [:year, :crspfundno, :lipperclassname, :lipperobjname])
    size(unique(funds[!,[:year,:crspfundno]]),1) == size(funds,1) || error(
      "non-unique fund-years found! $(size(unique(funds[!,[:year,:crspfundno]]),1) - size(funds,1))")


    #put an indicator for if the fun is an alt
    funds.alt = falses(size(funds,1))
    funds.lipperclassname .= (s::MString->ismissing(s) ? missing : lowercase(s)).(funds.lipperclassname)
    funds.lipperobjname .= (s::MString->ismissing(s) ? missing : lowercase(s)).(funds.lipperobjname)

    for r ∈ eachrow(funds)
      r.alt = (!ismissing(r.lipperclassname)) && sum((s->occursin(s,r.lipperclassname)).(altwords)) > 0
      r.alt && continue
      r.alt = (!ismissing(r.lipperobjname)) && sum((s->occursin(s,r.lipperobjname)).(altwords)) > 0
    end

    #get the assets data
    fundrets = CSV.read("$(prelimfundretspath).csv") |> DataFrame
    fundrets = initialclean(fundrets)
    #println(describe(fundrets))
    fundrets.assets = (s::MString->
      ismissing(s) ? missing : parseormissing(Float64, s)).(fundrets.mtna)
    fundrets = fundrets[(!ismissing).(fundrets.assets),:]
    fundrets = fundrets[fundrets.assets .> 0,:]
    select!(fundrets, [:year, :assets, :crspfundno])

    #take the average of assets for each fund-year
    aggfunc(v::AbstractVector{<:MFloat64}) = mean(skipmissing(v))
    fundretss = groupby(fundrets, [:year, :crspfundno])
    fundrets = combine(fundretss, valuecols(fundretss) .=> aggfunc)
    rename!(fundrets, (s::String-> s=>nofunction(s, "_aggfunc")).(names(fundrets)))
    fundrets.date = (y->Date(y,12,31)).(fundrets.year)


    #merge the assets with the descriptive data
    rename!(funds, :date=>:rawdate)
    funds = innerjoin(funds, fundrets, on=[:year, :crspfundno])

    #winsorize percentage data
    @inline winper(::Missing, any...) = missing
    @inline winper(per::Real, lo::Real = -100., hi::Real = 100.) = min(hi, max(per, lo))

    funds.l10assets = (log10).(funds.assets .* 1_000_000.)
    funds.wpercash = (winper).(funds.percash)
    funds.wpercom = (x->winper(x,-100,200)).(funds.percom)


    serialize("$(prelimfundspath).jls", funds)
  else
    funds = deserialize("$(prelimfundspath).jls")
  end


    #display(funds[funds.crspfundno .== 25619,:])
  return funds
end

#adds the assets across funds
function aggregatefunds(funds::AbstractDataFrame)

  local series::DataFrame

  aggfunc(v::AbstractVector{<:MFloat64}) = sum(skipmissing(v))
  seriess::GroupedDataFrame = groupby(view(funds, :, [:date, :assets]), :date)
  series = combine(seriess, valuecols(seriess) .=> aggfunc)

  #fix the names
  rename!(series, (s::String-> s=>nofunction(s, "_aggfunc")).(names(series)))
  rename!(series, :assets=>:total_assets)

  #repeat, but this time only for alt funds
  sfunds = view(funds, funds.alt, [:date, :assets, :alt])
  sfundss::GroupedDataFrame = groupby(view(sfunds, :, [:date, :assets]), :date)
  seriesalts = combine(sfundss, valuecols(sfundss) .=> aggfunc)
  rename!(seriesalts, (s::String-> s=>nofunction(s, "_aggfunc")).(names(seriesalts)))
  rename!(seriesalts, :assets=>:alt_assets)

  #join and melt the two aggregatations
  series = leftjoin(series, seriesalts, on=:date)
  series = stack(series, Not(:date))
  #println(describe(series))

  return series
end

#makes a scatterplot of fund assets
function graphscatterfunds(
  funds::AbstractDataFrame;
  outputfolder::String = PARAM[:figurepath],
  styleidx::Dict = error("styleidx is required"),
  Fx=:l10assets,
  Fy=:error("Fy is required"),
  Fgroup=:variable,
  graphname::String = error("graph name is required"),
  yaxisname::String = error("y-axis name is required"))

  #println(series)

  panel = @vlplot() +
    @vlplot( #first part contains the line and circle dots
      data = funds,
      width=styleidx[:widthscatter],
      height=styleidx[:heightscatter],
      title=nothing,
      mark={:circle,
      size=styleidx[:pointsize]},
      x={Fx,
        axis={title=styleidx[:xaxisscattername],
          grid=false,
          titleFontSize=styleidx[:textsize],
          #tickCount=11,
          titleFontWeight=:normal,
          labels=true},
        type=:quantitative,
        scale={zero=false}
        },
      y={Fy,
        axis={title=yaxisname,
            grid=false,
#            tickCount=4,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
        },
      )

  save("$outputfolder\\$graphname.png", panel)
end
