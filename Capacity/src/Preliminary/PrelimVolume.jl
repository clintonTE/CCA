#generate some preliminary graphs of beta

#-----> primary entry point for preliminary beta
function formpreliminarybeta()

  local beta::DataFrame
  local crsp::DataFrame

  #form the data to grph
  validpermnos = unique(beta.permno)
  crsp = prepcrsp(crsptype = :monthly,
    refreshcrsp=PARAM[:refreshpreliminarybeta],
    datapath=PARAM[:preliminarypath],
    validpermnos = validpermnos,
    prefix = "prelimvolume",
    crspcolumns = [:permno, :date, :price, :volume])::DataFrame

  beta = join(beta, crsp, on=[:date, :permno])

  series = aggregatebetadata(beta)
  PRELIM_STYLE_IDX::Dict = Dict(:textsize => 12.0,
      :linesize => 0.5,
      :shapesize => 0.25,
      #:xdomain => [1926, 2018],
      :width => 400,
      :height => 300,
      :xaxisname=>"Year",
      :legendtextsize => 9.0,
      :ydomain => [-1.0,2.0]#,
      #=:seriescolor => :black=#)

  #println(describe(series))
  graphbetadata(series,
    styleidx=PRELIM_STYLE_IDX,
    graphname="loadings",
    yaxisname = "Beta")

  graphbetadata(series,
    Fy=:valueabs,
    styleidx=PRELIM_STYLE_IDX,
    graphname="absloadings",
    yaxisname = "abs(Beta)")

end


function graphbetadata(
  series::AbstractDataFrame;
  outputfolder::String = FIGURE_PATH,
  styleidx::Dict = error("styleidx is required"),
  Fx=:date,
  Fy=:value,
  Fgroup=:variable,
  graphname::String = error("graph name is required"),
  yaxisname::String = error("y-axis name is required"))

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
        axis={axis=styleidx[:xaxisname],
          #formatType = "date",
          #format= "yyy",
          grid=false,
          titleFontSize=styleidx[:textsize],
          tickCount=11,
          titleFontWeight=:normal,
          labels=true},
          #=labelOpacity = {
            condition = {test = "datum.$(Fx) % 10 == 0",
              value = 0.0},
            value = 1.0}},=#
        type=:temporal,
        #scale={domain=styleidx[:xdomain]}
        },
      y={Fy,
        axis={title=yaxisname,
            grid=false,
            tickCount=4,
            titleFontSize=styleidx[:textsize],
            titleFontWeight=:normal},
          type=:quantitative,
          scale={zero=false,
            domain=styleidx[:ydomain]}
        },
      color={Fgroup,
        type=:nominal}
      )

  save("$outputfolder\\$graphname.pdf", panel)
end

#read in the appropriate file
function readprocessbetadata(;prelimbetapath::String = "$PRELIMINARY_PATH\\$PRELIM_BETA_DATA",
  refreshpreliminarybeta::Bool = error("refreshpreliminarybeta is required"))

  local beta::DataFrame
  local prelimbeta

  if refreshpreliminarybeta
    beta = CSV.read("$(prelimbetapath).csv") |> DataFrame
    names!(beta, (cleanname).(names(beta)))
    beta.date = (i::Int->Date("$i", WRDS_DATE_FORMAT)).(beta.date)
    serialize("$(prelimbetapath).jls", beta)
  else
    beta = deserialize("$(prelimbetapath).jls")
  end

  return beta
end

function aggregatebetadata(beta::AbstractDataFrame;
    aggfields::Vector{Symbol} = [:bmkt, :bsmb, :bhml, :bumd])
  local series::DataFrame
  local seriesabs::DataFrame

  #create the aggregated dataframe by hand
  sbeta::SubDataFrame = view(beta, :, [:date; aggfields; :mktcap])
  series = DataFrame(date = unique(beta.date))
  series.bmkt = Vector{MFloat64}(undef, size(series,1))
  series.bsmb = similar(series.bmkt)
  series.bhml = similar(series.bmkt)
  series.bumd = similar(series.bmkt)
  series.absbmkt = similar(series.bmkt)
  series.absbsmb = similar(series.bmkt)
  series.absbhml = similar(series.bmkt)
  series.absbumd = similar(series.bmkt)

  seriesrows::Dict = Dict(r.date=>r for r ∈ eachrow(series))

  ssbetas = groupby(sbeta, :date)
  @mpar for i::Int ∈ 1:length(ssbetas)
    ssbeta::SubDataFrame = ssbetas[i]
    r::DataFrameRow = seriesrows[ssbeta.date[1]]

    #aggregate each of the betas
    for f ∈ aggfields
      s3beta = view(ssbeta, completecases(ssbeta[!, [f, :mktcap]]), [f, :mktcap])
      (size(s3beta,1) == 0) && continue
      mktcap::Float64 = sum(s3beta.mktcap)
      r[f] = sum(s3beta.mktcap .* s3beta[!,f]) / mktcap
      r[Symbol(:abs, string(f))] = sum(s3beta.mktcap .* (abs).(s3beta[!,f])) / mktcap
    end
  end

  println(describe(series))
  #create a long series for each of the two graphs
  templong1 = melt(series[!,[:date; aggfields]], :date)

  #now create a long df for  the absolute values
  select!(series, Not(aggfields))
  rename!(series, (f::Symbol->Symbol(:abs, string(f))=>f).(aggfields))
  templong2 = melt(series[!,names(series)], :date, value_name=:valueabs)

  #join the two long dfs
  series = join(templong1, templong2, on=[:date, :variable])
  series = series[completecases(series),:]



  return series
end


#averages teh betas across permnos
function aggregatebetadataequalweighted(beta::AbstractDataFrame)
  local series::DataFrame
  local seriesabs::DataFrame

  sbeta::SubDataFrame = view(beta, :, [:date, :bmkt, :bsmb, :bhml, :bumd])
  println(describe(sbeta))
  aggfunc(v::AbstractVector{<:MFloat64}) = mean(skipmissing(v))
  series = aggregate(groupby(sbeta, :date), aggfunc)

  #fix the names
  rename!(series, (s::Symbol-> s=>nofunction(s, "_aggfunc")).(names(series)))

  #melt to graph
  series = melt(series, :date)

  #make the absolue value version
  aggfuncabs(v::AbstractVector{<:MFloat64}) = mean(skipmissing((abs).(v)))
  seriesabs = aggregate(groupby(sbeta, :date), aggfuncabs)
  rename!(seriesabs, (s::Symbol-> s=>nofunction(s, "_aggfuncabs")).(names(seriesabs)))
  seriesabs = melt(seriesabs, :date, value_name=:valueabs)

  series = join(series, seriesabs, on=[:date, :variable])


  return series
end
