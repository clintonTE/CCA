


#these functions build the symbols for each percentile ranking covering "period" years' return
function buildpercentilesym(period::Int, lbenchmark::NSymbol = nothing;
  lretsym::NSymbol = :lreturn)::Tuple{Symbol, Symbol}

  lretstring = string(something(lretsym, "zero"))
  lbenchstring::String = isnothing(lbenchmark) ? "" : string("net", lbenchmark)

  return (Symbol("p$lretstring", lbenchstring, period, "yr"),
    Symbol("L$(period)p$lretstring", lbenchstring, period, "yr"))
end

function buildpercentilesymtype(period::Int, lbenchmark::NSymbol = nothing;
  lretsym::NSymbol = :lreturn)::Tuple{Symbol, Symbol}

  lretstring = string(something(lretsym, "zero"))
  lbenchstring::String = isnothing(lbenchmark) ? "" : string("net", lbenchmark)

  return (Symbol("p$lretstring", lbenchstring, period, "yrtype"),
    Symbol("L$(period)p$lretstring", lbenchstring, period, "yrtype"))
end


#consider making a pooled version of this
function persistenceplot(data::NCCSData, period::Int, endyear::Int, lbenchmark::NSymbol;
  outputPath::String = GRAPH_PATH, typeplots::Bool = false, pooled::Bool=true)::Nothing

  #containers for the grpahs
  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()
  local subdf::SubDataFrame

  #declare and build the symbols used to get the fields
  local percentilesym::Symbol
  local percentilesymtype::Symbol
  local Lpercentilesym::Symbol
  local Lpercentilesymtype::Symbol

  local percentilebenchsym::Symbol
  local percentilebenchsymtype::Symbol
  local Lpercentilebenchsym::Symbol
  local Lpercentilebenchsymtype::Symbol


  #get the field symbols given the period and benchmark
  (percentilesym, Lpercentilesym) = buildpercentilesym(period, lbenchmark)
  (percentilesymtype, Lpercentilesymtype) = buildpercentilesymtype(period, lbenchmark)

  local dotsize::Vector = [0.4mm]

  if pooled
    subdf = view(data.df, data.df[:fisyr] .≤ endyear, :)
    subdf = view(data.df, ((endyear .- data.df[:fisyr]) .% period) .== 0, :)
  else
    subdf = view(data.df, data.df[:fisyr] .== endyear, :)
  end
  subdf = view(subdf, (!ismissing).(subdf[percentilesym]), :)
  subdf = view(subdf, (!ismissing).(subdf[Lpercentilesym]), :)


  push!(plotNames, "$(something(lbenchmark,""))$(period)yrPerf$(endyear)$(pooled ? "-pooled" : "")")
  #=push!(plots, plot(subdf, x=Lpercentilesym, y=percentilesym,
    #size=:lnetassets,
    Guide.xlabel("$(period)yr return in $(endyear-period)"),
    Guide.ylabel("$(period)yr return in $(endyear)"),
    Geom.point, style(highlight_width=0mm, point_size=0.3mm),
    Coord.cartesian(xmin = 0.0, xmax=1.0, ymin=0.0, ymax = 1.0)))=#
  push!(plots, plot(subdf, x=Lpercentilesym, y=percentilesym,
    #size=:lnetassets,
    Guide.xlabel("Lagged $period year cumulative return percentiles"),
    Guide.ylabel("$period year cumulative return percentiles"),
    Geom.hexbin(xbincount=100, ybincount=100),
    #style(highlight_width=0mm, point_size=0.3mm),
    Coord.cartesian(xmin = 0.0, xmax=1.0, ymin=0.0, ymax = 1.0)))

  if typeplots #create plots by type
    for t ∈ collect(keys(charitytype))
      ssubdf = view(subdf, subdf[:charitytype] .== t, :)


      push!(plotNames, "$(charitytype[t])$(period)yr$(lbenchmark)Perf$(endyear)$(pooled ? "-pooled" : "")")
      push!(plots, plot(subdf, x=Lpercentilesymtype, y=percentilesymtype,
        Guide.title(plotNames[end]),
        Guide.xlabel("$(period)yr return in $(endyear-period)"),
        Guide.ylabel("$(period)yr return in $(endyear)"),
        Geom.point, size=dotsize,
        Coord.cartesian(xmin = 0.0, xmax=1.0, ymin=0.0, ymax = 1.0)))
    end
  end

  for i ∈ 1:length(plots) #write the graphs
    draw(PDF("$outputPath\\$(plotNames[i]).pdf", 9inch, 7inch), plots[i])
    #savefig(plots[i], "$outputPath\\$(plotNames[i])$(graphsuffix).$ext")
    #println("Graph $(plotNames[i]) written.")
  end

  return nothing
end

#table is complex enough to require own type
struct PersistenceTable
  persistencecolumns::Vector{Vector{Union{MFloat64,Int}}}
  rowindex::Dict
  colindex::Dict

  typefields::Vector{Symbol}
  periods::Vector{Int}
  statfields::Vector{Symbol}
  pooled::Bool
end

function PersistenceTable(data::NCCSData, endyear::Int, lbenchmark::NSymbol;
  periods::Vector{Int}=PERSISTENCE_PERIODS,  pooled::Bool = true,
  cutoffs=PERSISTENCE_CUTOFFS)::PersistenceTable

  #variables for views on the data
  local subdf::SubDataFrame
  local ssubdf::SubDataFrame

  #variables to hold the field names
  local percentilesym::Symbol
  local Lpercentilesym::Symbol

  #it is helpful to structure the rows and columns via an index
  local typefields::Vector{Symbol} = [:pc, :pf, :co, :all]

  #form a dictionary with the column index by type then period
  local colindex::Dict =
    Dict((typefields[i],periods[j])=>j + (i-1)*length(periods)
      for i ∈ 1:length(typefields)
      for j ∈ 1:length(periods))
  local colindexkeys::Vector{Tuple{Symbol, Int}} = collect(keys(colindex))
  local ncols::Int = length(colindexkeys)

  #now build the row index
  local cutofffields::Vector{Symbol} = (f::Float64 -> Symbol("p",Int(round(f*100)))).(cutoffs)
  local statfields::Vector{Symbol} = [cutofffields; :ppersistent; :N]
  local rowindex::Dict = Dict(statfields[r]=>r for r ∈ 1:length(statfields))
  local nrows = length(statfields)

  #preallocate the columns
  persistencecolumns::Vector{Vector{Union{MFloat64,Int}}} =
    (i->Vector{Union{MFloat64,Int}}(undef, nrows)).(1:ncols)

  #holds the top and bottom cutoff for each interval
  cutoffintervals::Vector{Tuple{Float64,Float64}} =
    (i::Int->i==1 ? (0.0, cutoffs[i]) : (cutoffs[i-1], cutoffs[i])).(1:length(cutoffs))


  if !pooled #the not-pooled scenario is more precise but much less powerful
    subdf = view(data.df, data.df[:fisyr] .== endyear,:) #ensure same dataset for all periods
    periodyears = (y->[endyear]).(1:(length(periods))) # same year for all periods
    for p ∈ periods
      (percentilesym, Lpercentilesym) = buildpercentilesym(p, lbenchmark)

      subdf = view(subdf, (!ismissing).(subdf[percentilesym]), :)
      subdf = view(subdf, (!ismissing).(subdf[Lpercentilesym]), :)
    end
  end

  #filter down to a single type
  for p::Int ∈ periods
    if pooled
      subdf = view(data.df, data.df[:fisyr] .≤ endyear, :)
      subdf = view(subdf, ((endyear .- data.df[:fisyr]) .% p) .== 0, :)

      (percentilesym, Lpercentilesym) = buildpercentilesym(p, lbenchmark)
      subdf = view(subdf, (!ismissing).(subdf[percentilesym]), :)
      subdf = view(subdf, (!ismissing).(subdf[Lpercentilesym]), :)
    end

    for t::NSymbol ∈ typefields

      #build the symbols referencing the appropriate data column
      if t==:all
        ssubdf = subdf
        (percentilesym, Lpercentilesym) = buildpercentilesym(p, lbenchmark)
      else
        ssubdf = view(subdf, subdf[:charitytype] .== t, :)
        (percentilesym, Lpercentilesym) = buildpercentilesymtype(p, lbenchmark)
      end

      N::Int = size(ssubdf,1)

      #will hold the total number in each quantile
      NStarts::Vector{Int} = Vector{Int}(undef, length(cutoffs))

      #will hold the total number of perisstent in eachF quantile
      NPersistents::Vector{Int} = Vector{Int}(undef, length(cutoffs))

      for c::Int ∈ 1:length(cutoffs)

        #number of NPs in the first period
        NStarts[c] = sum((ssubdf[Lpercentilesym] .≥ cutoffintervals[c][1]) .&
          (ssubdf[Lpercentilesym] .< cutoffintervals[c][2]))

        #number of NPs in the first and seocnd period
        NPersistents[c] = sum((ssubdf[Lpercentilesym] .≥ cutoffintervals[c][1]) .&
          (ssubdf[Lpercentilesym] .< cutoffintervals[c][2]) .&
          (ssubdf[percentilesym] .≥ cutoffintervals[c][1]) .&
          (ssubdf[percentilesym] .< cutoffintervals[c][2])
          )


        #fill in the values
        persistencecolumns[colindex[t,p]][rowindex[cutofffields[c]]] = NPersistents[c]/NStarts[c]
      end

      #sum the persistence values
      persistencecolumns[colindex[t,p]][rowindex[:ppersistent]] = sum(NPersistents)/N
      persistencecolumns[colindex[t,p]][rowindex[:N]] = N

    end
  end

  return PersistenceTable(persistencecolumns,
    rowindex,
    colindex,
    typefields,
    periods,
    statfields,
    pooled)

end

function persistencetable(P::PersistenceTable, endyear::Int, lbenchmark::NSymbol;
  outputpath::String = TABLE_PATH,
  decimals::Int = DECIMALS_RETURN,
  )::Nothing

  local lbenchmarkstring::String = isnothing(lbenchmark) ? "" : "$lbenchmark"

  local ncols::Int = length(P.persistencecolumns) #number of columns excluding the year
  local nrows::Int = length(P.persistencecolumns[1])

  #collect the keys reverse the dicitonaries for notational convenience
  local colindexkeys::Vector{Tuple{Symbol,Int}} = collect(keys(P.colindex))
  local rowindexkeys::Vector{Symbol} = collect(keys(P.rowindex))

  local revcolindex::Dict = Dict(c=>P.colindex[colindexkeys[c]]  for c ∈ 1:ncols)
  local revrowindex::Dict = Dict(r=>P.rowindex[rowindexkeys[r]]  for r ∈ 1:nrows)

  #column headers
  local periodfields::Vector{Symbol} = (i::Int->Symbol("R", i,"y")).(P.periods)
  local typenames = (s->s==:all ? "all non-profits" : "$(charitytype[s])").(P.typefields)
  #typenames = (s->"$s")(typenames)

  #include make the column names
  colNames::Vector{Vector{String}} = Vector{Vector{String}}()
  push!(colNames, typenames)
  push!(colNames, [((s->string(s)).(periodfields) for t ∈ 1:length(P.typefields))...;])

  widthColNames::Vector{Vector{Int}} = Vector{Vector{Int}}()
  push!(widthColNames, [length(periodfields) for i ∈ 1:length(colNames[1])])
  push!(widthColNames, ones(Int, ncols))

  #make the row names
  local descRowNames::Vector{String} =
    (s->s == :ppersistent ? "\\% Persistant" : string(s)).(P.statfields)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = (
    (i::Int)->Vector{String}(undef, ncols)).(1:nrows)

  #NOTE: P has the data in COLUMNS while descContet absorbs data in ROWS (flipped indices)
  for r::Int ∈ 1:nrows
    for c::Int ∈ 1:ncols

      #Ns are ints
      if P.statfields[r] == :N
        #lookup the collumn number from the column index for the pth column symbol
        descContent[r][c] = num2Str(P.persistencecolumns[c][r], Ints=true)
      else
        descContent[r][c] = num2Str(P.persistencecolumns[c][r], decimals)
      end
    end
  end

  outtabletitle::String = (isnothing(lbenchmark) ? "Persistence of Non-Profit Returns $(endyear)"
    : "Persistence of Non-Profit Returns Net of $(lbenchmark) $(endyear)")
  outtable::String = texTable(titleCaption="Persistence of Non-Profit Returns $(endyear)",
      caption="""See Tex File""", #caption
      colNames= colNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true,
      widthColNames=widthColNames,
      alignmentstring=string("l ", join((w->join([ " | "; ("r" for i ∈ 1:w)...;])).(widthColNames[1])))
    )

  #periodstring::String = string((s->"-$s").(columnsymbols)...)
  pooledstring::String = P.pooled ? :"(pooled)" : ""
  writeNakedTable(outtable, path=outputpath,
    outName = "persistence$(lbenchmarkstring)$(endyear)$pooledstring-$(length(P.periods)).tex")
  return nothing
end



function persistencetables(data::NCCSData; persistenceyears::Vector{Int} = PERSISTENCE_YEARS,
  lbenchmarks::Vector{NSymbol} = LBENCHMARKS)::Nothing

  for s ∈ lbenchmarks
    for y ∈ persistenceyears

      P::PersistenceTable = PersistenceTable(data, y, s, pooled=false)
      persistencetable(P, y, s)

      P = PersistenceTable(data, y, s, pooled=true)
      persistencetable(P, y, s)
    end
  end

  return nothing
end

#creates all of the persistance plots
function persistenceplots(data::NCCSData,
  persistenceperiods::Vector{Int} = PERSISTENCE_PERIODS,
  persistenceyears::Vector{Int} = PERSISTENCE_YEARS,
  lbenchmarks::Vector{NSymbol} = LBENCHMARKS)

  for p::Int ∈ persistenceperiods, y::Int ∈ persistenceyears, s ∈ lbenchmarks
    persistenceplot(data, p, y, s)
  end

  println("completed persistence plots")
end

#writes out the data as a csv if desired
function writepercentiles(data::NCCSData;
  lbenchmark::NSymbol = nothing #=:bmFF3=#, period::Int = 4, endyear::Int = 2015,
  datapath::String = DATA_PATH, pooled::Bool=true)::Nothing

  local percentilesym::Symbol
  local Lpercentilesym::Symbol

  #get the percentile symbols
  (percentilesym, Lpercentilesym) = buildpercentilesym(period, lbenchmark)
  subdf = view(data.df, data.df[:fisyr] .≤ endyear, :)

  if pooled
    subdf = view(subdf, (((endyear .- subdf[:fisyr]) .% period) .== 0), :)
  else
    subdf = view(subdf, endyear .== subdf[:fisyr], :)
  end

  subdf = view(subdf, ((s::MFloat64, Ls::MFloat64)->
      (!ismissing(s)) && (!ismissing(Ls))).(
      subdf[percentilesym], subdf[Lpercentilesym]),:)

  #build the path and file name
  path =
    "$datapath\\percentilesout_$(percentilesym)$(pooled ? "-pooled" : "").csv"
  CSV.write(path, subdf)

  return nothing
end


function analyzepersistence(data::NCCSData; percentiles2files::Bool = true)
  #make the graphs
  persistenceplots(data)
  #make the table
  persistencetables(data)
  percentiles2files && writepercentiles(data::NCCSData)
end
