#TODO: Fix the cross-tab table
function createquantiles!(df::AbstractDataFrame, source::Symbol; target=Symbol(:q, source),
  groupfields::Union{Vector{Symbol}, Nothing} = [:fisyr])

  if target ∉ names(df)
    df[target] = Vector{Union{Missing,eltype(df[source])}}(undef, size(df,1))
  end

  if isnothing(groupfields) || length(groupfields)==0 #no grouping needed
    F = ecdf(collect(skipmissing(df[source])))

    if isfinite(F(0.0)) #we can shortcut this loop if there are no valid ecdf values
      df[target] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(df[source])
    end
  else #group the data using the provided fields
    for subdf::SubDataFrame ∈ groupby(df, :fisyr)
      F = ecdf(collect(skipmissing(subdf[source])))

      if isfinite(F(0.0)) #we can shortcut this loop if there are no valid ecdf values
        subdf[target] .= (f::MFloat64->ismissing(f) ? missing : F(f)).(subdf[source])
      end
  end
end

  return nothing
end

#creates a series of percentile intervals given the cutoffs
function createcutoffintervals(thresholds::Vector{Float64})
  intervals::Vector{Tuple{Float64,Float64}} = Vector{Tuple{Float64,Float64}}()

  sizehint!(intervals, length(thresholds))

  #make sure we are dealing with quantiles
  @assert minimum(thresholds) > 0.0
  @assert maximum(thresholds) == 1.0

  for i ∈ 1:length(thresholds)
    if i == 1
      push!(intervals, (0.0,thresholds[1]))
    else
      push!(intervals, (thresholds[i-1],thresholds[i]))
    end
  end

  return intervals
end


#holds all of the info to make a cross-tab table
struct CrossTabTable
  crosstabcolumns::Vector{Vector{Union{MFloat64,Int}}}

  rowshortnames::Vector{Symbol}
  #rowstatfields::Vector{Symbol}
  rowindex::Dict

  colshortnames::Vector{Symbol}
  #colstatfields::Vector{Symbol}
  colindex::Dict

  rowfield::Symbol
  colfield::Symbol
  tabfield::Symbol
  suffix::String
end

#a helper function to aggregate quantile columns for crosstab
function crosstabaggregate(target::AbstractVector; aggfunction::Function=median)::MFloat64
  local vcomplete::T where T <: Union{typeof(target), Vector{MFloat64}}

  if Missing <: eltype(target)
    vcomplete =  target[(!ismissing).(target)] #don't trust how missings are handled
  else
    vcomplete=target
  end

  if length(vcomplete) > 0
    return aggfunction(vcomplete)
  else
    return missing
  end
end

#makes fourblocker histograms of returns
function tabbedhistograms(df::AbstractDataFrame;
    outputpath::String = GRAPH_PATH,
    rowfield::Symbol=:qlreturn3yr,
    rowthresholds::Vector{Float64} = CROSSTAB_CUTOFFS,
    colfield = :qpprogramexpense,
    colthresholds::Vector{Float64} = CROSSTAB_CUTOFFS,
    tabfield::Symbol = :qpcontributions,
    suffix::String = "",
    shortnameindex::Dict = SHORT_NAME_INDEX,
    graphsuffix::String = GRAPH_SUFFIX)

  local plotnames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  rowcutoffintervals::Vector{Tuple{Float64,Float64}} = createcutoffintervals(rowthresholds)
  colcutoffintervals::Vector{Tuple{Float64,Float64}} = createcutoffintervals(colthresholds)

  graphindex::Dict = Dict{Tuple{Int,Int}, Int}()

  #iterate over the thresholds
  for r ∈ 1:length(rowcutoffintervals)
    for c ∈ 1:length(colcutoffintervals)
      dfsub::SubDataFrame =view(df, ((frow::MFloat64,fcol::MFloat64)->
        (!ismissing(frow)) && (!ismissing(fcol)) &&
        (colcutoffintervals[c][1] < fcol) && (fcol < colcutoffintervals[c][2]) &&
        (rowcutoffintervals[r][1] < frow) && (frow < rowcutoffintervals[r][2])
        ).(df[rowfield], df[colfield]), :)

      #index the location of the graph in the vector
      graphindex[(r,c)] = length(plots) + 1
      push!(plotnames,
        "$(shortnameindex[rowfield]): $(rowcutoffintervals[r][1])-$(rowcutoffintervals[r][2])
        $(shortnameindex[colfield]): $(colcutoffintervals[c][1])-$(colcutoffintervals[c][2])")

      #write the grpah
      push!(plots, plot(dfsub,
        x=tabfield,
        Guide.title(plotnames[end]),
        Guide.xlabel("contribution percentile"), Guide.ylabel("density"),
        Geom.line, Stat.histogram(bincount=50, density=true),
        Coord.cartesian(xmin=0.0, xmax=1.0)))

    end
  end

  finalplot::PlotContainer = title(gridstack(
    [plots[graphindex[(1,1)]] plots[graphindex[(1,2)]]
    plots[graphindex[(2,1)]] plots[graphindex[(2,2)]]]),
    "Normalized Contribution Percentile Histograms")

  outname::String = ("tabhist-$(shortnameindex[rowfield])X$(
    shortnameindex[colfield])BY$(shortnameindex[tabfield])$(suffix)_$(graphsuffix)")

  draw(PDF("$outputpath\\$(outname).pdf", 9inch, 7inch), finalplot)

  return nothing
end


#makes a histogram of betas
function plotβtab!(df::DataFrame;
  tabfield::Symbol = :b1lRet3yrCon,#:b1lRet3yrCon,
  suffix::String = "",
  outputpath::String = GRAPH_PATH,
  thresholds::Vector{Float64} = CROSSTAB_CUTOFFS,
  groupfield = :pprogramexpense_avg,
  shortnameindex::Dict = SHORT_NAME_INDEX,
  makequantiles::Bool = true,
  graphsuffix::String = GRAPH_SUFFIX
  )

  local plotnames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  cutoffintervals::Vector{Tuple{Float64,Float64}} = createcutoffintervals(thresholds)

  #make quintiles and classifers
  df = deepcopy(df[completecases(df[[groupfield, tabfield]]),:])
  makequantiles && createquantiles!(df, groupfield, groupfields = nothing)
  makequantiles && createquantiles!(df, tabfield, groupfields = nothing)
  tabfield = Symbol(:q, tabfield)
  groupfield = Symbol(:q, groupfield)
  df[:classification] = Vector{MSymbol}(undef, size(df,1))
  intervalnames = [:lowexp, :highexp]

  #given how much of a hack the below is...
  @assert length(cutoffintervals) == 2
  df[:classification] .= (f::Float64 ->
    f < thresholds[1] ? intervalnames[1] : intervalnames[2]).(df[groupfield])
  #=for c ∈ 1:length(cutoffintervals)
    dfsub::SubDataFrame =view(df, (f::Float64->
      (cutoffintervals[c][1] < f) && (f < cutoffintervals[c][2])
      ).(df[groupfield]), :)

    #add the classifier
    dfsub[:classification] .= intervalnames[c]

    #index the location of the graph in the vector
    #=push!(plotnames,
      "$(shortnameindex[groupfield]): $(cutoffintervals[c][1])-$(cutoffintervals[c][2])")

    #write the grpah
    push!(plots, plot(dfsub,
      x=tabfield,
      Guide.title(plotnames[end]),
      Guide.xlabel("contribution beta percentile"), Guide.ylabel("density"),
      Geom.line, Stat.histogram(bincount=25, density=false),
      Coord.cartesian(xmin=0.0, xmax=1.0)))=#

  end=#

  #=finalplot::PlotContainer = title(vstack(
    [plots[1],plots[2]]), "Normalized Contribution Percentile Histograms")=#


  finalplot = plot(df,
    x=tabfield, color=:classification,
    Guide.title("Histogram of beta percentile by governance category"),
    Guide.xlabel("contribution beta percentile"),
    Guide.ylabel("count"),
    Geom.line, Stat.histogram(bincount=25, density=true),
    Coord.cartesian(xmin=0.0, xmax=1.0))

  outname::String = (
    "tab-$(shortnameindex[groupfield])BY$(shortnameindex[tabfield])$(suffix)_$(graphsuffix)")

  draw(PDF("$outputpath\\$(outname).pdf", 9inch, 7inch), finalplot)

  return nothing
end


#creates the cross-tabulation table
function CrossTabTable(df::DataFrame;rowfield::Symbol=:qlreturn3yr,
  rowthresholds::Vector{Float64} = CROSSTAB_CUTOFFS,
  colfield = :qpprogramexpense,
  colthresholds::Vector{Float64} = CROSSTAB_CUTOFFS,
  yearrange::UnitRange = CROSSTAB_YEAR_RANGE,
  groupfield = :fisyr,
  tabfield::Symbol = :qpcontributions,
  aggfunction::Function = crosstabaggregate,
  suffix::String = "",
  shortnameindex::Dict = SHORT_NAME_INDEX)

  #build up indices and naming conventions
  local rowthresholdfields::Vector{Symbol} = (f::Float64 ->
    Symbol(rowfield, "p",Int(round(f*100)))).(rowthresholds)
  #local rowstatfields::Vector{Symbol} = [rowthresholdfields; :total; :N]
  local rowshortnames = (f::Float64 ->
    Symbol(shortnameindex[rowfield], "p",Int(round(f*100)))).(rowthresholds)
  rowshortnames = [rowshortnames; :total; :N]
  local nrows = length(rowshortnames)
  local rowindex::Dict = Dict(rowshortnames[r]=>r for r ∈ 1:length(rowshortnames))

  local colthresholdfields::Vector{Symbol} = (f::Float64 ->
    Symbol(colfield, "p",Int(round(f*100)))).(colthresholds)
  #local colstatfields::Vector{Symbol} = [colthresholdfields; :total; :N]
  local colshortnames = (f::Float64 ->
    Symbol(shortnameindex[colfield], "p",Int(round(f*100)))).(colthresholds)
  colshortnames = [colshortnames; :total; :N]
  local ncols::Int = length(colshortnames)
  local colindex::Dict = Dict(colshortnames[r]=>r for r ∈ 1:length(colshortnames))

  #this is what will actually hold the results
  crosstabcolumns::Vector{Vector{Union{MFloat64,Int}}} =
      (i->Vector{Union{MFloat64,Int}}(undef, nrows)).(1:ncols)

  #create a dataframe to work with
  dfcross::DataFrame = df[(y::MInt -> y ∈ yearrange).(
    df[:fisyr]), [rowfield, colfield, tabfield, groupfield]]

  #draw the graphs
  tabbedhistograms(dfcross,
    rowfield=rowfield,
    rowthresholds=rowthresholds,
    colfield=colfield,
    colthresholds=colthresholds,
    tabfield=tabfield,
    suffix=suffix,
    shortnameindex=shortnameindex)

  dfcross = dfcross[completecases(dfcross),:]
  #CSV.write("$WORKING_PATH\\crossdftest.csv", dfcross)

  #holds the top and bottom cutoff for each interval
  rowcutoffintervals::Vector{Tuple{Float64,Float64}} = createcutoffintervals(rowthresholds)
  colcutoffintervals::Vector{Tuple{Float64,Float64}} = createcutoffintervals(colthresholds)


  #now create the main table (We will need to make the row totals later)
  #NOTE: assumes that the cross tab columns come first!
  for c ∈ 1:length(colcutoffintervals)
    dfsub::SubDataFrame =view(dfcross, (f::MFloat64-> (!ismissing(f)) &&
      (colcutoffintervals[c][1] < f) && (f < colcutoffintervals[c][2])).(dfcross[colfield]), :)
    crosstabcolumns[c][rowindex[:total]] = aggfunction(dfsub[tabfield])
    crosstabcolumns[c][rowindex[:N]] = size(dfsub,1)

    for r ∈ 1:length(rowcutoffintervals)
      dfssub::SubDataFrame =view(dfsub, (f::MFloat64-> (!ismissing(f)) &&
        (rowcutoffintervals[r][1] < f) && (f < rowcutoffintervals[r][2])).(dfsub[rowfield]), :)
      crosstabcolumns[c][r] = aggfunction(dfssub[tabfield])
    end
  end

  #now make the rowwise subtotals
  for  r ∈ 1:length(rowcutoffintervals)
    dfsub::SubDataFrame =view(dfcross, (f::MFloat64-> (!ismissing(f)) &&
      (rowcutoffintervals[r][1] < f) && (f < rowcutoffintervals[r][2])).(dfcross[rowfield]), :)
    crosstabcolumns[colindex[:total]][r] = aggfunction(dfsub[tabfield])
    crosstabcolumns[colindex[:N]][r] = size(dfsub,1)
  end

  println(crosstabcolumns)

  #now make the grand total
  crosstabcolumns[colindex[:total]][rowindex[:total]] = aggfunction(dfcross[tabfield])
  crosstabcolumns[colindex[:N]][rowindex[:N]] = size(dfcross,1)

  return CrossTabTable(
    crosstabcolumns,

    rowshortnames,
    rowindex,

    colshortnames,
    colindex,

    rowfield,
    colfield,
    tabfield,
    suffix)
end


function writecrosstabtable(tab::CrossTabTable,
  outputpath::String = TABLE_PATH,
  decimals::Int = DECIMALS_RETURN,
  shortnameindex::Dict = SHORT_NAME_INDEX,
  graphsuffix::String = GRAPH_SUFFIX
  )

  local ncols::Int = length(tab.crosstabcolumns) #number of columns excluding the year
  local nrows::Int = length(tab.crosstabcolumns[1])

  #collect the keys reverse the dicitonaries for notational convenience
  local colindexkeys::Vector{Symbol} = collect(keys(tab.colindex))
  local rowindexkeys::Vector{Symbol} = collect(keys(tab.rowindex))

  #local revcolindex::Dict = Dict(c=>tab.colindex[colindexkeys[c]]  for c ∈ 1:ncols)
  #local revrowindex::Dict = Dict(r=>tab.rowindex[rowindexkeys[r]]  for r ∈ 1:nrows)

  #make the column names and sizes
  local colNames::Vector{Vector{String}} = Vector{Vector{String}}()
  push!(colNames, (s->string(s)).(tab.colshortnames))

  local widthColNames::Vector{Vector{Int}} = Vector{Vector{Int}}()
  push!(widthColNames, ones(Int, ncols))

  #make the row names
  local descRowNames::Vector{String} = (s->string(s)).(tab.rowshortnames)

  #need to print the descriptive rows
  descContent::Vector{Vector{String}} = (
    (i::Int)->Vector{String}(undef, ncols)).(1:nrows)

  #NOTE: tab has the data in COLUMNS while descContet absorbs data in ROWS (flipped indices)
  for r::Int ∈ 1:nrows
    for c::Int ∈ 1:ncols

      #Ns are ints
      if tab.rowshortnames[r] == :N || tab.colshortnames[c] == :N
        #lookup the collumn number from the column index for the pth column symbol
        descContent[r][c] = num2Str(tab.crosstabcolumns[c][r], Ints=true)
      else
        descContent[r][c] = num2Str(tab.crosstabcolumns[c][r], decimals)
      end
    end
  end

  outtabletitle::String = ("crosstab-$(shortnameindex[tab.rowfield])X$(
    shortnameindex[tab.colfield])BY$(shortnameindex[tab.tabfield])$(tab.suffix)")
  outtable::String = texTable(titleCaption=outtabletitle,
      caption="""See Tex File""", #caption
      colNames= colNames, #colNames
      descRowNames=descRowNames, #descRowNames
      descContent = descContent, #descContent
      nakedTable=true,
      widthColNames=widthColNames
    )

  writeNakedTable(outtable, path=outputpath, outName="$(outtabletitle)_$(graphsuffix).tex")

  return nothing
end
