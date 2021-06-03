


function maketablea2(dist::DataFrame;
    panela::Bool = true,
    panelb::Bool = true,
    )::Nothing

  #adjust the dataframe for the Harzmark Solomon stats
  hs::DataFrame = makehs(dist)

  #########################Panel a
  TA2_PANEL_A_ROWKEYS::Vector{Symbol} = [
    :ann,
    :interim,
    :ex,
    :annex,
    :ex40]

  TA2_PANEL_A_ROWLABEL::Dict = Dict(
    :ann=>"Announcement day",
    :interim=>"Interim",
    :ex=>"Ex-day",
    :annex=>"Ann-to-ex",
    :ex40=>"40 days after"
    )
  TA2_PANEL_A_COLNAMES::Vector{Vector{String}} = [["Return"]]

  #####################Panel b#######
  TA2_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :constant,
    :yield,
    :r2,
    :N,]

  TA2_PANEL_B_ROWLABEL::Dict = Dict(
    :constant=>"Constant",
    :yield=>"Div yield",
    :r2=>"\$R^2\$",
    :N=>"N",
    )

  TA2_PANEL_B_COLKEYS::Vector{Symbol} = [
    :ann,
    :interim,
    :ex,
    :annex,
    :ex40
  ]

  TA2_PANEL_B_COLNAMES::Vector{Vector{String}} = [[
    "Announcement day", "Interim", "Ex-day", "Ann-to-ex", "40 days after"]]

  TA2_PANEL_B_COLFILTER = Dict(
    :ann=>:isann,
    :interim=>:isinterim,
    :ex=>:isex,
    :annex=>:isannex,
    :ex40=>:isex40
  )

  TA2_PANEL_B_NAME::String = "ta2panelb-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  TA2_PANEL_B_Σ = clusteredΣ!

  panelb && tablea2panelb(hs,
    rowkeys=TA2_PANEL_B_ROWKEYS,
    rowlabel=TA2_PANEL_B_ROWLABEL,
    tablename=TA2_PANEL_B_NAME,
    colkeys = TA2_PANEL_B_COLKEYS,
    colfilter = TA2_PANEL_B_COLFILTER,
    colnames = TA2_PANEL_B_COLNAMES,
    Σfunc = TA2_PANEL_B_Σ
  )
end

function tablea2panelb(hs;
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Symbol} = error("colkeys is required"),
  colnames::Vector{Vector{String}} = error("colnames is required"),
  colfilter::Dict = error("colfilter is required"),
  tablename::String = error("table name is required"),
  Fret::Symbol = DIST_RET,
  Σfunc::Function = error("Σfunc is required"),
  tablepath::String = TABLE_PATH,
  decimals::Int = 3,)

  #the colmns vary by the views
  local specdfs::Vector{SubDataFrame} = Vector{SubDataFrame}()

  #only changing the dependent variable with the regressions
  specs::FMSpecs = FMSpecs(length(colkeys))
  xnames::Vector{Symbol} = [:intercept, :yield]
  clusterspec::Vector{Symbol} = [:permno, :exdate]
  xspec::Symbol=:yield
  push!(specs, yspec=:ann,  xspec=xspec,  xnames=xnames,  clusterspec=clusterspec)
  push!(specs, yspec=:interim,  xspec=xspec,  xnames=xnames,  clusterspec=clusterspec)
  push!(specs, yspec=:ex,  xspec=xspec,  xnames=xnames,  clusterspec=clusterspec)
  push!(specs, yspec=:annex,  xspec=xspec,  xnames=xnames,  clusterspec=clusterspec)
  push!(specs, yspec=:ex40,  xspec=xspec,  xnames=xnames,  clusterspec=clusterspec)

  #println(describe(hs))

  specs.clusterspecs[3] = :permno
  #set the filters
  #=push!(specdfs, view(hs, hs[!, colfilter[colkeys[1]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[2]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[3]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[4]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[5]]], [Fret; :yield; :day; :permno; :eid]))=#

  println(describe(hs))

  #generate the regression content (will display in Finometrics template)
  contentrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys[1:2])
  computeFMLMresults!(hs, specs#=, containsmissings=false=#)

  Ndescrows::Int = 2
  Ncols = length(colkeys)
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Ndescrows)
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys[3:4])

  for c ∈ 1:Ncols
    desccontent[1][c] = "$(num2str(R²(specs.results[c]), decimals))"
    desccontent[2][c] = "$(specs.results[c].N)"
  end

  #=for (i,lin) ∈ enumerate(specs.results)
    println("Iteration $i: ")
    display(Σfunc(lin))
  end=#

  tt::String = textable(specs.results, Σfunc, [:intercept, :yield],
    colnames=colnames,
    scaling=[1.,1.],
    descrownames=descrownames,
    desccontent=desccontent,
    decimaldigits=3)

  writetextable(tt, path=tablepath, outname="$tablename.tex")

  return nothing
end

#creates the neccessary fields to make the hs dataframe
function makehs(dist::DataFrame; Fret= DIST_RET)
  local hs::DataFrame

  sdist::SubDataFrame = view(dist, dist.oneex .* dist.match,:)

  #use this to aggregate over the days where appropriate
  @inline function aggfuncret(v::AbstractVector{<:MFloat64})
    (sum((ismissing).(v)) == length(v)) && return missing

    return prod(skipmissing(1.0 .+ v)) - 1.0
  end

  #form the root of the aggregated dataframe
  aggnames::Vector{Symbol} = [:permno, :eid, :yield, :exdate, :distcd, :dclrdate]
  sdisttoagg::SubDataFrame = view(sdist, (!ismissing).(sdist.yield),aggnames)
  hs = aggregate(groupby(sdisttoagg, :eid), v->v[1])

  #fix the names
  @inline nofunction(s::Symbol) = Symbol(replace(string(s), "_function"=>""))
  rename!(hs, (s::Symbol-> s=>nofunction(s)).(names(hs)))

  #these will be the outcome variables
  N::Int =size(hs,1)
  hs.ann = Vector{MFloat64}(undef, N)
  hs.interim = Vector{MFloat64}(undef, N)
  hs.ex = Vector{MFloat64}(undef, N)
  hs.annex = Vector{MFloat64}(undef, N)
  hs.ex40 = Vector{MFloat64}(undef, N)

  #prepare for the aggregation and create tags for deleting or keeping rows
  hs.tokeep = trues(N)
  hsindex::Dict = Dict(r.eid=>r for r::DataFrameRow ∈ eachrow(hs))

  hs.tokeep = trues(N)
  @assert issorted(sdist, [:eid, :day])

  ssdists::GroupedDataFrame = groupby(sdist, :eid)
  @mpar for i ∈ 1:length(ssdists)
    ssdist::SubDataFrame = ssdists[i]
    r::DataFrameRow = hsindex[ssdist.eid[1]]

    #find the index of the exdate if possible
    dvec::Vector{Int} = digits(ssdist.distcd[1])[4:-1:1]
    exidx::MInt = searchsortedfirst(ssdist.day, 0)
    local preexprice::MFloat64
    if ssdist.date[exidx] == r.exdate
      preexprice = ssdist.day[exidx-1] == -1 ? ssdist.price[exidx-1] : -1.0
    else
      exidx = missing
      preexprice = missing
    end

    dclridx::MInt = searchsortedfirst(ssdist.date, r.dclrdate)
    (ssdist.date[dclridx] ≠ r.dclrdate) && (dclridx = missing)

    #get the distribution codes and filter using the HS criteria
    #also some basic integrity checks
    if ((dvec[1] ≠ 1) ||
      (dvec[2] ≠ 2) ||
      (dvec[3] > 3) || #allow quarterly dividends only
      ismissing(exidx) ||
      ismissing(preexprice) ||
      (preexprice < 5.0) ||
      ismissing(dclridx))

      r.tokeep = false
      continue
    end

    r.ann = ssdist[dclridx, Fret]
    r.ex = ssdist[exidx, Fret]

    if exidx - dclridx > 1
      r.interim = aggfuncret(ssdist[(dclridx+1):(exidx-1), Fret])
    end

    r.annex = aggfuncret(ssdist[(dclridx):(exidx), Fret])
    r.ex40 = aggfuncret(ssdist[(d::Int->(d>0) && (d≤40)).(ssdist.day), Fret])

  end

  hs = hs[hs.tokeep,:]

  #WARNING delete later
  filter!(r::DataFrameRow->year(r.exdate) ≤ 2011, hs)
  #println(describe(hs))

  return hs
end
