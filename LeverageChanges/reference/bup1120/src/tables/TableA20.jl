


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

  #all the regression sets are the same, we are just sub-sampling the data
  specs::FMSpecs = FMSpecs(length(colkeys))
  (i->push!(specs,
    yspec=Fret,
    xspec=:yield,
    xnames=[:intercept, :yield],
    clusterspec=[:permno, :day])).(1:length(colkeys))

  #println(describe(hs))

  specs.clusterspecs[3] = :permno
  #set the filters
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[1]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[2]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[3]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[4]]], [Fret; :yield; :day; :permno; :eid]))
  push!(specdfs, view(hs, hs[!, colfilter[colkeys[5]]], [Fret; :yield; :day; :permno; :eid]))

  #generate the regression content (will display in Finometrics template)
  contentrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys[1:2])
  computeFMLMresults!(specdfs, specs, containsmissings=false)

  Ndescrows::Int = 2
  Ncols = length(colkeys)
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Ndescrows)
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys[3:4])

  for c ∈ 1:Ncols
    desccontent[1][c] = "$(num2str(R²(specs.results[c]), decimals))"
    desccontent[2][c] = "$(length(unique(specdfs[c].eid)))"
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

  hs=dist[dist.oneex .* dist.match,:]


#distcds::Vector{Vector{Int}} = (cd::Int->digits(cd)[4:-1:1]).(cumex.distcd) #(Reverse the order)
  hs = deleterows!(hs, (!).(completecases(view(hs, :, [:yield, Fret]))))

  #identify the announcement day
  N::Int = size(hs,1)
  hs.isann = Vector{Bool}(undef, N)
  hs.isinterim = Vector{Bool}(undef, N)
  hs.isex = Vector{Bool}(undef, N)
  hs.isannex = Vector{Bool}(undef, N)
  hs.isex40 = Vector{Bool}(undef, N)

  hs.tokeep = trues(N)
  @assert issorted(hs, [:eid, :day])

  hsss::GroupedDataFrame = groupby(hs, :eid)
  @mpar for i ∈ 1:length(hsss)
    shs::SubDataFrame = hsss[i]
    dvec::Vector{Int} = digits(shs[1,:distcd])[4:-1:1]
    idx::Int = searchsortedfirst(shs.day, -1)
    preexprice::Float64 = shs.day[idx] == -1 ? shs.price[idx] : -1.0

    #get the distribution codes and filter using the HS criteria
    if (dvec[1] ≠ 1) || (dvec[2] ≠ 2) || (dvec[3] > 5) || preexprice < 5.0
      shs.tokeep .= false
      continue
    end

    shs.isann .= shs.date .== shs.dclrdate
    shs.isex .= shs.day .== 0

    for r::DataFrameRow ∈ eachrow(shs)
      r.isinterim = (r.date > r.dclrdate) && (r.day < 0)
      r.isannex = r.isann || r.isinterim || r.isex
      r.isex40 = (r.day > 0) && (r.day ≤ 40)
    end
  end

  hs = hs[hs.tokeep,:]

  #println(size(hs))
  select!(hs, [:date, :isann, :isinterim, :isex, :isannex, :isex40,
    :permno, :day, Fret, :yield, :eid])

  #WARNING delete later
  #filter!(r::DataFrameRow->year(r.date) ≤ 2011, hs)
  #println(describe(hs))

  return hs
end
