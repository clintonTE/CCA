
#treat thsi table as a one-off since it has no common elements w/ other tables
function maketable8(dist::DataFrame;
    panela::Bool = true,
    panelb::Bool = true,
    panelc::Bool = true)


  ################common data
  #println("Size of file: $(length(unique(dist.eid)))")
  #println(describe(dist))
  #error("stop")

  ###################PANEL A############

  T8_PANEL_A_NAME::String = "t8panela-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  T8_PANEL_A_ROWKEYS::Vector{Symbol} = [
    :firstlast,
    :alldist,
    :ordinary,
    :match,
    :oneex,
    :dclrexok]

  T8_PANEL_A_ROWLABEL::Dict = Dict(
    :firstlast => "First and Last Dividends",
    :alldist => "All Distributions on CRSP",
    :ordinary => "Ordinary `CRSP 1232' Dividends",
    :match => "Matches to CRSP",
    :oneex => "Only 1 ex-date in window",
    :dclrexok => "\$\\ge$DCLR_EX_DAYS\$ Between Declaration and Ex-Date")

  panela && table8panela(dist,
      tablename= T8_PANEL_A_NAME,
      rowkeys=T8_PANEL_A_ROWKEYS,
      rowlabel=T8_PANEL_A_ROWLABEL,
      dclrexdays = DCLR_EX_DAYS,
      colnames = [[""]],
      widthcolnames = [[1]])

  ###################PANEL B############


  sdist::SubDataFrame = view(dist, dist.oneex .* dist.match .* dist.ordinary,:)

  T8_PANEL_B_ROWKEYS::Vector{Symbol} = [
    :dclrex,
    :cumulative,
    :Mcumulative,
    :after,
    :dclrex2,
    :cumulative2,
    :Mcumulative2,]

  T8_PANEL_B_ROWLABEL::Dict = Dict(
    :dclrex=>"Dclrtn to Ex",
    :cumulative=>"Cumulative",
    :Mcumulative=>"\$-\$Cumulative",
    :after=>"\\midrule After Selection",
    :dclrex2=>"\\midrule Dclrtn to Ex",
    :cumulative2=>"Cumulative",
    :Mcumulative2=>"\$-\$Cumulative",)

  T8_PANEL_B_COLKEYS::Vector{Int} = collect(0:9)
  T8_PANEL_B_COLNAMES::Vector{Vector{String}} = [["Days"; (string).(T8_PANEL_B_COLKEYS)]]

  T8_PANEL_B_NAME::String = "t8panelb-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  panelb && table8panelb(sdist,
    rowkeys=T8_PANEL_B_ROWKEYS,
    rowlabel=T8_PANEL_B_ROWLABEL,
    tablename=T8_PANEL_B_NAME,
    colkeys = T8_PANEL_B_COLKEYS,
    colnames = T8_PANEL_B_COLNAMES,
    dclrexdays = DCLR_EX_DAYS
  )

  ###################PANEL C############
  sdist = view(dist, dist.primary,:)
  T8_PANEL_C_ROWKEYS::Vector{Symbol} = [
    :return,
    :returnabs,
    :yield]

  T8_PANEL_C_ROWLABEL::Dict = Dict(
    :return=>"Net-of-Market Rates of\\\\Return, \$-12\$ to \$+12\$ days",
    :returnabs=>"\\dots in abs values",
    :yield=>"Div Yield \$\\delta, |\\ge 9\$ days")

  T8_PANEL_C_COLKEYS::Vector{Symbol} = [
    :min,
    :median,
    :max,
    :pos,
    :mean,
    :sd,
    :t,
    :N,
    :Nwinsorized]

  T8_PANEL_C_COLLABEL::Dict = Dict(
    :min=>"Min",
    :median=>"Median",
    :max=>"Max",
    :pos=>"\\%Pos",
    :mean=>"Mean",
    :sd=>"SD",
    :t=>"Trd-days",
    :N=>"\\#",
    :Nwinsorized=>"\\# Winsrzd")


  T8_PANEL_C_NAME::String = "t8panelc-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"

  panelc && table8panelc(sdist,
    rowkeys=T8_PANEL_C_ROWKEYS,
    rowlabel=T8_PANEL_C_ROWLABEL,
    tablename=T8_PANEL_C_NAME,
    colkeys = T8_PANEL_C_COLKEYS,
    collabel = T8_PANEL_C_COLLABEL,
    dclrexdays = DCLR_EX_DAYS
  )
end

function table8panelc(dist::AbstractDataFrame;
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Symbol} = error("colkeys is required"),
  collabel::Dict = error("collabel is required"),
  tablename::String = error("table name is required"),
  dclrexdays::Int = error("dclrexdays is required"),
  tablepath::String = TABLE_PATH,
  Fret::Symbol = DIST_RET,
  decimals::Int = 3,)

  local distagg::DataFrame

  colnames::Vector{Vector{String}} = [[(k::Symbol->string(collabel[k])).(colkeys)...;]]

  Nrows::Int = length(rowkeys)
  Ncols::Int = length(colnames[1])

  #pre-allocate space and index
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)

  #index the rows and cols for readability
  rowidx::Dict = Dict(s=>desccontent[i] for (i,s::Symbol) ∈ enumerate(rowkeys))
  colidx::Dict = Dict(s=>i for (i,s::Symbol) ∈ enumerate(colkeys))

  #does the formatting as required
  function plusnum(num::Real, d::Int=decimals)
    local s::String
    (abs(num) < 10^-10)  && (num = 0)
    if num > 10^-10
      s = "+$(num2str(num, d, Ints=true))"
    else
      s = "$(num2str(num, d, Ints=true))"
    end
    println("$s: $num")
    return s
  end

  #handles the common elements across rows
  function writepartialrow(r::Vector{String}, rets::Vector{Float64}, scaling::Int)
    r[colidx[:min]] = plusnum(minimum(rets) * scaling)
    r[colidx[:median]] = plusnum(median(rets) * scaling)
    r[colidx[:max]] = plusnum(maximum(rets) * scaling)
    r[colidx[:mean]] = num2str(mean(rets) * scaling, decimals)
    r[colidx[:sd]] = num2str(std(rets) * scaling, decimals)
  end

  rets::Vector{Float64} = collect(skipmissing(dist[!, Fret]))
  absrets::Vector{Float64} = (abs).(rets)

  #grab the dividend yield for each event
  sdist = view(dist, dist.dclrexok, :)
  distagg = aggregate(groupby(view(sdist, sdist.day.==0,
    [:eid, :yield, :winsorizedyield, :dclrexdays, :dclrextradingdays]), :eid), x->x[1])

  yields::Vector{Float64} = collect(skipmissing(distagg.yield_function))

  writepartialrow(rowidx[:return], rets, 100)
  writepartialrow(rowidx[:returnabs], absrets, 100)
  writepartialrow(rowidx[:yield], yields, 100)

  #now that we have done the common elements, compute the rest
  #start with return
  rowidx[:return][colidx[:pos]] = "$(plusnum(sum(rets .> 0)/length(rets) * 100,1))\\%"

  meandays::Float64 = mean(skipmissing(distagg.dclrextradingdays_function))
  rowidx[:return][colidx[:t]] = "$(plusnum(meandays, 1))"
  rowidx[:return][colidx[:N]] = "$(length(rets))"

  Nwinsorized::Int = sum(skipmissing(dist.winsorized))
  rowidx[:return][colidx[:Nwinsorized]] = "$Nwinsorized"

  #now compute additional rows for the divdend yield
  rowidx[:yield][colidx[:N]] = "$(length(yields))"

  Nwinsorizedyield::Int = sum(skipmissing(distagg.winsorizedyield_function))
  rowidx[:yield][colidx[:Nwinsorized]] = "$Nwinsorizedyield"

  #now compute the "stub" blanks
  rowidx[:returnabs][colidx[:pos]] = ""
  rowidx[:returnabs][colidx[:t]] = ""
  rowidx[:returnabs][colidx[:N]] = ""
  rowidx[:returnabs][colidx[:Nwinsorized]] = ""
  rowidx[:yield][colidx[:pos]] = ""
  rowidx[:yield][colidx[:t]] = ""

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent)

  writetextable(tt, path=tablepath, outname="$tablename.tex")

end


function table8panela(dist::DataFrame;
    tablename::String = error("tablename is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS,
    colnames::Vector{Vector{String}} = error("colnames is required"),
    rowkeys::Vector{Symbol} = error("rowkeys is required"),
    rowlabel::Dict = error("rowlabel is required"),
    dclrexdays::Int = error("dclrexdays is required"),
    widthcolnames::Vector{Vector{Int}} =  error("widthcolnames is required"))

  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)
  local sdist::SubDataFrame
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)

  #################################

  #preallocation and housekeeping
  Ncols::Int = sum(widthcolnames[1])
  Nrows::Int = length(descrownames)
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)
  alignmentstring::String = "lc"


  #make the content
  minyr, minmo = yearmonth(minimum(dist.exdate))
  maxyr, maxmo = yearmonth(maximum(dist.exdate))

  #index the rows for readability
  rowindex::Dict = Dict(s=>desccontent[i] for (i,s::Symbol) ∈ enumerate(rowkeys))


  rowindex[:firstlast][1] = "\\text{$minyr/$minmo to $maxyr/$maxmo}"

  rowindex[:alldist][1] = "\\text{N=}$(length(unique(dist.eid)))"

  sdist = view(dist, dist.ordinary, :)
  rowindex[:ordinary][1] = "\\text{N=}$(length(unique(sdist.eid)))"

  sdist = view(sdist, sdist.match, :)
  rowindex[:match][1] = "\\text{N=}$(length(unique(sdist.eid)))"

  sdist = view(sdist, sdist.oneex, :)
  rowindex[:oneex][1] = "\\text{N=}$(length(unique(sdist.eid)))"

  sdist = view(sdist, sdist.dclrexok, :)
  rowindex[:dclrexok][1] = "\\text{N=}$(length(unique(sdist.eid)))"


  #println(desccontent)


  #display(desccontent)
  #construct the tex table
  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    alignmentstring=alignmentstring,
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end

function table8panelb(dist::AbstractDataFrame;
  rowkeys::Vector{Symbol} = error("rowkeys is required"),
  rowlabel::Dict = error("rowlabel is required"),
  colkeys::Vector{Int} = error("colkeys is required"),
  colnames::Vector{Vector{String}} = error("collabel is required"),
  tablename::String = error("table name is required"),
  dclrexdays::Int = error("dclrexdays is required"),
  tablepath::String = TABLE_PATH,
  decimals::Int = DEFAULT_DECIMALS,)


  Nrows::Int = length(rowkeys)
  Ncols::Int = length(colnames[1]) - 1

  #pre-allocate space and index
  descrownames::Vector{String} = (s::Symbol->rowlabel[s]).(rowkeys)
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)
  desccondentidx::Dict = Dict(k=>desccontent[i] for (i::Int,k::Symbol) ∈ enumerate(rowkeys))

  #compute the number of events for each length of time between dclr and ex
  alldclrexdays::Vector{Int} = collect(0:(maximum(colkeys)+dclrexdays -1))
  dclrexdensity::Dict = Dict(Δ=>
    length(unique(view(dist, dist.dclrexdays .== Δ, :eid))) for Δ::Int ∈ alldclrexdays)

  #some important totals
  totevent::Int = length(unique(dist.eid))
  totselected::Int = length(unique(view(dist, dist.dclrexdays .≥ dclrexdays, :eid)))
  #construct the table
  for (i::Int,Δ::Int) ∈ enumerate(colkeys)

    dclrex::Int = dclrexdensity[Δ]
    desccondentidx[:dclrex][i] = num2str(dclrex, 0)

    #compute the cumulative sum
    cumulative::Int = sum((d->dclrexdensity[d]).(alldclrexdays[alldclrexdays .≤ Δ]))
    desccondentidx[:cumulative][i] = num2str(cumulative,0)

    #back out the negative cumulative sum
    Mcumulative::Int = totevent - cumulative
    desccondentidx[:Mcumulative][i] = num2str(Mcumulative,0)

    #second half table header
    Δ⁺::Int = Δ + 8
    desccondentidx[:after][i] = "$Δ⁺"

    if Δ⁺ == 8 #sepecial case to initialize the second part of the table
      desccondentidx[:dclrex2][i] = ""
      desccondentidx[:cumulative2][i] = ""
      desccondentidx[:Mcumulative2][i] = num2str(totselected,0)
    else
      dclrex2::Int = dclrexdensity[Δ⁺]
      desccondentidx[:dclrex2][i] = num2str(dclrex2, 0)

      #sums the density from the minimum selected period to Δ⁺
      cumulative2::Int = sum((d->
        dclrexdensity[d]).(alldclrexdays[(t->(t≥dclrexdays) && (t≤Δ⁺)).(alldclrexdays)]))

      desccondentidx[:cumulative2][i] = num2str(cumulative2,0)

      Mcumulative2::Int = totselected - cumulative2
      desccondentidx[:Mcumulative2][i] = num2str(Mcumulative2,0)
    end
  end

  #println(desccontent)

  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    rowlabelheader=true)

  writetextable(tt, path=tablepath, outname="$tablename.tex")

end
