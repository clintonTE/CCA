

plusnumber(t::Int) = t>0 ? "+$t" : "$t"

#treat thsi table as a one-off since it has no common elements w/ other tables
function maketable7(bb::DataFrame)
  T7_PANEL_A_COLKEYS = [:StdDev, :Mean, :Range, :Return, :AbsReturn]
  T7_PANEL_A_NAME::String = "t7panela-$(REPLICATION_TYPE[])-$(DIST_TYPE[])-$(OUT_SUFFIX[])"
  T7_PANEL_A_KEYROW_RANGES::Dict = Dict(-3=>-6:-3, 0=>-2:2, 3=>3:6)
  local days::Int = SEO_DAYS
  local Fret::Symbol = SEO_RET
  local minallowed::Int = SEO_MIN_ALLOWED

  requirefield!(bb, Fret, minallowed)
  deleterows!(bb, (!).(bb.oneex))

  table7panela(bb,
      days=days,
      tablename= T7_PANEL_A_NAME,
      colkeys = T7_PANEL_A_COLKEYS,
      keyrowranges = T7_PANEL_A_KEYROW_RANGES,
      Fret=Fret
      )
end

function requirefield!(df::DataFrame, Fret::Symbol, minallowed::Int; Fgroup=:permno)
  df.tokeep = trues(size(df,1))
  for sdf ∈ groupby(df, :permno)
    tokeep::Bool = sum((!ismissing).(sdf[!,Fret])) ≥ minallowed
    sdf.tokeep .= tokeep
  end

  deleterows!(df, (!).(df.tokeep))
  select!(df, Not(:tokeep))

  return nothing
end

function table7panela(bb::AbstractDataFrame;
    days = SEO_DAYS,
    tablename::String = error("tablename is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS,
    colkeys::Vector{Symbol} = error("colkeys is required"),
    colnames::Vector{Vector{String}} = [["", "Daily Cross-", "", ""],
      ["Event", "Sectional Returns", "", "Average"],
      ["Day"; (string).(colkeys);]],
    widthcolnames::Vector{Vector{Int}} = [[1,2,1,2], [1,2,1,2], ones(Int, length(colnames[3]))],
    descrownames::Vector{String} = (plusnumber).(-days:days),
    keyrowranges::Dict = error("keyrow ranges is required"),
    Fret::Symbol = error("Fret is required"),
    )

  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)
  #bb = bb[completecases[]]


  #################################
  #basic lookups
  i2day::Vector{Int} = collect(-days:days)
  sbbs::GroupedDataFrame = groupby(bb, :day)
  bbidx::Dict = Dict((i::Int->sbbs[i][1,:day]=>sbbs[i]).(1:length(sbbs)))

  #preallocation and housekeeping

  Ncols::Int = sum(widthcolnames[1])-1
  Nrows::Int = 2*days+1
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)

  #create an index to make this part of the code more readable
  colidx::Dict = Dict(colkey=>i for (i, colkey) ∈ enumerate(colkeys))

  #now populate the the table
  for (i,rowkey) ∈ enumerate(descrownames)
    day::Int = i2day[i] #get the period offset
    sbb::SubDataFrame = bbidx[day]

    #compute the following stats for each row
    daystd::Float64 = std(skipmissing(sbb[!, Fret]))*100
    daymean::Float64 = mean(skipmissing(sbb[!, Fret]))*100


    desccontent[i][colidx[:StdDev]] = num2str(daystd, decimals)
    desccontent[i][colidx[:Mean]] = num2str(daymean, decimals)

    if haskey(keyrowranges, day)

      dayₗ::Int = minimum(keyrowranges[day])
      dayₕ::Int = maximum(keyrowranges[day])
      desccontent[i][colidx[:Range]] = "\\text{$(plusnumber(dayₗ)) to $(plusnumber(dayₕ))}"

      ssbb = view(bb, (r::DataFrameRow->r.day ∈ keyrowranges[day]).(eachrow(bb)), [:permno, Fret])
      #=aggbb::DataFrame = aggregate(groupby(ssbb, :permno),
        v::AbstractVector{MFloat64}->prod(1 .+ skipmissing(v)))
      rets::Vector{Float64} = aggbb[!,Symbol(Fret, :_function)] .- 1

      retmean::Float64 = mean(skipmissing(rets))*100
      retabsmean::Float64 = mean(skipmissing((abs).(rets)))*100=#

      retmean::Float64 = mean(skipmissing(ssbb[!,Fret]))*100
      retabsmean::Float64 = mean(skipmissing((abs).(ssbb[!,Fret])))*100

      desccontent[i][colidx[:Return]] = num2str(retmean, decimals)
      desccontent[i][colidx[:AbsReturn]] = num2str(retabsmean, decimals)
    else
      desccontent[i][colidx[:Range]] = ""
      desccontent[i][colidx[:Return]] = ""
      desccontent[i][colidx[:AbsReturn]] = ""
    end
  end

  alignmentstring::String = "rrr c rr"

  #println(desccontent)


  #display(desccontent)
  #construct the tex table
  tt::String = textable(colnames=colnames,
    descrownames=descrownames,
    desccontent=desccontent,
    alignmentstring=alignmentstring,
    widthcolnames=widthcolnames,
    rowlabelheader=true)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end
