

plusnumber(t::Int) = t>0 ? "+$t" : "$t"

#treat thsi table as a one-off since it has no common elements w/ other tables
function maketable7(bb::DataFrame)
  T7_PANEL_A_COLKEYS = [:StdDev, :Mean, :Range, :Return, :AbsReturn]
  T7_PANEL_A_NAME::String = "t7panela-$(REPLICATION_TYPE[])"
  T7_PANEL_A_KEYROW_RANGES::Dict = Dict(-3=>-6:-3, 0=>-2:2, 3=>3:6)
  local days::Int = SEO_DAYS

  table7panela(bb,
      days=days,
      tablename= T7_PANEL_A_NAME,
      colkeys = T7_PANEL_A_COLKEYS,
      keyrowranges = T7_PANEL_A_KEYROW_RANGES,
      )
end

function table7panela(bb::DataFrame;
    days = SEO_DAYS,
    tablename::String = error("tablename is required"),
    tablepath::String = TABLE_PATH,
    decimals::Int = DEFAULT_DECIMALS,
    colkeys::Vector{Symbol} = error("colkeys is required"),
    colnames::Vector{Vector{String}} = [["Daily Cross-", "", ""],
      ["Sectional Returns", "", "Average"],
      (string).(colkeys)],
    widthcolnames::Vector{Vector{Int}} = [[2,1,2], [2,1,2], ones(Int, length(colnames[2]))],
    descrownames::Vector{String} = (plusnumber).(-days:days),
    keyrowranges::Dict = error("keyrow ranges is required"))

  local parallel::Bool = PARALLEL[] && (!DISABLE_FM_PARALLEL)


  #################################
  #basic lookups
  i2t::Vector{Int} = collect(-days:days)
  dayfields::Vector{Symbol} = getdayfields(days)
  Fₜ::Dict{Int, Symbol} = Dict(t=>dayfields[i] for (i,t) ∈ enumerate(i2t))

  #preallocation and housekeeping

  Ncols::Int = sum(widthcolnames[1])
  Nrows::Int = 2*days+1
  desccontent::Vector{Vector{String}} = (i->Vector{String}(undef, Ncols)).(1:Nrows)

  #create an index to make this part of the code more readable
  colidx::Dict = Dict(colkey=>i for (i, colkey) ∈ enumerate(colkeys))

  #now populate the the table
  for (i,rowkey) ∈ enumerate(descrownames)
    f::Symbol = dayfields[i]
    t::Int = i2t[i] #get the period offset

    #compute the following stats for each row
    daystd::Float64 = std(skipmissing(bb[!, f]))*100
    daymean::Float64 = mean(skipmissing(bb[!, f]))*100


    desccontent[i][colidx[:StdDev]] = num2str(daystd, decimals)
    desccontent[i][colidx[:Mean]] = num2str(daymean, decimals)

    if haskey(keyrowranges, t)

      tₗ::Int = minimum(keyrowranges[t])
      tₕ::Int = maximum(keyrowranges[t])
      rngfields::Vector{Symbol} = (i->Fₜ[i]).(keyrowranges[t])
      desccontent[i][colidx[:Range]] = "\\text{$(plusnumber(tₗ)) to $(plusnumber(tₕ))}"

      println(rngfields)
      rngmat::Matrix{Float64} = Matrix(bb[:, rngfields])
      rets::Vector{Float64} = vec(prod(rngmat .+ 1., dims=2)) .- 1
      retmean::Float64 = mean(rets)*100
      retabsmean::Float64 = mean((abs).(rets))*100

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
    widthcolnames=widthcolnames)

  writetextable(tt, path=tablepath, outname="$tablename.tex")
end
