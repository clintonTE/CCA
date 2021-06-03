
#MAIN ENTRY POINT
#final columns for the decompressed file

#TODO: Check the portfolio sorts and the factor returns
#(replicate on 2015/12/31)


#merges crsp and compustat
# note refresh merge will be automatically run if any of the other two are run
function formfactors(;
  refreshcomp::Bool = true,
  refreshcrsp::Bool = true,
  refreshmerge::Bool = true,
  refreshfactors::Bool = true,
  datapath::String=DATA_PATH,
  factname::String = FACT_NAME,
  injlsstream::Function = IN_JLS_STREAM,
  outjlsstream::Function = OUT_JLS_STREAM,
  validatemerged::Bool = true,
  parallel::Bool = false,
  outpath::String = OUT_PATH,
  writecsv::Bool = true)::Nothing

  local univ::DataFrame
  local fact::DataFrame

  #NOTE: Need to move this to the entry point if the entry point changes
  PARALLEL[] = parallel

  #refreshfactors must be true if any of the others are true
  refreshfactors = (refreshfactors || refreshmerge || refreshcrsp ||
    refreshcomp || validatemerged)

  #makes a serialization object so we don't have to keep parsing the CSV
  if refreshfactors || (!isfile("$datapath\\$factname.jls.lz4"))
    univ = makeuniv(refreshmerge=refreshmerge,  refreshcomp=refreshcomp,
      refreshcrsp=refreshcrsp, validatemerged=validatemerged) do univ::DataFrame
        formportfolios!(univ)
    end

    fact = formfactors(univ)

    #NOTE: may want to push the portfolio formation into univ once it is done
    #if so, the below will be just the factor portfolio
    outjlsstream("$datapath\\$factname.jls.lz4") do s
      serialize(s, fact)
    end

  else
    injlsstream("$datapath\\$factname.jls.lz4") do s
      fact = deserialize(s) #NOTE: DEBUG ONLY
    end
  end

  writecsv && (fact |> CSV.write("$outpath\\$factname.csv"))

  return nothing
end

#assigns the group for a portfolio
function assigngroups!(vals::AbstractVector{T},
  quantiles::AbstractVector{MFloat64},
  groups::AbstractVector{MInt},
  breakpoints::NTuple) where T

  #compute the quantiles
  qdist = StatsBase.ecdf((val::T->val).(vals))
  @inbounds @simd for i ∈ 1:length(vals)
    f::Float64 = vals[i]
    quantiles[i] = qdist(f)
  end

  #compute the groups
  groups .= (q::Float64-> findfirst(b::Float64->
      q≤b, breakpoints)).(quantiles)


  return nothing
end



function formportfolios!(univ::DataFrame;
    pcs::PortfolioConstructions = PortfolioConstructions(),
    weighton::Symbol = WEIGHT_ON, returnfield::Symbol=RETURN_FIELD)

  local sunivs::GroupedDataFrame

  N::Int = size(univ,1)
  #first pre-allocate the requisite fields
  for pc ∈ pcs
    univ[!, pc.Fquantile] = Vector{MFloat64}(undef, N)
    univ[!, pc.Fgroup] = Vector{MInt}(undef, N)
    univ[!, pc.Fweight] = Vector{MFloat64}(undef, N)
  end

  sunivs = groupby(univ, :date)
  for i ∈ 1:length(sunivs)
    suniv::SubDataFrame = sunivs[i]
    for pc ∈ pcs
      assigngroups!(suniv[!,pc.Fval],
        suniv[!, pc.Fquantile], suniv[!, pc.Fgroup], pc.breakpoints)
    end
  end

  println("size: ", N)
  #group field
  groupfields = (pc->pc.Fgroup).(pcs)
  univ.returndollars = Vector{MFloat64}(undef, N)
  for suniv ∈ groupby(univ, [:date; groupfields])
    for pc ∈ pcs
      suniv.returndollars .= suniv[!,weighton] .* suniv[!, returnfield]
    end
  end

  return univ
end

#this actually computes the factor returns
function computefactors(univ::DataFrame;
    pcs::PortfolioConstructions = PortfolioConstructions(),
    weighton::Symbol = WEIGHT_ON)

  fact::DataFrame = DataFrame(date=unique(univ.date))
  sort!(fact, :date)
  Nfact::Int = size(fact,1)
  fact.retid = collect(1:Nfact)

  #pre-allocate the return fields
  for pc ∈ pcs
     fact[!,pc.Fval] = zeros(Float64, Nfact)
  end

  #index the rows by date
  rfacts::Dict = Dict(r.date => r for r ∈ eachrow(fact))

  #now compute the returns
  groupfields::Vector{Symbol} = (pc->pc.Fgroup).(pcs)
  for suniv ∈ groupby(univ, :date)
    rfact = rfacts[suniv.date[1]]

    #split each date by triple sort
    for ssuniv ∈ groupby(suniv, groupfields)

      #this is the overall portfolio return
      portret::Float64 = sum(ssuniv.returndollars) / sum(ssuniv[!,weighton])
      #println("portret: $portret")

      for pc ∈ pcs
        #lookup the portfolio weight for the specific factor

        group::Int = ssuniv[1, pc.Fgroup]
        groupweight::Float64 = pc.groupweights[group]

        #println("groupweight: $groupweight")

        rfact[pc.Fval] += groupweight * portret
      end
    end
  end

  return fact
end
