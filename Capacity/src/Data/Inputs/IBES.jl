
#reads in the ibis data, th
function prepibes(;
  ibesname::String = PARAM[:ibesname],
  ibespath::String = PARAM[:ibespath],
  cilname::String = PARAM[:cilname],
  cilpath::String = PARAM[:cilpath],
  incsvstream::Function = IN_CSV_STREAM,
  csvextension::String = CSV_EXTENSION,
  ibescolumns::Union{Vector{Symbol}, Nothing} = PARAM[:ibescolumns])::DataFrame

  ibes::DataFrame = incsvstream("$ibespath\\$ibesname.$csvextension") |> CSV.File |> DataFrame
  rename!(ibes, (s::String->lowercase(s)).(names(ibes)))
  ibes.actdats = parseibes.(Date, ibes.actdats)
  ibes.anndats = parseibes.(Date, ibes.anndats)
  ibes.perioddate = parseibes.(Date, ibes.pends)

  #according to WRDS, announce date is used after 1/1/1999
  cutoffdate::Date = Date(1999,1,1)
  function pickdate(act,ann)
    #check if one is missing
    (ann === missing || act === missing) && return coalesce(ann, act, missing)

    ann ≥ cutoffdate && return ann
    return act
  end

  ibes.ibesadate = pickdate.(ibes.actdats, ibes.anndats)

  #read the cil linking file
  cil::DataFrame = CSV.read("$cilpath\\$cilname.csv", DataFrame)
  rename!(cil, (s::String->lowercase(s)).(names(cil)))
  cil = cil[cil.score .≤ 4, :] #only want likely matches- see WRDS documentation
  cil.edate = parsecil.(Date, cil.edate)
  cil.sdate = parsecil.(Date, cil.sdate)

  #NOTE: finish the linkage
  ibes = innerjoin(ibes, cil, on=:ticker)
  validlink(d, sdate, ::Missing) = true
  validlink(d, sdate, edate) = (d≥sdate) & (d≤edate)

  ibes = ibes[validlink.(ibes.ibesadate, ibes.sdate, ibes.edate), :]

  #we want the US currency row if there are multiple rows due to fx restatement
  ibes.uscurr = ibes.curr_act .=== "USD"
  sort!(ibes, [:permno, :perioddate, :uscurr])
  ibes = combine(last, groupby(ibes, [:permno, :perioddate]))

  #drop the suffix
  rename!(ibes, (n->n => (n[(end-4):end] === "_last" ? n[1:(end-5)] : n)).(names(ibes)))
  rename!(ibes, [:sdate=>:cilsdate, :edate=>:ciledate, :ticker=>:ibesticker])


  #now patch in a gvkey
  #onyl care that the link is valid on the proposed ibes adate
  ccm = prepccm(ccmcolumns=[:gvkey, :permno, :ccmeffdate, :ccmenddate])
  ibes = innerjoin(ibes, ccm, on=:permno)
  ibes = ibes[(ibes.ccmeffdate .≤ ibes.ibesadate) .& (ibes.ibesadate .≤ ibes.ccmenddate), :]

  #knock out dups by taking the first adate for any gvkey-period
  sort!(ibes, [:gvkey, :perioddate, :ibesadate])
  ibes = combine(first, groupby(ibes, [:gvkey, :perioddate]))

  #validate integrity of link
  @assert all((ibes.perioddate .=== lastdayofmonth.(ibes.perioddate)) .|
    (ibes.perioddate .=== missing))
  @assert ibes[nonunique(ibes,[:permno, :perioddate]),:] |> isempty
  @assert ibes[nonunique(ibes,[:gvkey, :perioddate]),:] |> isempty
  @assert ibes[nonunique(ibes,[:ibesticker, :perioddate]),:] |> isempty

  (!(ibescolumns === nothing)) && select!(ibes, ibescolumns)
  
  return ibes
end


parseibes(::Type{Date}, s::String, ibesdateformat::DateFormat=dateformat"yyyymmdd"
  ) = Dates.Date(s, ibesdateformat)
parseibes(::Type{<:Any}, ::Missing) = missing
parseibes(::Type{Date}, i::Int, ibesdateformat::DateFormat=dateformat"yyyymmdd") = parseibes(
  Date, "$i", ibesdateformat)

parsecil(::Type{<:Any}, ::Missing) = missing
parsecil(::Type{Date}, s::String, cildateformat::DateFormat=dateformat"dduuuyyyy"
  ) = Dates.Date(s, cildateformat)
