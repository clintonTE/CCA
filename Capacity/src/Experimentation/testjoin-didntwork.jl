using Revise

using Dates, DataFrames, BenchmarkTools, Random

import Base: ==

struct JoinableRange{Tuniformrange<:AbstractRange, Tval}
  uniformrange::Tuniformrange
  first::Tval
  last::Tval
end

JoinableRange(uniformrange::AbstractRange) = JoinableRange(
  uniformrange, minimum(uniformrange), maximum(uniformrange))

==(jr::JoinableRange, v) = ifelse((jr.first ≤ v) & (v ≤ jr.last), true, false)
==(jr::JoinableRange, ::Missing) = missing
==(v, jr::JoinableRange) = jr == v

function testfuzzydatejoin(N)
  df1 = DataFrame(grp=rand(collect(1:N), N))

  startdt = rand(collect(Date(1900,1,1):Day(1):Date(2019,12,31)), N)
  enddt = (startdt->rand(startdt:Day(1):Date(2020,12,31))).(startdt)

  rngs = ((startdt,enddt)->startdt:Day(1):enddt).(startdt,enddt)
  df1.jr = JoinableRange.(rngs)

  df2idx = shuffle!(collect(1:N))

  df2=DataFrame(grp=df1.grp[df2idx]|>deepcopy, dt = (jr->rand(jr.uniformrange)).(df1.jr[df2idx]))
  describe(df2) |> println

  df12 = innerjoin(df1,df2, on=[:grp, :jr=>:dt])

  #describe(df12) |> println

end

#@time testfuzzydatejoin(10^6)

using DataFrames, Dates

function mergecrspcompkernel(;crsp=error("crsp is required"), comp=error("comp is required"))
  ###Merge one permno at a time, linking each crsp row to a compid (comp row)
  scrsps = groupby(crsp, :permno)
  scomps = groupby(comp, :permno)
  comp.compid = 1:size(comp,1) |> collect
  crsp.compid = Vector{Union{Int, Missing}}(undef, size(crsp,1))

  #faster to work within the rows for a particular permno
  Threads.@threads for permno ∈ keys(scomps)
    scomp = scomps[permno]
    crsppermnokey = permno |> Tuple
    (!haskey(scrsps,crsppermnokey)) && continue
    #because the key originates from comp, cant use it directly on crsp- instead need its value
    scrsp = scrsps[crsppermnokey]

    #the idea here is to select all of the crsp rows between the effective date
    #and the end date of the comp row
    for r ∈ eachrow(scomp)
      targetrows = view(scrsp, (r.effdate .≤ scrsp.date) .& (scrsp.date .≤ r.enddate), :compid)
      if length(targetrows) > 0
        #under no circumstances should this be assigned, since date and permno should be unique
        (targetrows .=== missing) |> all || error(
          "multiple crsp rows found for a particular compid!!")
        targetrows .= r.compid
      end
    end
  end

  crsp = crsp[crsp.compid .!== missing, :]

  #now merge the linked rows
  rename!(comp, :permno=>:lpermno)
  panel = innerjoin(crsp, comp, on=:compid)

  #sanity checks
  @assert size(panel,1) == size(crsp,1)
  @assert (crsp.permno .== panel.permno) |> all
  @assert (crsp.date .== panel.date) |> all

  return panel
end
