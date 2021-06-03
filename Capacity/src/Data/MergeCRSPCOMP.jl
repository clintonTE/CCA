


function mergecrspcomp(;crsp=error("crsp is required"), comp=error("comp is required"),
    testexpensivemergetests::Bool = PARAM[:testexpensivemergetests])

  local crspold::DataFrame
  local compold::DataFrame
  local panel::DataFrame

  #this is just another layer of verification
  if testexpensivemergetests
    crspold = crsp |> deepcopy
    compold = comp |> deepcopy
  end

  #start with integrity checks
  (!issorted(crsp, [:permno, :date])) && sort!(crsp, [:permno, :date])
  (!issorted(comp, [:permno, :effdate])) && sort!(comp, [:permno, :fdate])

  @assert issorted(comp, [:permno, :fdate])
  @assert (nonunique(crsp, [:permno,:date]) .== false) |> all
  @assert (nonunique(comp, [:permno, :effdate]) .== false) |> all


  Ncrsp::Int = size(crsp,1)
  print("Beginning merge. Num crsp rows = $Ncrsp. . Time to complete: ")
  @time begin
    ###Merge one permno at a time, linking each crsp row to a compid (comp row)
    scrsps = groupby(crsp, :permno)
    scomps = groupby(comp, :permno)
    comp.compid = 1:size(comp,1) |> collect
    crsp.compid = Vector{MInt}(undef, size(crsp,1))

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
  end

  #integrity checks on the merge. Run these each time.
  @assert size(panel,1) == size(crsp,1)
  @assert (crsp.permno .== panel.permno) |> all
  @assert (crsp.date .== panel.date) |> all
  @assert (panel.permno .== panel.lpermno) |> all
  @assert (panel.effdate .≤ panel.date) |> all
  @assert (panel.date .≤ panel.enddate) |> all
  if !issorted(panel, [:permno, :date])
    sort!(panel, [:permno, :date])
  end

  @info "
  %%%%%%%%%%
    merge complete. Merged df dimensions = $(size(panel)),
    sizeof=$(sum(sizeof.(eachcol(panel)))),
    Panel rows = $(size(panel,1)),
    Match rate=$(round(size(panel,1)/Ncrsp*100,digits=4))%.
  %%%%%%%%%%%  "

  #print(describe(panel))

  testexpensivemergetests && expensivemergetests(crsp=crspold, comp=compold; panel)

  #after the integrity checks we no longer need the below columns
  select!(panel, Not([:lpermno, :effdate, :enddate]))

  #next expand out gvkey
end

#=function mergecrspcompkernel(;crsp=error("crsp is required"), comp=error("comp is required"))


  #integrity checks


  return panel
end=#

######################
###############Everything below this point relates to expensive merge comparison tests
######################


function expensivemergetests(; crsp=error("crsp is required"), comp=error("comp is required"),
  panel=error("panel is required for verification"))

  @warn "Running old merge code for verification.
  This may take a while. If this passes, no need to run this each time"

  oldmergecrspcomp(crsp=crsp|>deepcopy, comp=comp|>deepcopy; panel)
  noviewmergecrspcomp(; crsp, comp, verpanel=panel)


  @info "Completed comparison to crspcomp successfully- results fully validated"
end

#slower code (by ~20x) but clearer and more elegant
function noviewmergecrspcomp(;crsp=error("crsp is required"), comp=error("comp is required"),
  verpanel::DataFrame)


  local panel::DataFrame

  #start with integrity checks
  (!issorted(crsp, [:permno, :date])) && sort!(crsp, [:permno, :date])
  (!issorted(comp, [:permno, :effdate])) && sort!(comp, [:permno, :fdate])

  @assert issorted(comp, [:permno, :fdate])
  @assert (nonunique(crsp, [:permno,:date]) .== false) |> all
  @assert (nonunique(comp, [:permno, :effdate]) .== false) |> all


  Ncrsp::Int = size(crsp,1)
  comp.effrange = ((effdate, enddate)->effdate:Day(1):enddate).(comp.effdate, comp.enddate)
  print("Beginning NO-VIEW merge test. Num crsp rows = $Ncrsp. Time to complete: ")
  @time begin
  ###first, form an index for comp
  #build the index one column at a time
    comp.compid = 1:size(comp,1) |> collect
    fillrange(fr) = fill(fr[1], length(fr[2])) #repeats the values the appropriate number of times
    permno = mapreduce(fillrange,  append!, zip(comp.permno, comp.effrange), init=Vector{Int}())
    compdate = mapreduce(collect, append!, comp.effrange, init=Vector{Date}())
    compid = mapreduce(fillrange,  append!, zip(comp.compid, comp.effrange), init=Vector{Int}())
    #put the columns together
    compidx = DataFrame(permno=permno, date=compdate, compid=compid,  copycols=false)

    #integrity checks
    @assert issorted(compidx, [:permno, :date])
    @assert issorted(compidx, [:compid])
    @assert (nonunique(compidx, [:compid, :date]) .== false) |> all
    @assert (nonunique(compidx, [:permno, :date]) .== false) |> all

    #merge crsp and the comp idx
    crsp = innerjoin(crsp, compidx, on=[:permno, :date], validate=(true,true))

    #merge crsp and comp using compid
    rename!(comp, :permno=>:lpermno)
    panel = innerjoin(crsp, comp, on=:compid, validate=(false, true))
  end

  #integrity checks on the merge
  @assert size(panel,1) == size(crsp,1)
  @assert (crsp.permno .== panel.permno) |> all
  @assert (crsp.date .== panel.date) |> all
  @assert (panel.permno .== panel.lpermno) |> all
  @assert (panel.effdate .≤ panel.date) |> all
  @assert (panel.date .≤ panel.enddate) |> all
  if !issorted(panel, [:permno, :date])
    sort!(panel, [:permno, :date])
  end

  @assert (panel.permno .== verpanel.permno) |> all
  @assert (panel.date .== verpanel.date) |> all
  @assert (panel.effdate .== verpanel.effdate) |> all
  @assert (panel.enddate .== verpanel.enddate) |> all
  @assert (panel.gvkey .== verpanel.gvkey) |> all
  @assert (panel.compid .== verpanel.compid) |> all

  @info "NOVIEW merge complete. Match rate=$(round(size(panel,1)/Ncrsp*100,digits=4))%."

  return nothing
end


#some old merge code
#to be primarily used for occational verification
function oldmergecrspcomp(; crsp::DataFrame=error("original crsp df is required"),
    comp::DataFrame=error("original crsp df is required"),
    panel::DataFrame=error("verification panel required"))

  local univ::DataFrame
  local scomp::SubDataFrame
  local Ncrsp::Int = size(crsp,1)


  comp.lpermno = comp.permno |> deepcopy

  @assert issorted(crsp, [:permno, :date])
  @assert issorted(comp, [:lpermno, :effdate])
  leadwithin2!(comp, [:adate], :lpermno, date=:fdate)

  #create a record id to form the dictionary
  comp.compid = Vector{MInt}(1:(size(comp,1)))
  comp.validinterval = ((effdate,enddate)->
    (effdate):Day(1):(enddate)).(comp.effdate, comp.enddate)

  crsp.compid = Vector{MInt}(undef, Ncrsp)
  select!(comp, Not(:permno))


  #make the necessary allocations ahead of time
  print("matching...time to match (OLD CODE): ")
  @time begin

    scrsps::GroupedDataFrame = groupby(crsp, :permno)
    scomps::Dict = Dict([scomp.lpermno[1]=>scomp for scomp ∈ groupby(comp, :lpermno)])

    #this uses a similar matching algorithm to the main code, but uses older syntax from when
    #data frames were not as good.
    #so basically we first select a permno, match it to the comp permno, then
    #match each comp row to the crsp rows within the comp row's date range
    Threads.@threads for i ∈ 1:length(scrsps)
      local scrsp::SubDataFrame = scrsps[i]
      local permno::Int = scrsp.permno[1]

      if haskey(scomps, permno)
        local scomp = scomps[permno]
        for rscomp ∈ eachrow(scomp)
          effdate::Date =  rscomp[:effdate]
          enddate::Date = rscomp[:enddate]
          sscrsp::SubDataFrame = view(scrsp, (d->d≥effdate && d≤enddate).(scrsp[!, :date]),:)
          if size(sscrsp,1) ≥ 1
            (sum((!ismissing).(sscrsp.compid)) > 0) && error("multiple valid compids for one permno")
            sscrsp.compid .= rscomp[:compid]
          end
        end
      end
    end

    #println("sum crsp compid: $(sum(skipmissing(crsp.compid)))")
    println("OLD CODE matching complete: $(100*(Ncrsp - sum((ismissing).(crsp.compid)))/Ncrsp)% matched.")
    crsp = crsp[(!ismissing).(crsp.compid),:]
    univ = innerjoin(crsp, comp, on=:compid)
  end

  extracheckstestmerged(univ)
  println("OLD CODE Total records: $(size(univ))")

  if !issorted(univ, [:permno, :date])
    sort!(univ, [:permno, :date])
  end
  @assert size(univ,1) == size(panel,1)
  @assert (univ.permno .== panel.permno) |> all
  @assert (univ.adate .== panel.adate) |> all
  @assert (univ.gvkey .== panel.gvkey) |> all
  @assert (univ.date .== panel.date) |> all

  #=if TEST_OUTPUT
    d = Date(2015,12,31)
    CSV.write("output\\univ_mid_$d.csv", univ[(cd->cd==d).(univ[!, :date]), :])
  end=#

  return nothing
end

#do a bunch of integrity checks on the merged dataset
function extracheckstestmerged(univ::DataFrame;
  period2stale = PARAM[:compmaxadateinterval])

  Threads.@threads for r ∈ eachrow(univ)
    ismissing(r.gvkey) && continue
    #check date is valid
    d::Date = r[:date]
    (d ≥ r.ccmeffdate) || error("failed check: (d ≥ r.linkeffdate)")
    (d ≤ r.ccmenddate) || error("failed check: (d ≤ r.linkenddate)")
    (d ∈ r.validinterval) || error("failed check: (d ∈ r.validinterval)")
    (d ≥ r.adate) || error("failed check: (d > r.adate)")
    (r.Nadate === missing) || (d < r.Nadate) || error("failed check: (d ≤ r.Nadate)")

    (!(d ≤ lastdayofmonth(firstdayofmonth(r.adate) + period2stale))) && error(
      "lastdayofmonth(firstdayofmonth(r.adate) + period2stale)))
      $r")

    (r.lpermno == r.permno) || error("failed check: (r.lpermno == r.permno)")
  end

  #univ = sort(univ, [:permno, :adate, Fdate])

  sunivs::GroupedDataFrame = groupby(univ, :permno)
  Threads.@threads for i ∈ 1:length(sunivs)
    local suniv::SubDataFrame = sunivs[i]
    maxfiscal::Date= Date(1900,1,1)
    maxadate::Date = Date(1900,1,1)

    for r ∈ eachrow(suniv)
      ismissing(r.gvkey) && continue
      (r[:date] > maxadate) || error("failed check: r[Fdate] > maxadate
        permno=$(r[:permno]), r[Fdate]=$(r[Fdate]), maxadate=$maxadate\n\n$r")

      maxadate = r[:date]

      (maxfiscal ≤ r.fdate) ||  error("failed check: maxfiscal ≤ r.fdate")
      (maxfiscal ≠ r.fdate) && (maxfiscal=r.fdate)
    end

  end
end
