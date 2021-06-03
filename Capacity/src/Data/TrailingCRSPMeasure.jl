

function trailingmeasure!(F::TF, crsp::DataFrame, calendarinterval::DatePeriod,
  minfraction::Float64;
  Ftarget::Symbol = error("Ftarget is required"),
  Fmeasure::Symbol = error("Fmeasurename is required"),
  Fconditional::Union{Nothing, Symbol} = nothing,)::Nothing where TF

  applyconditional::Bool = (!isnothing(Fconditional))
  crspnames::Vector{String} = names(crsp)

  #check and make sure the target doesn't already have data in it (or somethign is probably wrong)
  (sum((ismissing).(crsp[!,Fmeasure])) == size(crsp,1)) ||  (error("$Fnew already contains data"))

  #make sure everything is sorted
  (issorted(crsp, [:permno, :date])) || error("dataframe not sorted by [:permno, :date]")

  #find the minimum date in the entire dataset
  alldates = unique(crsp.date) |> sort!
  initialminpoints = round(
    searchsortedfirst(alldates, alldates[1] + calendarinterval)*minfraction,RoundUp) |> Int

  #now compute the min number of points for each date
  function getminpoints(i)
    (alldates[i] ≤ alldates[1] + calendarinterval) && return initialminpoints
    lowerdateindex = searchsortedfirst(alldates, alldates[i] + Day(1) - calendarinterval)
    return round((i - lowerdateindex + 1) * minfraction,RoundUp) |> Int
  end
  minpointindex = Dict{Date, Int}(
    alldates[i] => getminpoints(i) for i ∈ initialminpoints:length(alldates))
  minstartdate = alldates[initialminpoints]
  @assert (alldates[alldates .≥ minstartdate] .|> (d)->haskey(minpointindex, d)) |> all

  #split and apply by permno
  requiredfields::Vector{Symbol} = [Ftarget; applyconditional ? Fconditional : []]
  completecrsp::SubDataFrame = view(crsp, completecases(crsp, requiredfields), :)
  scrsps::GroupedDataFrame = groupby(completecrsp, :permno)

  #when we are looking at eom values only
  @inline getj₀(d::Date, dates::AbstractVector{Date}) = searchsortedfirst(
    dates, d + Day(1) - calendarinterval)

  #the below is a low-effort performance booster- don't check arrays below this
  absminpoints::Int = minimum(values(minpointindex) |> collect)

  #now compute the trailing returns
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    nrows = size(scrsp, 1)
    (absminpoints < nrows) || continue #if there are enough  entries to compute the trailing returns

    #a hack to only compute calcualtions at end of month
    for j ∈ absminpoints:nrows
      applyconditional && (!scrsp[j,Fconditional]) && continue #if we are only computing end of month values
      currentdate = scrsp[j, :date]
      (currentdate < minstartdate) && continue
      j₀::Int = getj₀(currentdate, scrsp.date)
      N::Int = j-j₀ +1
      #if the returns are not stale and complete enough, record them
      (N < minpointindex[currentdate]) && continue
      scrsp[j, Fmeasure] = F(scrsp[j₀:j , Ftarget])
    end
  end

  return nothing
end

#applies a function to a vector of trailing variables
#this version uses a rolling calendar interval with a set number
#of minimum points
function trailingmeasure!(F::Function, crsp::DataFrame, calendarinterval::DatePeriod,
  minpoints::Int;
  Ftarget::Symbol = error("Ftarget is required"),
  Fmeasure::Symbol = error("Fmeasurename is required"),
  Fconditional::Union{Nothing, Symbol} = nothing)::Nothing

  applyconditional::Bool = (!isnothing(Fconditional))
  crspnames::Vector{String} = names(crsp)

  #check and make sure the target doesn't already have data in it (or somethign is probably wrong)
  (sum((ismissing).(crsp[!,Fmeasure])) == size(crsp,1)) ||  (error("$Fmeasure already contains data"))

  #make sure everything is sorted
  (issorted(crsp, [:permno, :date])) || error("dataframe not sorted by [:permno, :date]")

  #split and apply by permno
  requiredfields::Vector{Symbol} = [Ftarget; applyconditional ? Fconditional : []]
  completecrsp::SubDataFrame = view(crsp, completecases(crsp, requiredfields), :)
  scrsps::GroupedDataFrame = groupby(completecrsp, :permno)

  #when we are looking at eom values only
  @inline getj₀(d::Date, dates::AbstractVector{Date}) = searchsortedfirst(
    dates, d + Day(1) - calendarinterval)

  #now compute the trailing returns
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    nrows = size(scrsp, 1)
    (minpoints ≤ nrows) || continue #if there are enough  entries to compute the trailing returns

    #a hack to only compute calcualtions at end of month
    for j ∈ minpoints:nrows
      applyconditional && (!scrsp[j,Fconditional]) && continue #if we are only computing end of month values
      j₀::Int = getj₀(scrsp[j, :date], scrsp.date)
      N::Int = j-j₀ +1
      #if the returns are not stale and complete enough, record them
      (N < minpoints) && continue
      scrsp[j, Fmeasure] = F(scrsp[j₀:j , Ftarget])
    end
  end

  return nothing
end

#trailing index measure
function stricttrailingmeasure!(F::TF, crsp::DataFrame, calendarinterval::DatePeriod;
  Ftarget::Symbol = error("Ftarget is required"),
  Fmeasure::Symbol = error("Fmeasurename is required"),
  Fconditional::Union{Nothing, Symbol} = nothing,)::Nothing where TF

  applyconditional::Bool = (!isnothing(Fconditional))
  crspnames::Vector{String} = names(crsp)

  #check and make sure the target doesn't already have data in it (or somethign is probably wrong)
  (sum((ismissing).(crsp[!,Fmeasure])) == size(crsp,1)) ||  (error("$Fnew already contains data"))

  #make sure everything is sorted
  (issorted(crsp, [:permno, :date])) || error("dataframe not sorted by [:permno, :date]")

  #find the lagged dates in the entire dataset
  alldates = unique(crsp.date) |> sort!
  alldatesidx = Dict{Date, Int}(alldates[i] => i for i ∈ 1:length(alldates))
  initialminpoints = searchsortedfirst(alldates, alldates[1] + calendarinterval) |> Int

  #now compute the lagged dates
  getlaggedpoint(i) = alldates[searchsortedfirst(alldates, alldates[i] - calendarinterval)]

  laggeddateidx = Dict{Date, Date}(
    alldates[i] => getlaggedpoint(i) for i ∈ initialminpoints:length(alldates))
  minstartdate = alldates[initialminpoints]
  @assert (alldates[alldates .≥ minstartdate] .|> (d)->haskey(laggeddateidx, d)) |> all

  #low effort check to only look at rows after the minimum days across all intervals
  minintervaldays = minimum(((k,v)->alldatesidx[k]-alldatesidx[v]
    ).(keys(laggeddateidx), values(laggeddateidx)))

  #split and apply by permno
  requiredfields::Vector{Symbol} = [Ftarget; applyconditional ? Fconditional : []]
  completecrsp::SubDataFrame = view(crsp, completecases(crsp, requiredfields), :)
  scrsps::GroupedDataFrame = groupby(completecrsp, :permno)

  #now compute the trailing returns
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    nrows = length(scrsp.date)

    #if there aren't enough dates within a permno to create lagged dates, skip the permno
    (firstdate, lastdate) = (scrsp[1,:date], scrsp[end,:date])
    if ((nrows < minintervaldays) ||
        (lastdate < minstartdate) ||
        (laggeddateidx[lastdate] < firstdate))

      continue
    end

    #a hack to only compute calcualtions at end of month
    for j ∈ minintervaldays:length(scrsp.date)
      #if we are only computing end of month values, make sure we are at a month end
      applyconditional && (!scrsp[j,Fconditional]) && continue
      currentdate = scrsp[j, :date]

       #skip if currentdate is before the first entry in the lookup table
      currentdate < minstartdate && continue

      #now check if the lagged date is a valid index entry for this permno
      targetlaggeddate = laggeddateidx[currentdate]
      j₀ = searchsortedfirst(scrsp.date, targetlaggeddate)
      actuallaggeddate = scrsp[j₀, :date]
      (actuallaggeddate ≠ targetlaggeddate) && continue
      scrsp[j, Fmeasure] = F(scrsp[j₀:j , Ftarget])
    end
  end

  return nothing
end

#applies a function to a vector of trailing variables
#this version uses a rolling window computed by the number of entries as opposed to
#one computed by the dates
function trailingmeasure!(F::TF, crsp::DataFrame, windowsize::Int;
  Ftarget::Symbol = throw("Ftarget is required"),
  Fmeasure::Symbol = throw("Fmeasurename is required"),
  Fconditional::Union{Nothing, Symbol} = nothing,
  maxcalendarinterval::DatePeriod = throw("maxcalendarinterval is required")
  )::Nothing where TF

  applyconditional::Bool = (!isnothing(Fconditional))
  crspnames::Vector{Symbol} = names(crsp)

  #check and make sure the target doesn't already have data in it (or somethign is probably wrong)
  (sum((ismissing).(crsp[!,Fmeasure])) == size(crsp,1)) ||  (throw("$Fnew already contains data"))

  #make sure everything is sorted
  (issorted(crsp, [:permno, :date])) || throw("dataframe not sorted by [:permno, :date]")

  #split and apply by permno
  requiredfields::Vector{Symbol} = [Ftarget; applyconditional ? Fconditional : []]
  completecrsp::SubDataFrame = view(crsp, completecases(crsp, requiredfields), :)
  scrsps::GroupedDataFrame = groupby(completecrsp, :permno)

  #now compute the trailing returns
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    nrows = size(scrsp, 1)
    (nrows<windowsize) || continue #if there are enough  entries to compute the trailing returns

    for j ∈ windowsize:nrows
      applyconditional && (!scrsp[j,Fconditional]) && continue #if we are only computing end of month values
      staledate::Date = scrsp[j, :date] - maxcalendarinterval
      j₀::Int = j-windowsize+1

      (scrsp[j₀, :date] ≤ staledate)  && continue #elliminates stale windows
      scrsp[j, Fmeasurename] = F(scrsp[j₀:j, Ftarget])
    end
  end

  return nothing
end
