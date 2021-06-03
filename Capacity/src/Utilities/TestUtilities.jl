#applies a function to a vector of trailing variables
function trailingmeasuretest!(F::Function,
  crsp::DataFrame;
  windowsize::Int=error("window size required"),
  Ftarget::Symbol,
  Fmeasurename::Symbol,
  offset::Int=0,
  maxcalendarinterval::DatePeriod = error("max calendar size required"),
  sorted::Bool = false)::Nothing

  N::Int = size(crsp,1)

  #allocate spcae for the new fields
  crsp[!,Fmeasurename] = Vector{MFloat64}(undef, N)

  #need to sort
  (!sorted) && (sort!(crsp, [:permno, :date]))

  #avoid missing data in the returns for simplicity
  scrspnomissing::SubDataFrame = view(crsp, (!ismissing).(crsp[!,Ftarget]),:) #view(crsp, :, :)
  scrsps::GroupedDataFrame = groupby(scrspnomissing, :permno)

  #now compute the trailing returns
  Threads.@threads for i ∈ 1:length(scrsps)
    local scrsp = scrsps[i]

    Nscrsp = size(scrsp, 1)
    if windowsize ≤ Nscrsp #if there are enough  entries to compute the trailing returns
      for r ∈ windowsize:Nscrsp
        r₀::Int = r-windowsize-offset+1 #this si the first row of the period
        rₜ::Int = r - offset
        staledate::Date = scrsp[r,:date] - maxcalendarinterval
        mindate::Date = scrsp[r₀,:date]

        #if the returns are not stale, record them
        if staledate < mindate
          scrsp[r, Fmeasurename] = F(scrsp[r₀:rₜ , Ftarget])
        end
      end
    end
  end

  return nothing
end
#computes the trailing returns from the logged returns
function trailingsigmatest!(crsp::DataFrame; sorted::Bool = false)

  print("computing trailing vol using test routine. Time: ")
  trailingmeasuretest!(std, crsp,
    Ftarget=:ret,
    windowsize=PARAM[:sigmawindow],
    maxcalendarinterval=PARAM[:sigmamaxcalendarinterval],
    sorted=sorted,
    Fmeasurename=PARAM[:sigmaname])::Nothing
  println(" checksum: $(sum(skipmissing(crsp.sigma)))")

  print("computing trailing vol using primary routine. Time: ")
  trailingmeasure!(std, crsp,
    Ftarget=:ret,
    windowsize=PARAM[:sigmawindow],
    maxcalendarinterval=PARAM[:sigmamaxcalendarinterval],
    sorted=sorted,
    Fmeasurename=PARAM[:sigmaname])::Nothing
  println(" checksum: $(sum(skipmissing(crsp.sigma)))")
  return nothing
end
