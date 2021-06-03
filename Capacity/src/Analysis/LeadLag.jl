

function analyzeleadlag(;measure::DataFrame, measureinfo::NamedTuple,
  llname = PARAM[:llname])
  @unpack outputpath, resultid, returntype, datefrequency = measureinfo

  measure, llinfo = prepleadlag(; measure, measureinfo)
  llresults = analyzeleadlag(measure, measureinfo, llinfo)

  sort!(llresults, :betaxr)

  llresults |> CSV.write("$outputpath\\$(resultid)_$(llname)_results.csv")
  measure |> CSV.write("$outputpath\\$(resultid)_$(llname)_series.csv")

  return nothing
end

#this lags and leads the focal measure and return columns
function prepleadlag(;
  measure::DataFrame,
  measureinfo::NamedTuple,
  measureleadlags::Vector{Int} = PARAM[:llmeasureleadlags][measureinfo.datefrequency],
  returnwindows::Vector{Int} = PARAM[:llreturnwindows][measureinfo.datefrequency],
  additionalreturncolsidx = PARAM[:lladditionalreturncolsidx],
  additionalreturncols = additionalreturncolsidx |> keys |> collect,
  measurecolpatterns = PARAM[:llmeasurecolpatterns],
  detrendll = PARAM[:lldetrendll],
  additionallagcols = PARAM[:lladditionallagcols],
  controls= PARAM[:llcontrols])

  #NOTE: the below are a couple fast sanity checks
  #testshift(10)
  #testforwardreturnwindow(24)

  initialreturncols = [findcols(measure, :return); additionalreturncols]
  initialmeasurecols = [(p->findcols(measure, p)).(measurecolpatterns)...;]

  @assert !isempty(initialmeasurecols)
  @assert !isempty(initialreturncols)

  #the below will include shifted cols and variations on the return window
  measurecols = Symbol[]
  returncols = Symbol[]


  additionallagcols = [(s->findcols(measure, s)).(additionallagcols)...;]
  #fore each lead and lag and each measure col, apply the appropriate lagging
  for l ∈ measureleadlags
    prefix = Symbol(l>0 ? :N : :L, abs(l))
    for mcol ∈ initialmeasurecols
      shiftedmcol = Symbol(prefix, mcol)
      measure[!, shiftedmcol] = shift(measure[!, mcol], l)
      push!(measurecols, shiftedmcol) #record the newly created column
    end

    for lcol ∈ additionallagcols
      shiftedlcol = Symbol(prefix, lcol)
      measure[!, shiftedlcol] = shift(measure[!, lcol], l)
    end
  end

  #create the return windows
  for windowsize ∈ returnwindows
    prefix = Symbol(:FW, windowsize)
    for rcol ∈ initialreturncols
      returnwindowcol = Symbol(
        prefix, rcol ∈ additionalreturncols ? "_" : "", rcol)
      measure[!, returnwindowcol] = forwardreturnwindow(
        measure[!, rcol]; windowsize, logged=false)
      push!(returncols, returnwindowcol)
    end
  end

  #this will hold the lead lag pairs
  estimations = NamedTuple{(:rcol, :mcol, :ccols, :lag)}[]

  function extractlag(colname::String)
    prefix =  match(r"^L[0-9]+[A-Za-z]*_",colname)
    lagstring = match(r"[0-9]+", prefix)
    return parse(Int, lagstring)
  end
  extractlag(colname::Symbol) = extractlag(colname |> string)


  for rcol ∈ returncols, mcol ∈ measurecols
    local makepair
    strippedrcol = stripcolprefix(rcol)|>string
    strippedmcol = stripcolprefix(mcol)|>string

    #if the stripped versions match assume a pairing is desired
    if occursin(strippedrcol, strippedmcol)
      makepair=true
    #see if we explicitly want a pairing
    elseif (haskey(additionalreturncolsidx, strippedrcol |> Symbol) &&
        any((m->occursin(m|>string, strippedmcol)).(
          additionalreturncolsidx[strippedrcol |> Symbol])))
        makepair=true
    else
      makepair = false
    end
  #println("mcol: $mcol, rcol: $rcol,",
    #  " strippedmcol: $strippedmcol, strippedrcol: $strippedrcol",
    #  " makepair: $makepair")

    if makepair
      for control ∈ controls
        @unpack ccols, controlmatch = control
          if (controlmatch≡nothing) || (match(Regex(controlmatch), mcol |> string) !== nothing)
            if mcol ∈ ccols
              (length(ccols)==1) && continue #this means none of the controls are valid
              setdiff!(ccols, [mcol])
            end
            push!(estimations, (; rcol, mcol, ccols))
          end
        end
    end
  end

  #linearly detrends the results
  if detrendll
    @assert issorted(measure, [:date])

    originalestimations = estimations |> deepcopy
    for orig ∈ originalestimations
      @unpack rcol,mcol, ccols = orig

      #for now do not detrend the controls
      dtrcol, dtmcol = Symbol(:dtX,rcol), Symbol(:dtX,mcol)

      #do the actual detrending if necessary, being careful of missing values
      if dtrcol ∉ propertynames(measure)
        measure[!, dtrcol] = missings(Float64, nrow(measure))
        cleanmeasure = view(measure, completecases(measure, [rcol;]), :)
        cleanmeasure[:, dtrcol] .= detrend(cleanmeasure[:, rcol])
      end

      if dtmcol ∉ propertynames(measure)
        measure[!, dtmcol] = missings(Float64, nrow(measure))
        cleanmeasure = view(measure, completecases(measure, [mcol;]), :)
        cleanmeasure[:, dtmcol] .= detrend(cleanmeasure[:, mcol])
      end


      push!(estimations, (;rcol = dtrcol, mcol = dtmcol, ccols))
    end
  end







  llinfo = (; returncols, measurecols, measureleadlags, returnwindows, estimations)

  measure |> CSV.write(
    "$(PARAM[:testpath])\\$(PARAM[:analysismeasurefilename])_prepll.csv")
  return measure, llinfo
end

#function for simple linear detrending
function detrend(col; idx = 1:length(col) |> collect .|> Float64, idx1 = [idx ones(length(col))])
  @assert length(col) == length(idx) == size(idx1,1)

  m = cov(col, idx)/var(idx)
  b = mean(col) - m*mean(idx)

  chk = FMLM(idx1, col |> Vector{Float64})
  @assert m+10.0 ≈ chk.β[1]+10.0 "m=$m, chk.β[1]=$(chk.β[1])"
  @assert b+10.0 ≈ chk.β[2]+10.0

  return col .- (idx .* m .+ b)
end

function shift(v::AbstractVector{MFloat64}, i::Int)
  vout = circshift(v,-i)
  vout[ifelse(i>0, (end-i+1):end, 1:abs(i))] .= missing
  return vout
end

shift(v::AbstractVector, i::Int) = (i > 0 ?
  [v[(i+1):end]; missings(Float64,i);] :
  [missings(Float64, abs(i)); v[1:(end+i)];])

#averages returns over a forward window
function forwardreturnwindow(vin::AbstractVector{T}; windowsize, logged::Bool) where T
  v = deepcopy(vin)
  (windowsize==1) && return v #this is an identity scenario

  if !logged #compute the log returns
    @assert (vᵢ->vᵢ≡missing || vᵢ > -1.).(v) |> all
    v .= log.(v .+ 1.0)
  end

  #compute the window of lgoged returns
  cumlogret = v |> deepcopy |> Vector{Union{T,Missing}}

  for t ∈ 2:windowsize
    shiftedv = shift(v, t-1)
    @assert (shiftedv[(end  - t + 2) : end] .≡ missing) |> all
    cumlogret .+= shiftedv
  end

  logged && return cumlogret #unlog the returns
  return exp.(cumlogret) .- 1
end

function testshift(N=10)
  #throw("something is off, maybe with circshift")
  v = rand(N)
  @assert (shift(v,3) .≡ [v[4:end]; missings(3)]) |> all
  @assert (shift(v,-3) .≡ [missings(3); v[1:(end-3)]]) |> all
  @assert (shift(v,0) .≡ v) |> all

  vm = v |> Vector{MFloat64}
  @assert (shift(vm,3) .≡ [v[4:end]; missings(3)]) |> all "vm: $(shift(vm,3)); ver: $([v[4:end]; missings(3)])"
  @assert (shift(vm,-3) .≡ [missings(3); v[1:(end-3)]]) |> all
  @assert (shift(vm,0) .≡ v) |> all
end

function testforwardreturnwindow(N=24)
  v = rand(N) .- 0.5
  N1v = [v[2:end]; missing]
  N2v = [v[3:end]; missing; missing]

  approxormissing(x,y) = ((x≡y) | coalesce((x≈y),false))

  @assert (v .≈
    forwardreturnwindow(v,windowsize=1, logged=false)) |> all
  @assert approxormissing.((1 .+ v) .* (1 .+ N1v) .- 1,
    forwardreturnwindow(v,windowsize=2, logged=false)) |> all
  @assert approxormissing.((1 .+ v) .* (1 .+ N1v) .* (1 .+ N2v) .- 1,
    forwardreturnwindow(v,windowsize=3, logged=false)) |> all

  @assert (v .≈
    forwardreturnwindow(v,windowsize=1, logged=true)) |> all
  @assert approxormissing.((v .+ N1v),
    forwardreturnwindow(v,windowsize=2, logged=true)) |> all
  @assert approxormissing.((v .+ N1v .+ N2v),
    forwardreturnwindow(v,windowsize=3, logged=true)) |> all
end

#currently only analyzes the lag
function analyzeleadlag(measure, measureinfo, llinfo)

  #=function bsstatfunc(sim)
    r,m = sim[1], sim[2]

    ρ = cor(r,m)
    return Dict(
      :betaxr => ρ * std(r)/std(m),
      #:betamr => ρ * std(m)/std(r)
      )
  end

  nullar1statfunc(r,m) = bsstatfunc((r,m))=#
  function nullar1statfunc(Xy)
    X = @view Xy[:, 1:(end-1)]
    y = @view Xy[:, end]
    #println("X[1:3,:]", X[1:3,:])
    #println("r[1:3]", r[1:3])
    Dict(:betaxr => cholesky!(Symmetric(X'*X))\(X'*y),)
  end

  @unpack (neweylags, nwΣ,
    bsaggfunc, bsmeasurecolsregex, bsbenchmarkcolsregex, bssamples, bswidth,
      nullar1aggfunc, nullar1simulations) = measureinfo
  @unpack returncols, measurecols, measureleadlags, returnwindows, estimations = llinfo

  #create a container for the resutls
  results =   DataFrame(
    label=Symbol[], #label in case we concatenate this df with others
    controls=String[],
    measure=Symbol[], #measure field that we are benchmarking
    returnwindow=Symbol[], #benchmark field

    betaxr=MFloat64[],
      betaxrmw=MFloat64[],
      betaxrhac=MFloat64[],
      betaxrmbbp=MFloat64[],
      betaxrmbbz=MFloat64[],
      betaxrar1h0z=MFloat64[],
      betaxrar1h0p=MFloat64[],
      ibetaxrar1h0z=MFloat64[],
      ibetaxrar1h0p=MFloat64[],
      N=MInt[])

  #iterate through the cartesian product of column pairins
  for estimation ∈ estimations

    @unpack rcol, mcol, ccols = estimation
    teststr ="testing $rcol ~ $mcol + $ccols"

    occursin(r"^N", string(mcol)) && continue #an N would imply a lead column, not a lag column

    #make sure we have no missing entries
    jointnomissing = view(measure, completecases(measure, [rcol;mcol; ccols;]),
      [:date; rcol;mcol;ccols;])

    N=nrow(jointnomissing)
    r::Vector{Float64} = jointnomissing[:,rcol]
    #println([mcol, ccols])
    Xm::Matrix{Float64} = hcat(jointnomissing[:,[mcol; ccols]], ones(N)) |> Matrix{Float64}


    #pear = cor(r,m) #pearson correlation

    lmrm = FMLM(Xm, r)
    betaxr=lmrm.β[1]
    @assert ((lmrm.β[1] ≈ (qr(Xm)\r)[1]) || (lmrm.β[1] ≡ (qr(Xm)\r)[1]))
    betaxrmw = betaxr/modifiedwhiteΣ!(lmrm)[1,1]^0.5
    betaxrhac = betaxr/nwΣ(lmrm)[1,1]^0.5

    #to save time only use the bootstrap and simulation some of the time
    if occursin(bsmeasurecolsregex, string(mcol))
      bs = regmbb(Xm, r; bssamples, bswidth)

      #check dimensions of results line up
      @assert length(bs.beta.se) == length(bs.beta.p) == length(ccols) + 2 "
        bs.beta: $(bs.beta) ccols: $ccols"

      betaxrmbbz =  betaxr/bs.beta.se[1]
      betaxrmbbp =  bs.beta.p[1]

      println(teststr)
      ar1h0teststats = Dict(:betaxr => betaxr#=, :betamr => betamr=#)
      ar1h0 = comparear1undernull(Xm,r, nullar1statfunc,
        (M,s)->nullar1aggfunc(M,s,ar1h0teststats),
        Nsimulations=nullar1simulations)
      #println("$teststr: ar1 complete")

      @assert length(ar1h0.betaxr.se) == length(ar1h0.betaxr.p) == length(ccols) + 2

      betaxrar1h0z =  betaxr/ar1h0.betaxr.se[1]
      betaxrar1h0p =  ar1h0.betaxr.p[1]

      iar1h0 = comparearima11undernull(Xm, r, nullar1statfunc,
        (M,s)->nullar1aggfunc(M,s,ar1h0teststats),
        Nsimulations=nullar1simulations)
      #println("$teststr: iar1 complete")
      @assert length(iar1h0.betaxr.se) == length(iar1h0.betaxr.p) == length(ccols) + 2

      ibetaxrar1h0z =  betaxr/iar1h0.betaxr.se[1]
      ibetaxrar1h0p =  iar1h0.betaxr.p[1]

    else
      betaxrmbbp=missing;
      betaxrmbbz=missing;

      betaxrar1h0p=missing;
      betaxrar1h0z=missing;

      ibetaxrar1h0p=missing;
      ibetaxrar1h0z=missing;
    end

    push!(results,(;
      label=:predict,
      measure=mcol,
      returnwindow=rcol,
      controls = join(string.(ccols), ", "),
      betaxr, betaxrmw, betaxrhac, betaxrmbbz,  betaxrmbbp,
      betaxrar1h0z,  betaxrar1h0p, ibetaxrar1h0z,  ibetaxrar1h0p,
      N))
  end


  return  results


  #form the prediction pairs
  #first match on LP

end
