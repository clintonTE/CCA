

iterresultsname(λ⁺, ::Val{:ols}) = "ans_$(PARAM[:crspdatalabel])" *
  boundrangestr() *
  (PARAM[:placebovol] ? "_SHUFFLED" : "") *
  "_$(PARAM[:measuretype])" *
  "_$(PARAM[:itersolvetype])" *
  "_$(PARAM[:iterresultssuffix])" *
  "_$(Int(round(λ⁺)))" *
  "_$(Dates.format(now(),"yyyymmdd_HHMM"))"
iterresultsname(::Val{:ols}) = "ans_$(PARAM[:crspdatalabel])" *
  boundrangestr() *
  (PARAM[:placebovol] ? "_SHUFFLED" : "") *
  "_$(PARAM[:measuretype])" *
  "_$(PARAM[:itersolvetype])" *
  "_$(PARAM[:iterresultssuffix])" *
  "_$(Dates.format(now(),"yyyymmdd_HHMM"))"

iterresultsname(λ⁺) = iterresultsname(λ⁺, Val(PARAM[:itermodel]))
iterresultsname() = iterresultsname(Val(PARAM[:itermodel]))

iterresultsname(::Val{:bootstrap}) = iterresultsname(Val(:ols))
iterresultsname(λ⁺, ::Val{:bootstrap}) = iterresultsname(λ⁺, Val(:ols))

#optresultsname() = optresultsname(Val(PARAM[:itermodel]))
#primary function for estimating the measure using optimizations
#note that this still uses the same volume parts construct as flux
function formiterregression(panel::DataFrame, ms::MeasureSpec,
  IterModel::Val = Val(PARAM[:itermodel]);
  returns::DataFrame, funds,
  estimationmindate= PARAM[:globalmindate],
  estimationmaxdate= PARAM[:globalmaxdate])

  #issorted(panel, [:permno, :date]) || error("panel must be sorted")
  sort!(panel, :date)

  #bound the dates
  #the lag of the first date is generally expected to be missing
  #if there is no cutoff then it shouldn't matter, but if the dates are constrained then it will
  measuredatekey = Symbol(PARAM[:crspfrequencysuffix],PARAM[:crsptypesuffix])
  if haskey(estimationmindate, measuredatekey)
    mindate=estimationmindate[measuredatekey]
    @info "overriding mindate with mindate=$mindate"
    panel = panel[panel.date .≥ mindate,:]
    for FRLw ∈ Fweights(ms, :FRLw)
      panel[panel.date .== minimum(panel.date), FRLw] .= missing
    end
  end
  if haskey(estimationmaxdate, measuredatekey)
    panel = panel[panel.date .≤ estimationmaxdate,:]
  end

  #in contrast to other methods, drop any column where we don't have
  #the lag weight (improves performance I think)

  #conditionforturing!(panel::SubDataFrame, ms::MeasureSpec)

  zs::ZSpec = condition!(panel, ms)
  #=
  for p ∈ propertynames(zs)
    (p===:originalms) && continue
    prop = getproperty(zs, p)
    if length(propertynames(prop)) == 0
      println("$p: $prop")
    else
      for subp ∈ propertynames(prop)
        subprop = getproperty(prop, subp)
        println("$subp: $subprop")
      end
    end
  end

  throw("stop")=#

  #testing
  if PARAM[:testmeasure]
    testiter(panel, zs)
  end

  if PARAM[:refreshmeasure]
    λ⁺, results::DataFrame = itermodel(panel, zs, IterModel; funds)
    destandardize!(results, zs, intercept=false)
    results = computemeasurestatistics(panel, results, ms, returns)

    #this will be the name of the results file
    resultname::String = iterresultsname(λ⁺)

    results |> CSV.write("$(PARAM[:outputpath])\\$(resultname).csv")

  end

  return nothing
end

function computemeasurestatistics(panel, results, ms::MeasureSpec, returns;
    regressiontype::Symbol = PARAM[:iterregressiontype],
    volscale::Float64 = PARAM[:volscale])


  Fw::Vector{Symbol} = Fweights(ms, :Fw)
  FRtw::Vector{Symbol} = Fweights(ms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(ms, :FRLw)
  FW::Vector{Symbol} = ms.FW
  FWgroup = Fgroupedcontrols(ms)

  #construct the model parts using cpu only data types and no controls
  p = modelparts(panel, Float64, Vector{Float64}, Matrix{Float64},
    weightdata = (;Fw, FRtw, FRLw), Fvol=ms.Fvol,
    ;FW, FWgroup, ms)

  dims = p.dims
  ws = p.dat[:Fw]
  Rtws = p.dat[:FRtw]
  RLws = p.dat[:FRLw]
  W = p.dat[:FW]
  Wgroup = p.dat[:FWgroup]
  ts::Vector{Int} = p.ts
  v = p.dat[:Fvol]

  #construct simplest possible structures for computing predicted volumes
  @info "Computing raw controls for output file"
  Xvraw = XYIterControl(ws, Rtws, RLws, v, ts) #raw implies no controls or transformations
  @info "Computing actual controls for output file"
  Xv = XYIterControl(ws, Rtws, RLws, v, ts, W=W, Wgroup=Wgroup)
  Θ = VolumePartsIter(dims.T, dims.K, ts, Float64, Matrix{Float64}, Vector{Float64})

  #read in the estimations from the results file
  Fcoefs::Vector{Symbol} = findcols(results, :asset)
  Fzcoefs::Vector{Symbol} = findcols(results, :z)
  Fgcoefs::Vector{Symbol} = findcols(results, :growth)
  Fvcoefs::Vector{Symbol} = (s->Symbol(replace(string(s), r"^A_"=>"V_"))).(Fcoefs)
  Fvscoefs::Vector{Symbol} = (s->Symbol(replace(string(s), r"^A_"=>"Vs_"))).(Fcoefs)
  A = Matrix{Float64}(results[!, Fcoefs])
  Z = Matrix{Float64}(results[!, Fzcoefs])
  G = Matrix{Float64}(results[2:end, Fgcoefs])

  #now follow the sequence for estimating the volume, minus the actual estimation
  someones = ones(Float64, dims.K, 1)
  prodmatG = cumprodbyrow(G' |> Matrix{Float64})
  Ã = hcat(someones, prodmatG) |> Matrix{Float64}

  @assert Xvraw.xsection.xM === nothing
  RHSraw::Matrix{Float64} = projectbyt(Ã, Xvraw.xsection) #no controls or transofmrations
  RHS::Matrix{Float64} = hcat(projectbyt(Ã, Xv.xsection), Xv.W̃)

  #now read in the coefficients from the results file
  #we can also compare to the A in the results file as a check
  local β::Vector{Float64}=Vector{Float64}(undef, dims.K)
  if regressiontype === :none
    noregressionbase::Float64 = PARAM[:iternoregressionbase]
    β .= ones(Float64, dims.K) .* noregressionbase .* A[1,:] ./ Z[1,:] .* 2
  else
    β .= reg(RHS, Xv.ṽ, Val(regressiontype))[1:dims.K] .* 2
  end
  all(β .≈ A[1,:]) || throw(
    "Fcoefs: $Fcoefs
    A[1,:]: $(A[1,:])
    β: $β")
  @assert all(β' .* Ã' .≈ A)


  @info "Computing the volume for each strategy/stock/date"
  #now compute the volume for each strategy x stock x date
  volmat::Matrix{Float64} = RHSraw .* β'

  RHSpurchases = projectpurchasesbyt(Ã, Xvraw.xsection) #no controls or transofmrations
  purchasemat = RHSpurchases .* β'
  @assert abs.(purchasemat) ≈ volmat

  @assert size(volmat) === (length(ts), dims.K)
  @assert issorted(results, :date)
  alldates = (t->results.date[t]).(ts)
  volbystock::DataFrame = hcat(
    DataFrame(date=alldates),
    DataFrame(volmat, Fvcoefs),
    DataFrame(purchasemat, Fvscoefs),
    )

  #aggregate over all stocks
  strategyvol = combine(groupby(volbystock, :date),
    Fvcoefs .=> sum,
    Fvscoefs .=> sum, renamecols=false)

  aggvol = combine(groupby(panel, :date),
    :date => last => :date,
    :Wdvol => ((v)->sum(skipmissing(v))) => :desc_Wvolume,
    :dvolscaled => ((v)->sum(skipmissing(v))) => :desc_volume,
    :dvolscaled => ((v)->std(skipmissing(v))) => :desc_vol_std
    )

  stats = leftjoin(strategyvol, aggvol, on=:date)
  sort!(stats, :date)

  #compute the sum of all market caps in each quantile, differentiating long and short
  @info "Computing market cap statistics"
  Fmcgrosss = (s->Symbol(:mcgross_, s)).(Fw)
  Fmclongs = (s->Symbol(:mclong_, s)).(Fw)
  Fmcshorts = (s->Symbol(:mcshort_, s)).(Fw)
  Fmctotal = :mctotal
  aggmc = combine(groupby(panel,:date),
    :date => last => :date,
    :mc => ((x)->sum(x) |> (x)->x/volscale) => :mctotal,
    zip(Fw, [:mc for i ∈ 1:dims.K]) .|> collect .=>
      ((w, mc) -> (w .> zero(Float64)) .* mc |> sum |> (x)->x/volscale) .=>
      Fmclongs,
    zip(Fw, [:mc for i ∈ 1:dims.K]) .|> collect .=>
      ((w, mc) -> (w .< zero(Float64)) .* mc |> sum |> (x)->x/volscale) .=>
      Fmcshorts,)
  transform!(aggmc, zip(Fmclongs, Fmcshorts) .|> collect .=>
    ((l,s)-> l + s) .=>
    Fmcgrosss)

  stats = outerjoin(stats, aggmc, on=:date, validate=(true,true))
  sort!(stats, :date)

  #started working on showing the controls but no time
  #create the control names
  dtstrings = (d->Dates.format(d, "_yy_m_d")).(alldates|> unique |> sort!)
  FWgroupnames::Vector{Vector{Symbol}} = (dtstr->(Symbol.(FWgroup,dtstr))).(dtstrings)
  #FWgroupnames::Vector{Vector{Symbol}} = (F->(Symbol.(F,dtstrings))).(FWgroup)
  #FWgroupnames::Vector{Vector{Symbol}} = dtstrings .|> (dtstr->(F->(Symbol(F,dtstr)).(FWgroup)))

  #create the expanded full control matrix
  #NOTE: below commented code is useful for verificaiton, but it doesn't check for singular
  #columns so I don't use it anymore
  #=@info "Creating expanded control matrix"
  @time begin
    Wdf = ((xM,Fnames)->DataFrame(xM.W,Fnames)).(Xv.xsection.xM, FWgroupnames)
    Wgroupfull = foldl((df1,df2)->vcat(df1,df2,cols=:union), Wdf)
    for c ∈ eachcol(Wgroupfull)
      c .= c .|> (x)->coalesce(x, 0.0)
    end
    expanded::Matrix{Float64} = Wgroupfull |> Matrix{Float64}
  end

  @info "Verifying expanded control matrix"
  @time @assert (expanded .≈ expandcontrolsbyt(Wgroup, ts)) |> all=#

  @info "Checking groupped regression results match full regression results"
  expanded::Matrix{Float64}, controlsused = expandcontrolsbyt(Wgroup, ts) |> (nt
    )->(nt.expanded, nt.Wgroupcolsused)

  #run the full regression to get the controls
  #while we are at it, we can check that the regression ran as expected
  if regressiontype === :none
    RHSfull = hcat(expanded, Xv.W)
    βfull = svd(Symmetric(RHSfull' * RHSfull))\(RHSfull'*v) |> vec
    βgroup = βfull[1:(end-size(Xv.W,2))]
    fullerr = abs.(sum(RHS[:,1:dims.K],dims=2) .+ RHSfull*βfull .- Xv.v)
    βW̃ = (size(Xv.W̃,2) > 0) ? reg(Xv.W̃, Xv.ṽ, Val{:svd}()) : TV()
    βproj = [ones(dims.K); βW̃]
  else
    RHSfull = hcat(RHSraw,expanded, Xv.W)
    βfull =  cholesky!(Symmetric(RHSfull' * RHSfull))\(RHSfull'*v) |> vec
    βgroup = βfull[(size(RHSraw,2) + 1):(end-size(Xv.W,2))]
    fullerr = abs.(RHSfull*βfull .- Xv.v)
    βproj = reg(RHS, Xv.ṽ, Val(regressiontype))
  end


  #the below check verifies that the residuals are projected as expected
  projerr = abs.(RHS*βproj .- Xv.ṽ)
  parterr = abs.(sum(RHSraw,dims=2))
  @assert (fullerr ≈ projerr) "
    fullerr[1:5]: $(fullerr[1:5])
    projerr[1:5]: $(projerr[1:5])
    parterr[1:5]: $(parterr[1:5])
    Median err: Full: $(median(fullerr)) proj: $(median(projerr))
    Mean err: Full: $(mean(fullerr)) proj: $(mean(projerr))
    Median dev: Full: $(median(fullerr .- projerr))
    mean dev: Full: $(mean(fullerr .- projerr))
    Max abs dif: $(maximum(fullerr .- projerr))
    SSE err: Full: $(sum(fullerr.^2)) Proj: $(sum(projerr.^2))"

  @info "verificaiton of projection matrix passed"
  #reshape the control coefficients into time series with as much care as possible
  #@assert length(βgroup) == length(unique(ts)) * size(Wgroup,2)
  βgroupconformed = Vector{MFloat64}(undef, length(unique(ts)) * size(Wgroup,2))
  βgroupconformed[controlsused] .= βgroup
  βgroupmat = reshape(βgroupconformed, size(Wgroup,2), dims.T-1)'
  @assert (vec(βgroupmat') .≡ βgroupconformed) |> all

  rhsdf::DataFrame = DataFrame(date=alldates|>deepcopy, resid=projerr)
  controldf = combine(groupby(rhsdf,:date),
      :date => last => :date,
      :resid=> ((v)->std(skipmissing(v)))=>:se_resid)
  @assert issorted(controldf, :date)
  for (c,F) ∈ enumerate(FWgroup)
    #controldf[!, F] = [missing; βgroupmat[:,c]]
    controldf[!, F] = βgroupmat[:,c]
  end

  @info "combining controls, extra stats, and returns with the results dataframe"

  stats = outerjoin(stats, controldf, on=:date, validate=(true,true))
  sort!(stats, :date)
  results = leftjoin(results, stats, on=:date, validate=(true,true))

  results = leftjoin(results, returns, on=:date, validate=(true,true))
  sort!(results, :date)


  #patch in data from the return df
  #@assert all(completecases(results[2:end, [:volume; :Wvolume]], [:volume; :Wvolume]))
  if !all(completecases(results[2:end, [:volume; :Wvolume]], [:volume; :Wvolume]))
    results |> CSV.write("$(PARAM[:testpath])\\results_error.csv")
    @assert all(completecases(results[2:end, [:volume; :Wvolume]], [:volume; :Wvolume]))
  end

  resultsdescvolumes = view(results, completecases(results, [:desc_volume, :desc_Wvolume]), :)
  @assert all(resultsdescvolumes[:, :volume] .≈ resultsdescvolumes[:, :desc_volume] )
  @assert all(resultsdescvolumes[:, :Wvolume] .≈ resultsdescvolumes[:, :desc_Wvolume])

  #the data from the returns file are more compelte
  select!(results, Not([:volume, :Wvolume]))

  for ξ ∈ ms.ξs
    @assert all(completecases(results[2:end, [ξ.Fξ]]))
    rename!(results, Dict(ξ.Fξ=>Symbol(:R_, ξ.Fξ)))
  end

  #####now merge in the returns data

  results = ivoltests(panel, results, ms, returns)

  return  results
end
