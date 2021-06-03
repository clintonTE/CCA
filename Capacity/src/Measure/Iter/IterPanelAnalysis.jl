
function paformat(β;scalefactor, decimals, prefix, suffix)
  valstring = num2str(β, decimals; scalefactor)
  return "$(prefix)$(valstring)$(suffix)"
end

pabetaformat =(; scalefactor=1.0, decimals=2, prefix="", suffix="\\:\\:")
pabetax2format =(; scalefactor=2.0, decimals=2, prefix="", suffix="\\:\\:")
papercentageformat =(; scalefactor=100.0, decimals=2, prefix="", suffix="\\:\\:")
paseformat = (;scalefactor=1.0, decimals=2, prefix="\\mathit{(", suffix=")}")
pasex2format = (;scalefactor=1.0, decimals=2, prefix="\\mathit{(", suffix=")}")
paNformat = (;scalefactor=1.0, decimals=0, prefix="", suffix="")
par2format = (;scalefactor=100.0, decimals=2, prefix="", suffix="")

#put all of the field information here
pafields = Dict{Symbol, Any}(
  :rhs=>Dict(
    :A_WXmomentum121s=>(;
      id=:A_WXmomentum121s,
      label="\$A_0\$ (estimated) L+S",
      format=pabetax2format,
      se=(;
        id=:A_WXmomentum121s_se,
        label="",
        format=pasex2format)),
    :A_WXmomentum121s0=>(;
      id=:A_WXmomentum121s0,
      label="\$A_0\$ (\$G=1 \\: \\forall t\$) L+S",
      format=pabetax2format,
      se=(;
        id=:A_WXmomentum121s0_se,
        label="",
        format=pasex2format)),
    :intercept=>(;
      id=:intercept,
      label="intercept",
      format=pabetaformat,
      se=(;
        id=:intercept_se,
        label="",
        format=paseformat)),
    ),
    :description => Dict(
      :permno=>(;
        id=:permno,
        label="firm FE"),
      :date=>(;
        id=:date,
        label="time FE"),
      :N=>(;
        id=:N,
        label="N",
        format=paNformat),
      :R2 => (
        ;
        id=:R2,
        label = "\\midrule \$R^2\$",
        format = par2format),
      :R2delta => (
        ;
        id=:R2delta,
        label = "\$R^2 - R^2(mean)\$",
        format = par2format),
    )
  )

  #start of a table spec
function patablefull()
  Fmeasures=[
    [:intercept],
    #[],
    [:A_WXmomentum121s, :intercept],
    #[:A_WXmomentum121s],
    #[:A_WXmomentum121s],
    #[:A_WXmomentum121s],
    [:A_WXmomentum121s0, :intercept],
    #[:A_WXmomentum121s0],
    #[:A_WXmomentum121s0],
    #[:A_WXmomentum121s0],
    ]
  tablename="$(PARAM[:crspdatalabel])_pafull"
  Ffes=[
    Symbol[],
    #[:date],
    Symbol[],
    #[:date],
    #[:permno],
    #[:permno, :date],
    Symbol[],
    #[:date],
    #[:permno],
    #[:permno, :date]
    ]

  ses = [()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    ()->Vcov.cluster(:date,:permno),
    ()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    #()->Vcov.cluster(:date,:permno),
    ]
  Fvols = [:dvolscaled for i ∈ 1:length(Fmeasures)]
  Fse=:se

  Fdescriptions = [:R2, :R2delta, :N, #=:date, :permno=#]
  Frowheaders = [:A_WXmomentum121s; :A_WXmomentum121s0; :intercept; Fdescriptions;]
  alignmentstring="l | rrr" #"l | r | rrr | rrr"

  #verify integrity here for the sake of debugging
  allmeasures = unique([Fmeasures...;])
  @assert length(Fmeasures) == length(Ffes)
  @assert length(Fvols) == length(Ffes)
  @assert setdiff(allmeasures, Frowheaders) |> isempty

  #build the row list
  Frows = Symbol[]
  for Frowheader ∈ Frowheaders
    push!(Frows, Frowheader)
    @assert haskey(pafields[:rhs], Frowheader) || haskey(pafields[:description], Frowheader
      ) "field not found in pafield: Frowheader=$Frowheader"

    if haskey(pafields[:rhs], Frowheader)
      push!(Frows, Symbol(Frowheader, :_, Fse))
      @assert pafields[:rhs][Frowheader].se.id ≡ Symbol(Frowheader, :_, Fse) "
        could not find $(Symbol(Frowheader, :_, Fse)) ∈ Frowheader $(pafields[:rhs][Frowheader])"
    end
  end


  return (; Fmeasures, Frows, tablename, Ffes, ses, Fvols, Fdescriptions, Fse, alignmentstring)
end

#run it with the Gs and without, cross-sectionally as well
#include standard errors- maybe fixed effects?
function ivoltests(panel, results, ms::MeasureSpec, returns)

  #get the field names
  Fξs = ms.Fξs
  Fw::Vector{Symbol} = Fweights(ms, :Fw)
  FRtw::Vector{Symbol} = Fweights(ms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(ms, :FRLw)
  FW::Vector{Symbol} = ms.FW
  FWgroup = Fgroupedcontrols(ms)

  Fcoefs::Vector{Symbol} = findcols(results, :asset)
  Fgcoefs::Vector{Symbol} = findcols(results, :growth)
  Fcoefs0 = Fcoefs .|> (Fcoef->Symbol(Fcoef, 0))


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
  permno = p.dat[:permno]
  dates = p.dat[:date]
  cleanpanel = p.cleanpanel

  @assert (dates .== (t->results.date[t]).(ts)) |> all

  #construct simplest possible structures for computing predicted volumes
  @info "Computing raw controls for output file"
  Xvraw = XYIterControl(ws, Rtws, RLws, v, ts) #raw implies no controls or transformations
  Θ = VolumePartsIter(dims.T, dims.K, ts, Float64, Matrix{Float64}, Vector{Float64})

  A = Matrix{Float64}(results[!, Fcoefs])
  G = Matrix{Float64}(results[2:end, Fgcoefs])
  #now follow the sequence for estimating the volume, minus the actual estimation
  someones = ones(Float64, dims.K, 1)
  prodmatG = cumprodbyrow(G' |> Matrix{Float64})
  Ã = hcat(someones, prodmatG) |> Matrix{Float64}

  RHSraw::Matrix{Float64} = projectbyt(Ã, Xvraw.xsection) #no controls or transofmrations
  RHSraw0::Matrix{Float64} = projectbyt(ones(Float64, size(Ã)), Xvraw.xsection)

  #form the sample for regressions
  expanded::Matrix{Float64}, controlsused = expandcontrolsbyt(Wgroup, ts) |> (nt
    )->(nt.expanded, nt.Wgroupcolsused)
  focal = hcat(
    DataFrame(date=dates, permno=permno),
    DataFrame(RHSraw, Fcoefs),
    DataFrame(RHSraw0, Fcoefs0))
  focal[!, PARAM[:Fvol]] = v

  #construct intercepts by hand
  focal.intercept = ones(length(v))

  ############################
  patable(focal; patablefull()...)
  pax = paxsectional(focal; Fcoefs,Fcoefs0)
  results = leftjoin(results, pax, on=:date, validate=(true,true))
  results = sort!(results, :date)

  return results
end

function patable(focal;
    tablename,
    Fmeasures,
    Fvols::TFvols,
    Fse,
    ses,
    Ffes,
    Frows,
    Fdescriptions,
    verbose = PARAM[:tabpaverbose],
    tablepath = PARAM[:tablepath],
    colnames=[1:length(Fmeasures) .|> i->"($i)",],
    widthcolnames= [ones(Int,length(Fmeasures)),],
    alignmentcolnames = [
      (rc-> "c").(colnamerow) for (r,colnamerow) ∈ enumerate(colnames)],
    alignmentstring = string("l",["r" for i ∈ 1:length(Fmeasures)]...),

    #construct the info for each row
    #first extract lightly constructed field objects IN ORDER
    datafields = reduce(vcat, (reduce(vcat, Fmeasures) |> unique) .|> (
      (s)->[pafields[:rhs][s], pafields[:rhs][s][Fse]])),
    descriptionfields = Fdescriptions .|> ((s)->pafields[:description][s])

    ) where {TFvols}

  #create an index of row information. order is not preserved.
  Frowindex = Dict(f.id=>f for f ∈ [datafields; descriptionfields])

  #double check alignment
  @assert (id->Frowindex[id].id≡id).(Frows) |> all

  #extract the row labels, using the order preserved in Frowid
  descrownames = (id->Frowindex[id].label).(Frows)

  #form the regression formulai
  formulai = FormulaTerm[]
  fmlmrhs = Any[]
  for (Fmeasure, Ffe, Fvol) ∈ zip(Fmeasures, Ffes, Fvols)

    termrhs = term.([Fmeasure; fe.(Ffe)]) |> sum
    push!(fmlmrhs, Meta.parse(join((string).([Fmeasure; Ffe]),"+") * "+ 0"))

    if :intercept ∈ Fmeasure
      @assert Ffe |> isempty "Intercept detected in fixed effects spec!!"
      termrhs += term(0) # do not want an itnercept if we are doing the regressions by hand
    else
      @assert Ffe |> !isempty "No intercept detected without FE!"
    end

    push!(formulai, term(Fvol) ~ termrhs)
  end

  reg(f,err) = FixedEffectModels.reg(focal, f, err())

  #run the full regressions
  femlms = reg.(formulai, ses)

  verbose && println.(femlms)


  allFfes = reduce(vcat, Ffes) |> unique
  regressioncols = DataFrame([id=>String[] for id ∈ Frows])

  #make sure that the first column si the baseline
  @assert (length(Fmeasures[1]) == 1) && (Fmeasures[1][1] ≡ :intercept)
  r2base = r2(femlms[1])
  Ndates = length(unique(focal.date))
  for (c, femlm) ∈ enumerate(femlms)
    #get the main coefficients and SEs
    coefindex = Dict(Symbol(n) => A for (n,A) ∈ zip(coefnames(femlm), femlm.coef))
    tindex = Dict(Symbol(n, :_, Fse)=>A/σ for (n, A, σ) ∈ zip(
      coefnames(femlm), femlm.coef, vcov(femlm) |> (Σ)->diag(Σ).^0.5))
    regressioncolvals = merge(coefindex, tindex)

    #check the coefficients
    Ncoef = Fmeasures[c] |> length

    #= TOO SLOWWWWWWW
    fmlm = FMLM(focal, fmlmrhs[c], Fvols[c],
      #=withinsym=length(Ffes[c]) > 0 ? Ffes[c][1] : nothing,=# clustersyms=[:permno, :date])
    fmlmdf = DataFrame(coef=Fmeasures[c])
    fmlmdf.beta = fmlm.β[1:Ncoef]
    fmlmdf.femlmbeta = (F -> coefindex[F]).(Fmeasures[c])
    fmlmdf.femlmse = (F -> coefindex[F]/tindex[Symbol(F, :_, Fse)]).(Fmeasures[c])
    fmlmdf.robust = diag(modifiedwhiteΣ!(fmlm))[1:Ncoef].^0.5
    fmlmdf.clustered = diag(clusteredΣ!(fmlm))[1:Ncoef].^0.5
    fmlmdf.hac2way = diag(neweywestΣ!(fmlm, nrow(focal)^0.33333 |> ceil |> Int))[1:Ncoef].^0.5

    fmlm1way = FMLM(focal, fmlmrhs[c], Fvols[c],
      #=withinsym=length(Ffes[c]) > 0 ? Ffes[c][1] : nothing,=# clustersyms=[:permno])
    fmlmdf.panelhac1way = diag(neweywestpanelΣ!(fmlm1way, Ndates^0.33333 |> ceil |> Int))[1:Ncoef].^0.5

    @assert fmlmdf.beta ≈ fmlmdf.femlmbeta
    @info "$(Fvols[c]) ~ $(fmlmrhs[c])"
    println(fmlmdf)=#


    #format the values
    regressioncol = Dict(k =>paformat(regressioncolvals[k]; Frowindex[k].format...)
      for k ∈ keys(regressioncolvals))

    Fempties = setdiff(Frows,
      [Fdescriptions; Fmeasures[c]; (n->Symbol(n, :_, Fse)).(Fmeasures[c])])
    regressioncol = merge(regressioncol, Dict(F=>"" for F ∈ Fempties))

    #now do fixed effects
    regressioncol = merge(regressioncol,Dict(F=>F ∈ Ffes[c] ? "X" : "" for F ∈ allFfes))

    regressioncol[:R2] = paformat(r2(femlm); Frowindex[:R2].format...)
    regressioncol[:R2delta] = paformat(r2(femlm) - r2base; Frowindex[:R2delta].format...)
    regressioncol[:N] = paformat(
      sum(completecases(focal, [Fmeasures[c]; Ffes[c];])); Frowindex[:N].format...)

    #@info "propertynames: $(propertynames(regressioncols))"
    #@info "row keys: $(keys(regressioncol))"
    push!(regressioncols, regressioncol)
  end

  @assert (propertynames(regressioncols) .≡ Frows) |> all
  desccontent = eachcol(regressioncols) |> collect .|> Vector


  t = textable(;colnames, descrownames, desccontent,
    widthcolnames, alignmentstring, alignmentcolnames)

  writetextable(t, path=tablepath, outname="$tablename.tex")
end

function paxsectional(focal;
    Fcoefs,Fcoefs0, Fvol=PARAM[:Fvol],
    )


  Fcoefsx = Fcoefs .|> (Fcoef->Symbol(:x, Fcoef))
  Fcoefsxse = Fcoefs .|> (Fcoef->Symbol(:xse, Fcoef))
  Fcoefsx0 = Fcoefs0 .|> (Fcoef->Symbol(:x, Fcoef))
  Fcoefsx0se = Fcoefs0 .|> (Fcoef->Symbol(:xse, Fcoef))
  Fcoefsxd = Fcoefs .|> (Fcoef->Symbol(:xt, Fcoef))

  #this will hold the results
  xs = DataFrame(date=Date[])
  for F ∈ [Fcoefsx; Fcoefsxse; Fcoefsx0; Fcoefsxse]
    xs[!, F] = MFloat64[]
  end

  #some housekeeping
  intercepts = [:intercept for i ∈ 1:length(Fcoefs)]
  focal.intercept0 = deepcopy(focal.intercept)
  intercepts0 = [:intercept0 for i ∈ 1:length(Fcoefs)]
  Fvols = [Fvol for i ∈ 1:length(Fcoefs)]


  #run xsectional regressions
  function runxsectional(t)
    f = term(Fvol) ~ sum(term.(Fcoefs)) + term(:intercept) + term(0)
    femlm = FixedEffectModels.reg(t, f, Vcov.robust())

    #run the regressions
    coefindex= Dict(Symbol(:x, n, ) => 2*A for (n,A) ∈ zip(coefnames(femlm), femlm.coef))
    seindex = Dict(Symbol(:xse, n, )=> 2*σ for (n, σ) ∈ zip(
      coefnames(femlm), vcov(femlm) |> (Σ)->diag(Σ).^0.5))

    #run the 0 regressions (note the 0 is already added)
    f0 = term(Fvol) ~ sum(term.(Fcoefs0)) + term(:intercept0) + term(0)
    femlm0 = FixedEffectModels.reg(t, f0, Vcov.robust())
    coefindex0= Dict(Symbol(:x, n, ) => 2*A for (n,A) ∈ zip(coefnames(femlm0), femlm0.coef))
    seindex0 = Dict(Symbol(:xse, n, )=> 2*σ for (n, σ) ∈ zip(
      coefnames(femlm0), vcov(femlm0) |> (Σ)->diag(Σ).^0.5))

    return Ref(merge(coefindex, coefindex0, seindex, seindex0))
  end


  xfocal = groupby(focal, :date)
  xresults = combine(xfocal, zip(Fcoefs, Fcoefs0, Fvols, intercepts, intercepts0) .|> collect .|> AsTable
    .=> runxsectional .=> Fcoefsxd)

  #splat all the gathered regression results
  for i ∈ 1:length(Fcoefs)
    transform!(xresults,
      Fcoefsxd => ((d)-> (di->di[Fcoefsx[i]]).(d)) => Fcoefsx[i],
      Fcoefsxd => ((d)-> ((di->di[Fcoefsxse[i]]).(d))) => Fcoefsxse[i],
      Fcoefsxd => ((d)-> ((di->di[Fcoefsx0[i]]).(d))) => Fcoefsx0[i],
      Fcoefsxd => ((d)-> ((di->di[Fcoefsx0se[i]]).(d))) => Fcoefsx0se[i])
  end
  transform!(xresults,
     Fcoefsxd => ((d)->(di->di[:xintercept]).(d))=>:xintercept,
     Fcoefsxd => ((d)-> (di->di[:xseintercept]).(d))=>:xseintercept,
     Fcoefsxd => ((d)-> (di->di[:xintercept0]).(d))=>:xintercept0,
     Fcoefsxd => ((d)-> (di->di[:xseintercept0]).(d))=>:xseintercept0)

  select!(xresults, Not(Fcoefsxd))
  return xresults
end
