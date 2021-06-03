
function generateweeklycomomentum(;
    minprice=PARAM[:plminprice],
    minmarketcap=PARAM[:plminmarketcap],
    minnysequantile=PARAM[:plminnysequantile],
    Fplmomentum::Symbol=PARAM[:plmomentum],
    momentum = PARAM[:momentumstrategies][Fplmomentum],
    lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
    Fplret::Symbol=PARAM[:plret])

  #this check helps sync up the measure output frequency (weekly/monthly) with that of
  #the comomentum
  PARAM[:crspfrequencysuffix] == "w" || throw("Weekly comomentum should only be
    created when crspfrequencysuffix=\"w\"")

  #housekeeping

  selectsummarydata_pl(df; kwargs...) = selectsummarydata(df; minprice, minmarketcap, minnysequantile, kwargs...)


  crspw::DataFrame = makeweekly(selectsummarydata_pl,
    binprefix="wpl_",
    mainpanel=false,
    computeivol=false)
  validatedates(crspw.date |> unique |> sort!, frequency=:week)
  @assert issorted(crspw, [:permno, :date])
  conditionmomentumstrategies!(crspw)

  crspw = crspw[completecases(crspw, Fplret), :]

  #first do the weekly data
  comomw = generatecomomentum(deepcopy(crspw))
  simpledif(v) = [missing; v[2:end] .- v[1:(end-1)]]

  #comomw = lpfilter(comomw;
  #  focalcols = [:comom10, :comom1, :comom, :lGcomom10, :lGcomom1, :lGcomom], lpfilterperiods)
  comomw.Gcomom =[missing; comomw.comom[2:end] .- comomw.comom[1:(end-1)];]
  return comomw
end


function generatemonthlycomomentum(;
    minprice=PARAM[:plminprice],
    minmarketcap=PARAM[:plminmarketcap],
    minnysequantile=PARAM[:plminnysequantile],
    Fplmomentum::Symbol=PARAM[:plmomentum],
    momentum = PARAM[:momentumstrategies][Fplmomentum],
    lpfilterperiods::Vector{<:Union{Day, Nothing}} = PARAM[:analysislpfilterperiods],
    patchmonthly::Bool = PARAM[:plpatchmonthly],
    Fplret::Symbol=PARAM[:plret])

  PARAM[:crspfrequencysuffix] == "m" || throw("Monthly comomentum should only be
    created when crspfrequencysuffix=\"m\"")


  #housekeeping
  selectsummarydata_pl(df; kwargs...) = selectsummarydata(
    df; minprice, minmarketcap, minnysequantile, kwargs...)

  crspw::DataFrame = makeweekly(selectsummarydata_pl,
    binprefix="wpl_",
    mainpanel=false,
    computeivol=false)
  validatedates(crspw.date |> unique |> sort!, frequency=:week)
  @assert issorted(crspw, [:permno, :date])

  crspw = crspw[completecases(crspw, Fplret), :]

  if patchmonthly #this implies using the monthly momentum
    (PARAM[:refreshweekly] ≠ PARAM[:refreshmonthly]) && @warn(
      "PARAM[:refreshweekly]=$(PARAM[:refreshweekly]) while
      PARAM[:refreshmonthly]=$(PARAM[:refreshmonthly])!!! This could imply
      out of sync crsp files as used to compute comomentum.")
    crspm::DataFrame = makemonthly(selectsummarydata_pl,
      binprefix="mpl_",
      mainpanel=false,
      computeivol=false)
    validatedates(crspm.date |> unique |> sort!, frequency=:month)
    @assert issorted(crspm, [:permno, :date])


    crspm = crspm[completecases(crspm, Fplret), :]

    #for the monthly data, we first need to splice in the monthly deciles into the weekly data
    crspw.yearmonth = year.(crspw.date) .+ month.(crspw.date) ./ 100
    crspm.yearmonth = year.(crspm.date) .+ month.(crspm.date) ./ 100
    Nrowold = nrow(crspw)
    select!(crspw, Not(Fplmomentum))
    conditionmomentumstrategies!(crspm)
    crspw = leftjoin(crspw, crspm[!, [:permno, :yearmonth, Fplmomentum]], on=[:permno, :yearmonth],
      validate=(false, true))
    if !issorted(crspw, [:permno, :yearmonth])
      sort!(crspw, [:permno, :yearmonth])
    end
    @assert Nrowold == nrow(crspw)
  else
    conditionmomentumstrategies!(crspw)
  end


  comomm = generatecomomentum(crspw)

  #drop the extraneous rows
  @assert issorted(comomm, :date)
  @assert allunique(comomm.date)
  comomm.yearmonth = year.(comomm.date) .+ month.(comomm.date) ./ 100
  comomm = combine(last, groupby(comomm, :yearmonth), renamecols=false)
  select!(comomm, Not(:yearmonth))

  #run through an LP filter
  #comomm = lpfilter(comomm;
  #  focalcols = [:comom10, :comom1, :comom, :lGcomom10, :lGcomom1, :lGcomom], lpfilterperiods)
  comomm.Gcomom =[missing; comomm.comom[2:end] .- comomm.comom[1:(end-1)];]

  return comomm
end

function generatecomomentum(crsp::DataFrame;
  testdate::Date = PARAM[:testdate],
  minwindownonmissing = PARAM[:plminwindownonmissing],
  Fplret::Symbol=PARAM[:plret],
  Fplmomentum::Symbol=PARAM[:plmomentum]
  )

  crspprefix = Dict("m"=>PARAM[:crspmonthly], "w"=>PARAM[:crspweekly])[PARAM[:crspfrequencysuffix]]

  #merge the monthly momentum deciles with the weekly data

  #merge in ff
  crsp = mergeff(crsp, datefrequency=:week, Fgroup=:permno)
  if !issorted(crsp, [:permno, :date])
    sort!(crsp, [:permno, :date])
  end

  #need a symbol id for each column
  crsp.spermno = crsp.permno .|> (i)->Symbol(:p, i)

  #these are conditioned on when computing the correlations
  crsp.some1s = ones(nrow(crsp))
  Fcontrols = [:F_ff3m_mktrf, :F_ff3m_smb, :F_ff3m_hml]
  Fwindowcols = [:date; :spermno; Fplret; Fcontrols]
  crsp = crsp[completecases(crsp, Fwindowcols), :]

  #compute number missing in past year
  #this will be useful when we compute the cross-sectional correlations
  crsp.plretnmissing1y = missings(Float64, nrow(crsp))
  trailingmeasure!(length, crsp, Year(1), 1,
    Ftarget=Fplret,
    Fmeasure=:plretnmissing1y,)

  function deciles(v::AbstractVector{Float64})::Vector{Int} where T
    q = ecdf(v)
    #the extra min is due to the top value only
    d = (q.(v) .* 10 .+ 1) .|> x->min(floor(Int, x),10)
    @assert (1,10) ≡ extrema(d) "extrema(d): $(extrema(d))"
    return d
  end


  #create a version for iterating on deciles
  crspdeciles = crsp[completecases(crsp, [Fplmomentum]), [:date, :spermno, Fplmomentum, :plretnmissing1y, :mc]]
  crspdeciles.plmomentum = Vector{Float64}(crspdeciles[!, Fplmomentum])
  transform!(groupby(crspdeciles, :date),
    :plmomentum => deciles => :d_plmomentum,
    :plretnmissing1y => maximum => :maxnonmissing)

  #this is useful to examine the completeness
  completeness = combine(groupby(crspdeciles, :date),
    nrow => :totalpermno,
    [:plretnmissing1y,:maxnonmissing] =>
      ((plretnmissing1y,maxnonmissing) ->
        sum(plretnmissing1y .≥ minwindownonmissing * maxnonmissing)) =>
      :completepermno,
    :maxnonmissing => last => :maxnonmissing,
    [:plretnmissing1y,:maxnonmissing, :mc] =>
      ((plretnmissing1y,maxnonmissing, mc) ->
        sum((plretnmissing1y .≥ minwindownonmissing * maxnonmissing) .* mc)) =>
      :completemc,
    :mc => sum => :totalmc)
  completeness.perkept = completeness.completepermno ./ completeness.totalpermno
  completeness.permckept = completeness.completemc ./ completeness.totalmc
  #println("complete by date\n", completeness)

  @info "Average comomentum $crspprefix completeness: $(mean(completeness.perkept))
    Value-weighted: $(mean(completeness.permckept))"

  #drop unused deciles
  crspdeciles = crspdeciles[(crspdeciles.d_plmomentum .== 1) .| (crspdeciles.d_plmomentum .== 10), :]
  crspdeciles = crspdeciles[
    crspdeciles.plretnmissing1y .≥ minwindownonmissing .* crspdeciles.maxnonmissing, :]
  select!(crspdeciles, Not([:maxnonmissing, :mc]))

  #now compute the cross-correlations
  crspdeciles.xcor = missings(Float64, nrow(crspdeciles))

  #force controls to exclude missing values
  controls = sort!(unique(crsp[!, [:date; Fcontrols;]]), :date)
  for Fcontrol ∈ Fcontrols
    controls[!, Fcontrol] = Vector{Float64}(controls[!, Fcontrol])
  end

  #calls the partial correlation function after accounting for missing
  function mpartialcor(x::AbstractVector{<:MFloat64},
      y::AbstractVector{<:MFloat64},
      Z::AbstractMatrix{Float64})
    pcorinds = (x .!== missing) .& (y .!== missing)

    return partialcor(x[pcorinds], y[pcorinds], Z[pcorinds, :])
  end

  #Approach below is as follows:
  #1) For each date, select a lagging 1 year window
  #2) For each of the relevant deciles, identify the permnos with complete enough series over the window
  #3) Unstack these into a TxN matrix. Compute the EW portoflio returns after accounting for missings
  #4) Compute the partial correlation between each column and the EW portfolio returns

  gcrspdeciles = groupby(crspdeciles, [:date])
  @assert allunique(controls.date)
  Threads.@threads for i ∈ 1:length(gcrspdeciles)
    scrspdeciles = gcrspdeciles[i]
    dt = scrspdeciles[1,:date]
    dtM1y=dt - Year(1)

    #lagging one year window of returns
    crspwindow = view(crsp, (dtM1y .< crsp.date) .& (crsp.date .≤ dt),Fwindowcols)

    #compute the control matrix- only need to do this once per decile
    controlswindow = view(controls, (dtM1y .< controls.date) .& (controls.date .≤ dt),Fcontrols)


    controlmat = Matrix{Float64}(controlswindow) #control matrix for correlations
    for dec ∈ groupby(scrspdeciles, :d_plmomentum)
      #select only complete permno from the particular decile
      decwindow = innerjoin(crspwindow, dec[!, [:spermno]], on=:spermno)
      if !issorted(decwindow, [:spermno])
        sort!(decwindow, [:spermno])
      end
      @assert dec[1, :plretnmissing1y] .≥ length(unique(crspwindow.date)) .* minwindownonmissing

      #create the TxN matrix and compute the correlations
      crspbypermno = unstack(decwindow, :date, :spermno, Fplret)

      if (dt == testdate) && (dec[1, :d_plmomentum] == 10)
        crspwindow |> CSV.write(
          "$(PARAM[:testpath])\\$(crspprefix)_comom-crspwindowd-10.csv")
        controlswindow |> CSV.write(
          "$(PARAM[:testpath])\\$(crspprefix)_comom-ctrlwindowd-10.csv")
        crspbypermno  |> CSV.write(
          "$(PARAM[:testpath])\\$(crspprefix)_comom-crsppermnod-10.csv")
      end

      decilemat = Matrix{MFloat64}(crspbypermno[!, Not(:date)])

      #compute the total returns while taking care with missing values
      deciletotret = sum((x->coalesce(x,0.0)).(decilemat), dims=2) |> vec
      Ndecile = sum(decilemat .!== missing, dims=2) |> vec

      #NOTE- far better way to do the below, probably...
      colnames = setdiff(propertynames(crspbypermno), [:date])
      for (i, ri) ∈ enumerate(eachcol(decilemat))
        try
          dec[i, :xcor] = mpartialcor(ri, (deciletotret .- ri) ./ (Ndecile .- 1), controlmat)
        catch err
          dec[i, :xcor] = missing
          @warn "Failed to compute partial cor $err
            i=$i, colname: $(colnames[i]), dt: $dt"
        end
      end
    end
  end

  #print out some data for testing purposes
  testview = view(crspdeciles,
    ((s->parse(Int,string(s)[2:end])).(crspdeciles[!, :spermno])
      .% PARAM[:testpermnomult]) .== 0, :)
  testview |> CSV.write(
    "$(PARAM[:testpath])\\$(crspprefix)_comomdeciles.csv")

  #average the cross-sectional correlations, first within decile, then across top/bottom deciles
  crspdeciles = crspdeciles[crspdeciles.xcor .!== missing, :]
  crspdeciles = crspdeciles[isfinite.(crspdeciles.xcor), :]
  meanif(v, cond) = sum(v .* cond) ./ sum(cond)
  crspdeciles.d1 = crspdeciles.d_plmomentum .== 1
  crspdeciles.d10 = crspdeciles.d_plmomentum .== 10
  comom = combine(groupby(crspdeciles, :date),
    [:xcor, :d10] => meanif => :comom10,
    [:xcor, :d1] => meanif => :comom1)
  comom.comom = (comom.comom10 .+ comom.comom1) ./ 2
  sort!(comom, :date)

  #create log versions
  #comom.lGcomom = [missing; log.(comom.comom[2:end] ./ comom.comom[1:(end-1)])]
  #comom.lGcomom10 = [missing; log.(comom.comom10[2:end] ./ comom.comom10[1:(end-1)])]
  #comom.lGcomom1 = [missing; log.(comom.comom1[2:end] ./ comom.comom1[1:(end-1)])]

  return comom
end


function adjustindustry(crsp;
  industryfilename = PARAM[:plindustryfilename],
  industryfilepath = PARAM[:plindustryfilepath],
  sicfilename = PARAM[:plsicfilename],
  sicfilepath = PARAM[:plsicfilepath],
  missingunmatched=PARAM[:plindustrymissingunmatched],
  plretmin::Float64=PARAM[:plretmin],
  Fmc=:mc)

  lowercasesymbol(s) = s |> string |> lowercase |> Symbol

  #stack each sic code and with an ff code
  sic::DataFrame = CSV.File("$sicfilepath\\$(sicfilename).csv") |> DataFrame
  sort!(sic, [:ffcode, :beginsic])
  sic.ffdesc = sic.ffdesc .|> cleanpropertyname

  #create a lookup table of ffcd vs ffdesc
  fflookup = unique(sic[!, [:ffcode, :ffdesc]])
  @assert all(fflookup.ffcode .== collect(1:30))
  ffdesc = fflookup.ffdesc

  sic = reduce(append!!,
    (r->DataFrame(ffcode=r.ffcode, ffdesc=r.ffdesc,
      siccd=r.beginsic:r.endsic)).(eachrow(sic)))

  #validate
  @assert sum(nonunique(sic, :siccd)) == 0
  crsp.siccd = crsp.siccd .|> (i)->parsecrsp(Int, i)
  unknownsic = setdiff(unique(crsp.siccd), sic.siccd)

  #join the siccds
  crsp = leftjoin(crsp, sic, on=:siccd, matchmissing=:equal, validate=(false, true))
  if !issorted(crsp, [:permno, :date])
    sort!(crsp, [:permno, :date])
  end
  @assert sum(nonunique(crsp, [:permno, :date])) == 0

  #check what we are missing
  avgmissingsicmc = mean(combine(groupby(crsp, :date), [:ffcode, Fmc] =>
    ((ffcode, mc)->sum((ffcode .=== missing) .* mc)) => :missingsicmc).missingsicmc)
  avgtotalmc =  mean(combine(groupby(crsp, :date), [Fmc] => sum =>  :totmc).totmc)
  #@info "unmatched siccds: $unknownsic\n average mc missing: $avgmissingsicmc vs $avgtotalmc"

  #need to handle the mismatches- probably just link it to missing
  industry = loadff(industryfilename,
    ffpath=industryfilepath,
    fieldprefix = "",
    ffscale=true)

  Frets = ffdesc |> deepcopy
  if (!missingunmatched)
    crsp.ffcode[(ismissing).(crsp.ffcode)] .= 31 #new special code
    industry.unknown = zeros(nrow(industry))
    push!(Frets, :unknown)
  end

  if !issorted(crsp, [:permno, :date])
    sort!(crsp, [:permno, :date])
  end


  #now merge the returns
  crspindustry = mergereturns(crsp[!, [:permno, :date, :ret, :siccd]], industry,
    Fgroup=:permno,
    Frets=Frets,
    datefrequency=:day,
    testlabel="industry")

  industrymatrix = Matrix(crspindustry[!, Frets])

  #acquire the industry returns by using ffcode as the column#
  crspindustry.industryret = ((i,j)->
    j≡missing ? missing : industrymatrix[i,j]).(1:nrow(crsp), crsp.ffcode)

  #the below basically just drops the first row for which no industry returns are found
  nfound = combine(groupby(crspindustry, :date), [:industryret]=>
    ((iret)->(iret .!== missing)|>sum) => :numnonmissing)
  crspindustry = innerjoin(crspindustry, nfound, on=:date)
  sort!(crspindustry, :date)

  #optionally set fields with missing industry ret to 0.0
  if !missingunmatched
    crspindustry.industryret[ismissing.(crspindustry.industryret) .&
      (crspindustry.numnonmissing .> 0)] .= 0.0
  end

  crspindustry.plret = crsp.ret .- crspindustry.industryret

  crspindustry.plret .= crspindustry.plret .|> finiteormissing
  crspindustry.plret .=  ((x)->max(x,plretmin)).(crspindustry.plret)

  testpermno::Int = PARAM[:testpermnomult]
  testview = view(crspindustry, (crspindustry[!, :permno] .% testpermno) .== 0, :)
  testview |> CSV.write(
    "$(PARAM[:testpath])\\$(PARAM[:crspdatalabel])_industrypl.csv")

  return crspindustry.plret
end

#widget to create the cumulative return indices for plret
#=function computecumulativeindustryreturns(crsp)
  returns = crsp[!, [:permno, :date, :plret]] |> deepcopy

  #compute the cumulative adjusted returns
  returns.plcumret = missings(Float64, nrow(returns))
  returnsplretcomplete = view(returns, crsp.plret .!== missing, :)
  returnsplretcomplete.plcumret .= computecrspcumulativereturns(
    returnsplretcomplete|>DataFrame, Fret = :plret)

  return returns.plcumret
end=#
