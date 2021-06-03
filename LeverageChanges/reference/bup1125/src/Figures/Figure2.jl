function makefigure2(dist::DataFrame;
  Fret::Symbol = DIST_RET)

  local sdist::SubDataFrame

  generateagg(df::AbstractDataFrame) = aggregatebyday(x->mean(x)*100, df, selectcols=[:day, Fret])

  ydomain::Vector{Float64} = [-0.2,0.2]
  yaxistext::String = "Mean (daily cross- in %)"
  figurename="figure2"

  sdist = view(dist, (r->r.primary && (!ismissing(r[Fret]))).(eachrow(dist)), :)

  makedividendeventplots(generateagg, sdist,
    Fret=Fret,
    ydomain=ydomain,
    yaxistext=yaxistext,
    figurename=figurename
    )
end
