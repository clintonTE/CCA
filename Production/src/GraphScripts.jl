


function analyzescenarios(Ψ::Simulation,
  τscenarios::Vector{Vector{Float64}},
  τscenarionames::Vector{Symbol};
  prefix::String = "",
  outpaper::String=OUT_PAPER,
  extension::String = "pdf")

  local x::UnitRange = 1:Ψ.nsims
  local y::Matrix{Float64}

  gr()
  nscenarios::Int = length(τscenarios)

  highproductionplots = Vector()
  lowproductionplots = Vector()
  investmentplots = Vector()
  capitalplots = Vector()
  valueplots = Vector()

  for i ∈ 1:nscenarios
    Random.seed!(11)
    updatefirmτ!(Ψ.Φ, τscenarios[i])
    Simulate!(Ψ)

    legend::Union{Bool, Symbol} = i==1 ? :topright : false
    xlab::String = i==nscenarios ? "period" : ""

    p = plot(x, Ψ.Yh, lab = "high state", legend=legend)
    plot!(p, x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    #plot!(p, x, Ψ.Y)
    xaxis!(p, xlab)#"period")
    yaxis!(p, "production")
    plot!(p, title = "High State Allocations:  $(τscenarionames[i])")
    push!(highproductionplots, p)

    p = plot(x, Ψ.Yl, lab = "low state", legend=legend)
    plot!(p, x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    #plot!(p, x, Ψ.Y)
    xaxis!(p, xlab)
    yaxis!(p, "production")
    plot!(p, title = "Low State Allocations:  $(τscenarionames[i])")
    push!(lowproductionplots, p)

    p = plot(x, Ψ.I, lab = "investment", legend=legend)
    plot!(p, x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    #plot!(p, x, Ψ.Y)
    xaxis!(p, xlab)
    yaxis!(p, "investment")
    plot!(p, title = "Firm Investment:  $(τscenarionames[i])")
    push!(investmentplots, p)

    p = plot(x, Ψ.K, lab = "capital", legend=legend)
    plot!(p, x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    #plot!(p, x, Ψ.Y)
    xaxis!(p, xlab)
    yaxis!(p, "capital")
    plot!(p, title = "Firm Capital:  $(τscenarionames[i])")
    push!(capitalplots, p)

    p = plot(x, Ψ.V, lab = "firm value", legend=legend)
    #plot!(twinx(p), x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    plot!(p, x, Ψ.S , lab = "economic state (H/L)", line=:dot)
    #plot!(p, x, Ψ.Y)
    xaxis!(p, xlab)
    yaxis!(p, "firm value")
    plot!(p, title = "Firm Value:  $(τscenarionames[i])")
    push!(valueplots, p)
  end

  statespervar::Int = Int(round(sqrt(Ψ.Φ.N)))

  outproduction = plot(highproductionplots..., layout= (nscenarios,1), size=(1000, 900))
  savefig(outproduction, "$outpaper\\$(prefix)highscen_s$(statespervar)_n$(Ψ.nsims).$extension")

  outproduction = plot(lowproductionplots..., layout= (nscenarios,1), size=(1000, 900))
  savefig(outproduction, "$outpaper\\$(prefix)lowscen_s$(statespervar)_n$(Ψ.nsims).$extension")

  outproduction = plot(investmentplots..., layout= (nscenarios,1), size=(1000, 900))
  savefig(outproduction, "$outpaper\\$(prefix)investment_s$(statespervar)_n$(Ψ.nsims).$extension")

  outproduction = plot(capitalplots..., layout= (nscenarios,1), size=(1000, 900))
  savefig(outproduction, "$outpaper\\$(prefix)capital_s$(statespervar)_n$(Ψ.nsims).$extension")

  outproduction = plot(valueplots..., layout= (nscenarios,1), size=(1000, 900))
  savefig(outproduction, "$outpaper\\$(prefix)value_s$(statespervar)_n$(Ψ.nsims).$extension")
end

mutable struct τmeans
  τvals::Vector{Float64}
  highproductionmean::Vector{Float64}
  lowproductionmean::Vector{Float64}
  investmentmean::Vector{Float64}
  capitalmean::Vector{Float64}
  valuemean::Vector{Float64}
end

#this method does the tax analysis
function analyzetaxes(Ψ::Simulation,
  τvals::Vector{Float64};
  outpath = OUT_PATH, prefix::String = "")


  #these will hold the means
  nvals::Int = length(τvals)
  highproductionmean = Vector{Float64}(undef, nvals)
  lowproductionmean = Vector{Float64}(undef, nvals)
  investmentmean = Vector{Float64}(undef, nvals)
  capitalmean = Vector{Float64}(undef, nvals)
  valuemean = Vector{Float64}(undef, nvals)

  statespervar::Int = Int(round(sqrt(Ψ.Φ.N)))
  for i ∈ 1:nvals
    Random.seed!(11)
    updatefirmτ!(Ψ.Φ, [τvals[i], 0.0])
    Simulate!(Ψ)

    highproductionmean[i] = mean(Ψ.Yh)
    lowproductionmean[i] = mean(Ψ.Yl)
    investmentmean[i]  = mean(Ψ.I)
    capitalmean[i]  = mean(Ψ.K)
    valuemean[i]  = mean(Ψ.V)
  end

  writeto::String = "$outpath\\$(prefix)highstatetax_s$(statespervar)_n$(Ψ.nsims).jls"
  T::τmeans = τmeans(τvals, highproductionmean, lowproductionmean, investmentmean, capitalmean, valuemean)
  serialize(writeto, T)
  println("Wrote tax summary data to $writeto")

  for i ∈ 1:nvals
    Random.seed!(11)
    updatefirmτ!(Ψ.Φ, [0.0, τvals[i]])
    Simulate!(Ψ)

    highproductionmean[i] = mean(Ψ.Yh)
    lowproductionmean[i] = mean(Ψ.Yl)
    investmentmean[i]  = mean(Ψ.I)
    capitalmean[i]  = mean(Ψ.K)
    valuemean[i]  = mean(Ψ.V)
  end


  writeto = "$outpath\\$(prefix)lowstatetax_s$(statespervar)_n$(Ψ.nsims).jls"
  T= τmeans(τvals, highproductionmean, lowproductionmean, investmentmean, capitalmean, valuemean)
  serialize(writeto, T)
  println("Wrote tax summary data to $writeto")
end

function graphtaxanalysis(Ψ::Simulation;
  outpath = OUT_PATH,
  prefix::String = "",
  statespervar::Int = Int(round(sqrt(Ψ.Φ.N))),
  nsims::Int = Ψ.nsims,
  outpaper::String=OUT_PAPER,
  extension::String = "pdf"
  )

  gr()

  local plotholder::Vector = Vector()
  local plotnames::Vector{String} = Vector{String}()

  writefrom::String = "$outpath\\$(prefix)highstatetax_s$(statespervar)_n$(Ψ.nsims).jls"
  Thigh::τmeans = deserialize(writefrom)
  writefrom = "$outpath\\$(prefix)lowstatetax_s$(statespervar)_n$(Ψ.nsims).jls"
  Tlow = deserialize(writefrom)


  Ts::Vector{τmeans} = [Thigh, Tlow]
  taxstatenames::Vector{String} = ["high_state_taxes", "low_state_taxes"]

  for i ∈ 1:length(Ts)
    legend::Union{Bool, Symbol} = false
    xlab = "Tax Rates"

    p = plot(Ts[i].τvals, Ts[i].highproductionmean, lab = "high production alloc.", legend=:topright, size=(1000, 900))
    plot!(p, Ts[i].τvals, Ts[i].lowproductionmean , lab = "low production alloc.")
    xaxis!(p, xlab)#
    yaxis!(p, "mean production")
    plot!(p, title = "$(taxstatenames[i]): High and Low State Production")
    push!(plotholder, p)
    push!(plotnames, "production_$(taxstatenames[i])" )


    p = plot(Ts[i].τvals, Ts[i].capitalmean, lab = "capital", legend=legend, size=(1000, 900))
    xaxis!(p, xlab)#
    yaxis!(p, "mean capital")
    plot!(p, title = "$(taxstatenames[i]): Capital")
    push!(plotholder, p)
    push!(plotnames, "capital_$(taxstatenames[i])" )

    p = plot(Ts[i].τvals, Ts[i].investmentmean, lab = "investment", legend=legend, size=(1000, 900))
    xaxis!(p, xlab)#
    yaxis!(p, "mean investment")
    plot!(p, title = "$(taxstatenames[i]): Investment")
    push!(plotholder, p)
    push!(plotnames, "capital_$(taxstatenames[i])" )

    p = plot(Ts[i].τvals, Ts[i].valuemean, lab = "value", legend=legend, size=(1000, 900))
    xaxis!(p, xlab)#
    yaxis!(p, "mean value")
    plot!(p, title = "$(taxstatenames[i]): Value")
    push!(plotholder, p)
    push!(plotnames, "value_$(taxstatenames[i])" )

    #=p = plot(Ts[i].τvals, Ts[i].lowproductionmean, lab = "capital", legend=legend, size=(1000, 900))
    xaxis!(p, xlab)#
    yaxis!(p, "mean production")
    plot!(p, title = "$(taxstatenames[i]): Capital")
    push!(plotholder, p)
    push!(plotnames, "value_$(taxstatenames[i])" )=#
  end

  outplotholder = Vector()

  push!(outplotholder, plot(plotholder[1:4]...,  layout= (2,2)))
  push!(outplotholder, plot(plotholder[5:8]...,  layout= (2,2)))
  for i ∈ 1:length(outplotholder)
    savefig(outplotholder[i], "$outpaper\\$(prefix)$(taxstatenames[i])_s$(statespervar)_n$(Ψ.nsims).$extension")
  end

end

function graphscript!(Ψ::Simulation;
  τvals = Vector{Float64}(),
  τscenarios::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),
  τscenarionames::Vector{Symbol} = Vector{Symbol}(),
  scenariographs::Bool = true,
  prefix::String = "", taxanalysis::Bool=true, taxgraphs::Bool= true)

  scenariographs && analyzescenarios(Ψ, τscenarios, τscenarionames,prefix=prefix)
  taxanalysis && analyzetaxes(Ψ::Simulation, τvals, prefix=prefix)
  taxgraphs && graphtaxanalysis(Ψ::Simulation, prefix=prefix)
end
