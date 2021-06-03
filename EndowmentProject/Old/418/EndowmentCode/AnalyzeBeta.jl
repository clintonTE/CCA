
function βdf(data::NCCSData; minPointsForβ::Int = MIN_POINTS_FOR_BETA,
  path::String=WORKING_PATH, charitytypes::Symbol = :all)
  local βdf::DataFrame = DataFrame(ein = unique(data.df[:ein]))
  local Nβ::Int = size(βdf,1)
  local reg::FMLM

  local plotNames::Vector{String} = Vector{String}()
  local plots::Vector{PlotContainer} = Vector{PlotContainer}()

  #NOTE: DELETE THE BELOW


  (charitytypes ≠ :all) && (βdf = βdf[βdf[:charitytype] .== charitytypes, :])

  #initialize the rest of the dataframe
  βdf[:b0return] = Vector{MFloat64}(missing,Nβ)
  βdf[:b1return] = similar(βdf[:b0return])

  βdf[:b0contributions] = similar(βdf[:b0return])
  βdf[:b1contributions] = similar(βdf[:b0return])

  βdf[:b0contribonret] = similar(βdf[:b0return])
  βdf[:b1contribonret] = similar(βdf[:b0return])

  βdf[:b0expenses] = similar(βdf[:b0return])
  βdf[:b1expenses] = similar(βdf[:b0return])

  βdf[:avgwealth] = similar(βdf[:b0return])
  βdf[:N] = Vector{MInt}(missing, Nβ)
  βdf[:charitytype] = Vector{MSymbol}(missing, Nβ)

  #index the dataframe
  βdfindex::Dict = Dict(βdf[i,:ein] => βdf[i,:] for i::Int  ∈ 1:Nβ)

  local YSpecs::Vector{Symbol} = Vector{Symbol}() # holds the LHS
  local XSpecs::Vector{FMExpr} = Vector{FMExpr}() #holds the RHS
  local XNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() # holds the names of the x variables


  #regress excess returns on excess s&p
  push!(XSpecs, Meta.parse("sp500ret"))
  push!(XNames, [:intercept; :sp500ret])
  push!(YSpecs, :excess)

  #regress contributions on excess s&p
  push!(XSpecs, Meta.parse("sp500ret"))
  push!(XNames, [:intercept; :lagsp500ret])
  push!(YSpecs, :pcontributions)

  #regress contributions on excess returns
  push!(XSpecs, Meta.parse("excess"))
  push!(XNames, [:intercept; :excess])
  push!(YSpecs, :pcontributions)

  #regress contributions on excess returns
  push!(XSpecs, Meta.parse("sp500ret"))
  push!(XNames, [:intercept; :sp500ret])
  push!(YSpecs, :pexpenses)

  println("beginning beta computations")
  #print(data.df[50000:50100, [:ein, :datefiscal, :return, :sp500ret]])
  for subdf ∈ groupby(data.df, :ein) #compute the beta for each firm
    N::Int = sum(completecases(subdf[[:excess, :sp500ret, :pcontributions, :pexpenses]]))
    βrow = βdfindex[subdf[1,:ein]]
    βrow[:N] = N
    βrow[:charitytype] = subdf[1,:charitytype]


    if N ≥ minPointsForβ
      local models::Vector{FMLM} = Vector{FMLM}() #holds the regression models

      #####run the regressions
      try
        for i ∈ 1:length(XSpecs)
          push!(models, FMLM(subdf, XSpecs[i], YSpecs[i],
            XNames = XNames[i], YName=YSpecs[i]))
        end
      catch err
        println(subdf[[:ein, :datefiscal, :excess, :sp500ret]])
        error("Regression failed. $err")
      end

      #store the values
      βrow[:b0return] = models[1].β[1]
      βrow[:b1return] = models[1].β[2]
      βrow[:b0contributions] = models[2].β[1]
      βrow[:b1contributions] = models[2].β[2]
      βrow[:b0contribonret] = models[3].β[1]
      βrow[:b1contribonret] = models[3].β[2]
      βrow[:b0expenses] = models[4].β[1]
      βrow[:b1expenses] = models[4].β[2]
      βrow[:avgwealth] =
        mean(skipmissing((ismissing).(subdf[:excess]) .* subdf[:adjnetassets]))
    end


  end


  #histogram for betas of returns
  βdf = βdf[completecases(βdf[[:b1return, :b0return]]),:]
  push!(plotNames, "return_beta_hist_$charitytypes")
  push!(plots, plot(βdf[(f->abs(f)<1.1).(βdf[:b1return]),:], x=:b1return,
    Guide.title("histogram of return betas"), Geom.histogram(bincount=25)))


  #histogram for betas of the contributions
  push!(plotNames, "return_beta_contrib_$charitytypes")
  βdf = βdf[completecases(βdf[[:b1contributions, :b0contributions]]),:]
  push!(plots, plot(βdf[(f->abs(f)<1.1).(βdf[:b1contributions]),:], x=:b1contributions,
    Guide.title("histogram of contribution betas"), Geom.histogram(bincount=25)))

  #histogram for betas of the contribonret
  βdf = βdf[completecases(βdf[[:b1contribonret, :b0contribonret]]),:]
  push!(plotNames, "return_beta_ret_on_contrib_$charitytypes")
  push!(plots,  plot(βdf[(f->abs(f)<1.1).(βdf[:b1contribonret]),:], x=:b1contribonret,
    Guide.title("histogram of return on contribution betas"), Geom.histogram(bincount=25)))

    #histogram for betas of expenses
  βdf = βdf[completecases(βdf[[:b1expenses, :b0contribonret]]),:]
  push!(plotNames, "return_beta_expense_$charitytypes")
  push!(plots, plot(βdf[(f->abs(f)<1.1).(βdf[:b1expenses]),:], x=:b1expenses,
      Guide.title("histogram of expense betas"), Geom.histogram(bincount=25)))

  println("return beta computations complete with median $(median(βdf[:b1return])),
    mean $(mean(βdf[:b1return])),
    weighted mean $(mean(βdf[:b1return], Weights(Vector{Float64}(βdf[:avgwealth])))),
    sd $(std(βdf[:b1return])), and se $(std(βdf[:b1return])/size(βdf,1)^0.5)")

  println("contribution on return beta computations complete with median $(median(βdf[:b1contributions])),
    mean $(mean(βdf[:b1contributions])),
    weighted mean $(mean(βdf[:b1contributions], Weights(Vector{Float64}(βdf[:avgwealth])))),
    sd $(std(βdf[:b1contributions])), and se $(std(βdf[:b1contributions])/size(βdf,1)^0.5)")


  println("contribution on return beta computations complete with median $(median(βdf[:b1contribonret])),
    mean $(mean(βdf[:b1contribonret])),
    weighted mean $(mean(βdf[:b1contribonret], Weights(Vector{Float64}(βdf[:avgwealth])))),
    sd $(std(βdf[:b1contribonret])), and se $(std(βdf[:b1contribonret])/size(βdf,1)^0.5)")

  println("expense beta computations complete with median $(median(βdf[:b1expenses])),
    mean $(mean(βdf[:b1expenses])),
    weighted mean $(mean(βdf[:b1expenses], Weights(Vector{Float64}(βdf[:avgwealth])))),
    sd $(std(βdf[:b1expenses])), and se $(std(βdf[:b1expenses])/size(βdf,1)^0.5)")

  bronbc = FMLM(βdf, Meta.parse("b1return"), :b1contributions,
    XNames = [:b1return], YName=:b1contributions)
  println("beta contributions on beta returns is $(bronbc.β)")
  #NOTE: THinking of a table which has something like the above stats with a pooled beta stat too

  for i ∈ 1:length(plots) #write the graphs
    try
      draw(SVG("$path\\$(plotNames[i]).svg", 9inch, 7inch), plots[i])
    catch err
      @error "ERROR writing $(plotNames[i]). Error message $err"
      println("Size: $(size(subdf))")
    end
    #println("Graph $(plotNames[i]) written.")
  end

  return βdf
end


function processβ(data::NCCSData)::NCCSData
  βdf(data)

  return data
end
