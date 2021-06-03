



#runs a regression on every firm
function regressbynonprofit!(df::DataFrame;
  YField::Symbol = :lreturn,
  bmfield::Symbol = Symbol(:pred_, YField),
  XFields::Vector{Symbol} = [:lsp500ret],
  groupbyfield::Symbol = :ein,
  suppressintercept::Bool = false,
  XNames::Vector{Symbol} = [:intercept; XFields],
  YName::Symbol = YField,
  βNames::Vector{Symbol} = (s::Symbol->Symbol(bmfield, "b", s)).(XNames),
  minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA,
  predictwithresults::Bool = true,
  regstring::String = string(suppressintercept ? "0 + " : "", join((string).(XFields), " + ")),
  checkunique::Bool = false
  )::Nothing


  #make sure we have the right number of names
  local N::Int = size(df,1)

  @assert length(XNames) == length(βNames)
  @assert length(XFields) == length(βNames) - 1 + suppressintercept

  local XSpec::FMExpr = Meta.parse(regstring)

  #preallcoate the arrays for the coefficients and prediction
  #println(skipmissing(df[1:1000, :pcontributions]))
  for f::Symbol ∈ βNames
    df[f] = Vector{MFloat64}(undef, N)
  end

  predictwithresults && (df[bmfield] = Vector{MFloat64}(undef, N))
  #split/apply/combine
  ctr=0

  #cache the fields we will use to test if there are enough data points
  checkfields::Vector{Symbol} = [YField; XFields;]

  #loop over all unique eins
  for subdf ∈ groupby(df, groupbyfield)

    local subN::Int
    if checkunique
      ssubdf = view(subdf, completecases(subdf[checkfields]), checkfields)
      subN = size(unique(ssubdf[XFields]),1)
    else
      subN = sum(completecases(subdf[checkfields]))
    end

    #make sure we have enough points for the regression
    if subN ≥ minpointspluskforβ + length(XNames)
      local lm::FMLM

      try
        lm = FMLM(subdf, XSpec, YField, XNames = XNames, YName=YName)
      catch err
        println(subdf[[:ein; :fisyr; checkfields]])
        error("ERROR: $err")
      end

      #predict the resutls and write down the coefficients
      predictwithresults && (subdf[bmfield] .= suppressintercept ? 0.0 : lm.β[1])
      for i ∈ 1:length(βNames)
        subdf[βNames[i]] .=  lm.β[i]
        if predictwithresults && (XNames[i] ≠ :intercept) #include the factor influence
          subdf[bmfield] .+= lm.β[i] .* subdf[XNames[i]]
        end
      end
    end
  end

  return nothing
end

#helper function signature which takes in the data object
regressbynonprofit!(data::NCCSData;
  YField::Symbol = :lreturn,
  bmfield::Symbol = Symbol(:pred_, YField),
  XFields::Vector{Symbol} = [:lsp500ret],
  groupbyfield::Symbol = :ein,
  suppressintercept::Bool = false,
  XNames::Vector{Symbol} = [:intercept; XFields],
  YName::Symbol = YField,
  βNames::Vector{Symbol} = (s::Symbol->Symbol(bmfield, "b", s)).(XNames),
  minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA) = (
    regressbynonprofit!(data.df,
      bmfield=bmfield,
      YField = YField,
      XFields = XFields,
      groupbyfield = groupbyfield,
      suppressintercept = suppressintercept,
      XNames = XNames,
      YName = YName,
      βNames = βNames,
      minpointspluskforβ = minpointspluskforβ))

function βdf(data::NCCSData; minpointspluskforβ::Int = MIN_POINTS_PLUS_K_FOR_BETA,
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
  #βdf[:N] = Vector{MInt}(missing, Nβ)
  βdf[:charitytype] = Vector{MSymbol}(missing, Nβ)

  local YFields::Vector{Symbol} = Vector{Symbol}() # holds the LHS
  local XFields::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() #holds the RHS
  local βNames::Vector{Vector{Symbol}} = Vector{Vector{Symbol}}() # holds the names of the x variables


  #regress excess returns on excess s&p
  push!(XFields, [:sp500ret])
  push!(βNames, [:b0return, :b1return])
  #push!(XNames, [:intercept; :sp500ret])
  push!(YFields, :excess)

  #regress contributions on excess s&p
  push!(XFields, [:sp500ret])
  push!(βNames, [:b0contributions, :b1contributions])
  #push!(XNames, [:intercept; :lagsp500ret])
  push!(YFields, :pcontributions)

  #regress contributions on excess returns
  push!(XFields, [:excess])
  push!(βNames, [:b0contribonret, :b1contribonret])
  #push!(XNames, [:intercept; :excess])
  push!(YFields, :pcontributions)

  #regress contributions on excess returns
  push!(XFields, [:sp500ret])
  push!(βNames, [:b0expenses, :b1expenses])
  #push!(XNames, [:intercept; :sp500ret])
  push!(YFields, :pexpenses)

  println("beginning beta computations")
  for i ∈ 1:length(YFields)

    regressbynonprofit!(data.df, YField=YFields[i], XFields=XFields[i],
      βNames=βNames[i], minpointspluskforβ=minpointspluskforβ)
  end

  allβ::Vector{Symbol} = [βNames...;]

  #println(data.df[100:110,[:ein; allβ;]])
  #println("Complete cases of b0return: $(sum(!ismissing(data.df[:b0return])))")

  #βdf = join(βdf, data.df[[:ein; :charitytype; allβ;]], on = :ein, kind=:left)
  #index the dataframe
  βdfindex::Dict = Dict(βdf[i,:ein] => βdf[i,:] for i::Int  ∈ 1:Nβ)
  for subdf ∈ groupby(data.df, :ein) #compute the average wealth for each firm
    subN::Int = sum(completecases(subdf[[:excess, :sp500ret, :pcontributions, :pexpenses]]))
    βrow = βdfindex[subdf[1,:ein]]
    βrow[:charitytype] = subdf[1,:charitytype]


    if subN ≥ minpointspluskforβ
      for f ∈ allβ
        βrow[f] = subdf[1,f]
      end

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
