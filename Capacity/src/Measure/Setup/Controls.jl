#TODO: dispatch on each type of control and grouped control to create the fields
#Maybe start with the dispatch to get everything working

#form an empty parametric type for dispatch
struct Control{ControlType} end
Control(ControlType::Symbol) = Control{ControlType}()

function formcontrols!(panel::DataFrame, ms::MeasureSpec;)

  formglobalcontrols!(panel, ms)
  formcontrolsbyt!(panel, ms)
  formcontrolsbykt!(panel, ms)
  checkcontrols(panel, ms)

  @view(panel[(size(panel,1)÷2-1000):(size(panel,1)÷2),:]
    ) |> CSV.write("$(PARAM[:testpath])\\testctrl$(PARAM[:crspdaily])_panel.csv")

end

#compute the weights of the market portfolio
function formmarketweights!(panel)
  panel.wmc = Vector{MFloat64}(undef, size(panel,1))

  spanels::GroupedDataFrame = groupby(panel, :date)
  Threads.@threads for t ∈ 1:length(spanels)
    spanel = spanels[t]
    spanel.wmc .= spanel.mc ./ sum(spanel.mc)
  end
end

#####################
#GLOBAL CONTROLS
######################
function formglobalcontrols!(panel::DataFrame, ms::MeasureSpec,
    controlsglobal::Vector{<:Control}=PARAM[:controlsglobal] .|> Control |> Vector{Control})

  for FW ∈ ms.FW
    #integrity check
    (FW ∈ propertynames(panel)) && error("$(FW) already exists in dataframe")
    panel[!,FW] = Vector{MFloat64}(undef, size(panel,1))
  end

  Threads.@threads for W ∈ controlsglobal
    formglobalcontrol!(panel, W)
  end
end

function formglobalcontrol!(panel, ::Control{:intercept};
  prefix::Symbol = PARAM[:controlsprefix])
  FW = Symbol(prefix, :intercept)
  @assert FW ∈ propertynames(panel)
  @assert all(ismissing.(panel[!,FW]))

  panel[!, FW] .= 1.0
end

function formglobalcontrol!(panel, ::Control{C};
  prefix::Symbol = PARAM[:controlsprefix]) where C
  FW = Symbol(prefix, C)
  @assert FW ∈ propertynames(panel)
  @assert all(ismissing.(panel[!,FW]))

  panel[!, FW] .= panel[!, C]
end


#####################
#Controls by t
######################
function formcontrolsbyt!(panel::DataFrame, ms::MeasureSpec,
    controlsbyt::Vector{<:Control}=PARAM[:controlsbyt] .|> Control |> Vector{Control})

  for FW ∈ ms.FWbyt
    #integrity check
    (FW ∈ propertynames(panel)) && error("$(FW) already exists in dataframe")
    panel[!,FW] = Vector{MFloat64}(undef, size(panel,1))
  end

  spanels::GroupedDataFrame = groupby(panel, :date)
  for W ∈ controlsbyt
    #Threads.@threads for t ∈ 1:length(spanels)
      #spanel = spanels[t]
      formcontrolbyt!(panel, W)
    #end
  end
end

function formcontrolbyt!(panel::AbstractDataFrame, ctrl::Control{:fixedt};
  prefix::Symbol = PARAM[:controlsprefix])
  FW = Symbol(prefix, :fixedt)
  @assert all(ismissing.(panel[!,FW]))

  panel[:, FW] .= 1.0
end

function formcontrolbyt!(panel, ::Control{C};
  prefix::Symbol = PARAM[:controlsprefix]) where C
  FW = Symbol(prefix, C)
  @assert FW ∈ propertynames(panel)
  @assert all(ismissing.(panel[!,FW]))

  panel[!, FW] .= panel[!, C]
end


#####################
#Controls by k and t
######################

function formcontrolsbykt!(panel::DataFrame, ms::MeasureSpec,
    controlsbykt::Vector{<:Control}=Control.(PARAM[:controlsbykt]) |> Vector{Control})


  for Fξ ∈ ms.Fξs
    for FW ∈ ms.FWbykt[Fξ]
      #integrity check
      (FW ∈ propertynames(panel)) && error("$(FW) already exists in dataframe!!")
      panel[!,FW] = Vector{MFloat64}(undef, size(panel,1)) #allocate space for the field
    end
  end

  #@info propertynames(panel)
  ##spanels::GroupedDataFrame = groupby(panel, :date)
  for W ∈ controlsbykt
    #Threads.@threads for t in 1:length(spanels)
    #spanel = spanels[t]
    formcontrolbykt!(panel, ms::MeasureSpec, W)
    #end
  end
end


#controls for interinvestor flows by |Lwk*R_ik| + |wk*R_k|+ |Lwk*R_ik-wk*R_k)
# ∝ |Lwk*R_ik/R_k| + |wk|+ |Lwk*R_ik/R_k-wk|) with adjustment (see Thorne)
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:adjustedinterinvestor};
  prefix::Symbol = PARAM[:controlsprefix])


  #this is iterating over K
  for (ξ, Fw, FRLw, FLw) ∈ zip(ms.ξs, Fweights(ms, :Fw), Fweights(ms, :FRLw), Fweights(ms, :FLw))
    FW = Symbol(prefix, :adjustedinterinvestor, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))

    spanels = groupby(panel, :date)
    Threads.@threads for k ∈ keys(spanels)
      spanel = spanels[k]
      all(ismissing.(spanel[!, FLw])) && continue
      @assert (spanel[1,:date] .== spanel[!,:date]) |> all
      Rtk = 1.0+sum(skipmissing(spanel[!, FRLw]))
      @assert Rtk ≈ 1+ sum(skipmissing(spanel[!, FLw] .* (1 .+ spanel[!, :ret])))
      spanel[:, FW] .= (abs.(spanel[!, FRLw] ./ Rtk) .+ abs.(spanel[!, Fw])
        .- abs.(spanel[!, FRLw] ./ Rtk .- spanel[!, Fw]))
    end
    #=if !(Rtk ≈ 1 + sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
      @info "!(Rtk ≈ sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret])))) !!!"
      error("
        1+sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))): $(
          sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
        Rtk: $Rtk")
    end=#


  end

  return nothing
end

#controls for interinvestor flows by |wk|.+|RLwk|
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:interinvestorold};
  prefix::Symbol = PARAM[:controlsprefix])


  #this is iterating over K
  for (ξ, Fw, FRLw, FLw) ∈ zip(ms.ξs, Fweights(ms, :Fw), Fweights(ms, :FRLw), Fweights(ms, :FLw))
    FW = Symbol(prefix, :interinvestorold, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))

    spanels = groupby(panel, :date)
    Threads.@threads for k ∈ keys(spanels)
      spanel = spanels[k]
      all(ismissing.(spanel[!, FLw])) && continue
      @assert (spanel[1,:date] .== spanel[!,:date]) |> all
      Rtk = 1.0+sum(skipmissing(spanel[!, FRLw]))
      @assert Rtk ≈ 1+ sum(skipmissing(spanel[!, FLw] .* (1 .+ spanel[!, :ret])))
      spanel[:, FW] .= abs.(spanel[!, FRLw] ./ Rtk) .+ abs.(spanel[!, Fw])
    end
    #=if !(Rtk ≈ 1 + sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
      @info "!(Rtk ≈ sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret])))) !!!"
      error("
        1+sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))): $(
          sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
        Rtk: $Rtk")
    end=#


  end

  return nothing
end

#controls for interinvestor flows by |RLwk|+|wk|
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:interinvestor};
  prefix::Symbol = PARAM[:controlsprefix])


  #this is iterating over K
  for (ξ, Fw#=, FRLw, FLw=#) ∈ zip(ms.ξs, Fweights(ms, :Fw)#=, Fweights(ms, :FRLw), Fweights(ms, :FLw)=#)
    FW = Symbol(prefix, :interinvestor, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))

    spanels = groupby(panel, :date)
    Threads.@threads for k ∈ keys(spanels)
      spanel = spanels[k]
      #all(ismissing.(spanel[!, FLw])) && continue
      @assert (spanel[1,:date] .== spanel[!,:date]) |> all
      #Rtk = 1.0+sum(skipmissing(spanel[!, FRLw]))
      #@assert Rtk ≈ 1+ sum(skipmissing(spanel[!, FLw] .* (1 .+ spanel[!, :ret])))
      spanel[:, FW] .= abs.(spanel[!, Fw])
    end
    #=if !(Rtk ≈ 1 + sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
      @info "!(Rtk ≈ sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret])))) !!!"
      error("
        1+sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))): $(
          sum(skipmissing(panel[!, FLw] .* (1 .+ panel[!, :ret]))))
        Rtk: $Rtk")
    end=#


  end

  return nothing
end

#controls for inflow from market by |wk-wm|-|wk|
#(the last term ensures we count the allocation in the focal variable)
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:marketin};
  prefix::Symbol = PARAM[:controlsprefix])

  for (ξ, Fw) ∈ zip(ms.ξs, Fweights(ms, :Fw))
    FW = Symbol(prefix, :marketin, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))
    panel[:,FW] .= abs.(panel[!, Fw] .- panel.wmc) .- abs.(panel[!, Fw])
  end

  return nothing
end

#controls for inflow from market by |wk-wm|-|wk|
#(the last term ensures we count the allocation in the focal variable)
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:marketinshort};
  prefix::Symbol = PARAM[:controlsprefix])

  for (ξ, Fw) ∈ zip(ms.ξs, Fweights(ms, :Fw))
    FW = Symbol(prefix, :marketinshort, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))
    panel[:,FW] .= abs.(panel[!, Fw] .+ panel.wmc) .- abs.(panel[!, Fw])
  end

  return nothing
end

#controls for outflow from market by |RLwk-wm|-|RLwk|
#(the last term ensures we count the allocation in the focal variable)
function formcontrolbykt!(panel, ms::MeasureSpec, ::Control{:marketout};
  prefix::Symbol = PARAM[:controlsprefix])

  for (ξ, FRLw) ∈ zip(ms.ξs, Fweights(ms, :FRLw))
    FW = Symbol(prefix, :marketout, ξ.Fξ)
    @assert all(ismissing.(panel[!, FW]))
    panel[:,FW] .= abs.(panel[!, FRLw] .- panel.wmc) .- abs.(panel[!, FRLw])
  end

  return nothing
end


#####CHECK Controls
#every control must have a check, so this funciton is kinda long
function checkcontrols(panel::DataFrame, ms;
  controls = [PARAM[:controlsglobal]; PARAM[:controlsbyt]; PARAM[:controlsbykt];] .|> Control,
  prefix::Symbol = PARAM[:controlsprefix])

  Fcontrols = Vector{Symbol}()

  #############List of control check functions starts here
  #fallback method
  function checkcontrol(::Control{C}) where C
    Fcontrol = Symbol(prefix,C)
    try
      @assert (panel[!,Fcontrol] .=== panel[!,C]) |> all
    catch err
      panel[1:20_000,:] |> CSV.write("$(PARAM[:testpath])\\checkcontroldump$C.csv")

      @error "Error $err for control $C" exception=(err, catch_backtrace())
      throw(err)
    end
    push!(Fcontrols, Fcontrol)
  end

  #global intercept method
  function checkcontrol(::Control{:intercept})
    Fcontrol = Symbol(prefix,:intercept)
    @assert (panel[!,Fcontrol] .== 1.0) |> all
    push!(Fcontrols, Fcontrol)
  end


  function checkcontrol(::Control{:fixedt})
    Fcontrol = Symbol(prefix,:fixedt)
    @assert (panel[!,Fcontrol] .== 1.0) |> all
    push!(Fcontrols, Fcontrol)
  end


  checkcontrol(W::T) where T<:Union{Control{:interinvestor},
    Control{:adjustedinterinvestor}, Control{:interinvestorold}} = (ξ->checkcontrol(ξ,W)).(ms.ξs)
  function checkcontrol(ξ, ::Control{:adjustedinterinvestor})
    Fcontrol = Symbol(prefix,:adjustedinterinvestor, ξ.Fξ)

    #check missings were handled correctly
    completerows = completecases(panel, [ξ.X[:Fw], ξ.X[:FRLw]])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    #now check the contents

    for spanel ∈ groupby(view(panel, completerows, :), :date)
      Radj::Float64 = 1.0 + sum(spanel[!, ξ.X[:FRLw]])
      @assert (spanel[!,Fcontrol] .≈
        abs.(spanel[!, ξ.X[:Fw]]) .+ abs.(spanel[!, ξ.X[:FRLw]]./Radj)
        .- abs.(spanel[!, ξ.X[:FRLw]]./Radj .- spanel[!, ξ.X[:Fw]])) |> all
    end
    push!(Fcontrols, Fcontrol)
  end

  function checkcontrol(ξ, ::Control{:interinvestorold})
    Fcontrol = Symbol(prefix,:interinvestorold, ξ.Fξ)

    #check missings were handled correctly
    completerows = completecases(panel, [ξ.X[:Fw], ξ.X[:FRLw]])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    #now check the contents
    for spanel ∈ groupby(view(panel, completerows, :), :date)
      Radj::Float64 = 1.0 + sum(spanel[!, ξ.X[:FRLw]])
      @assert (spanel[!,Fcontrol] .≈
        abs.(spanel[!, ξ.X[:Fw]]) .+ abs.(spanel[!, ξ.X[:FRLw]]./Radj)) |> all
    end

    push!(Fcontrols, Fcontrol)
  end

  function checkcontrol(ξ, ::Control{:interinvestor})
    Fcontrol = Symbol(prefix,:interinvestor, ξ.Fξ)

    #check missings were handled correctly
    completerows = completecases(panel, [ξ.X[:Fw]#=, ξ.X[:FRLw]=#])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    panel[!,Fcontrol] .≈ #=abs.(panel[!, ξ.X[:FRLw]])+=#abs.(panel[!, ξ.X[:Fw]])

    push!(Fcontrols, Fcontrol)
  end

  checkcontrol(W::Control{:marketin}) = (ξ->checkcontrol(ξ,W)).(ms.ξs)
  function checkcontrol(ξ, ::Control{:marketin})
    Fcontrol = Symbol(prefix,:marketin, ξ.Fξ)

    completerows = completecases(panel, [ξ.X[:Fw]])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    @assert (panel[!,Fcontrol] .≈
      abs.(panel.wmc .- panel[!, ξ.X[:Fw]]) .- abs.(panel[!, ξ.X[:Fw]])) |> all
    push!(Fcontrols, Fcontrol)
  end
  checkcontrol(W::Control{:marketinshort}) = (ξ->checkcontrol(ξ,W)).(ms.ξs)
  function checkcontrol(ξ, ::Control{:marketinshort})
    Fcontrol = Symbol(prefix,:marketinshort, ξ.Fξ)

    completerows = completecases(panel, [ξ.X[:Fw]])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    @assert (panel[!,Fcontrol] .≈
      abs.(panel.wmc .+ panel[!, ξ.X[:Fw]]) .- abs.(panel[!, ξ.X[:Fw]])) |> all
    push!(Fcontrols, Fcontrol)
  end
  checkcontrol(W::Control{:marketout}) = (ξ->checkcontrol(ξ,W)).(ms.ξs)
  function checkcontrol(ξ, ::Control{:marketout})
    Fcontrol = Symbol(prefix,:marketout, ξ.Fξ)

    completerows = completecases(panel, [ξ.X[:FRLw]])
    @assert @views (!ismissing).(panel[completerows, Fcontrol]) |> all
    @assert @views (ismissing).(panel[(!).(completerows), Fcontrol]) |> all

    Fcontrol = Symbol(prefix,:marketout, ξ.Fξ)
    spanel = view(panel, completerows, :)
    @assert (spanel[!,Fcontrol] .≈
      abs.(spanel.wmc .- spanel[!, ξ.X[:FRLw]]) .- abs.(spanel[!, ξ.X[:FRLw]])) |> all
    push!(Fcontrols, Fcontrol)
  end
  ######end recalculation checks

  for W ∈ controls
    checkcontrol(W)
  end

  #check that all controls are accounted for
  Fcontrolscheck::Vector{Symbol} = [ms.FW; ms.FWbyt]
  for Fξ ∈ ms.Fξs
    for FW ∈ ms.FWbykt[Fξ]
      push!(Fcontrolscheck, FW)
    end
  end
  @assert issetequal(Fcontrols, Fcontrolscheck)
end
