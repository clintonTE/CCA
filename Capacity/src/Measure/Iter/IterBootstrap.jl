
#conjugate gradient solver
#Note on line searches- HZ much better than backtracking
function solveA!(loss, Θ::AbstractVolumePartsIter{TM,TV,T},
    Xv::AbstractXYIter{TM, TV, T},
    SolveType::Val{:bootstrap};
    RegressionType::Val = Val(PARAM[:iterbootstrapregressiontype])) where {TM, TV, T}

  dims = (K=size(Θ.A,1), T=size(Θ.A,2))

  if (Θ.G[1,1] .≈ Θ.G) |> all
    throw("All G values are the same!!! Implies G was not sourced
    bootstrapped correctly from the source")
  end

  #compute the loss and the values from A
  λ = updateA₁!(loss, Θ, Xv, RegressionType)
  updateAfromGAₜ!(Θ,Xv, 1)
  @assert λ ≈ lossfromA(Θ, Xv, RegressionType)

  #lower::Vector{T} = fill(lowerbound, (dims.T-1) * dims.K)
  #upper::Vector{T} = fill(upperbound, (dims.T-1) * dims.K)

  b = IOBuffer()
  write(b, "\n\n***Boostrap additional simulation parameters:\n" *
    "SolveType: $SolveType" *
    "measuretype: $(PARAM[:measuretype])\n" *
    "regressionmethod: $(PARAM[:regressionmethod])\n" *
    "controlsbykt: $(PARAM[:controlsbykt])\n" *
    "controlsbyt: $(PARAM[:controlsbyt])\n" *
    "controlsglobal: $(PARAM[:controlsglobal])\n" *
    "growth type: N/A (level since bootstrap)\n" *
    "momentum weight spec: $(PARAM[:momentumweightspec])\n" *
    "sorted default thresholds: $(PARAM[:sorteddefaultthresholds])\n" *
    "regression method (bootstrap): $(RegressionType)\n" *
    "Types (Matrix, Vector, T): $TM, $TV, $T\n" *
    "bounds: disabled\n"
    )
  write(b, "\n$(formatparameters())\n")

  iterinfo = String(take!(b))
  #println(iterinfo)
  open("$(PARAM[:outputpath])\\$(iterresultsname(λ))_info(bs).txt", "w+") do f
    write(f, iterinfo)
  end
  return λ
end
