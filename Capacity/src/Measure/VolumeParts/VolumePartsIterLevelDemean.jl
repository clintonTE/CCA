#mostly legacy code used in testing


#the below funciton is now only for testing purposes
function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterLevel,
  RegressionType::Val) where {TM<:Matrix,TV<:Vector,T}
  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = initializeÃandfactor(Θ,Xv)

  #run the regression for A₁
  try
    A[:,1] .= reg(factored, Xv.v, RegressionType)
  catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage > 1000) ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || error("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])")

    A[:,1] .= reg(factored, Xv.v, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallbackreg fallback val=$(A[:,1])"
  end

  #println("A[:,1]: $(A[:,1])")

  #absμ .= vec(sum(absμₖ,dims=2))
  absμ = factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end

#primary prediction alogrithm for the demeaned version
function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterDemean,  RegressionType::Val
    ) where {TM,TV,T}
  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  prevA₁ = A[:,1]
  factored::TM = Xv.M * initializeÃandfactor(Θ,Xv)


  #run the regression for A₁
  A[:,1] .= reg(factored, Xv.ṽ, RegressionType)
  #@eval Main demeanbeta = $(A[:,1] |> Vector)

  absμ = factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end


(Θ::VolumePartsIter)(Xv::AbstractXYIter{TM,TV,T},
    ::Val{:leveldemean}) where {TM<:CuMatrix,TV<:CuVector,T} = error(
    "Do not use CUDA for the demeaning version")


function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterLevel, RegressionType::Val
  ) where {TM<:CuMatrix,TV<:CuVector,T}

  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = initializeÃandfactor(Θ,Xv)

  #run the regression for A₁
  try
    A[:,1] .= reg(factored, Xv.v, RegressionType)
    #@info "culevel A[:,1]: $(A[:,1])"
    #@info "culevel factored[20,:]: $(factored[20,:])"
    #@info "culevel A[:,20]: $(A[:,20])"
  catch err #NOTE: below may be broken, could try to fix if needed
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage) > 1000 ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || @warn("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])
      attempting correction")
    factored = (x->min(x,LARGE_VAL)).(factored)

    A[:,1] .= reg(factored, Xv.v, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallback val=$(A[:,1])"
  end

  #println("A[:,1]: $(A[:,1])")

  absμ = factored * (x->Float64(abs(x))).(A[:,1])


  return absμ
end

#NOTE: this is only for testing purposes, but is it even needed for that???
#handles the intercept case
function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterLevel,
  ::Union{Val{:interceptzygote}, Val{:intercept}},

  ) where {TM<:CuMatrix,TV<:CuVector,T}

  @unpack G,A = Θ
  factored::TM = hcat(fill!(similar(Xv.v), T(1.0)), initializeÃandfactor(Θ,Xv))

  #run the regression for A₁
  A₁ = reg(factored, Xv.v, Val{:choleskyzygote}())
  absμ = factored * A₁


  return absμ
end

#not sure why I need to write this again...but it throws a method error otherwise
function (Θ::AbstractVolumePartsIter{TM,TV,T})(Xv::XYIterLevel,
  ::Union{Val{:interceptzygote},Val{:intercept}},
  ) where {TM<:Matrix,TV<:Vector,T}

  @unpack G,A = Θ
  factored::TM = hcat(fill!(similar(Xv.v), T(1.0)), initializeÃandfactor(Θ,Xv))
  #run the regression for A₁
  A₁ = reg(factored, Xv.v, Val{:choleskyzygote}())
  absμ = factored * A₁

  return absμ
end



function (Θ::VolumePartsIter)(Xv::XYIterLevel{TM,TV,T}, ::Val{:intercept}
  ) where {TM<:CuMatrix,TV<:CuVector,T}

  @unpack G,A = Θ

  #initialize A₁.==1 and all other As to the amount implied by G
  #print("prod(Θ.G): f0=$(prod(Θ.G|>Matrix)) x64=$(prod(Float64.(Θ.G |> Matrix) )) |")
  prevA₁ = A[:,1]
  factored::TM = hcat(Xv.tv1s, initializeÃandfactor(Θ,Xv))
  local b::TV
  #expectations = mean(factored, dims=1)
  #v̄ = mean(Xv.v)

  #run the regression for A₁
  try
    b = reg(factored, Xv.v, Val{:cholesky}())
    #@info "culevel A[:,1]: $(A[:,1])"
    #@info "culevel factored[20,:]: $(factored[20,:])"
    #@info "culevel A[:,20]: $(A[:,20])"
  catch err
    (err == InterruptException()) && throw(err)
    errmessage = "$err"
    errmessage = length(errmessage)  > 1000 ? errmessage[1:1000] : errmessage
    all((isfinite.(factored))) || @warn("Factored includes nonfinite values!
      prevA₁: $(prevA₁)
      ****************\nG: $G
      ****************\nfactored[1:1000,:]: $(factored[1:1000,:])
      attempting correction")
    factored .= (x->min(x,LARGE_VAL)).(factored)
    #expectations .= mean(factored, dims=1)

    b = reg(factored, Xv.v, Val{:fallbackreg}())
    @warn "regression method $RegressionType failed due to $errmessage. Using fallback val=$(A[:,1])"
  end

  #intercept = v̄ - dot(vec(expectations), A[:,1])

  #println("A[:,1]: $(A[:,1])")
  A[:,1] .= b[2:end]
  absμ = factored * b


  return absμ
end



function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T}, Xv::XYIterLevel,
  RegressionType::Val) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), Nexpanded=size(Xv.ws,1))

  #special version for zygote
  #expand a matrix for A
  (RegressionType === Val{:intercept}()) && error("Don't use the pre-compute method with an intercept regression")

  Ã = hcat(ones(T, dims.K, 1), cumprodbyrow(G))
  Ãexpanded = view(Ã, Θ.expand.index.Ã...)
  LÃexpanded = view(Ã, Θ.expand.index.LÃ...)
  factored = genabsμₖ.(Ãexpanded, LÃexpanded, Xv.ws', Xv.Rtws', Xv.RLws')' #0.4ms


  absA₁ = (x->Float64(abs(x))).(reg(factored, Xv.v, RegressionType))#0.5ms
  #println("CPU absAₜ: $absAₜ")

  #sum across signals
  absμ = factored * absA₁
  #absμ = sum(Aexpanded, dims=1) |> vec
  return absμ
end

#used in testing
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T}, Xv::XYIterDemean,
  RegressionType::Val) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), Nexpanded=size(Xv.ws,1))

  #special version for zygote
  #expand a matrix for A
  (RegressionType === Val{:intercept}()) && error("Don't use the pre-compute method with an intercept regression")

  Ã = hcat(ones(T, dims.K, 1), cumprodbyrow(G))
  Ãexpanded = view(Ã, Θ.expand.index.Ã...)
  LÃexpanded = view(Ã, Θ.expand.index.LÃ...)
  factored = Xv.M * genabsμₖ.(Ãexpanded', LÃexpanded', Xv.ws, Xv.Rtws, Xv.RLws) #0.4ms


  absA₁ = (x->Float64(abs(x))).(reg(factored, Xv.ṽ, Val{:choleskyzygote}()))#0.5ms
  #@info ("absA₁ (:zygoteleveldemean): $(absA₁)")

  #sum across signals
  absμ = factored * absA₁
  #absμ = sum(Aexpanded, dims=1) |> vec
  return absμ
end



function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T}, Xv::XYIterLevel,
  RegressionType::Val) where {TM<:CuMatrix,TV<:CuVector,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), Nexpanded=size(Xv.ws,1))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G)
  Ã = hcat(someones, prodmatG)

  Ãexpanded = CuMatrix{T}(Ã[Θ.expand.index.Ã...])
  LÃexpanded = CuMatrix{T}(Ã[Θ.expand.index.LÃ...])


  factored = genabsμₖ.(Ãexpanded, LÃexpanded, Xv.ws', Xv.Rtws', Xv.RLws')' #0.4ms

  #run the regression (use a spacial non-inplace version)

  #println(sum(cholesky(factored'*factored)\(factored'*Xv.v)))
  #println("Cu absAₜT: $absAₜT")

  absA₁ = reg(factored, Xv.v, RegressionType)


  absμ = factored * absA₁

  #=return absμ=#

  return absμ
end

#verion for regressions with intercept
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T}, Xv::XYIterLevel,
  RegressionType::Val{:interceptzygote}) where {TM<:CuMatrix,TV<:CuVector,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), Nexpanded=size(Xv.ws,1))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G)
  Ã = hcat(someones, prodmatG)

  Ãexpanded = CuMatrix{T}(Ã[Θ.expand.index.Ã...])
  LÃexpanded = CuMatrix{T}(Ã[Θ.expand.index.LÃ...])
  factored = CUDA.hcat(
    Xv.tv1s, genabsμₖ.(Ãexpanded, LÃexpanded, Xv.ws', Xv.Rtws', Xv.RLws')') #0.4ms

  absA₁ =reg(factored, Xv.v, Val{:choleskyzygote}())
  absμ = factored * absA₁
  #absμ = sum(Aexpanded, dims=1) |> vec
  return absμ
end
#cpu version of above
function (Θ::AbstractVolumePartsIter{TM,TV,T})(G::Matrix{T}, Xv::XYIterLevel,
  RegressionType::Val{:interceptzygote}) where {TM<:Matrix,TV<:Vector,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), Nexpanded=size(Xv.ws,1))

  someones = ones(T, dims.K, 1)
  prodmatG = cumprodbyrow(G)
  Ã = hcat(someones, prodmatG)

  Ãexpanded = view(Ã, Θ.expand.index.Ã...)
  LÃexpanded = view(Ã, Θ.expand.index.LÃ...)
  factored = hcat(Xv.tv1s, genabsμₖ.(Ãexpanded, LÃexpanded, Xv.ws', Xv.Rtws', Xv.RLws')') #0.4ms

  absA₁ =reg(factored, Xv.v, Val{:choleskyzygote}())
  absμ = factored * absA₁
  #absμ = sum(Aexpanded, dims=1) |> vec
  return absμ
end



function lossforzygote(G::Matrix{T}, Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel;
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}()) where {TM,TV,T}

  #ε = factored * Aₜ .- Xv.ṽ
  Σε2 = sum((Θ(G, Xv, RegressionType) .- Xv.v).^2) #0.3ms

  return Σε2

end

function lossforzygote(G::Matrix{T}, Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean;
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}()) where {TM,TV,T}

  #ε = factored * Aₜ .- Xv.ṽ
  Σε2 = sum((Θ(G, Xv, RegressionType) .- Xv.ṽ).^2) #0.3ms

  return Σε2
end


#computs the predicted values for each strategy k and time t and demeans the result
#via the projection matrix
function projectvolumefromA(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean) where {TM,TV,T}

  @unpack A, LA = Θ.expand
  @unpack ws, Rtws, RLws, M = Xv
  #first map each coefficient to its appropriate

  absμₖ = genabsμₖ.(A', LA', ws, Rtws, RLws)

  multiply!(M, absμₖ)
  return absμₖ
end

#same as the above without demeaning
function projectvolumefromA(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel) where {TM,TV,T}

  @unpack A, LA = Θ.expand
  @unpack ws, Rtws, RLws= Xv
  #first map each coefficient to its appropriate

  absμₖ =genabsμₖ.(A', LA', ws, Rtws, RLws)

  return absμₖ
end

function projectvolumefromA(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel
    ) where {TM<:CuMatrix, TV<:CuVector, T}

  @unpack A, LA = Θ.expand
  @unpack ws, Rtws, RLws = Xv
  #first map each coefficient to its appropriate

  dims = (K=size(Θ.A,1), T = size(Θ.A,2))

  genabsμₖ.(A', LA', ws, Rtws, RLws)
  #absμₖ .= cudagenabsμₖ.(A', LA', ws, RLws)

  return absμₖ
end


##############Iter loss methods
#again, mostly legacy

function iterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel, RegressionType::Val;
  ) where {TM<:CuMatrix, TV<:CuVector, T}

  Δ = Θ(Xv, RegressionType) .- Xv.v
  λ = CUDA.sum(Δ .* Δ)# / Ncleanpanel

  return λ
end


iterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel, RegressionType::Val;
  ) where {TM, TV, T} = sum(
    (Θ(Xv, RegressionType) .- Xv.v).^2)


iterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean, RegressionType::Val;
  ) where {TM, TV, T} = sum(
    (Θ(Xv, RegressionType) .- Xv.ṽ).^2)


function currentloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean, absμₖ::TM) where {TM, TV, T}

  someones::Vector{T} = ones(T, size(Θ.A,1))
  return sum((absμₖ*someones .- Xv.ṽ).^2)
end

function currentloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean,
   absμₖ::TM) where {TM<:CuMatrix, TV<:CuVector, T}

  someones::CuVector{T} = CUDA.ones(T, size(Θ.A,1))
  return sum((absμₖ*someones .- Xv.ṽ).^2)
end

function currentloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel,
   absμₖ::TM) where {TM<:CuMatrix, TV<:CuVector, T}

  someones::CuVector{T} = CUDA.ones(T, size(Θ.A,1))
  return sum((absμₖ*someones .- Xv.v).^2)
end


function testiterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel;
    RegressionType::Val=Val(PARAM[:iterregressiontype])) where {TM, TV, T}


  vhattest = iterloss(deepcopy(Θ), Xv, RegressionType)

  vhatver = sum((Θ(Xv, RegressionType) .- Xv.v) .^ 2)

  @assert vhatver ≈ vhattest
end

function testiterloss(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean;
  RegressionType::Val=Val(:choleskyzygote)) where {TM, TV, T}

  vhattest = iterloss(deepcopy(Θ), Xv, RegressionType)

  vhatver = sum((Θ(Xv, RegressionType) .- Xv.ṽ) .^ 2)

  @assert vhatver ≈ vhattest
end
