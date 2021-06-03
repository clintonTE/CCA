

##analytical gradient- now used mostly for test (too slow for regular use)
DgXki(ti, Xik, wLAsignXik, t, Gkt) = (ti>t)*(Xik/Gkt) + (t==ti)*wLAsignXik


function dfactor!(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T}, DgX; dinfo...) where {TM,TV,T}
  @unpack Gkt, k, t, wLAsignX, dtsM1, dims, X= dinfo

  DgX[:, k] .= DgXki.(dtsM1, X[:,k], wLAsignX[:,k], t, Gkt)

  return nothing
end

#try to compute the analytical gradient
function ∇lossag(G::Matrix{T}, Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlLevel;
  RegressionType::Val = Val(PARAM[:iterregressiontype])) where {TM<:CuMatrix,TV<:CuVector,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2),KW=size(Θ.A,1) + size(Xv.W,2),
    Nexpanded=size(Xv.ws,1), expandedinds = collect(1:size(Xv.ws,1)) |> CuVector{Int})
  @unpack v = Xv


  #pre-compute some common values
  X::TM = hcat(initializeÃandfactor(Θ, Xv), Xv.W)
  cuI = TM(I, dims.KW, dims.KW) #need an identity matrix
  XtXinv = cholesky!(X'*X)\cuI
  XXtXinv = X*XtXinv
  XtXinvXtv = XtXinv * (X'*v)
  εt = (X*XtXinvXtv-Xv.v)'


  wLAsignX = signabsμₖ.(Θ.expand.Ã', Θ.expand.LÃ', Xv.ws, Xv.Rtws, Xv.RLws) .* Xv.ws .* abs.(Θ.expand.LÃ')

  #preallocate this carefully high performance
  DgX::TM = CUDA.zeros(T, dims.Nexpanded, dims.KW)

  DG = similar(G)
  Pv = similar(Xv.v)
  dtsM1 = Θ.tsM1 |> CuVector{Int}
  for k ∈ 1:dims.K
    for t ∈ 1:(dims.T-1)


      dfactor!(Θ, Xv, DgX, Gkt=G[k,t]; X, wLAsignX, k, t, dims, dtsM1)
      copyto!(Pv, DgX*XtXinvXtv -
        XXtXinv*(DgX'*X+X'*DgX)*XtXinvXtv +
        XXtXinv*(DgX'*v))
      DG[k,t] = 2.0 * (εt * Pv)

    end
    DgX[:,k] .= T(0.0)
  end

  return DG
end


#try to compute the analytical gradient
function ∇lossag(G::Matrix{T}, Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlLevel
  ) where {TM,TV,T}

  dims = (K=size(Θ.A,1), T = size(Θ.A,2), KW=size(Θ.A,1) + size(Xv.W,2),
    Nexpanded=size(Xv.ws,1), expandedinds = collect(1:size(Xv.ws,1)) |> CuVector{Int})
  @unpack v = Xv


  #pre-compute some common values
  #note the reversed intercept location
  X::TM = hcat(initializeÃandfactor(Θ, Xv), Xv.W)
  cuI = TM(I, dims.KW, dims.KW) #need an identity matrix
  XtXinv = cholesky!(X'*X)\cuI
  XXtXinv = X*XtXinv
  XtXinvXtv = XtXinv * (X'*v)
  εt = (X*XtXinvXtv-Xv.v)'


  wLAsignX = signabsμₖ.(Θ.expand.Ã', Θ.expand.LÃ', Xv.ws, Xv.RLws) .* Xv.ws .* abs.(Θ.expand.LÃ')

  #preallocate this carefully high performance
  DgXs::Vector{TM} = [zeros(T, dims.Nexpanded, dims.KW) for i ∈ 1:Threads.nthreads()]
  Pvs = [similar(Xv.v) for i ∈ 1:Threads.nthreads()]

  DG = similar(G)

  dtsM1 = Θ.tsM1
  for k ∈ 1:dims.K
    Threads.@threads for t ∈ 1:(dims.T-1)
      DgX = DgXs[Threads.threadid()]
      Pv = Pvs[Threads.threadid()]

      dfactor!(Θ, Xv, DgX, Gkt=G[k,t]; X, wLAsignX, k, t, dims, dtsM1)
      copyto!(Pv, DgX*XtXinvXtv -
        XXtXinv*(DgX'*X+X'*DgX)*XtXinvXtv +
        XXtXinv*(DgX'*v))
      DG[k,t] = 2.0 * (εt * Pv)

    end
    for DgX ∈ DgXs
      DgX[:,k] .= T(0.0)
    end
  end

  return DG
end

function testag(Θ, Xv::AbstractXYIter{TM,TV,T},
  PredictionType::Val=Val{:level}(),
  GrowthType::Val = Val(PARAM[:itergrowthtype])) where {TM,TV,T}
  CUDA.allowscalar(false)

  @info "Beginning X_g accuracy test. PredictionType=$PredictionType, GrowthType=$GrowthType"

  dims = (K=size(Θ.A,1), T = size(Θ.A,2),
    Nexpanded=size(Xv.ws,1), expandedinds = collect(1:size(Xv.ws,1)) |> CuVector{Int})

  Θtest = deepcopy(Θ)
  Θtest.G .= exp.(Θtest.G)

  local Ãffactor, Ãag
  #first check accuracy of d/dG(X)
  function ffactor(G)
    someones = ones(T, dims.K, 1)
    prodmatG = cumprodbyrow(G)
    Ã = hcat(someones, prodmatG)
    Ãffactor = deepcopy(Ã)
    Ãexpanded = CuMatrix{T}(Ã[Θtest.expand.index.Ã...])
    LÃexpanded = CuMatrix{T}(Ã[Θtest.expand.index.LÃ...])
    factored = genabsμₖ.(Ãexpanded, LÃexpanded, Xv.ws', Xv.Rtws', Xv.RLws')' #0.4ms

    return sum(factored)
  end

  function ∇Xag(G)
    @assert (G|>Matrix) == (Θtest.G|>Matrix)
    X = initializeÃandfactor(Θtest, Xv) |> deepcopy
    DgX = zeros(T, dims.Nexpanded, dims.K) |> TM
    wLAsignX = signabsμₖ.(Θtest.expand.Ã', Θtest.expand.LÃ', Xv.ws, Xv.Rtws, Xv.RLws) .* Xv.ws .* abs.(Θtest.expand.LÃ')
    DG = similar(G)
    Ãag = deepcopy(Θtest.Ã)
    dtsM1 = Θ.tsM1 |> CuVector{Int}
    for k ∈ 1:dims.K
      for t ∈ 1:(dims.T-1)
        dfactor!(Θtest, Xv, DgX, Gkt=G[k,t]; wLAsignX, k, t, dims, X, dtsM1)
        DG[k,t] = sum(DgX)
      end
      DgX[:,k] .= 0.0
    end

    return DG
  end

  ∇X = gradient(ffactor, Θtest.G |> Matrix)
  ∇Xans = ∇X[1]
  ∇agXans = ∇Xag(Θtest.G |> Matrix)
  @assert (Ãag |> Matrix) ≈ (Ãffactor |> Matrix) #check product matrix
  if !(∇agXans ≈ ∇Xans)
    error("!(∇Xag(Gtest) ≈ ∇X(Gtest))!!
      ∇Xag: $(∇agXans)
      ∇X: $(∇Xans)")
  end
  @info "ag X accuracy test passed"

  #acquire gradient functions
  ∇zygote = ∇lossforzygote(Θ, Xv, GrowthType, GradType=Val{:zygote}())
  ∇ag = ∇lossforzygote(Θ, Xv, GrowthType, GradType=Val{:ag}())

  #create the inputs and space for the outputs
  g = Θ.g |> TV |> deepcopy
  ∇gzygote = g |> TV |> deepcopy |> Vector
  ∇gag = g |> TV |> deepcopy |> Vector

  #check the accuracy
  @info "beginning ag accuracy test"

  ∇zygote(∇gzygote, g)
  ∇ag(∇gag, g)
  if !(∇gag ≈ ∇gzygote)
    error("!(∇gag ≈ ∇gzygote)!!
      ∇gag: $∇gag
      ∇gzygote: $∇gzygote")
  end
  @info "ag accuracy test passed"

  #check performance
  if PARAM[:benchgrad]
    print("Zygote performance: ")
    @btime CUDA.@sync $∇zygote($∇gzygote, $g)
    print("Analytic performance: ")
    @btime CUDA.@sync $∇ag($∇gag, $g)
  end

  return nothing
end


##
#this used to be part of VolumePartsIter, now it is here for testing purposes
function volumepartsxsection(Θ::VolumePartsIter)
  @unpack A,G,ts = Θ

  dims = (T=size(A,2), K=size(A,1))
  xsections = Vector{NamedTuple}(undef, dims.T-1)
  for t ∈ 2:dims.T
    #count the number of stocks in each window (t and t-1)
    Nt = sum(ts.==t)

    #repeat t for each stock
    Ainds = repeat([t],Nt)
    LA_inds = repeat([t-1], Nt)
    xsections[t-1] = (A = view(A, :, Ainds),
      LA = view(A, :, LA_inds),
      G = view(G, :, LA_inds))
  end
  xsection::NamedTuple = (
    A = (x->x.A).(xsections),
    LA = (x->x.LA).(xsections),
    G = (x->x.G).(xsections),
  )

  return xsection
end

###These are legacy code routines moved from VolumePartsIter, still sometimes used in testing
#creates a factored object and associated cross-secitons for building up the
#regression rhs
genfactoredLA(wsᵢ, Rtwsᵢ, RLwsᵢ,Gᵢt)= abs(Rtwsᵢ  - RLwsᵢ) .+ abs(Rtwsᵢ  - Gᵢt * wsᵢ)
#factors the |G*LA*w-LA*RLw| expression(pulls out the LA)
function factorLA(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  ) where {TM,TV,T}
#use the pre-allocation
  factored =  genfactoredLA.(Xv.ws, Xv.Rtws, Xv.RLws, Θ.expand.G')

  return factored
end

function unfactorLA(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  factoredLA::TM) where {TM, TV, T}
  absμₖ = factoredLA .* (abs).(Θ.expand.LA)'

  return absμₖ
end

#factors the |G*LA*w-LA*RLw| expression(pulls out the LA) for each crosssection
function xfactorLA(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},xsection::NamedTuple,
  ) where {TM,TV,T}

  @assert !((TM <: CuArray) || TV <: CuArray)
  Θxsection = volumepartsxsection(Θ)
  #use the pre-allocation
  factoredLA::Vector{TM} = deepcopy(xsection.ws)
  Threads.@threads for i ∈ 1:(size(Θ.A,2)-1)
    factoredLA[i] .= genfactoredLA.(xsection.ws[i],xsection.Rtws[i],xsection.RLws[i], Θxsection.G[i]')
  end

  return factoredLA
end

#gpu version of above
#primarily for testing since the xsection is no longer pre-allocated
function xfactorLA(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},xsection::NamedTuple,
  ) where {TM<:CuMatrix,TV<:CuVector,T}

  Θxsection = volumepartsxsection(Θ)

  #form the crosssection
  factoredLA::Vector{TM} = CUDA.deepcopy(xsection.ws)
  for i ∈ 1:(size(Θ.A,2)-1)
    factoredLA[i] .= genfactoredLA.(xsection.ws[i],xsection.Rtws[i],xsection.RLws[i], Θxsection.G[i]')
  end

  return factoredLA
end


function unfactorA(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  factoredA::TM) where {TM, TV, T}
  absμₖ = factoredA .* abs.(Θ.expand.A)'

  return absμₖ
end


#copies relevant data into global variables for testing purposes
function testglobal(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterLevel,
    ) where {TM,TV,T}
  global thetag, Xvg = deepcopy(Θ), deepcopy(Xv)
  global factoredg = factorLA(Θ, Xv) |> Matrix{T}
  global df = DataFrame(factoredg)
  df.ts = (i->"t$i").(thetag.ts)
  df.v = Xvg.v |> Vector{T}

  global v = df.v |> TV
  global X = factoredg |> TM
  return nothing
end

#copies relevant data into global variables for testing purposes
function testglobal(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean,
    ) where {TM,TV,T}
  global thetag, Xvg = deepcopy(Θ), deepcopy(Xv)
  global factoredg = factorLA(Θ, Xv) |> Matrix{T}
  global df = DataFrame(factoredg)
  df.ts = (i->"t$i").(thetag.ts)
  df.v = Xvg.v |> Vector{T}
  #f = @formula(0 ~ x1 + x2 + ts + 0)
  global MM = hcat(factoredg, Xv.M.F|> collect |> Matrix{T})
  global F = MM[:,(1+size(Xv.ws,2)):end] |> TM
  global v = df.v |> TV
  global X = factoredg |> TM
  return nothing
end

function testfixed(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
    PredictionType::Val,
    RegressionType::Val=Val{PARAM[:iterregressiontype]}(),
    benchmeasure::Bool = true) where {TM,TV,T}

  CUDA.allowscalar(false)
  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  factoredg = factorLA(Θ, Xv) |> Matrix{T}

  #F = Xv.M.F |> TM
  v = Xv.v |> TV
  X = factoredg |> TM

  @assert Xv.M*v ≈ Xv.ṽ

  #now check the demeaned predictions by strategy
  absμₖdemean::Matrix{Float64} = similar(Xv.ws)
  D = (ws=Xv.ws, Rtws=Xv.Rtws, RLws=Xv.RLws)
  @unpack ws, Rtws, RLws = dataxsection((ws=Xv.ws, Rtws=Xv.Rtws, RLws=Xv.RLws),Θ.ts,TM,TV,T)
  @unpack xabsμₖ = dataxsection((xabsμₖ=absμₖdemean,),Θ.ts,Matrix{Float64},Vector{Float64},Float64)

  Θxsection=volumepartsxsection(Θ)
  for t ∈ 1:(dims.T-1)
    xabsμₖ[t] .= genabsμₖ.(Θxsection.A[t]', Θxsection.LA[t]', ws[t], Rtws[t], RLws[t])
    xabsμₖ[t] .-= mean(xabsμₖ[t], dims=1) #demeaning step
  end
  absμₖdemean .= reduce(vcat, xabsμₖ)
  absμₖproject::Matrix{Float64} = projectstrategyvolume(Θ, Xv)

  absμₖrmserror = sqrt(mean((absμₖdemean .- absμₖproject).^2))
  if T.(absμₖdemean) ≈ T.(absμₖproject)
    @info "testfixed absμₖ test passed with rmserror=$absμₖrmserror"
  else
    println("absμₖ test failed!!! absμₖdemean[1:10,:]")
    printmln(absμₖdemean[1:10,:])
    println("absμₖproject[1:10,:]")
    printmln(absμₖproject[1:10,:])
    error("absμₖ test failed!!! rmserror=$absμₖrmserror")
  end

  vhatdemean = sum(absμₖdemean, dims=2) |> vec
  vhatproject = sum(absμₖproject, dims=2) |> vec

  vrmserror = sqrt(maximum((vhatdemean .- vhatproject).^2))
  if T.(vhatdemean)  ≈ T.(vhatproject)
    @info "testfixed v test passed with rmserror=$vrmserror"
  else
    error("objective test failed!!! rmserror=$vrmserror")
  end


end

#compares the zygote loss/prediction against the standard loss/prediction
function testpredictionzygote(Θ::AbstractVolumePartsIter, Xv::AbstractXYIter,
  ::Type{T}=PARAM[:itertype];
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}()) where T


  #now get the zygote prediction
  G::Matrix = Θ.G |> Matrix
  vhatzygotever = Θ(Xv, RegressionType)
  vhatzygote = Θ(G, Xv, RegressionType)  |> Vector
  (T.(vhatzygotever |> Vector) ≈ T.(vhatzygote)) || error("vhatzygote≠vhattest
    vhattest[1:10]: $(vhatzygote[1:10])
    vhatzygote[1:10]: $(vhatzygotever[1:10])")

  λzygote = lossforzygote(G, Θ, Xv)

  Δver = vhatzygotever .- Xv.ṽ
  λzygotever = sum(Δver .* Δver)

  #!(λzygote ≈ λzygotever) && @warn "λzygote !≈ λzygotever. If zygote fixes cumprod,
  #  replace the zygote cumprod and we can get much better with this."
  (λzygote ≈ λzygotever) || error("λzygote ≠ λzygotever
    λzygote: $λzygote
    λzygotever: $λzygotever
    vhattest[1:10]: $(vhatzygote[1:10])
    vhatzygote[1:10]: $(vhatzygotever[1:10])")
end

#compares the zygote loss/prediction against the standard loss/prediction
function testcudazygote(Θ::VolumePartsIter, Xv::XYIterControl,
  Θcpu::VolumePartsIter, Xvcpu::XYIterControl,
  ::Type{T}=PARAM[:itertype];
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}(),
  GrowthType::Val = Val{PARAM[:itergrowthtype]}(),
  benchmeasure::Bool = PARAM[:benchitermeasure]) where T

  ###first compare standard predictions

  vhatcpu = Θcpu(Xvcpu, RegressionType) |> Vector
  vhat = Θ(Xv, RegressionType)  |> Vector
  (T.(vhatcpu) ≈ T.(vhat)) || error("vhat≠vhatcpu (cuda zygote test)
    vhat[1:10]: $(vhat[1:10])
    vhatcpu[1:10]: $(vhatcpu[1:10])")

  λ = iterloss(Θ, Xv, RegressionType)
  λcpu = iterloss(Θcpu, Xvcpu, RegressionType)

  @assert λ ≈ λcpu

  ###compare cpu and cuda zygote prediction
  vhatzygotecpu = Θcpu(Θcpu.G, Xvcpu, RegressionType) |> Vector
  vhatzygote = Θ(Θcpu.G, Xv, RegressionType)  |> Vector
  (T.(vhatzygotecpu |> Vector) ≈ T.(vhatzygote)) || error("vhatzygote≠vhattest (cuda zygote test)
    vhatzygote[1:10]: $(vhatzygote[1:10])
    vhatzygotecpu[1:10]: $(vhatzygotecpu[1:10])")

  λzygote = lossforzygote(Θcpu.G, Θ, Xv,)
  λzygotecpu = lossforzygote(Θcpu.G, Θcpu, Xvcpu,)

  Δcpu = vhatzygotecpu .- Xvcpu.ṽ
  λzygotecpu = sum(Δcpu .* Δcpu)
  @assert (λzygote ≈ λzygotecpu) || error("λzygote ≈/ λzygotecpu!!!
    λzygote: $λzygote ; λzygotecpu: $λzygotecpu")

  ###now compare the gradients
  x = deepcopy(Θcpu.g |> Vector)

  if typeof(GrowthType) <: Val{:log}
    x .= log.(abs.(x))
    Θcpu.g .= abs.(Θcpu.g)
    Θ.g .= abs.(Θ.g)
  end

  outcpu = similar(x)
  out = similar(x)

  ∇Θcpu = ∇lossforzygote(Θcpu, Xvcpu, GrowthType)
  ∇Θ = ∇lossforzygote(Θ, Xv, GrowthType)

  ∇Θcpu(outcpu, x)
  ∇Θ(out, x)

  @assert Θcpu.g ≈ (Θ.g |> Vector)
  @assert Θcpu.G[Θcpu.expand.index.G...] ≈ (Θ.G[Θ.expand.index.G...]  |> Matrix)

  (out ≈ outcpu) || error("(out ≈/ outcpu)!!!  out: $out outcpu: $outcpu")


  if benchmeasure
    print("Loss performance cpu:")
    @btime iterloss($Θ, $Xv, $RegressionType)
    print("Loss performance gpu:")
    @btime iterloss($Θcpu, $Xvcpu, $RegressionType)


    print("Gradient performance cpu:")
    @btime $∇Θcpu($outcpu, $x)
    print("Gradient performance gpu:")
    @btime $∇Θ($out, $x)
  end
end

#compares the zygote loss/prediction against the standard loss/prediction
function testpredictionzygote(Θ::VolumePartsIter, Xv::XYIterDemean,
  ::Type{T}=PARAM[:itertype];
  RegressionType::Val=Val{PARAM[:iterregressiontypezygote]}()) where T

  #generate the standard pre-alloc prediction
  factored = initializeÃandfactor(Θ, Xv)
  Θ.A[:, 1] .= reg(Xv.M * factored, Xv.ṽ,Val{:glm}())
  #@info "Θ.A[:, 1] (testpredictionzygote): $(Θ.A[:, 1])"

  #now get the zygote prediction
  vhatzygotever =  factored * (abs).(Θ.A[:,1])
  vhatzygote = Θ(Θ.G, Xv, Val{:choleskyzygote}())  |> Vector
  (vhatzygotever≈vhatzygote) || error("vhatzygote≠vhattest
    vhattest[1:10]: $(vhatzygote[1:10])
    vhatzygote[1:10]: $(vhatzygotever[1:10])")

  @assert lossforzygote(Θ.G, Θ, Xv) ≈ sum((Xv.ṽ .- (Float64).(vhatzygotever)).^2)
end


#some validation tests for xsection integrity
function validatexsection(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
    p::NamedTuple, xsection::NamedTuple) where {TM,TV,T}
  @unpack dims, tidx, cleanpanel = p

  #@unpack xsection = Xv
  Θxsection=volumepartsxsection(Θ)

  for i ∈ 1:(dims.T-1)
    xsectionsize = size(Θxsection.A[i]')
    (xsectionsize == size(Θxsection.LA[i]')) || error("unpected size of Θxsection.LA[i]'")
    (xsectionsize == size(Θxsection.G[i]')) || error("unpected length of Θxsection.G[i]'")
    (xsectionsize == size(xsection.ws[i]) == size(xsection.RLws[i])) || error(
      "Unexpected array sizing for xsection $i: xsectionsize=$(xsectionsize)
      xsection.ws[i]=$(size(xsection.ws[i]))
      xsection.RLws[i]=$(size(xsection.RLws[i]))")
    @assert xsectionsize[1] == length(xsection.ṽ[i])
    xsectiondate = cleanpanel.date[Θ.ts .== i+1]
    @assert sum(xsectiondate .== xsectiondate[1]) == length(xsectiondate) #all the same date
    @assert i+1 == tidx[xsectiondate[1]]
  end

  return nothing
end

validatexsection(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  p::NamedTuple) where {TM,TV,T} = validatexsection(
    Θ, Xv, p, dataxsection((ws=Xv.ws, RLws=Xv.RLws, ṽ=Xv.ṽ, v=Xv.v), Θ.ts, TM, TV, T))


function testG(Θ::VolumePartsIter, Xv::AbstractXYIter{TM,TV,T},
  xsection::NamedTuple) where {TM,TV,T}

  #@unpack xsection = Xv
  Θxsection=volumepartsxsection(Θ)
  #first, generate predicatable test data and run the update fuinction
  dims = (K=size(Θ.A,1), T=size(Θ.A,2))
  testrowA::TM = exp.([isodd(i) ? log(T(1.1)) : log(T(0.9)) for i ∈ 1:dims.T])' |> TM
  testrowA .*= (i->T(-1)^T(i+1)).(1:dims.T)' |> TM
  testrowA .= cumprod(testrowA, dims=2)#, dims=2)
  testrowG::TM = testrowA[:,2:end] ./ testrowA[:,1:(end-1)]

  #exp.(collect(1:dims.T) .* (i->isodd(i)*log(T(1.1)) + iseven(i))).(1:dims.T)'
  # testrowG
  Θ.G .= testrowG
  Θ.A .= hcat(ones(T, dims.K,1) |> TM, cumprodbyrow(Θ.G))
  Θ.G .= Θ.A[:, 2:end] ./ Θ.A[:, 1:(end-1)]

  #test updateAfromGAₜ!
  Θcpy=deepcopy(Θ)
  Θcpy.A[:,25] .= T(25)
  updateAfromGAₜ!(Θcpy, Xv,25)
  ΠG = hcat(ones(T, dims.K,1) |> TM, cumprodbyrow(Θ.G))
  @assert Θcpy.A ≈ ΠG .* T(25) ./  ΠG[:, 25]

  #create test factored flows
  factoredLA = deepcopy(factorLA(Θ, Xv))
  xfactoredLA = deepcopy(xfactorLA(Θ, Xv, xsection))

  @assert factoredLA == reduce(vcat, xfactoredLA)

  #make sure the values of G check out
  testG::TM = reduce(vcat, [testrowG for i ∈1:dims.K])
  testA =hcat(ones(T,dims.K,1)|>TM, cumprodbyrow(testG))
  testG .= testA[:,2:end] ./ testA[:,1:(end-1)]

  (Θ.G |> Matrix! ≈ testG |> Matrix!) || error(
    "Differences found between G and test
    G[:,1:5]: $(Θ.G[:, 1:5]) expected $(testG[:, 1:5])")

  #verify the overall data structure is as expected
  expected = Matrix(abs.(Θ.expand.LA' .* Xv.Rtws .- Θ.expand.LA' .* Xv.RLws) .+
    abs.(Θ.expand.LA' .* Xv.Rtws .- Θ.expand.A' .* Xv.ws))
  @assert expected ≈ Matrix!(factoredLA .* abs.(Θ.expand.LA)')
  @assert expected ≈ Matrix!(unfactorLA(Θ,Xv,factoredLA))
  @assert expected ≈ Matrix!(genabsμₖ.(Θ.expand.A', Θ.expand.LA', Xv.ws, Xv.Rtws, Xv.RLws))

  #verify each crossection is as expected
  expectedprojectt = Vector{Matrix{Float64}}()
  for i ∈ 1:(dims.T-1)
    Nᵢ::Int = size(xsection.ws[i],1)
    expected = Matrix!(abs.(Θxsection.LA[i]' .* xsection.Rtws[i] .- Θxsection.LA[i]' .* xsection.RLws[i]
      ) .+ abs.(Θxsection.LA[i]' .* xsection.Rtws[i] .- Θxsection.A[i]' .* xsection.ws[i]))
    xpredLA = Matrix!(xfactoredLA[i]) .* Matrix!(abs.(Θxsection.LA[i]))'

    @assert (expected ≈ xpredLA) "
      size expected:  $(size(expected))
      size Matrix(xpredLA):  $(size(xpredLA))
      expected[:,1:5]: $(expected[1:5,:])
      predLA[:,1:5]: $(xpredLA[1:5,:])"

    #the below is a simple component test of the projectt methods and
    #genabsμk
    if Xv.xsection.xM ≡ nothing
      expected = Matrix!(abs.(Θxsection.LA[i]' .* Xv.xsection.xRtws[i] .- Θxsection.LA[i]' .* Xv.xsection.xRLws[i]
        ) .+ abs.(Θxsection.LA[i]' .* Xv.xsection.xRtws[i] .- Θxsection.A[i]' .* Xv.xsection.xws[i]))
      @assert expected ≈ Matrix!(projectt(Θxsection.A[i], Θxsection.LA[i],
        Xv.xsection.xws[i], Xv.xsection.xRtws[i], Xv.xsection.xRLws[i], nothing))
      push!(expectedprojectt, expected)
    else
      expected = Matrix!(Xv.xsection.xM[i] * (abs.(Θxsection.LA[i]' .* Xv.xsection.xRtws[i] .- Θxsection.LA[i]' .* Xv.xsection.xRLws[i]
        ) .+ abs.(Θxsection.LA[i]' .* Xv.xsection.xRtws[i] .- Θxsection.A[i]' .* Xv.xsection.xws[i])))
      @assert expected ≈ Matrix!(projectt(Θxsection.A[i], Θxsection.LA[i],
        Xv.xsection.xws[i], Xv.xsection.xRtws[i], Xv.xsection.xRLws[i], Xv.xsection.xM[i]))
      push!(expectedprojectt, expected)
    end
  end

  combinedexpectedproj = reduce(vcat, expectedprojectt)
  @assert combinedexpectedproj ≈ Matrix!(projectbyt(testA, Xv.xsection))
  @info "completed test of G"


  return nothing
end

#verify that standardized predictions are equivelent to the unstandardized
#(after the scaling adjustment)
function teststandardization(panel::AbstractDataFrame,  zs::ZSpec,
  ::Type{TM}, ::Type{TV}, ::Type{T}) where {TM, TV, T}
  @unpack Zms = zs
  ms = zs.originalms

  ZFw::Vector{Symbol} = Fweights(Zms, :Fw)
  ZFRtw::Vector{Symbol} = Fweights(Zms, :FRtw)
  ZFRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  Fw::Vector{Symbol} = Fweights(ms, :Fw)
  FRtw::Vector{Symbol} = Fweights(ms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(ms, :FRLw)

  p = modelparts(panel, T, TV, TM,
    weightdata = (;Fw, FRtw, FRLw),
    Fvol=ms.Fvol, zs=zs,)
  Zp = modelparts(panel, T, TV, TM,
    weightdata = (;Fw=ZFw, FRtw=ZFRtw, FRLw=ZFRLw),
    Fvol=Zms.Fvol, ms=Zms, zs=zs)

  dims = p.dims
  ws::TM = p.dat[:Fw]
  Rtws::TM = p.dat[:FRtw]
  RLws::TM = p.dat[:FRLw]
  v::TV = p.dat[:Fvol]
  ts::Vector{Int} = p.ts
  Nexpanded::Int = length(ts)


  #use the standardized verions
  Xv::XYIterLevel = XYIterLevel(ws,Rtws, RLws, v, ts)
  Θ = VolumePartsIter(dims.T, dims.K, ts, T)

  Zws::TM = Zp.dat[:Fw]
  ZRtws::TM = Zp.dat[:FRtw]
  ZRLws::TM = Zp.dat[:FRLw]
  Zv::TV = Zp.dat[:Fvol]

  #now test the level versions
  ZXv::XYIterLevel = XYIterLevel(Zws,ZRtws,ZRLws, Zv, ts)
  ZΘ = VolumePartsIter(dims.T, dims.K, ts, T)

  #placeholder values for G
  Θ.G .= (T(0.75) .+ T(0.5) .* rand(Uniform(),size(Θ.G))) .* rand((-1,1,1,1), size(Θ.G)) |> TM
  ZΘ.G .= Θ.G

  #now make the predictions
  iterloss(Θ, Xv, Val(:cholesky))
  iterloss(ZΘ, ZXv, Val(:cholesky))

  #sync up A
  updateAfromGAₜ!(Θ,Xv, 1)
  updateAfromGAₜ!(ZΘ,ZXv, 1)

  ZA = ZΘ.A |> deepcopy |> Matrix
  A = Θ.A |> deepcopy |> Matrix

  Fcoefs::Vector{Symbol} = (ξ->ξ.X[:Fw]).(ms.ξs)
  FZcoefs::Vector{Symbol} = (ξ->ξ.X[:Fw]).(Zms.ξs)

  #de-scale ZA
  ZAscale = (s->zs.scale[s]).(FZcoefs)
  σᵥ::Float64 = zs.scale[zs.Zms.Fvol]
  println("mean scaling ratio (pre-descaling): $(mean(A ./ ZA))")
  AfromZA = (ZA .* σᵥ) ./ ZAscale

  if !( AfromZA ≈ A)
    println("σᵥ: $σᵥ")
    println("ZAscale: $ZAscale")

    println("AfromZA: ")
    printmln(AfromZA[:,1:10])
    println("A:")
    printmln(A[:,1:10])
    println("AfromZA/A:")
    printmln(AfromZA[:,1:10] ./ A[:,1:10])
    error("!( AfromZA ≈ A) !!!")
  end
  @info "Passed standardization test"

  return nothing
end

#this fucntion tests the equivelency of level, intercept, demean specs against the control spec
#note we may want to get rid of this or at least put it in a seperate file
#along with all the legacy specs
function testcontrolequivelency(panel::AbstractDataFrame,  zs::ZSpec,
  ::Type{TM}, ::Type{TV}, ::Type{T};) where {TM,TV,T}


  @info "Beginning control equivlency test"

  CUDA.allowscalar(false)
  Fw::Vector{Symbol} = Fweights(zs.Zms, :Fw)
  FRtw::Vector{Symbol} = Fweights(zs.Zms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(zs.Zms, :FRLw)
  p = modelparts(panel, T, TV, TM,
    weightdata = (;Fw, FRtw, FRLw),
    Fvol=zs.Zms.Fvol, zs=zs)


  dims = p.dims
  ws::TM = p.dat[:Fw]
  Rtws::TM = p.dat[:FRtw]
  RLws::TM = p.dat[:FRLw]
  v::TV = p.dat[:Fvol]
  ts::Vector{Int} = p.ts
  Nexpanded::Int = length(ts)


  Xvverlevel::XYIterLevel = XYIterLevel(ws, Rtws, RLws, v, ts)
  Xvverintercept::XYIterLevel = XYIterLevel(ws ,Rtws, RLws, v, ts)
  Xvverdemean::XYIterDemean = XYIterDemean(ws|>Array ,Rtws|>Array, RLws|>Array, v|>Array, ts) #no gpu version here
  Xvtestlevel::XYIterControl= XYIterControl(ws ,Rtws, RLws, v, ts)
  Xvtestintercept::XYIterControl= XYIterControl(ws ,Rtws, RLws, v, ts, W=ones(Nexpanded,1)|>TM)
  Xvtestdemean::XYIterControl= XYIterControl(ws ,Rtws, RLws, v, ts, Wgroup=ones(Nexpanded,1)|>TM)

  #placeholder values for G


  Θverlevel = VolumePartsIter(dims.T, dims.K, ts, T)
  G = (T(0.75) .+ T(0.5) .* rand(Uniform(),size(Θverlevel.G))) |> TM
  Θverlevel.G .= G
  g = vec(log.(G)) |> Array #use the log form since this is what I use in practice

  #form the other prediction objects
  Θverintercept= Θverlevel |> deepcopy
  Θverdemean = VolumePartsIter(dims.T, dims.K, ts, T, Matrix{T}, Vector{T})
  Θverdemean.G .= G |> Array
  Θtestlevel = Θverlevel |> deepcopy
  Θtestintercept = Θverlevel |> deepcopy
  Θtestdemean = Θverlevel|> deepcopy

  #compute losses
  lossverlevel = iterloss(Θverlevel, Xvverlevel, Val(:cholesky))
  lossverintercept = iterloss(Θverintercept, Xvverintercept, Val(:intercept))
  lossverdemean = iterloss(Θverdemean, Xvverdemean, Val(:cholesky))
  losstestlevel = iterloss(Θtestlevel, Xvtestlevel, Val(:cholesky))
  losstestintercept = iterloss(Θtestintercept, Xvtestintercept, Val(:cholesky))
  losstestdemean = iterloss(Θtestdemean, Xvtestdemean, Val(:cholesky))

  #computepredictions
  μverlevel = Θverlevel(Xvverlevel, Val{:cholesky}()) |> Vector
  μverintercept = Θverintercept(Xvverintercept, Val{:intercept}()) |> Vector
  μverdemean = Θverdemean(Xvverdemean, Val{:cholesky}()) |> Vector
  μtestlevel = Θtestlevel(Xvtestlevel, Val{:cholesky}()) |> Vector
  μtestintercept = Θtestintercept(Xvtestintercept, Val{:cholesky}()) |> Vector
  μtestdemean = Θtestdemean(Xvtestdemean, Val{:cholesky}()) |> Vector

  #check the pre-conditioned volume
  @assert (Xvtestdemean.ṽ |> Array )≈ (Xvverdemean.ṽ |> Array)

  #check the predictions
  @assert μverlevel ≈ μtestlevel
  @assert μverintercept ≈ μtestintercept
  #@assert μverdemean ≈ μtestdemean
  (μverdemean ≈ μtestdemean) || error("!(μverdemean ≈ μtestdemean)!!
    μverdemean[1:10] $(μverdemean[1:10])
    μtestdemean[1:10] $(μtestdemean[1:10])")

  #check the aggregate loss
  @assert lossverlevel ≈ losstestlevel
  @assert lossverintercept ≈ losstestintercept
  @assert lossverdemean ≈ losstestdemean


  #acquire gradient functions
  ∇verlevel! = ∇lossforzygote(Θverlevel, Xvverlevel,
    Val(:log), RegressionType=Val(:choleskyzygote))
  ∇verintercept! = ∇lossforzygote(Θverintercept, Xvverintercept,
    Val(:log), RegressionType=Val(:interceptzygote))
  ∇verdemean! = ∇lossforzygote(Θverdemean, Xvverdemean,
    Val(:log), RegressionType=Val(:choleskyzygote))
  ∇testlevel! = ∇lossforzygote(Θtestlevel, Xvtestlevel,
    Val(:log), RegressionType=Val(:choleskyzygote))
  ∇testintercept! = ∇lossforzygote(Θtestintercept, Xvtestintercept,
    Val(:log), RegressionType=Val(:choleskyzygote))
  ∇testdemean! = ∇lossforzygote(Θtestdemean, Xvtestdemean,
    Val(:log), RegressionType=Val(:choleskyzygote))

  #these will hold the gradient output
  gverlevel  = similar(g)
  gverintercept = similar(g)
  gverdemean = similar(g)
  gtestlevel = similar(g)
  gtestintercept = similar(g)
  gtestdemean = similar(g)

  global g = g
  global gverlevel = gverlevel
  global verlevel = ∇verlevel!

  ∇verlevel!(gverlevel, g)
  ∇verintercept!(gverintercept, g)
  ∇verdemean!(gverdemean, g)
  ∇testlevel!(gtestlevel, g)
  ∇testintercept!(gtestintercept, g)
  ∇testdemean!(gtestdemean, g)

  if PARAM[:benchcontrol]
    print("Benchmark ∇verlevel!(gverlevel, g): ")
    @btime CUDA.@sync $∇verlevel!($gverlevel, $g)
    print("Benchmark ∇verintercept!(gverintercept, g): ")
    @btime CUDA.@sync $∇verintercept!($gverintercept, $g)
    print("Benchmark ∇verdemean!(gverdemean, g): ")
    @btime CUDA.@sync $∇verdemean!($gverdemean, $g)
    print("Benchmark ∇testlevel!(gtestlevel, g): ")
    @btime CUDA.@sync $∇testlevel!($gtestlevel, $g)
    print("Benchmark ∇testintercept!(gtestintercept, g): ")
    @btime CUDA.@sync $∇testintercept!($gtestintercept, $g)
    print("Benchmark ∇testdemean!(gtestdemean, g): ")
    @btime CUDA.@sync $∇testdemean!($gtestdemean, $g)
  end
  @assert gverlevel ≈ gtestlevel
  @assert gverintercept ≈ gtestintercept
  @assert gverdemean ≈ gtestdemean

  @info "Test equivelancy completed successfully."
end

####main entry point for testing
function testiter(panel::AbstractDataFrame,  zs::ZSpec,
    ::Type{T}=PARAM[:itertype],
    ::Type{TV}=PARAM[:itergpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:itergpu] ? CUDA.CuMatrix{T} : Matrix{T};
    tol::Float64 = PARAM[:testtol]
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  @unpack Zms = zs
  #NOTE: disabled as experiment if PARAM[:fluxmodel] ≠ :levellsq
  #  error("PARAM[:fluxmodeltype]≠:levellsq is not coded for use with Jacobian")
  #end


  testcumprodbyrow()
  choltest()
  testXYIterControl()
  teststandardization(panel, zs, TM, TV, T)

  global gtestcontrolequivelency = ()->testcontrolequivelency(panel, zs, TM, TV, T)
  testcontrolequivelency(panel, zs, TM, TV, T)

  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRtw::Vector{Symbol} = Fweights(Zms, :FRtw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)
  essentialfields = [Fw; FRtw; FRLw; Zms.Fvol; ]

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, T, TV, TM,
    weightdata = (;Fw, FRtw, FRLw),
    Fvol=Zms.Fvol, zs=zs)

  @warn "test gpu commented out, should work after dependencies of packages are sorted"
  #testgpu(p)

  dims = p.dims
  ws::TM = p.dat[:Fw]
  Rtws::TM = p.dat[:FRtw]
  RLws::TM = p.dat[:FRLw]
  v::TV = p.dat[:Fvol]
  ts::Vector{Int} = p.ts
  Nexpanded::Int = length(ts)
  #this creates the factored objects- which are accurate given that initially G=1
  #absΔ(wsᵢ::CuArray,RLwsᵢ)=CUDA.abs(wsᵢ-RLwsᵢ)
  absΔ(wsᵢ,Rtwsᵢ,RLwsᵢ)=abs(Rtwsᵢ-RLwsᵢ)


  #now test the level versions
  Xvlevel = XYIterControl(ws , Rtws,RLws, v, ts, W=ones(T, Nexpanded, 1) |> TM) #includes intercept
  Xvdemean = XYIterControl(ws , Rtws,RLws, v, ts, Wgroup=ones(T, Nexpanded, 1) |> TM)

  Θ = VolumePartsIter(dims.T, dims.K, ts, T)

  Dlevel = (ws=Xvlevel.ws, Rtws=Xvlevel.Rtws, RLws=Xvlevel.RLws,)
  xsectionlevel = dataxsection(Dlevel,Θ.ts,TM, TV, T)

  testG(Θ, Xvlevel, xsectionlevel)
  @assert Θ.G |> Matrix! ≈ Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] |> Matrix!

  #WARNING WARNING WARNING uncomment below!!!
  testpredictionzygote(deepcopy(Θ),Xvlevel)
  testiterloss(deepcopy(Θ),Xvlevel)

  Ddemean = (ws=Xvdemean.ws, Rtws=Xvdemean.Rtws, RLws=Xvdemean.RLws, ṽ = Xvdemean.ṽ,)
  xsectiondemean = dataxsection(Ddemean,Θ.ts,TM, TV, T)
  validatexsection(Θ, Xvdemean, p, xsectiondemean)

  testG(Θ, Xvdemean, xsectiondemean)
  @assert Θ.G |> Matrix! ≈ Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] |> Matrix!

  testXv(Θ, Xvlevel, Xvdemean)

  testpredictionzygote(deepcopy(Θ),Xvdemean,)
  testiterloss(deepcopy(Θ),Xvdemean,)


  if (TM<:CuArray)
    pcpu = modelparts(panel, T, Vector{T}, Matrix{T},
      weightdata = (;Fw, FRtw, FRLw),
      Fvol=Zms.Fvol,
      zs=zs)

    wscpu::Matrix{T} = pcpu.dat[:Fw]
    Rtwscpu::Matrix{T} = pcpu.dat[:FRtw]
    RLwscpu::Matrix{T} = pcpu.dat[:FRLw]
    vcpu::Vector{T} = pcpu.dat[:Fvol]

    Xvlevelcpu::XYIterControl = XYIterControl(wscpu, Rtwscpu, RLwscpu, vcpu, ts, W= Xvlevel.W|>Matrix)
    Θcpu = VolumePartsIter(dims.T, dims.K, ts, T, Matrix{T}, Vector{T})
    Θcpu.G .= Θ.G |> Matrix
    Θcpu.A .= Θ.A |> Matrix

    Xvdemeancpu::XYIterControl = XYIterControl(wscpu, Rtwscpu, RLwscpu, vcpu, ts, Wgroup= ones(T, Nexpanded, 1))

    @info "beginning tests on Xvlevel gpu v cpu"
    testcudazygote(deepcopy(Θ),Xvlevel,Θcpu|>deepcopy, Xvlevelcpu)

    @info "beginning tests on Xvdemean gpu v cpu"
    testcudazygote(deepcopy(Θ),Xvdemean,Θcpu|>deepcopy, Xvdemeancpu)
    testag(Θ, Xvlevel)
  end


  return nothing
end
