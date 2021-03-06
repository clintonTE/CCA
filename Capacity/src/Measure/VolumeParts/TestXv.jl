#A bunch of mostly legacy methods
#keep them mostly for testing purposes

###################Additional XYIter types############################
#construct the fixed effect matrix by hand so we have the correct coding
function fixedeffectmatrix(v::TV, ::Type{elTV}=eltype(TV),
  ::Type{T}=PARAM[:itertype]) where {TV, elTV, T}

  #code each unique value
  vunique::TV = unique(v)
  vcoding::Dict{elTV, Int} = Dict{elTV, Int}(k=>i for (i,k) â enumerate(vunique))

  F::Matrix{T} = zeros(T, length(v), length(vcoding))

  Threads.@threads for i â 1:length(v)
    F[i,vcoding[v[i]]] = T(1.0)
  end

  return F, vcoding
end



abstract type AbstractXYIterLevel{
    TM<:AbstractMatrix, TV<:AbstractVector, T<:Real} <: AbstractXYIter{TM, TV, T} end
abstract type AbstractXYIterDemean{
    TM<:AbstractMatrix, TV<:AbstractVector, T<:Real} <: AbstractXYIter{TM, TV, T} end
#holds all unchanging data as well as helper pre-allocations
struct XYIterLevel{TM<:AbstractMatrix,
    TV<:AbstractVector,
    T<:Real} <: AbstractXYIterLevel{TM,TV,T}
  ws::TM
  Rtws::TM
  RLws::TM
  v::TV
  tv1s::TV

  #constructor from minimum componeents
  XYIterLevel(ws::TM,Rtws::TM, RLws::TM,  v::TV, ::Tts,
    ::Type{T}=eltype(TM)) where {TM, TV, Tts, T} = (
    new{TM, TV, T}(ws, Rtws, RLws, v, ones(T, length(v)) |> TV))
end


#holds all unchanging data as well as helper pre-allocations
#this version also includes the projection matrix for demeaning
struct XYIterDemean{TM<:AbstractMatrix,
    TV<:AbstractVector,
    TProjection<:FixedAnnihilator,
    TFcoding<:Dict,
    T<:Real} <: AbstractXYIter{TM,TV,T}
  ws::TM
  Rtws::TM
  RLws::TM
  v::TV
  vĖ::TV

  M::TProjection
  Fcoding::TFcoding


  #constructor from minimum componeents
  function XYIterDemean(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts,
      ::Type{T}=eltype(TM)) where {TM,TV, Tts, T}

    F, Fcoding = fixedeffectmatrix(ts)


    #construct a compact form of the annihilator matrix
    #only pre-multiplication of the matrix on a matrix or vector is supported
    M = FixedAnnihilator(F, TM)

    #demean/project v (this won't change, so we might as well cache it)
    vĖ = M * v

    #simple integrity check
    D = (vĖ=v, )
    xsection = dataxsection(D, ts, TM, TV, T)
    xsection.vĖ .= (vĖâ-> vĖâ .- mean(vĖâ)).(xsection.vĖ) #demean

    #NOTE: perhaps uncomment the below for safety
    @assert foldl(vcat, xsection.vĖ) â vĖ

    return new{TM,TV, typeof(M), typeof(Fcoding), T}(ws, Rtws, RLws, v, vĖ, M, Fcoding)
  end
end


###############################################
#compare the within results with the fixed effects results
function testwithin(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean;
  factored::AbstractMatrix=error("factored is required"),
  vĖ::AbstractVector=error("vĖ is required"),
  undemeaned::NamedTuple=error("undemeanedfactoredv is required"),
  tsâ=Î.ts)  where {TM, TV,T}


  @assert (length(tsâ) == length(vĖ) == size(factored,1)
    == length(undemeaned.v) == size(undemeaned.factored,1))

  @info "point 0wi"
  bfocal = (Float64).(Vector(reg(factored, vĖ, Val{:cholesky}()) |> vec))
  @info "point5wi"
  #now make a fixed effects version of the within specification
  dfver = DataFrame(Matrix{Float64}(undemeaned.factored), :auto)
  #rhs = Meta.parse(join(names(dfver),"+"))
  rhs = Meta.parse("0 + $(join(names(dfver),"+")) + t")
  form = @eval @formula v ~ $rhs
  dfver.t = (i->string(:t,i)).(tsâ)
  dfver.v = Vector{Float64}(undemeaned.v)

  local bver::Vector{Float64}
  #try
  bver = FMLM(dfver,  rhs, :v).Îē[1:(length(bfocal))]

  (Float32.(bver) â Float32.(bfocal)) || error(" bver â  bfocal\nbver:$bver\nbfocal: $bfocal
    lm(form, dfver): $(lm(form, dfver))")

  @info "within test complete"
  return nothing
end

function testwithin(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlDemean;
  factored::AbstractMatrix=error("factored is required"),
  vĖ::AbstractVector=error("vĖ is required"),
  undemeaned::NamedTuple=error("undemeanedfactoredv is required"),
  tsâ=Î.ts)  where {TM, TV,T}


  @assert (length(tsâ) == length(vĖ) == size(factored,1)
    == length(undemeaned.v) == size(undemeaned.factored,1))

  @info "point 0wi"
  bfocal = (Float64).(Vector(reg(factored, vĖ, Val{:cholesky}()) |> vec))
  @info "point5wi"
  #now make a fixed effects version of the within specification
  dfver = DataFrame(Matrix{Float64}(undemeaned.factored), :auto)
  #rhs = Meta.parse(join(names(dfver),"+"))
  rhs = Meta.parse("0 + $(join(names(dfver),"+")) + t")
  form = @eval @formula v ~ $rhs
  dfver.t = (i->string(:t,i)).(tsâ)
  dfver.v = Vector{Float64}(undemeaned.v)

  local bver::Vector{Float64}
  #try
  bver = FMLM(dfver,  rhs, :v).Îē[1:(length(bfocal))]

  (Float32.(bver) â Float32.(bfocal)) || error(" bver â  bfocal\nbver:$bver\nbfocal: $bfocal
    lm(form, dfver): $(lm(form, dfver))")

  @info "within test complete"
  return nothing
end


#pulls the expected volume for each strategy, stock, and time, then factors out Aâ
#dated- probably should only be called for testing purposes
function factorX(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::Union{XYIterLevel, XYIterControlLevel},
  t::Int) where {TM, TV,T}

  @unpack A, LA = Î.expand
  @unpack ws, Rtws, RLws= Xv
  #first map each coefficient to its appropriate


  factored = genabsÎžâ.(A', LA', ws, Rtws, RLws)

  factored ./= abs.(Î.A[:,t])'

  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  factoredver::TM = initializeAĖandfactor(Î,Xv) .* abs.(Î.A[:,1])' ./ abs.(Î.A[:,t])'


  @assert factoredver â factored

  return factored
end

function factorX(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean,
  t::Int) where {TM, TV,T}

  @unpack A, LA, AĖ, LAĖ = Î.expand
  @unpack ws, Rtws, RLws, M = Xv
  #first map each coefficient to its appropriate

  Î G = similar(Î.A)
  Î G[:,1] .= T(1.0)
  Î G[:,2:end] .= cumprodbyrow(Î.G)
  #Aâ = T(1.0) ./ Î G[:,dims.T-5] .* Î.A[:,dims.T - 5]
  @assert Î.A â Î G .* Î.A[:,1]

  factored = genabsÎžâ.(A', LA', ws, Rtws, RLws)


  factored ./= abs.(Î.A[:,t])'
  multiply!(M, factored)



  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  factoredver::TM = (Xv.M * initializeAĖandfactor(Î,Xv)) .* abs.(Î.A[:,1])' ./ abs.(Î.A[:,t])'

  @assert Î.A â Î G .* Î.A[:,1]
  @assert Î.AĖ â Î G
  #println(size(Î.A))
  #println(size(Î.AĖ))

  if !(factoredver â factored)
    print("G")
    printmln(Î.G[:,1:10])
    print("At/LAt")
    printmln((Î.A[:,2:end]./Î.A[:,(1:(end-1))])[:,1:10])
    print("AĖt/LAĖt")
    printmln((Î.AĖ[:,2:end]./Î.AĖ[:,(1:(end-1))])[:,1:10])
    print("factoredver:")
    printmln(factoredver[1:10,:])
    print("factored:")
    printmln(factored[1:10,:])
    error("!(factoredver â factored)!!")
  end

  return factored
end


function factorX(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlDemean,
  t::Int) where {TM, TV,T}

  @unpack A, LA, AĖ, LAĖ = Î.expand
  @unpack ws, Rtws, RLws  = Xv
  #first map each coefficient to its appropriate

  Î G = similar(Î.A)
  Î G[:,1] .= T(1.0)
  Î G[:,2:end] .= cumprodbyrow(Î.G)
  #Aâ = T(1.0) ./ Î G[:,dims.T-5] .* Î.A[:,dims.T - 5]
  @assert Î.A â Î G .* Î.A[:,1]

  factored = projectbyt(Î.A, Xv.xsection)
  factored ./= abs.(Î.A[:,t])'


  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  initializeÃandG!(Î, Xv)
  factoredver::TM = projectbyt(Î.AĖ, Xv.xsection) .* abs.(Î.A[:,1])' ./ abs.(Î.A[:,t])'

  @assert Î.A â Î G .* Î.A[:,1]
  @assert Î.AĖ â Î G
  #println(size(Î.A))
  #println(size(Î.AĖ))

  if !(factoredver â factored)
    print("G")
    printmln(Î.G[:,1:10])
    print("At/LAt")
    printmln((Î.A[:,2:end]./Î.A[:,(1:(end-1))])[:,1:10])
    print("AĖt/LAĖt")
    printmln((Î.AĖ[:,2:end]./Î.AĖ[:,(1:(end-1))])[:,1:10])
    print("factoredver:")
    printmln(factoredver[1:10,:])
    print("factored:")
    printmln(factored[1:10,:])
    error("!(factoredver â factored)!!")
  end

  return factored
end



#we can recycle the testing code for the msot part since
function testXv_(buildXv_, Î::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel, Xvdemean::XYIterControlDemean;
  t=3, runtestwithin::Bool = false) where {TM,TV,T}

  dims = (T=size(Î.A,2), K=size(Î.A,1))

  D = (ws=Xvdemean.ws, Rtws=Xvdemean.Rtws, RLws=Xvdemean.RLws, v=Xvdemean.v, vĖ=Xvdemean.vĖ)
  xsection = dataxsection(D, Î.ts, TM, TV, T)

  Drep = (factored=similar(Xvlevel.ws), undemeaned=similar(Xvlevel.ws),
    v= deepcopy(Xvdemean.v), vĖ= deepcopy(Xvdemean.vĖ))
  xrep = dataxsection(Drep, Î.ts, TM,TV,T)
  Îxsection = volumepartsxsection(Î)

  factored, vĖ = buildXv_(Î,Xvdemean,t), Xvdemean.vĖ
  @unpack ts = Î
  Aâ = Î.A[:,t]

  #check each cross-section
  tsidxs = collect(1:length(ts))
  for s â 1:(dims.T-1)
    idxs = tsidxs[ts .== s + 1]
    (xsection.ws[s] â Xvdemean.ws[idxs, :]) || error("Unexpected crosssection
      printmln(Xv.xsection.ws[s]):\n $(printmln(Matrix(Xvdemean.xsection.ws[s])))
      printmln(Xv.ws[idxs, :]):\n $(printmln(Matrix(Xvdemean.ws[idxs, :])))")
    @assert xsection.RLws[s] â Xvdemean.RLws[idxs,:]
    @assert xsection.Rtws[s] â Xvdemean.Rtws[idxs,:]
    @assert xsection.vĖ[s] â Xvdemean.vĖ[idxs]
    expected = abs.((xsection.Rtws[s]  .- xsection.RLws[s]).* Îxsection.LA[s]'
      ) .+ abs.(xsection.ws[s] .* Îxsection.A[s]' .- xsection.Rtws[s] .* Îxsection.LA[s]')

    #now do the undemeaning
    xrep.undemeaned[s] = factored[idxs, :]
    xrep.undemeaned[s] .+= mean(expected, dims=1) ./  transpose(Aâ)

    #the expected factored value against what is actually there
    expected .-= mean(expected, dims=1)
    @assert expected â  transpose(Aâ) .* factored[idxs, :]
    xrep.factored[s] = expected

    expectedv = Xvdemean.v[idxs]
    xrep.v[s] = vĖ[idxs]
    xrep.v[s] .+= mean(expectedv)
    expectedv .-= mean(expectedv)
    @assert expectedv â xsection.vĖ[s]
    @assert expectedv â vĖ[idxs]
    xrep.vĖ[s] = expectedv
  end
  #An important check- the factored X should have the same prediction
  X = abs.((Xvdemean.Rtws  .- Xvdemean.RLws) .* Î.expand.LA'
    ) .+ abs.(Xvdemean.ws .* Î.expand.A' .- Xvdemean.Rtws .* Î.expand.LA')
  undemeaned = (factored=reduce(vcat, xrep.undemeaned),
    v=reduce(vcat, xrep.v))
  repundemeanedAâ = (undemeaned.factored) .* transpose(Aâ)
  @assert repundemeanedAâ â X
  @assert undemeaned.v â Xvdemean.v

  #check the resulting concatenation
  @assert factored .*  transpose(Aâ) â reduce(vcat, xrep.factored)
  @assert vĖ â reduce(vcat, xrep.vĖ)

  factoredlevel = buildXv_(Î,Xvlevel, t)
  @assert undemeaned.factored â factoredlevel "
    undemeaned.factored: $(undemeaned.factored[1:10])
    factoredlevel: $(factoredlevel[1:10])"

  runtestwithin && testwithin(Î,Xvdemean;factored, vĖ, undemeaned,tsâ=ts)
  return nothing
end

#simpler level only version
function testXv_(buildXv_, Î::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel;
  t=3) where {TM,TV,T}

  dims = (T=size(Î.A,2), K=size(Î.A,1))

  factored= buildXv_(Î,Xv,t, Val(:level))
  @unpack ts = Î
  Aâ = Î.A[:,t]

  #check each cross-section
  tsidxs = collect(1:length(ts))
  expected = abs.((Xv.Rtws  .- Î.expand.LA' .* Xv.RLws) .* Î.expand.LA'
    ) .+ abs.(Î.expand.LA' .* Xv.Rtws .- Î.expand.A' .* Xv.ws)

  #An important check- the factored X should have the same prediction

  unfactored = factored .* transpose(Aâ)
  @assert unfactored â factored

  return nothing
end

#compares two regression methods to make sure they have the same results
#the verify method must be cpu powered
function testregressionequiv(Î::AbstractVolumePartsIter{TM, TV,T},
  Xv::AbstractXYIter,
  PredictionType::Val,
  RegressionTypeTest::Val=Val(PARAM[:iterregressiontype]),
  RegressionTypeVerify::Val=Val{:glm}()) where {TM, TV,T}

  Xtest=factorX( Î, Xv, 7)

  #@info "point 0"
  btest = (Float64).(Vector(reg(Xtest,Xv.v, RegressionTypeTest)[1:size(Xtest,2)] |> vec))
  #@info "point 5"
  bverify = (Float64).(reg(Xtest |> Matrix, Xv.v |> Vector, RegressionTypeVerify))

  @assert btest â bverify

  return nothing
end

testregressionequiv(Î::AbstractVolumePartsIter{TM, TV,T},
  Xv::AbstractXYIter,
  PredictionType::Val,
  RegressionTypeTest::Val{:intercept},
  RegressionTypeVerify::Val=Val{:glm}()) where {TM, TV,T} = testregressionequiv(
    Î, Xv, PredictionType, Val{:cholesky}(), RegressionTypeVerify)



#tests the factorX routine and methodology
function testfactoredX(Î::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel, Xvdemean::XYIterControlDemean
  ) where {TM, TV,T}

  testXv_(factorX, Î, Xvlevel, Xvdemean, t=3, runtestwithin=true)
  testXv_(factorX, Î, Xvlevel, Xvdemean, t=1, runtestwithin=false)

  testregressionequiv(Î, Xvdemean,  Val(:leveldemean))

  @info "completed tests of factoredX, including factoring equivelence and a regression check"
end

#level only factorX test routine
function testfactoredX(Î::AbstractVolumePartsIter{TM,TV,T}, Xvlevel::XYIterControlLevel) where {TM, TV,T}

  testXv_(factorX, Î, Xvlevel, t=3)
  testXv_(factorX, Î, Xvlevel, t=1)

  testregressionequiv(Î, Xvlevel,  Val(:level))

  @info "completed tests of factoredX, including factoring equivelence and a regression check"
end

function testinitialized(Î::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlLevel) where {TM, TV,T}

  #make sure our factoring works as expected
  Îinit = deepcopy(Î)
  Îinit.A .= zero(T)
  initializeAfromG!(Îinit, Xv)
  @assert Îinit.A[:, 2:end] ./ Îinit.A[:, 1:(end-1)] |> Matrix! â Îinit.G|> Matrix!
  local factored::TM = factorinitializedX(Îinit,Xv)
  local expected::TM = abs.((Xv.Rtws .- Xv.RLws) .* Î.expand.LA'
    ) .+ abs.(Xv.Rtws .*Î.expand.LA' .- (Xv.ws .*Î.expand.A'))

  predicted = factored .* Î.A[:,1]' |> Matrix!
  (expected |> Matrix! â predicted) || error(
    "unexpected predicted values
    predicted[1:5,:]: $(predicted[1:5,:])
    expected[1:5,:]: $(expected[1:5,:])"
  )
  @assert expected |> Matrix! â unfactorX(Îinit, Xv, Float64.(Î.A[:,1]),factored) |> Matrix!

  #stub function to make the below compatable with the rest of testing
  function factorinitializedX_(Î, Xv, t::Int, ::Any)
    @assert t==1
    factorinitializedX(Î, Xv)
  end


  return nothing
end

function testinitialized(Î::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel,
  Xvdemean::XYIterControlDemean) where {TM, TV,T}

  dims = (T=size(Î.A,2), K=size(Î.A,1))

  #make sure our factoring works as expected
  Îinit = deepcopy(Î)
  Îinit.AĖ .= zero(T)
  initializeÃandG!(Îinit, Xvdemean)
  factored = projectbyt(Îinit.AĖ, Xvdemean.xsection)
  @assert Îinit.AĖ[:, 2:end] ./ Îinit.AĖ[:, 1:(end-1)] |> Matrix â Îinit.G|> Matrix
  initializeAĖandG!(Î, Xvlevel)
  local expected::TM = abs.(Î.expand.LAĖ' .* Xvlevel.Rtws .- Î.expand.LAĖ' .* Xvlevel.RLws
    ) .+ abs.(Î.expand.AĖ' .*Xvlevel.ws .-  Î.expand.LAĖ' .* Xvlevel.Rtws)
  @assert expected â broadcast(genabsÎžâ,  Îinit.expand.AĖ', Îinit.expand.LAĖ', Xvdemean.ws, Xvdemean.Rtws, Xvdemean.RLws)
  expected .*= Î.A[:,1]'
  #need to demean
  for t â 1:(dims.T-1)
    inds = collect(1:length(Î.ts))[Î.ts .== t+1]
    expected[inds, :] .-= mean(expected[inds, :], dims=1)
  end


  predicted = factored .* Îinit.A[:,1]' |> Matrix!
  (expected |> Matrix! â predicted) || error(
    "unexpected predicted values
    predicted[1:5,:]: $(predicted[1:5,:])
    expected[1:5,:]: $(expected[1:5,:])"
  )

  #stub function to make the below compatable with the rest of testing
  function factorinitializedX_(Î, Xvdemean::XYIterControlDemean, t::Int)
    @assert t==1
    initializeÃandG!(Î, Xvdemean)
    factored = projectbyt(Î.AĖ, Xvdemean.xsection)
    return factored

  end
  function factorinitializedX_(Î, Xv::AbstractXYIter, t::Int)
    @assert t==1
    factored = abs.((Xvdemean.Rtws .-  Xvdemean.RLws) .* Î.expand.LAĖ'
      ) .+ abs.(Î.expand.AĖ' .* Xvdemean.ws .- Xvdemean.Rtws .* Î.expand.LAĖ')
    return factored

  end
  testXv_(factorinitializedX_, Îinit, Xvlevel, Xvdemean, t=1)
end

function testXv(Îin::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel,
  Xvdemean::XYIterControlDemean) where {TM, TV,T}

  @assert Xvlevel.xsection.xM === nothing #no cross-sectional projection for the level variant

  Î::VolumePartsIter = deepcopy(Îin)

  #needed since xsection is no longer bundled with Xv
  D = (ws=Xvdemean.ws, Rtws=Xvdemean.Rtws, RLws=Xvdemean.RLws, v=Xvdemean.v, vĖ=Xvdemean.vĖ)
  xsection = dataxsection(D, Î.ts, TM, TV, T)

  #make sure we have test data
  dims = (T=size(Î.A,2), K=size(Î.A,1))
  Î.G::TM .= T(0.75) .+ TM(rand(T, size(Î.G))) ./ T(2)

  #sync the test data with G
  Î.A[:,dims.T - 5] .= 1.0
  updateAfromGAâ!(Î, Xvdemean, dims.T - 5)
  Î G = similar(Î.A)
  Î G[:,1] .= T(1.0)
  Î G[:,2:end] .= cumprod(Î.G, dims=2)
  Aâ = T(1.0) ./ Î G[:,dims.T-5] .* Î.A[:,dims.T - 5]
  @assert Î.A â Î G .* Aâ


  testfactoredX(deepcopy(Î), Xvlevel, Xvdemean)
  testinitialized(Î, Xvlevel, Xvdemean)

  return nothing
end

function testXv(Îin::AbstractVolumePartsIter{TM, TV,T}, Xvlevel::XYIterControlLevel) where {TM, TV,T}
  Î::VolumePartsIter = deepcopy(Îin)

  #make sure we have test data
  dims = (T=size(Î.A,2), K=size(Î.A,1))
  Î.G::TM .= T(0.75) .+ TM(rand(T, size(Î.G))) ./ T(2)

  #sync the test data with G
  Î.A[:,dims.T - 5] .= 1.0
  updateAfromGAâ!(Î, Xvlevel, dims.T - 5)
  Î G = similar(Î.A)
  Î G[:,1] .= T(1.0)
  Î G[:,2:end] .= cumprod(Î.G, dims=2)
  Aâ = T(1.0) ./ Î G[:,dims.T-5] .* Î.A[:,dims.T - 5]
  @assert Î.A â Î G .* Aâ


  testfactoredX(deepcopy(Î), Xvlevel)
  testinitialized(Î, Xvlevel)

  return nothing
end



#expand each group to a full length vector and pad with zeros
#assume each group is already consecutively numbered as an integer starting at 1
function expandsections(Ws::AbstractVector{TM}, groups::Vector{Int},::Type{T}=eltype(TM)
    ) where {TM<:AbstractMatrix,T}
  uniquegroups = unique(groups) |> sort!
  @assert all(collect(1:length(Ws)) == uniquegroups)
  dims=(; Nexpanded=length(groups), K=size(Ws[1],2), Ngroups=length(uniquegroups))

  expanded = zeros(T, length(groups), dims.K*dims.Ngroups) |> TM
  expanded .= T(0)

  for g â uniquegroups
    inds = (1:dims.Nexpanded)[groups .== g]
    startcol = (g-1)*dims.K+1
    sexpanded = expanded[inds, startcol:(startcol + dims.K - 1)]

    #NOTE: these checks are low performance, and should probably be disabled if this
    #function is used in any performance sensitive context
    @assert all((sexpanded |> Matrix) .== T(0))
    @assert size(sexpanded) === size(Ws[g])
    @view(expanded[inds, startcol:(startcol + dims.K - 1)]) .= Ws[g]
  end

  return expanded
end



#set up a placeholder, run the full regression and make sure its equivelent
function testXYIterControl(::Type{TM}, ::Type{TV}, ::Type{T}=eltype(TM), dims=(;N=100, K=3, T=5)
  ) where {TM, TV, T}
  #make placeholder data
  ws = rand(T, dims.N,dims.K) |> TM
  Rtws = rand(T, dims.N,dims.K) |> TM
  RLws = rand(T, dims.N,dims.K) |> TM
  ts = rand(2:dims.T, dims.N) |> sort!
  Wgroup = rand(T, dims.N, 2*dims.K) |> TM
  W = rand(T, dims.N, dims.K) |> TM
  v = rand(T, dims.N) |> TV

  Xv = XYIterControl(ws, Rtws, RLws, v, ts; W, Wgroup)

  testgen(wsáĩĒ, RtwsáĩĒ, RLwsáĩĒ) = abs(RtwsáĩĒ - RLwsáĩĒ) .+ abs(RtwsáĩĒ)

  #form the full matrix
  Wsfull = expandsections((xMt->xMt.W).(Xv.xsection.xM), ts .- 1)
  X = testgen.(ws, Rtws, RLws) #not exactly the same as what is used, but serves as a placeholder
  Xver = hcat(X, W, Wsfull)
  bver = cholesky(Symmetric(Xver'*Xver))\(Xver'*Xv.v)

  #check the projected version is equivelent
  XĖ = reduce(vcat, ((xwsâ, xRtwsâ, xRLwsâ, xMâ)->xMâ * (testgen.(xwsâ, xRtwsâ, xRLwsâ))).(
    Xv.xsection.xws, Xv.xsection.xRtws, Xv.xsection.xRLws, Xv.xsection.xM))
  Xtest = hcat(XĖ, Xv.WĖ)
  btest = cholesky(Xtest'*Xtest)\(Xtest'*Xv.vĖ)

  @info "size(Xver): $(size(Xver))"
  @info "size(Xtest): $(size(Xtest))"

  @assert bver[1:(2dims.K)] |> Vector â btest[1:end] |> Vector
  @assert (Xv.vĖ .- Xtest*btest)|>Vector â (Xv.v .- Xver*bver)|>Vector
  @info "testXYIterControl passed"

end

function testXYIterControl()
  testXYIterControl(Matrix{Float64}, Vector{Float64})
  testXYIterControl(CuMatrix{Float64}, CuVector{Float64})
end
#testXYIterControl()
