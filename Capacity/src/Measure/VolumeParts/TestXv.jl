#A bunch of mostly legacy methods
#keep them mostly for testing purposes

###################Additional XYIter types############################
#construct the fixed effect matrix by hand so we have the correct coding
function fixedeffectmatrix(v::TV, ::Type{elTV}=eltype(TV),
  ::Type{T}=PARAM[:itertype]) where {TV, elTV, T}

  #code each unique value
  vunique::TV = unique(v)
  vcoding::Dict{elTV, Int} = Dict{elTV, Int}(k=>i for (i,k) ∈ enumerate(vunique))

  F::Matrix{T} = zeros(T, length(v), length(vcoding))

  Threads.@threads for i ∈ 1:length(v)
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
  ṽ::TV

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
    ṽ = M * v

    #simple integrity check
    D = (ṽ=v, )
    xsection = dataxsection(D, ts, TM, TV, T)
    xsection.ṽ .= (ṽₜ-> ṽₜ .- mean(ṽₜ)).(xsection.ṽ) #demean

    #NOTE: perhaps uncomment the below for safety
    @assert foldl(vcat, xsection.ṽ) ≈ ṽ

    return new{TM,TV, typeof(M), typeof(Fcoding), T}(ws, Rtws, RLws, v, ṽ, M, Fcoding)
  end
end


###############################################
#compare the within results with the fixed effects results
function testwithin(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean;
  factored::AbstractMatrix=error("factored is required"),
  ṽ::AbstractVector=error("ṽ is required"),
  undemeaned::NamedTuple=error("undemeanedfactoredv is required"),
  tsₛ=Θ.ts)  where {TM, TV,T}


  @assert (length(tsₛ) == length(ṽ) == size(factored,1)
    == length(undemeaned.v) == size(undemeaned.factored,1))

  @info "point 0wi"
  bfocal = (Float64).(Vector(reg(factored, ṽ, Val{:cholesky}()) |> vec))
  @info "point5wi"
  #now make a fixed effects version of the within specification
  dfver = DataFrame(Matrix{Float64}(undemeaned.factored), :auto)
  #rhs = Meta.parse(join(names(dfver),"+"))
  rhs = Meta.parse("0 + $(join(names(dfver),"+")) + t")
  form = @eval @formula v ~ $rhs
  dfver.t = (i->string(:t,i)).(tsₛ)
  dfver.v = Vector{Float64}(undemeaned.v)

  local bver::Vector{Float64}
  #try
  bver = FMLM(dfver,  rhs, :v).β[1:(length(bfocal))]

  (Float32.(bver) ≈ Float32.(bfocal)) || error(" bver ≠ bfocal\nbver:$bver\nbfocal: $bfocal
    lm(form, dfver): $(lm(form, dfver))")

  @info "within test complete"
  return nothing
end

function testwithin(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlDemean;
  factored::AbstractMatrix=error("factored is required"),
  ṽ::AbstractVector=error("ṽ is required"),
  undemeaned::NamedTuple=error("undemeanedfactoredv is required"),
  tsₛ=Θ.ts)  where {TM, TV,T}


  @assert (length(tsₛ) == length(ṽ) == size(factored,1)
    == length(undemeaned.v) == size(undemeaned.factored,1))

  @info "point 0wi"
  bfocal = (Float64).(Vector(reg(factored, ṽ, Val{:cholesky}()) |> vec))
  @info "point5wi"
  #now make a fixed effects version of the within specification
  dfver = DataFrame(Matrix{Float64}(undemeaned.factored), :auto)
  #rhs = Meta.parse(join(names(dfver),"+"))
  rhs = Meta.parse("0 + $(join(names(dfver),"+")) + t")
  form = @eval @formula v ~ $rhs
  dfver.t = (i->string(:t,i)).(tsₛ)
  dfver.v = Vector{Float64}(undemeaned.v)

  local bver::Vector{Float64}
  #try
  bver = FMLM(dfver,  rhs, :v).β[1:(length(bfocal))]

  (Float32.(bver) ≈ Float32.(bfocal)) || error(" bver ≠ bfocal\nbver:$bver\nbfocal: $bfocal
    lm(form, dfver): $(lm(form, dfver))")

  @info "within test complete"
  return nothing
end


#pulls the expected volume for each strategy, stock, and time, then factors out Aₜ
#dated- probably should only be called for testing purposes
function factorX(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::Union{XYIterLevel, XYIterControlLevel},
  t::Int) where {TM, TV,T}

  @unpack A, LA = Θ.expand
  @unpack ws, Rtws, RLws= Xv
  #first map each coefficient to its appropriate


  factored = genabsμₖ.(A', LA', ws, Rtws, RLws)

  factored ./= abs.(Θ.A[:,t])'

  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  factoredver::TM = initializeÃandfactor(Θ,Xv) .* abs.(Θ.A[:,1])' ./ abs.(Θ.A[:,t])'


  @assert factoredver ≈ factored

  return factored
end

function factorX(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterDemean,
  t::Int) where {TM, TV,T}

  @unpack A, LA, Ã, LÃ = Θ.expand
  @unpack ws, Rtws, RLws, M = Xv
  #first map each coefficient to its appropriate

  ΠG = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G)
  #A₀ = T(1.0) ./ ΠG[:,dims.T-5] .* Θ.A[:,dims.T - 5]
  @assert Θ.A ≈ ΠG .* Θ.A[:,1]

  factored = genabsμₖ.(A', LA', ws, Rtws, RLws)


  factored ./= abs.(Θ.A[:,t])'
  multiply!(M, factored)



  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  factoredver::TM = (Xv.M * initializeÃandfactor(Θ,Xv)) .* abs.(Θ.A[:,1])' ./ abs.(Θ.A[:,t])'

  @assert Θ.A ≈ ΠG .* Θ.A[:,1]
  @assert Θ.Ã ≈ ΠG
  #println(size(Θ.A))
  #println(size(Θ.Ã))

  if !(factoredver ≈ factored)
    print("G")
    printmln(Θ.G[:,1:10])
    print("At/LAt")
    printmln((Θ.A[:,2:end]./Θ.A[:,(1:(end-1))])[:,1:10])
    print("Ãt/LÃt")
    printmln((Θ.Ã[:,2:end]./Θ.Ã[:,(1:(end-1))])[:,1:10])
    print("factoredver:")
    printmln(factoredver[1:10,:])
    print("factored:")
    printmln(factored[1:10,:])
    error("!(factoredver ≈ factored)!!")
  end

  return factored
end


function factorX(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlDemean,
  t::Int) where {TM, TV,T}

  @unpack A, LA, Ã, LÃ = Θ.expand
  @unpack ws, Rtws, RLws  = Xv
  #first map each coefficient to its appropriate

  ΠG = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprodbyrow(Θ.G)
  #A₀ = T(1.0) ./ ΠG[:,dims.T-5] .* Θ.A[:,dims.T - 5]
  @assert Θ.A ≈ ΠG .* Θ.A[:,1]

  factored = projectbyt(Θ.A, Xv.xsection)
  factored ./= abs.(Θ.A[:,t])'


  #WARNING- below is fort testing only, should be removed if this
  #function is used in anything high performance
  initializeÃandG!(Θ, Xv)
  factoredver::TM = projectbyt(Θ.Ã, Xv.xsection) .* abs.(Θ.A[:,1])' ./ abs.(Θ.A[:,t])'

  @assert Θ.A ≈ ΠG .* Θ.A[:,1]
  @assert Θ.Ã ≈ ΠG
  #println(size(Θ.A))
  #println(size(Θ.Ã))

  if !(factoredver ≈ factored)
    print("G")
    printmln(Θ.G[:,1:10])
    print("At/LAt")
    printmln((Θ.A[:,2:end]./Θ.A[:,(1:(end-1))])[:,1:10])
    print("Ãt/LÃt")
    printmln((Θ.Ã[:,2:end]./Θ.Ã[:,(1:(end-1))])[:,1:10])
    print("factoredver:")
    printmln(factoredver[1:10,:])
    print("factored:")
    printmln(factored[1:10,:])
    error("!(factoredver ≈ factored)!!")
  end

  return factored
end



#we can recycle the testing code for the msot part since
function testXv_(buildXv_, Θ::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel, Xvdemean::XYIterControlDemean;
  t=3, runtestwithin::Bool = false) where {TM,TV,T}

  dims = (T=size(Θ.A,2), K=size(Θ.A,1))

  D = (ws=Xvdemean.ws, Rtws=Xvdemean.Rtws, RLws=Xvdemean.RLws, v=Xvdemean.v, ṽ=Xvdemean.ṽ)
  xsection = dataxsection(D, Θ.ts, TM, TV, T)

  Drep = (factored=similar(Xvlevel.ws), undemeaned=similar(Xvlevel.ws),
    v= deepcopy(Xvdemean.v), ṽ= deepcopy(Xvdemean.ṽ))
  xrep = dataxsection(Drep, Θ.ts, TM,TV,T)
  Θxsection = volumepartsxsection(Θ)

  factored, ṽ = buildXv_(Θ,Xvdemean,t), Xvdemean.ṽ
  @unpack ts = Θ
  Aₜ = Θ.A[:,t]

  #check each cross-section
  tsidxs = collect(1:length(ts))
  for s ∈ 1:(dims.T-1)
    idxs = tsidxs[ts .== s + 1]
    (xsection.ws[s] ≈ Xvdemean.ws[idxs, :]) || error("Unexpected crosssection
      printmln(Xv.xsection.ws[s]):\n $(printmln(Matrix(Xvdemean.xsection.ws[s])))
      printmln(Xv.ws[idxs, :]):\n $(printmln(Matrix(Xvdemean.ws[idxs, :])))")
    @assert xsection.RLws[s] ≈ Xvdemean.RLws[idxs,:]
    @assert xsection.Rtws[s] ≈ Xvdemean.Rtws[idxs,:]
    @assert xsection.ṽ[s] ≈ Xvdemean.ṽ[idxs]
    expected = abs.((xsection.Rtws[s]  .- xsection.RLws[s]).* Θxsection.LA[s]'
      ) .+ abs.(xsection.ws[s] .* Θxsection.A[s]' .- xsection.Rtws[s] .* Θxsection.LA[s]')

    #now do the undemeaning
    xrep.undemeaned[s] = factored[idxs, :]
    xrep.undemeaned[s] .+= mean(expected, dims=1) ./  transpose(Aₜ)

    #the expected factored value against what is actually there
    expected .-= mean(expected, dims=1)
    @assert expected ≈  transpose(Aₜ) .* factored[idxs, :]
    xrep.factored[s] = expected

    expectedv = Xvdemean.v[idxs]
    xrep.v[s] = ṽ[idxs]
    xrep.v[s] .+= mean(expectedv)
    expectedv .-= mean(expectedv)
    @assert expectedv ≈ xsection.ṽ[s]
    @assert expectedv ≈ ṽ[idxs]
    xrep.ṽ[s] = expectedv
  end
  #An important check- the factored X should have the same prediction
  X = abs.((Xvdemean.Rtws  .- Xvdemean.RLws) .* Θ.expand.LA'
    ) .+ abs.(Xvdemean.ws .* Θ.expand.A' .- Xvdemean.Rtws .* Θ.expand.LA')
  undemeaned = (factored=reduce(vcat, xrep.undemeaned),
    v=reduce(vcat, xrep.v))
  repundemeanedAₜ = (undemeaned.factored) .* transpose(Aₜ)
  @assert repundemeanedAₜ ≈ X
  @assert undemeaned.v ≈ Xvdemean.v

  #check the resulting concatenation
  @assert factored .*  transpose(Aₜ) ≈ reduce(vcat, xrep.factored)
  @assert ṽ ≈ reduce(vcat, xrep.ṽ)

  factoredlevel = buildXv_(Θ,Xvlevel, t)
  @assert undemeaned.factored ≈ factoredlevel "
    undemeaned.factored: $(undemeaned.factored[1:10])
    factoredlevel: $(factoredlevel[1:10])"

  runtestwithin && testwithin(Θ,Xvdemean;factored, ṽ, undemeaned,tsₛ=ts)
  return nothing
end

#simpler level only version
function testXv_(buildXv_, Θ::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel;
  t=3) where {TM,TV,T}

  dims = (T=size(Θ.A,2), K=size(Θ.A,1))

  factored= buildXv_(Θ,Xv,t, Val(:level))
  @unpack ts = Θ
  Aₜ = Θ.A[:,t]

  #check each cross-section
  tsidxs = collect(1:length(ts))
  expected = abs.((Xv.Rtws  .- Θ.expand.LA' .* Xv.RLws) .* Θ.expand.LA'
    ) .+ abs.(Θ.expand.LA' .* Xv.Rtws .- Θ.expand.A' .* Xv.ws)

  #An important check- the factored X should have the same prediction

  unfactored = factored .* transpose(Aₜ)
  @assert unfactored ≈ factored

  return nothing
end

#compares two regression methods to make sure they have the same results
#the verify method must be cpu powered
function testregressionequiv(Θ::AbstractVolumePartsIter{TM, TV,T},
  Xv::AbstractXYIter,
  PredictionType::Val,
  RegressionTypeTest::Val=Val(PARAM[:iterregressiontype]),
  RegressionTypeVerify::Val=Val{:glm}()) where {TM, TV,T}

  Xtest=factorX( Θ, Xv, 7)

  #@info "point 0"
  btest = (Float64).(Vector(reg(Xtest,Xv.v, RegressionTypeTest)[1:size(Xtest,2)] |> vec))
  #@info "point 5"
  bverify = (Float64).(reg(Xtest |> Matrix, Xv.v |> Vector, RegressionTypeVerify))

  @assert btest ≈ bverify

  return nothing
end

testregressionequiv(Θ::AbstractVolumePartsIter{TM, TV,T},
  Xv::AbstractXYIter,
  PredictionType::Val,
  RegressionTypeTest::Val{:intercept},
  RegressionTypeVerify::Val=Val{:glm}()) where {TM, TV,T} = testregressionequiv(
    Θ, Xv, PredictionType, Val{:cholesky}(), RegressionTypeVerify)



#tests the factorX routine and methodology
function testfactoredX(Θ::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel, Xvdemean::XYIterControlDemean
  ) where {TM, TV,T}

  testXv_(factorX, Θ, Xvlevel, Xvdemean, t=3, runtestwithin=true)
  testXv_(factorX, Θ, Xvlevel, Xvdemean, t=1, runtestwithin=false)

  testregressionequiv(Θ, Xvdemean,  Val(:leveldemean))

  @info "completed tests of factoredX, including factoring equivelence and a regression check"
end

#level only factorX test routine
function testfactoredX(Θ::AbstractVolumePartsIter{TM,TV,T}, Xvlevel::XYIterControlLevel) where {TM, TV,T}

  testXv_(factorX, Θ, Xvlevel, t=3)
  testXv_(factorX, Θ, Xvlevel, t=1)

  testregressionequiv(Θ, Xvlevel,  Val(:level))

  @info "completed tests of factoredX, including factoring equivelence and a regression check"
end

function testinitialized(Θ::AbstractVolumePartsIter{TM,TV,T}, Xv::XYIterControlLevel) where {TM, TV,T}

  #make sure our factoring works as expected
  Θinit = deepcopy(Θ)
  Θinit.A .= zero(T)
  initializeAfromG!(Θinit, Xv)
  @assert Θinit.A[:, 2:end] ./ Θinit.A[:, 1:(end-1)] |> Matrix! ≈ Θinit.G|> Matrix!
  local factored::TM = factorinitializedX(Θinit,Xv)
  local expected::TM = abs.((Xv.Rtws .- Xv.RLws) .* Θ.expand.LA'
    ) .+ abs.(Xv.Rtws .*Θ.expand.LA' .- (Xv.ws .*Θ.expand.A'))

  predicted = factored .* Θ.A[:,1]' |> Matrix!
  (expected |> Matrix! ≈ predicted) || error(
    "unexpected predicted values
    predicted[1:5,:]: $(predicted[1:5,:])
    expected[1:5,:]: $(expected[1:5,:])"
  )
  @assert expected |> Matrix! ≈ unfactorX(Θinit, Xv, Float64.(Θ.A[:,1]),factored) |> Matrix!

  #stub function to make the below compatable with the rest of testing
  function factorinitializedX_(Θ, Xv, t::Int, ::Any)
    @assert t==1
    factorinitializedX(Θ, Xv)
  end


  return nothing
end

function testinitialized(Θ::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel,
  Xvdemean::XYIterControlDemean) where {TM, TV,T}

  dims = (T=size(Θ.A,2), K=size(Θ.A,1))

  #make sure our factoring works as expected
  Θinit = deepcopy(Θ)
  Θinit.Ã .= zero(T)
  initializeÃandG!(Θinit, Xvdemean)
  factored = projectbyt(Θinit.Ã, Xvdemean.xsection)
  @assert Θinit.Ã[:, 2:end] ./ Θinit.Ã[:, 1:(end-1)] |> Matrix ≈ Θinit.G|> Matrix
  initializeÃandG!(Θ, Xvlevel)
  local expected::TM = abs.(Θ.expand.LÃ' .* Xvlevel.Rtws .- Θ.expand.LÃ' .* Xvlevel.RLws
    ) .+ abs.(Θ.expand.Ã' .*Xvlevel.ws .-  Θ.expand.LÃ' .* Xvlevel.Rtws)
  @assert expected ≈ broadcast(genabsμₖ,  Θinit.expand.Ã', Θinit.expand.LÃ', Xvdemean.ws, Xvdemean.Rtws, Xvdemean.RLws)
  expected .*= Θ.A[:,1]'
  #need to demean
  for t ∈ 1:(dims.T-1)
    inds = collect(1:length(Θ.ts))[Θ.ts .== t+1]
    expected[inds, :] .-= mean(expected[inds, :], dims=1)
  end


  predicted = factored .* Θinit.A[:,1]' |> Matrix!
  (expected |> Matrix! ≈ predicted) || error(
    "unexpected predicted values
    predicted[1:5,:]: $(predicted[1:5,:])
    expected[1:5,:]: $(expected[1:5,:])"
  )

  #stub function to make the below compatable with the rest of testing
  function factorinitializedX_(Θ, Xvdemean::XYIterControlDemean, t::Int)
    @assert t==1
    initializeÃandG!(Θ, Xvdemean)
    factored = projectbyt(Θ.Ã, Xvdemean.xsection)
    return factored

  end
  function factorinitializedX_(Θ, Xv::AbstractXYIter, t::Int)
    @assert t==1
    factored = abs.((Xvdemean.Rtws .-  Xvdemean.RLws) .* Θ.expand.LÃ'
      ) .+ abs.(Θ.expand.Ã' .* Xvdemean.ws .- Xvdemean.Rtws .* Θ.expand.LÃ')
    return factored

  end
  testXv_(factorinitializedX_, Θinit, Xvlevel, Xvdemean, t=1)
end

function testXv(Θin::AbstractVolumePartsIter{TM,TV,T},
  Xvlevel::XYIterControlLevel,
  Xvdemean::XYIterControlDemean) where {TM, TV,T}

  @assert Xvlevel.xsection.xM === nothing #no cross-sectional projection for the level variant

  Θ::VolumePartsIter = deepcopy(Θin)

  #needed since xsection is no longer bundled with Xv
  D = (ws=Xvdemean.ws, Rtws=Xvdemean.Rtws, RLws=Xvdemean.RLws, v=Xvdemean.v, ṽ=Xvdemean.ṽ)
  xsection = dataxsection(D, Θ.ts, TM, TV, T)

  #make sure we have test data
  dims = (T=size(Θ.A,2), K=size(Θ.A,1))
  Θ.G::TM .= T(0.75) .+ TM(rand(T, size(Θ.G))) ./ T(2)

  #sync the test data with G
  Θ.A[:,dims.T - 5] .= 1.0
  updateAfromGAₜ!(Θ, Xvdemean, dims.T - 5)
  ΠG = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprod(Θ.G, dims=2)
  A₀ = T(1.0) ./ ΠG[:,dims.T-5] .* Θ.A[:,dims.T - 5]
  @assert Θ.A ≈ ΠG .* A₀


  testfactoredX(deepcopy(Θ), Xvlevel, Xvdemean)
  testinitialized(Θ, Xvlevel, Xvdemean)

  return nothing
end

function testXv(Θin::AbstractVolumePartsIter{TM, TV,T}, Xvlevel::XYIterControlLevel) where {TM, TV,T}
  Θ::VolumePartsIter = deepcopy(Θin)

  #make sure we have test data
  dims = (T=size(Θ.A,2), K=size(Θ.A,1))
  Θ.G::TM .= T(0.75) .+ TM(rand(T, size(Θ.G))) ./ T(2)

  #sync the test data with G
  Θ.A[:,dims.T - 5] .= 1.0
  updateAfromGAₜ!(Θ, Xvlevel, dims.T - 5)
  ΠG = similar(Θ.A)
  ΠG[:,1] .= T(1.0)
  ΠG[:,2:end] .= cumprod(Θ.G, dims=2)
  A₀ = T(1.0) ./ ΠG[:,dims.T-5] .* Θ.A[:,dims.T - 5]
  @assert Θ.A ≈ ΠG .* A₀


  testfactoredX(deepcopy(Θ), Xvlevel)
  testinitialized(Θ, Xvlevel)

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

  for g ∈ uniquegroups
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

  testgen(wsᵢ, Rtwsᵢ, RLwsᵢ) = abs(Rtwsᵢ - RLwsᵢ) .+ abs(Rtwsᵢ)

  #form the full matrix
  Wsfull = expandsections((xMt->xMt.W).(Xv.xsection.xM), ts .- 1)
  X = testgen.(ws, Rtws, RLws) #not exactly the same as what is used, but serves as a placeholder
  Xver = hcat(X, W, Wsfull)
  bver = cholesky(Symmetric(Xver'*Xver))\(Xver'*Xv.v)

  #check the projected version is equivelent
  X̃ = reduce(vcat, ((xwsₜ, xRtwsₜ, xRLwsₜ, xMₜ)->xMₜ * (testgen.(xwsₜ, xRtwsₜ, xRLwsₜ))).(
    Xv.xsection.xws, Xv.xsection.xRtws, Xv.xsection.xRLws, Xv.xsection.xM))
  Xtest = hcat(X̃, Xv.W̃)
  btest = cholesky(Xtest'*Xtest)\(Xtest'*Xv.ṽ)

  @info "size(Xver): $(size(Xver))"
  @info "size(Xtest): $(size(Xtest))"

  @assert bver[1:(2dims.K)] |> Vector ≈ btest[1:end] |> Vector
  @assert (Xv.ṽ .- Xtest*btest)|>Vector ≈ (Xv.v .- Xver*bver)|>Vector
  @info "testXYIterControl passed"

end

function testXYIterControl()
  testXYIterControl(Matrix{Float64}, Vector{Float64})
  testXYIterControl(CuMatrix{Float64}, CuVector{Float64})
end
#testXYIterControl()
