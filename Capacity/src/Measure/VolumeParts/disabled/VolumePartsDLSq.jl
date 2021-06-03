

#holds the time-indexed parameters and a pre-allocation for the row-indexed versions
struct VolumePartsDLSq{
    TΠ<:AbstractVector,
    TA<:AbstractMatrix,
    TV₀<:AbstractVector,
    Tts<:AbstractVector,
    Texpand<:NamedTuple,
    elTΠ<:Real} <: AbstractVolumePartsLSq{TΠ, TA, Tts, Texpand, elTΠ}

  Π::TΠ
  A::TA
  V₀::TV₀

  #not to be accessed directly- indices for the expand function
  ts::Tts
  expand::Texpand

  #inner constructor needed to identify elTΠ
  function VolumePartsDLSq(Π::TΠ, A::TA, V₀::TV₀, ts::Tts, expand::Texpand) where {
      TΠ, TA, TV₀, Tts, Texpand}

    @assert (Π == reduce(vcat, [vec(A), vec(V₀)]))

    return new{TΠ, TA, TV₀, Tts, Texpand, eltype(TΠ)}(Π,A,V₀,ts,expand)
  end
end

#pre-allocates for the vpexpand
function VolumePartsDLSq(A::TA, V₀::TV₀, ts::Tts) where {
    TA<:AbstractMatrix,
    TV₀<:AbstractVector,
    Tts<:AbstractVector}

  #create vector versions for the parameter vector
  vecA = vec(A)
  vecV₀ = vec(V₀)

  #this will hold the parameter data
  Π = reduce(vcat, [vecA, vecV₀])
  start=1
  linkedA = less_unsafe_wrap(TA, source=Π, start=start, dims=size(A))
  start+=length(vecA)
  linkedV₀ = less_unsafe_wrap(TV₀, source=Π, start=start, dims=size(V₀))
  start+=length(vecV₀)
  @assert start - 1 == length(Π)

  #final checks on the linked dataframes
  @assert sum(linkedA .== A) + sum(linkedV₀ .== V₀) == length(Π)

  expand = (A = view(linkedA, :, ts),
    LA = view(linkedA, :, ts .- 1),
    V₀ = view(linkedV₀, ts .- 1))

  Θ = VolumePartsDLSq(Π, linkedA, linkedV₀, ts, expand)
  return Θ
end


#basic constructor from dimension arguments
function VolumePartsDLSq(T::Int, K::Int, ts::Vector{Int}, ::Type{Tvp}=PARAM[:opttype],
  ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{Tvp} : Vector{Tvp},
  ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{Tvp} : Matrix{Tvp},
    ) where {Tvp<:Real, TV<:AbstractVector{Tvp}, TM<:AbstractMatrix{Tvp}}

  #note- includes a work-around until Flux.gpu works
  A = TM(abs.(randn(Tvp, K,T)))
  V₀ = TV(abs.(randn(Tvp, T-1)))

  return VolumePartsDLSq(A, V₀, ts)
end


#creates new VolumeParts with copies of the trainable portion and pass-by-reference
# of the non-trainable indices
function copytrainable(Θ::VolumePartsDLSq)
  #copy trainable portions by value
  A  = deepcopy(Θ.A)
  V₀  = deepcopy(Θ.V₀)

  #copy the non-trainable portions by reference
  return VolumePartsDLSq(A, V₀, Θ.ts)
end

#needed due to the segemented array structure
function Base.deepcopy(Θ::VolumePartsDLSq)
  A  = deepcopy(Θ.A)
  V₀  = deepcopy(Θ.V₀)
  ts = deepcopy(Θ.ts)

  return VolumePartsDLSq(A,V₀,ts)
end


function (Θ::VolumePartsDLSq)(X::NamedTuple, ::AbstractVector, ::Val{:level})
  @unpack ws, RLws = X
  #first map each coefficient to its appropriate

  absμₖ = abs.(Θ.expand.A .* ws' .- Θ.expand.LA .* RLws')
  absμ = vec(sum(absμₖ,dims=1)) .+ exp.(Θ.expand.V₀)

  return absμ
end

function (Θ::VolumePartsDLSq)(X::NamedTuple, ::CuVector{T}, ::Val{:level}) where T
  @unpack ws, RLws = X
  @unpack A, LA, V₀ = Θ.expand
  #first map each coefficient to its appropriate

  absμₖ = cudaabsμₖ.(A, LA, ws', RLws')
  absμ = CUDA.vec(CUDA.sum(absμₖ,dims=1)) .+ CUDA.exp.(V₀)


  return absμ
end

function (Θ::VolumePartsDLSq)(X::NamedTuple, ::AbstractVector, ::Val{:levellog})
  @unpack ws, RLws = X
  #first map each coefficient to its appropriate

  absμₖ = abs.((exp).(Θ.expand.A) .* ws' .- (exp).(Θ.expand.LA) .* RLws')
  absμ = vec(sum(absμₖ,dims=1)) .+ exp.(Θ.expand.V₀)

  return absμ
end


Flux.@functor VolumePartsDLSq
Flux.params(Θ::VolumePartsDLSq) = Flux.params(Θ.A, Θ.V₀,)





struct JacobianCache{TVolumeParts<:AbstractVolumeParts,
    TJcolmap<:NamedTuple,
    Tindex<:NamedTuple,
    TbaseJ<:AbstractMatrix}

  Jcolmap::TJcolmap
  index::Tindex
  baseJ::TbaseJ

  #inner constructor captures the type of VolumeParts
  function JacobianCache(::TVolumeParts, Jcolmap::TJcolmap, index::Tindex, baseJ::TbaseJ) where {
      TVolumeParts,TJcolmap, Tindex, TbaseJ}
    return new{TVolumeParts, TJcolmap, Tindex, TbaseJ}(Jcolmap, index, baseJ)
  end
end

function jacobiancolmap(Θ::VolumePartsDLSq)

  #map each entry in the matrix to a corresponding entry in A
  start::Int = 1
  start, Amap = map2linear(start, size(Θ.A))
  start, V₀map = map2linear(start, size(Θ.V₀))

  #checks that should always be true
  @assert start == length(Θ) + 1
  @assert Θ.A == reshape(Θ.Π[vec(Amap)], size(Θ.A))
  @assert Θ.V₀ == Θ.Π[V₀map]

  for (i,Πidx) ∈ enumerate(Amap)
    @assert Θ.A[i] == Θ[Πidx]
  end

  for (i,Πidx) ∈ enumerate(V₀map)
    @assert Θ.V₀[i] == Θ[Πidx]
  end

  return (A=Amap, V₀=V₀map)
end

#NOTE: maybe try a list of cartesian coordinates to speed things up?

#pre-allocates the sparsity array
#a nice feature here is the magnitudes are constant for A, just the signs change
function baseleveljacobian(Θ::VolumePartsDLSq, Jcolmap::NamedTuple, X::NamedTuple,
  ::Ty, ::Type{T}) where {Ty<:AbstractVector, T}
  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "PARAM[:predictiontype] $(PARAM[:predictiontype]) not supported for manual jacobians!")
  @unpack ws, RLws = X

  #useful quantities for indexing
  N::Int = size(Θ.A,2)
  K::Int = size(Θ.A,1)
  Nexpand::Int = length(Θ.ts)
  Nnonsparse = 2*Nexpand*K + Nexpand #two expanded A matrices + V0 worth of coefficients

  #stores the sparse matrix vectors in ijv format
  baseJvec = (
    I=Vector{Int}(undef, Nnonsparse),
    J=Vector{Int}(undef, Nnonsparse),
    V=Vector{T}(undef, Nnonsparse))#spzeros(T, Nexpand, length(Θ.Π))

  #coordA is a vector of linear indices that maps, for each k and row i,
  #the relevant output row (i) and input col (j) in the sparse J matrix
  coordA::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef, Nexpand, K)
  coordLA::Matrix{CartesianIndex} = Matrix{CartesianIndex}(undef, Nexpand, K)
  coordV₀ = Vector{CartesianIndex}(undef, Nexpand)

  #first set the coordinates
  coordctr::Int = 1
  for k ∈ 1:K
    for i ∈ 1:length(Θ.ts)
      t=Θ.ts[i]

      Aj = Jcolmap.A[k,t]
      LAj = Jcolmap.A[k,t-1]

      coordA[i,k] = CartesianIndex(i,Aj)
      coordLA[i,k] = CartesianIndex(i,LAj)

      #baseJ[coordA[i,k]] = ws[i,k]
      #baseJ[coordLA[i,k]] = -RLws[i,k]

      #add the row, column, and value of the Jacobian
      baseJvec.I[coordctr], baseJvec.J[coordctr],  baseJvec.V[coordctr] = i, Aj, ws[i,k]
      coordctr+=1
      baseJvec.I[coordctr], baseJvec.J[coordctr],  baseJvec.V[coordctr] = i, LAj, -RLws[i,k]
      coordctr+=1
    end
  end
  for i ∈ 1:length(Θ.ts)
    t=Θ.ts[i]
    V₀j = Jcolmap.V₀[t-1]
    coordV₀[i] =CartesianIndex(i,V₀j)
    #baseJ[coordV₀[i]] = T(1)

    baseJvec.I[coordctr], baseJvec.J[coordctr],  baseJvec.V[coordctr] = i, V₀j, T(1)
    coordctr+=1
  end


  baseJ = sparse(baseJvec.I, baseJvec.J, baseJvec.V)

  #now find the linear index among values within the sparse array
  nzs = findnz(baseJ)
  coords::Vector{Tuple} = collect(zip(nzs[1],nzs[2]))
  vals::Vector{T} = nzs[3]


  #maps the cartesian coordinates to the value row
  indexA::Matrix{Int} = Matrix{Int}(undef, Nexpand,K)
  indexLA::Matrix{Int} = Matrix{Int}(undef, Nexpand,K)
  indexV₀::Vector{Int} = Vector{Int}(undef, Nexpand)

  coordrow::Dict = Dict(zip(coords, 1:length(coords)))
  findcoord(ci) = coordrow[(ci[1],ci[2])]

  for i ∈ 1:Nexpand
    for k ∈ 1:K
      try
        indexA[i,k] = findcoord(coordA[i,k])
        indexLA[i,k] = findcoord(coordLA[i,k])
      catch err
        @info "i=$i, k=$k, coordA=$(coordA[i,k]), coordLA=$(coordA[i,k]), ws=$(ws[i,k])"
        error(err)
      end
    end
    indexV₀[i] = findcoord(coordV₀[i])
  end

  #verify the values are as expected
  print("Verifying baseJ...")

  @time begin
    for k ∈ 1:K
      @assert ((i,j)->baseJ[i,j]).(1:Nexpand, Jcolmap.A[k, Θ.ts]) == vals[indexA[:,k]]
      @assert ((i,j)->baseJ[i,j]).(1:Nexpand, Jcolmap.A[k, Θ.ts .- 1]) == vals[indexLA[:,k]]
    end
    @assert ws == vals[indexA]
    @assert -1 .* RLws == vals[indexLA]
    @assert ((i,j)->baseJ[i,j]).(1:Nexpand, Jcolmap.V₀[Θ.ts .- 1]) == vals[indexV₀]
  end

  index = (A=indexA, LA=indexLA, V₀=indexV₀)


  return index, baseJ
end

#front-loads and caches some of the computation
function JacobianCache(Θ::VolumePartsDLSq, X::NamedTuple, y::AbstractVector,::Type{T}) where T<:Real
  Jcolmap = jacobiancolmap(Θ)
  index, baseJ = baseleveljacobian(Θ,Jcolmap,X,y,T)

  return JacobianCache(Θ, Jcolmap, index, baseJ)
end

cudasignᵢ(Aᵢ, LAᵢ, wsᵢ, RLwsᵢ) = CUDA.sign(Aᵢ * wsᵢ - LAᵢ * RLwsᵢ)
function jacobian!(J, Jcache::JacobianCache, Θ::VolumePartsDLSq, X::NamedTuple, ::AbstractVector,
    ::Val{:level})::Nothing
  size(J) == (length(Θ.ts), length(Θ)) || error("incorrect dimensions $(size(out)) for J template")
  @unpack ws, RLws = X
  N::Int = size(Θ.A,2)
  K::Int = size(Θ.A,1)
  Nexpand::Int = length(Θ.ts)

  Jvals = nonzeros(J)
  baseJvals = nonzeros(Jcache.baseJ)

  #get the expand versions of all entries
  #A = Θ.expand.A
  #LA = Θ.expand.LA
  #V₀ = (exp).(Θ.expand.V₀)

  #get the signing of each entry
  signs = sign.(Θ.expand.A' .* ws .- Θ.expand.LA' .* RLws)
  Jvals[Jcache.index.A] .= baseJvals[Jcache.index.A] .* signs
  Jvals[Jcache.index.LA] .= baseJvals[Jcache.index.LA] .* signs
  Jvals[Jcache.index.V₀] .= exp.(Θ.expand.V₀)

  return nothing

end

function jacobian!(J, Jcache::JacobianCache, Θ::VolumePartsDLSq, X::NamedTuple, ::AbstractVector,
    ::Val{:levellog})::Nothing
  size(J) == (length(Θ.ts), length(Θ)) || error("incorrect dimensions $(size(out)) for J template")

  @unpack ws, RLws = X
  N::Int = size(Θ.A,2)
  K::Int = size(Θ.A,1)
  Nexpand::Int = length(Θ.ts)

  Jvals = nonzeros(J)
  baseJvals = nonzeros(Jcache.baseJ)

  A = exp.(Θ.expand.A')
  LA = exp.(Θ.expand.LA')

  #get the signing of each entry, and compute the derivitives
  signs = sign.(A .* ws .- LA .* RLws)
  Jvals[Jcache.index.A] .= baseJvals[Jcache.index.A] .* A .* signs
  Jvals[Jcache.index.LA] .= baseJvals[Jcache.index.LA] .* LA .* signs
  Jvals[Jcache.index.V₀] .= exp.(Θ.expand.V₀)

  return nothing

end

#returns an operating function which makes the jacobian, as well as a template of the matrix
function JVolumePartsDLSq(Θ::VolumePartsDLSq, X::NamedTuple, v::AbstractVector,::Type{T},
    PredictionType::Val;
    Jcache = JacobianCache(Θ, X, v, T)) where T<:Real


  #check the expand views
  J = deepcopy(Jcache.baseJ)

  #generate the function to be used
  J!(J, Θ) = jacobian!(J, Jcache, Θ, X, v, PredictionType)

  return J, J!
end


#this is a somewhat lighter weight version of the above for more regular testing
function testjacobian(Θtest::VolumePartsDLSq, Xtest::NamedTuple, vtest::AbstractVector,
  ::Type{T}=PARAM[:fluxtype],
  ::Type{TV}=PARAM[:fluxgpu] ? CUDA.CuVector{T} : Vector{T},
  ::Type{TM}=PARAM[:fluxgpu] ? CUDA.CuMatrix{T} : Matrix{T},
  PredictionType::Val = Val(PARAM[:predictiontype]);
  tol::Float64 = PARAM[:testtol],
  J::Union{AbstractMatrix} = error("J prototype matrix is required]"),
  J! = error("J! is required"),
  dims::NamedTuple = error("dims is required!"),
  )where {
    T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  @info "***Verifying Jacobian procedure..."

  if PARAM[:fluxmodel] ≠ :levellsq
    error("PARAM[:fluxmodeltype]≠:levellsq is not coded for use with Jacobian")
    return nothing
  end
  PARAM[:predictiontype] ∈ [:level, :levellog] || error(
    "invalid predictiontype $PredictionType for model type $(PARAM[:fluxmodeltype])")

  CUDA.allowscalar(false)

  #create the verification object
  Xver = (ws=fluxgpu(Xtest.ws,T), RLws=fluxgpu(Xtest.RLws,T))
  vver::TV = fluxgpu(vtest,T)
  Θver = VolumePartsLSq(dims.T, dims.K, Θtest.ts, T)
  Θver .= Θtest

  hatver = Vector{T}(Θver(Xver, vver, PredictionType))
  hattest = Θtest(Xtest, vtest, PredictionType)
  sum(abs.(hatver .- hattest)) < tol || error("prediction deviation (not jacobian) exceeds tolerance")

  @assert Θtest.Π == Θver.Π
  @assert Θtest.A == Θver.A
  @assert Θtest.V₀ == Θver.V₀
  @assert Θtest.expand.A ≈ view(Θver.A, Θver.expand.Aidx...)
  @assert Θtest.expand.LA ≈ view(Θver.A, Θver.expand.LAidx...)
  @assert Θtest.expand.V₀ ≈ view(Θver.V₀, Θver.expand.V₀idx)

  Π=Flux.params(Θver)
  length(Π) == 0 && error("No parameters noticed by Flux!")

  loss(X::NamedTuple,y::AbstractVector) = sum(Θver(Xver, vver, PredictionType))

  #now get the gradient for verification
  ∇Θ! = ∇VolumeParts!(loss, Xver,vver, Θver)
  ∇ver = similar(Θver.Π)
  ∇Θ!(∇ver, Θver)

  J!(J, Θtest)

  ∇test = vec(sum(J, dims=1))

  Δ = sum(abs.(∇test .- ∇ver))
  if Δ > tol
    #@info "∇test: $∇test"
    #@info "∇ver: $∇ver"
    @info "listing differences..."
    ctr::Int = 0
    maxctr::Int = 2
    for i ∈ 1:length(∇test)
      if !(∇test[i] ≈ ∇ver[i]) || (abs(∇test[i] - ∇ver[i]) > tol)
        @info "!: ∇test[$i], ∇ver[$i]=($(∇test[i]), $(∇ver[i])); Δ=$(∇test[i] - ∇ver[i])"
        ctr+=1;
      end
      if ctr ≥ maxctr
        @info "max $ctr reached"
        break
      end
    end
    error("FAILED jacobian test with Δ=$Δ")
  else
    @info "passed jacobian test with Δ=$Δ"
  end

end


function testvolumepartsindicesdirectlsq(::Type{Tvp}=PARAM[:opttype],
  ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{Tvp} : Vector{Tvp},
  ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{Tvp} : Matrix{Tvp}
    ) where {Tvp<:Real, TV<:AbstractVector{Tvp}, TM<:AbstractMatrix{Tvp}}
  Θ = VolumePartsDLSq(200,5, rand(2:200, 10^5))
  Θcopy = deepcopy(Θ)

  #first test getidx
  (Θ.A == reshape(Θ[1:1000], :, 200)) || error("inconsistent getidx dimensions Θ.A")
  (Θ.V₀ == Θ[1001:1199]) || error("inconsistent getidx dimensions Θ.V₀")

  #now test setidx
  testvals = rand(PARAM[:fluxtype], 1199)
  testidx = Random.randperm(1199)

  @assert allunique(testidx)
  Θ[testidx] .= testvals

  testcombined = hcat(testidx, testvals)
  testcombined = testcombined[sortperm(testcombined[:,1]),:]
  orderedtestvals = Vector(vec(testcombined[:,2]))
  @assert sum(orderedtestvals .== Vector(Θ.Π)) == length(Θ.Π)

  Θcopy.A .= TM(reshape(orderedtestvals[1:1000], size(Θcopy.A)))
  Θcopy.V₀ .= TV(orderedtestvals[1001:1199])

  (Θcopy.A ≈ Θ.A) || error("inconsistent setidx! dimensions Θ.A")
  (Θcopy.V₀ ≈ Θ.V₀) || error("inconsistent setidx! dimensions Θ.V₀")


  @info "test vol parts indices completed successfully"
end

#convenience constructor which creates a verificaiton object prior to performing the test
function testprediction(Θtest::VolumePartsDLSq,
  X::NamedTuple, y::AbstractVector,
  PredictionType::Val,
  ::Type{T}=PARAM[:opttype]) where T


  NK, NT = size(Θtest.A)
  Θver = VolumePartsLSq(NT, NK, Θtest.ts, T)
  Θver.Π .= Vector(Θtest.Π)


  yhattest = Θtest(X,y,PredictionType) |> Vector
  yhatver = Θver(X,y,PredictionType) |> Vector
  rmserror = sqrt(mean((yhattest .- yhatver).^2))

  if PARAM[:testmeasure]
    @info "Performance of DLSq objective"
    @btime $Θtest($X,$y,$PredictionType)

    @info "Performance of baseline LSq objective"
    @btime $Θver($X,$y,$PredictionType)
  end

  if yhattest ≈ yhatver
    @info "objective test passed with rmserror=$rmserror"
  else
    error("objective test failed!!! rmserror=$rmserror")
  end

  return nothing
end

#NOTE- test make sure this settup works
function testdlsq(panel::AbstractDataFrame,  ms::MeasureSpec,
    ::Type{T}=PARAM[:opttype],
    ::Type{TV}=PARAM[:optgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:optgpu] ? CUDA.CuMatrix{T} : Matrix{T};
    PredictionType::Val = Val(PARAM[:predictiontype]),
    tol::Float64 = PARAM[:testtol]
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  if PARAM[:fluxmodel] ≠ :levellsq
    error("PARAM[:fluxmodeltype]≠:levellsq is not coded for use with Jacobian")
  end

  testvolumepartsindicesdirectlsq()

  Fw::Vector{Symbol} = Fweights(ms, :Fw)
  FRLw::Vector{Symbol} = Fweights(ms, :FRLw)
  essentialfields = [Fw; FRLw; ms.Fvol; ]

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, ms,T, TV, TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=ms.Fvol)
  dims = p.dims
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  ts::Vector{Int} = p.ts
  Xtest = (ws=ws, RLws=RLws)
  vtest::TV = p.dat[:Fvol]

  Θtest = VolumePartsDLSq(dims.T, dims.K, ts, T)

  if PARAM[:optgpu]
    @info "gpu usage on- derivitives will not work with DLSq+GPU. Skipping Jacobian test"
  else
    J, J! = JVolumePartsDLSq(Θtest, Xtest, vtest, T, PredictionType)
    testjacobian(Θtest, Xtest, vtest, J=J, J! = J!, dims=dims)
  end

  return testprediction(Θtest, Xtest, vtest, PredictionType)
end
