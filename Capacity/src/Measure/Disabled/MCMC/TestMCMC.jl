
#builds mcmc xsections, currently only used in testing
function xsectionmcmc(Θ::VolumePartsMCMC,::Type{TM}, ::Type{TV}, ::Type{T}) where {TM,TV,T}
  @unpack Xv,ts = Θ

  #=D = (ws=Xv.ws', RLws=Xv.RLws', Ã=Θ.expand[:Ã])
  xsection = dataxsection(D,Θ.ts,AbstractMatrix{T}, AbstractVector{T}, T)
  LD = G=Θ.expand[:G],σ²=Θ.expand[:σ²], ω=Θ.expand[:ω]=#

  dims = (T=size(Θ.Ã,2), K=size(Θ.Ã,1))
  xsections = Vector{NamedTuple}(undef, dims.T-1)
  for t ∈ 2:dims.T
    #count the number of stocks in each window (t and t-1)
    Nt = sum(ts.==t)

    #repeat t for each stock
    #Ainds = repeat([t],Nt)
    #LAinds = repeat([t-1], Nt)
    xinds = ts .== t
    xsections[t-1] = (
      ws = view(Θ.Xv.ws, xinds, :),
      RLws = view(Θ.Xv.RLws, xinds, :),
      v = view(Θ.Xv.v, xinds),
      #X̃ = view(Θ.X̃, xinds, :),
      #ṽ  = view(Θ.ṽ, xinds),
      Ã = view(Θ.expand[:Ã], :, xinds),
      LÃ = view(Θ.expand[:LÃ], :, xinds),
      G = view(Θ.expand[:G], :, xinds),
      σ²   = view(Θ.expand[:σ²], xinds),
      ω = view(Θ.expand[:ω], xinds),
      )
  end
  xsection::NamedTuple = (
    ws = (x->x.ws).(xsections),
    RLws = (x->x.RLws).(xsections),
    v = (x->x.v).(xsections),
    #X̃ = (x->x.X̃).(xsections),
    #ṽ = (x->x.ṽ).(xsections),

    Ã = (x->x.Ã).(xsections),
    LÃ = (x->x.LÃ).(xsections),
    G = (x->x.G).(xsections),
    σ²   = (x->x.σ²).(xsections),
    ω = (x->x.ω).(xsections),
  )

  return xsection
end


function testcapacitychains()
  hyper=(A = (a=1.1,), V=(v1=2.1, v2=2.2), M=(m1=rand(2,2), m2=3.2))
  testparams = (A = 2.3, V=rand(5), M=rand(2,5))
  cc = CapacityChains(Float64, Nburn=25, Ndraws=10, paramstemplate=testparams, hyper=hyper)
  Nparams = 16
  @assert size(cc) === (10,16)

  ###test the indices
  #prior to looping
  for i ∈ 1:10
    @assert i == cc.idx
    @assert i == cc.i
    increment!(cc)
  end

  #after the first loop of the buffer
  for i ∈ 1:15
    @assert (i-1)%10 + 1 == cc.idx
    @assert i + 10 == cc.i
    increment!(cc)
  end

  #after the first reset
  @assert cc.idx == 1
  @assert cc.i == 26
  @assert cc.previdx ==10

  ###test the state vector
  Ψ::CCState = CCState(cc)
  v = rand(Nparams)
  Ψ .= v
  @assert all(v .== Ψ)
  updateandincrement!(cc, Ψ)
  @assert all(cc.value[cc.previdx,:,1] .== Ψ)
  @assert cc.i == 27
  @assert cc.idx == 2
  @assert cc.previdx == 1
  for (n,p) ∈ pairs(Ψ.params) #test each params
    @assert all(reshape(cc.params[n][cc.previdx,:,1],size(p)) .== p)

    #=(all(reshape(cc.params[n][cc.previdx, :, 1],size(p)) .== p)) || error(
      "mismatch!!! cc.idx: $(cc.idx), previdx: $(cc.previdx)
      reshape(cc.params[$n][$(cc.previdx), :, 1],size(p)): $(reshape(cc.params[n][:,cc.previdx,1],size(p)))
      p: $(p)")=#
  end

  #test the params individually
  Ψ.A[] = testparams.A
  Ψ.V .= testparams.V
  Ψ.M .= testparams.M
  @assert all([(v->[v...;]).(values(testparams))...;] .== Ψ)
  updateandincrement!(cc,Ψ)
  @assert (cc.i, cc.idx, cc.previdx) === (28,3,2)
  for (n,p) ∈ pairs(Ψ.params) #test each params
    if n ≡ :A
      @assert reshape(cc.params[:A][cc.previdx,:,1],size(p))[] === testparams.A
    elseif n ≡ :V
      @assert all(reshape(cc.params[:V][cc.previdx,:,1],size(p)) .=== testparams.V)
    elseif n ≡ :M
      @assert all(reshape(cc.params[:M][cc.previdx,:,1],size(p)) .=== testparams.M)
    else
      @assert false
    end
  end

end

function testA₀(Θ::AbstractVolumePartsMCMC{TM,TV,T},
  Ntest::Int=PARAM[:mcmcntest]) where {TM,TV,T}

  @unpack Ψ, Xv, xsection, Ã = deepcopy(Θ)
  @unpack A₀, G, ω, σ² = Ψ
  #@info "hyper keys=$(keys(Θ.cc.hyper))"
  @unpack μA₀, VA₀inv = Θ.cc.hyper[:A₀]

  dims = (K=size(Ã, 1), T=size(Ã, 2), )

  A₀ver = Matrix(undef, Ntest, dims.K)
  A₀test = similar(A₀ver)

  function verA₀!()
    Ã .= hcat(ones(T, dims.K), cumprod(G, dims=2))

    #store the cross-sectional matrices
    #Λₜ = Vector{TM}(undef, dims.T-1)
    #Āₜ = Vector{TV}(undef, dims.T-1)
    VA₀invmat = Matrix(VA₀inv, dims.K, dims.K)
    μA₀vec = fill(μA₀, dims.K)
    Λₜ = [Matrix{T}(undef, dims.K, dims.K) for t ∈ 1:(dims.T-1)]
    ΛĀₜ = [Vector{T}(undef, dims.K) for t ∈ 1:(dims.T-1)]
    τₜ = 1 ./ σ²
    Threads.@threads for t ∈ 1:(dims.T-1)
      X = abs.(xsection[:ws,t] .* xsection[:Ã,t]' .- xsection[:RLws,t] .* xsection[:LÃ,t]')
      XtX = X' * X
      ṽ = xsection[:v,t] .- xsection[:ω,t]
      Â = cholesky(XtX)\(X'*ṽ)
      ΛĀₜ[t] .= τₜ[t] .* (XtX*Â .+ VA₀inv*μA₀vec)
      Λₜ[t] .= τₜ[t] .* (XtX .+ VA₀invmat)

      #=if !(Λₜ[t] ≈ Λₜ[t]') #fixes some roundoff error
        println("XtX: ")
        printmln(XtX)
        println("Λₜ[t]: ")
        printmln(Λₜ[t])
        error("$t: something is wrong with Λₜ[t]: !(Λₜ[t] ≈ Λₜ[t]')")
      end=#
    end

    ΛĀ = sum(ΛĀₜ)
    Λ = sum(Λₜ)

    Σ = Λ\I

    if Σ ≈ Σ' #fixes some roundoff error
      Σ .= (Σ .+ Σ') ./ 2
    else
      printmln(Σ)
      error("something is wrong with Σ: !(Σ ≈ Σ')")
    end

    Ā = Σ*ΛĀ
    #@info sum(abs.(Σ .- Σ'))
    d = MvNormal(Ā, Σ)
    return rand(d)

    #=for k ∈ shuffle(1:dims.K)
      notk = k .≠ 1:dims.K
      Σkkinv = cholesky(Σ[notk,notk])\I

      μk = Ā[k] + Σ[k:k, notk] * Σkkinv * (A₀[notk] .- Ā[notk])=#

  end

  #verA₀!()
  #@assert false

  if PARAM[:mcmcbenchtest]
    print("time updateA₀!: ")
    @btime updateA₀!($Θ)
    print("time verA₀!: ")
    @btime $verA₀!()
  end

  for i ∈ 1:Ntest
    A₀test[i,:] .= updateA₀!(Θ)
    A₀ver[i,:] .= verA₀!()
  end

  #halfNtest = Ntest ÷ 2
  #A₀vermean = mean(A₀ver[halfNtest:end, :], dims=1)
  #A₀testmean = mean(A₀test[halfNtest:end, :], dims=1)

  A₀vermean = mean(A₀ver, dims=1)
  A₀testmean = mean(A₀test, dims=1)

  Δ = abs.(A₀vermean .- A₀testmean)
  @info "
    A₀vermean: $A₀vermean
    A₀testmean: $A₀testmean
    Δ: $Δ"

end

function testposteriors(Θ::AbstractVolumePartsMCMC)
  testA₀(Θ)
end

#we use a lot of views in VolumeParts- this tests them
function testmcmcexpand(Θ::AbstractVolumePartsMCMC{TM,TV,T}, xsectionver) where {TM, TV, T}
  dims = (K=size(Θ.Ã,1), T=size(Θ.Ã,2))
  @unpack ts=Θ

  #these sort the views as though they were first divided into xsections
  function sorttvert(M::AbstractMatrix, ts)
    allts = sort!(unique(ts))
    reduce(vcat, (t->view(M, ts .== t, :)).(allts))
  end
  function sorttvert(M::AbstractVector, ts)
    allts = sort!(unique(ts))
    reduce(vcat, (t->view(M, ts .== t)).(allts))
  end

  function sortthoriz(M::AbstractMatrix, ts)
    tst = ts'
    allts = sort!(unique(ts))
    reduce(hcat, (t->view(M, :, vec(tst .== t))).(allts))
  end

  #concatenate the cross-sections in xsection to verify the Θ.expand views
  vdata = [:ws, :RLws, :v, :σ², :ω, ]
  hdata = [:Ã, :LÃ, :G, ]
  for s ∈ [vdata; hdata]
    local catfunc
    local sortfunc
    if s ∈ vdata #check if we arranged the data vertically or horizontally
      catfunc = vcat
      sortfunc = sorttvert
    elseif s ∈ hdata
      catfunc = hcat
      sortfunc = sortthoriz
    else
      @assert false
    end

    #check the complete expanded view
    vermatrix = reduce(catfunc, xsectionver[s])
    testmatrix = sortfunc(Θ.expand[s], ts)
    testxmatrix = reduce(catfunc, (t->Θ.xsection[s,t]).(1:(dims.T-1)))
    #@assert size(vermatrix) === size(testmatrix) === size(testxmatrix)
    all(vermatrix .== testmatrix .== testxmatrix)||error(
      "!all(vermatrix .== testmatrix .== testxmatrix)
      vermatrix[1:3,1:3]: $(vermatrix[1:3,1:3])
      testmatrix[1:3,1:3]: $(testmatrix[1:3,1:3])
      testxmatrix[1:3,1:3]: $(testxmatrix[1:3,1:3])
      "
    )



  #concatenate the cross-sections in xsection to verify the Θ.expand views (now for each t)
    for t ∈ 1:(dims.T-1)

      partts = @view ts[ts .≥ t + 1]
      vermatrix = reduce(catfunc, xsectionver[s][t:(dims.T-1)])
      testmatrix = sortfunc(Θ.expand[s,t], partts)
      testxmatrix = reduce(catfunc, (t->Θ.xsection[s,t]).(t:(dims.T-1)))
      @assert size(vermatrix) === size(testmatrix) === size(testxmatrix)
      @assert all(vermatrix .== testmatrix .== testxmatrix)

      @assert all(Θ.xsection[s,t] .== xsectionver[s][t])
    end
  end
return nothing
end

function testmcmc(panel::AbstractDataFrame,  zs::ZSpec,
    ::Type{T}=PARAM[:mcmctype],
    ::Type{TV}=PARAM[:mcmcgpu] ? CUDA.CuVector{T} : Vector{T},
    ::Type{TM}=PARAM[:mcmcgpu] ? CUDA.CuMatrix{T} : Matrix{T};
    prior::Symbol=PARAM[:mcmcprior],
    Nburn::Int=PARAM[:mcmcnburn],
    Ndraws::Int=PARAM[:mcmcndraws]
    ) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  testcapacitychains()

  @unpack Zms = zs
  Fw::Vector{Symbol} = Fweights(Zms, :Fw)
  FRLw::Vector{Symbol} = Fweights(Zms, :FRLw)

  #form the data into a more usable format and expand the arrays
  p = modelparts(panel, zs, T, TV, TM,
    weightdata = (Fw=Fw, FRLw=FRLw),
    Fvol=Zms.Fvol)

  dims = p.dims
  ws::TM = p.dat[:Fw]
  RLws::TM = p.dat[:FRLw]
  v::TV = p.dat[:Fvol]
  ts::Vector{Int} = p.ts
  Ascale::Vector{T} = p.Ascale
  Nexpanded::Int = length(ts)
  #this creates the factored objects- which are accurate given that initially G=1
  absΔ(wsᵢ::CuArray,RLwsᵢ)=CUDA.abs(wsᵢ-RLwsᵢ)
  absΔ(wsᵢ,RLwsᵢ)=abs(wsᵢ-RLwsᵢ)


  #now test the level versions
  Xv = XYIter(ws ,RLws, v, ts, Ascale)
  hyper=PARAM[prior]
  initial=(
    A₀ = TV(rand(T, dims.K)),
    G = TM(ones(T, dims.K, dims.T-1)),
    σ² = rand(T, dims.T-1),
    ω = rand(T, dims.T-1)
    )
  cc = CapacityChains(T, Nburn=Nburn, Ndraws=Ndraws, paramstemplate=initial, hyper=hyper)

  Θ = VolumePartsMCMC(TM, TV, T; initial, cc, Xv, ts, )

  xsection = xsectionmcmc(Θ, TM, TV, T)

  testmcmcexpand(Θ,xsection)

  testposteriors(Θ)

#=  testG(Θ, Xvlevel, _Xv, xsectionlevel)
  @assert Θ.G |> Matrix! ≈ Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] |> Matrix!

  #WARNING WARNING WARNING uncomment below!!!
  testprediction(Θ, Xvlevel, deepcopy(_Xv))
  testpredictionzygote(deepcopy(Θ),Xvlevel,deepcopy(_Xv))
  testiterloss(deepcopy(Θ),Xvlevel,deepcopy(_Xv))


  if !(TM<:CuArray)

    Xvdemean::XYIterDemean = XYIterDemean(ws ,RLws, v, ts, Ascale)
    Ddemean = (ws=Xvdemean.ws, RLws=Xvdemean.RLws, ṽ = Xvdemean.ṽ,)
    xsectiondemean = dataxsection(Ddemean,Θ.ts,TM, TV, T)
    validatexsection(Θ, Xvdemean, _Xv, p, xsectiondemean)

    testG(Θ, Xvdemean, _Xv, xsectiondemean)
    @assert Θ.G |> Matrix! ≈ Θ.A[:,2:end] ./ Θ.A[:,1:(end-1)] |> Matrix!

    testXv(Θ, Xvlevel, Xvdemean, _Xv)

    testprediction(Θ, Xvdemean, deepcopy(_Xv))
    testpredictionzygote(deepcopy(Θ),Xvdemean,deepcopy(_Xv))
    testiterloss(deepcopy(Θ),Xvdemean,deepcopy(_Xv))
  else
    #create cpu versions
    pcpu = modelparts(panel, zs,T, Vector{T}, Matrix{T},
      weightdata = (Fw=Fw, FRLw=FRLw),
      Fvol=Zms.Fvol)

    wscpu::Matrix{T} = pcpu.dat[:Fw]
    RLwscpu::Matrix{T} = pcpu.dat[:FRLw]
    vcpu::Vector{T} = pcpu.dat[:Fvol]

    Xvlevelcpu::XYIterLevel = XYIterLevel(wscpu ,RLwscpu, vcpu, ts, Ascale)
    Θcpu = VolumePartsIter(dims.T, dims.K, ts, T, Matrix{T}, Vector{T})
    Θcpu.G .= Θ.G |> Matrix
    Θcpu.A .= Θ.A |> Matrix
    _Xvcpu = XYIterAlloc(Θ, Xvlevelcpu)

    testcudazygote(deepcopy(Θ),Xvlevel,deepcopy(_Xv), Θcpu, Xvlevelcpu,  _Xvcpu)
  end=#


  return nothing
end
