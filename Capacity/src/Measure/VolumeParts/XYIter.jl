#=using Revise
using LinearAlgebra, CUDA, BenchmarkTools

##WARNING- comment out when live
include("../Annihilator.jl")=#



#creates crosssections of the data and working arrays
#note this is now only used in testing- use the annihilator matrix instead of the below
function dataxsection(D::NamedTuple, ts::Vector{Int},
  ::Type{TM}, ::Type{TV}, ::Type{T}=eltype(TM),
  ) where {
    T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}

  #technically we don't need this, but if its not sorted the process if much harder to reverse
  @assert issorted(ts)
  @assert sum(ts.==1) == 0
  allts = sort!(unique([1; ts...;]))
  dims= (T = length(allts), Nexpanded=length(ts))
  @assert collect(1:dims.T) == allts #validate integrity over time

  Dkeys = keys(D)
  Dtypes = (typeof).(values(D))

  xsections = Vector{NamedTuple}(undef, dims.T-1)
  allinds = collect(1:dims.Nexpanded)
  for t ∈ 1:(dims.T-1)
    #we build the named tuple by first creating a dictionary, then converting to a named tuple
    xinds = allinds[ts.==t+1] #cache this since we will call it a lot

    xsections[t] = (; broadcast(zip(Dkeys,Dtypes)) do (k,DT) #form a named tuple for each time t
      if DT <: TM
        k=>TM(D[k][xinds,:])
      elseif DT <: TV
        k=>TV(D[k][xinds])
      else
        @assert false
      end
    end...)
  end

  #now make a tuple of the data fields, where each field is a vector for all t
  xsection = (; [k=>(x->x[k]).(xsections) for k ∈ Dkeys]...)

  return xsection
end


abstract type AbstractXYIter{TM<:AbstractMatrix, TV<:AbstractVector, T<:Real} end
abstract type AbstractXYIterControl{TM<:AbstractMatrix, TV<:AbstractVector,
   T<:Real, TxM<:Union{Nothing, Vector{<:AbstractAnnhilator}}} <: AbstractXYIter{TM,TV,T} end

const XYIterControlLevel = AbstractXYIterControl{
  <:AbstractMatrix, <:AbstractVector, <:Real, <:Nothing}
const XYIterControlDemean= AbstractXYIterControl{
  <:AbstractMatrix, <:AbstractVector, <:Real, <:AbstractVector}

function findnonsingularcols(Mraw::AbstractMatrix)
  dims = (rows = size(Mraw,1), cols = size(Mraw,2))
  rank(Mraw) == dims.cols && throw("No singular cols to drop!")

  Mcols = [true; falses(dims.cols-1)]
  for c ∈ 2:(dims.cols)
    Mcols[c] = true
    sM = view(Mraw, :, Mcols)
    Mcols[c] = rank(sM) == sum(Mcols)
  end
  sM = view(Mraw, :, Mcols)
  @assert rank(sM) == size(sM,2) == rank(Mraw)
  return Mcols
end

rankrt(M::AbstractMatrix{T}) where T = rank(M, eps(T)^0.5)

function findnonsingularcols2(Mraw::AbstractMatrix)
  dims = (rows = size(Mraw,1), cols = size(Mraw,2))
  isposdef!(Symmetric(Mraw' * Mraw)) &&  (rankrt(Mraw) == dims.cols) && throw(
    "No singular cols to drop!")

  Mcols = [true; falses(dims.cols-1)]
  for c ∈ 2:(dims.cols)
    Mcols[c] = true
    M = Mraw[:, Mcols]
    Mcols[c] = (rankrt(M) == sum(Mcols)) && isposdef!(Symmetric(M'*M))
  end
  M = Mraw[:, Mcols]
  @assert rankrt(M) == size(M,2)
  @assert isposdef!(Symmetric(M'*M))
  return Mcols
end

#this function groups the controls as well as checking the matrices for degeneracy
function groupcontrolsbyt(Wgroup::TM, ts, ::Type{T}=eltype(TM))where {TM, T}
  xWgroupraw = dataxsection((xWgroup=Wgroup |> Matrix,), ts, Matrix{Float64}, Vector{Float64}, T)

  #checking for control independence
  xWgroupcpu = similar(xWgroupraw.xWgroup)
  KW::Int = size(Wgroup,2)
  droppedcols::Vector{Int} = zeros(Int, KW)
  print("Checking for control column independence t∈2:$(length(xWgroupcpu)).
    Singular cols found, if any, at t=[")
  Wgroupcolsused = trues((maximum(ts)-1)*KW) #this is an index of the control cols
  for (tM1, Wraw) ∈ enumerate(xWgroupraw.xWgroup)
    if (rankrt(Wraw) == size(Wraw,2)) && isposdef!(Symmetric(Wraw' * Wraw))
      xWgroupcpu[tM1] = Wraw
      @assert isposdef!(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1])) "
        matrix: $(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1]))"
    else
      nonsingular = view(Wgroupcolsused,((tM1 - 1)*KW+1):(tM1*KW))
      @assert (nonsingular .=== true) |> all
      nonsingular .= findnonsingularcols2(Wraw)
      print("$(tM1+1), ")
      droppedcols .+= (!).(nonsingular)
      #@info "Singular control columns found in group at t=$(tM1+1)." *
      #  "Dropped cols $((1:(length(nonsingular)))[(!).(nonsingular)])"
      xWgroupcpu[tM1] = Wraw[:, nonsingular]
      @assert isposdef!(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1])) "
        nonsingular: $nonsingular
        matrix: $(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1]))"

    end

    @assert isposdef!(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1])) "
      matrix: $(Symmetric(xWgroupcpu[tM1]' * xWgroupcpu[tM1]))"
  end
  println("]", sum(droppedcols) > 0 ?
    reduce(*, ["\ncol c dropped (ntimes): ";
      ((c)->"$(c[1]) ($(c[2])), ").(enumerate(droppedcols));]) : "")

  #convert the type if necessary
  xWgroup = (xWgroup = TM === Matrix{Float64} ? xWgroupcpu : TM.(xWgroupcpu),)

  return (;xWgroup, Wgroupcolsused)
end

#this creates the expanded control matrix
function expandcontrolsbyt(Wgroup::TM, ts, ::Type{T}=eltype(TM)) where {TM,T}

  xWgroup, Wgroupcolsused = groupcontrolsbyt(Wgroup, ts) |> (nt)->(nt.xWgroup, nt.Wgroupcolsused)
  K = size(Wgroup,2)
  N = length(ts)

  expanded = similar(Wgroup, N, K * length(ts |> unique))
  #expanded = similar(Wgroup, N, sum(Wgroupcolsused))
  expanded .= 0.0

  rctr = 1
  for (tM1, xWgroupt) ∈ enumerate(xWgroup.xWgroup)
    Nt = size(xWgroupt,1)

    #be careful here to not capture non-singular columns
    #proceed in two steps
    submatrix = view(expanded, rctr:(rctr + Nt-1), (K*(tM1 - 1)+1):(tM1*K))
    subcolumnsused = view(Wgroupcolsused, (K*(tM1 - 1)+1):(tM1*K))
    submatrix = view(submatrix, :, subcolumnsused)

    @assert length(subcolumnsused) == K
    @assert all(submatrix .== 0.0)
    submatrix .= xWgroupt
    rctr += Nt
  end

  #now run some checks
  @assert rctr == N+1
  #@assert (sum(expanded, dims=2) .≈ sum(Wgroup, dims=2)) |> all
  nonsingularexpanded = expanded[:, Wgroupcolsused]
  checkcolsums = reduce(vcat, (xWgroupt -> sum(xWgroupt,dims=1) |> vec).(xWgroup.xWgroup))
  @assert (checkcolsums .≈ vec(sum(nonsingularexpanded, dims=1))) |> all
  checkrowsums = reduce(vcat, (xWgroupt -> sum(xWgroupt,dims=2) |> vec).(xWgroup.xWgroup))
  @assert (checkrowsums .≈ vec(sum(nonsingularexpanded, dims=2))) |> all

  return (;expanded=nonsingularexpanded, Wgroupcolsused)
end


#holds all unchanging data as well as helper pre-allocations
#this version also includes the projection matrix for demeaning
struct XYIterControl{TM<:AbstractMatrix,
    TV<:AbstractVector,
    T<:Real,
    TxM,
    Txsection} <: AbstractXYIterControl{TM,TV,T, TxM}
  ws::TM
  Rtws::TM
  RLws::TM
  v::TV
  ṽ::TV #ṽ = xM * W

  W::TM
  W̃::TM #W̃ = xM * W

  xsection::Txsection

  #constructor from minimum componeents assuming a Wgroup cross-secitonal component
  function XYIterControl(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts, W::TM,  Wgroup::TM,
      ::Type{T}=eltype(TM)) where {TM<:AbstractMatrix,TV, Tts<:AbstractVector, T}

    issorted(ts) || error("ts needs to be sorted, otherwise row orderwing will be misalligned
      due to the reduce->vcat")

    #no need for anything fancy if there are no cross-secitonal controls
    if size(Wgroup,2) == 0
      smallxsection = (;xM=nothing, dataxsection((;xws=ws, xRtws=Rtws,xRLws=RLws), ts, TM, TV, T)...)
      return new{TM,TV, T, Nothing, typeof(smallxsection)}(ws, Rtws, RLws, v, v, W, W, smallxsection)
    end

    #(size(Wgroup,2) == 0) && return new{TM,TV, T, Nothing}(ws, RLws, v, v, W, W,nothing)

    local W̃::TM

    #construct a compact form of the annihilator matrix
    #only pre-multiplication of the matrix on a matrix or vector is supported
    #start by sectioning off components. These will be discarded

    #Form the grouped projection matrices
    xv = dataxsection((xv=v,), ts, TM, TV, T)
    xWgroup = groupcontrolsbyt(Wgroup, ts).xWgroup
    xM = (xWgroupt->DenseAnnihilator(xWgroupt)).(xWgroup.xWgroup)

    #project W onto the columnspace of Ws to form W̃. Similarly project v onto Ws
    ṽ = reduce(vcat, ((vt,Mt)->Mt*vt).(xv.xv, xM))

    #if size(W,2) > 0
    xW = dataxsection((xW=W,), ts, TM, TV, T)
    W̃ = reduce(vcat, ((Wt,Mt)->Mt*Wt).(xW.xW, xM))
    #=else
      W̃ = W
    end=#

    #now bundle the relevant components together for easy group calculations
    D = (xws=ws, xRtws=Rtws, xRLws=RLws)
    xsection = (;xM=xM, dataxsection(D, ts, TM, TV, T)...)

    #a quick integrity check
    xsectioncheck= dataxsection((;xṽ=ṽ, xW̃=W̃), ts, TM, TV, T)
    for (xvt, xṽt, xWt, xW̃t, xMt) ∈ zip(
        xv.xv, xsectioncheck.xṽ, xW.xW, xsectioncheck.xW̃, xsection.xM)


      #=if !((xṽt |> Vector) ≈ (
        xvt - xMt.W * (cholesky(Symmetric(xMt.W' * xMt.W))\(xMt.W' * xvt))) |> Vector)
        println("ERROR !((xṽt |> Vector) ≈ (
          xvt - xMt.W * (cholesky(Symmetric(xMt.W' * xMt.W))\\(xMt.W' * xvt))) |> Vector)")
        println("xvt:")
        printmln(xṽt)
        println("xvt - (...* xvt)")
        printmln((
          xvt - xMt.W * (cholesky(Symmetric(xMt.W' * xMt.W))\(xMt.W' * xvt))) |> Vector)
        printmln("xvt: ")
        printmln(xvt)
        error("stop")
      end=#

      heuristicatol = mean(abs.(xvt)*1e-4)
      heuristicisapprox(x,y) = isapprox(x,y,atol=heuristicatol)
      @assert heuristicisapprox(xṽt |> Vector,
        (xvt - xMt.W * (cholesky(Symmetric(xMt.W' * xMt.W))\(xMt.W' * xvt))) |> Vector)# |> all
      @assert heuristicisapprox(xW̃t |> Matrix,
        (xWt - xMt.W * (cholesky(Symmetric(xMt.W' * xMt.W))\(xMt.W' * xWt))) |> Matrix)# |> all
    end

    return new{TM,TV, T, typeof(xM), typeof(xsection)}(ws, Rtws, RLws, v, ṽ, W, W̃, xsection)
  end
end

XYIterControl(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts, ::Type{T}=eltype(TM);
    W::TM = similar(ws, size(ws,1),0),
    Wgroup::TM = similar(ws, size(ws,1),0), #empty matrix
    ) where {TM,TV, Tts, T} = XYIterControl(ws, Rtws, RLws,  v, ts, W, Wgroup, T)


#entry point for default constructor of XYIter object
XYIter(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts,;
  PredictionType::Val=Val(PARAM[:predictiontype])
  ) where {TM, TV, Tts} = XYIter(ws, Rtws, RLws, v, ts, PredictionType)

XYIter(ws::TM, Rtws::TM, RLws::TM,  v::TV, Z::TM, ts::Tts, ;
  PredictionType::Val=Val(PARAM[:predictiontype])
  ) where {TM, TV, Tts} = XYIter(ws, Rtws, RLws, v, Z, ts, PredictionType)

XYIter(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts, ::Val{:leveldemean}
  ) where {TM, TV, Tts} = XYIterDemean(ws, Rtws, RLws, v, ts)

XYIter(ws::TM, Rtws::TM, RLws::TM,  v::TV, ts::Tts, ::Val{:level}
  ) where {TM, TV, Tts} = XYIterLevel(ws, Rtws, RLws, v, ts)
