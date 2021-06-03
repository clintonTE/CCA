

#general function to validate dates in a time series
#dates should already be unique and sorted
function validatedates(dates::Vector{Date};
    frequency::NSymbol=throw("frequency is required"), validateexactfrequency::Bool=false)
  @assert allunique(dates)
  @assert issorted(dates)

  (frequency === nothing) && return true

  #not much to do with daily data, we already checked the data was unique
  if Symbol(lowercase(string(frequency))) == :day
    (!validateexactfrequency) || throw("Exact validation of daily frequency not supported")
    return true

  #weekly data should have a minimum frequency of 7 days
  elseif Symbol(lowercase(string(frequency))) == :week
    dateseow = (lastdayofweek).(dates)
    #=datesyw = broadcast(dates) do d
      wd = week(d)
      if (wd==1) && (year(d) ≠ year(lastdayofweek(d)))
        return year(lastdayofweek(d)) .+ wd ./ 100
      elseif (wd ≥ 52) && (year(d) ≠ year(firstdayofweek(d)))
        return year(firstdayofweek(d)) .+ wd ./ 100
      end
      return year(d) .+ wd ./ 100
    end=#
    if !allunique(dateseow)
      nonuniquedateindices = nonunique(DataFrame(datesyw=dateseow))
      nonuniquedates = dates[nonuniquedateindices]
      throw("weeks and years of weekly data are not unique. Nonunique dates:
        $(collect(zip(dates, dateseow))[(d->d∈nonuniquedates).(dates)])")
    end
    if validateexactfrequency #make sure there are no missing or misaligned weeks
      validdates= firstdayofweek(minimum(dates)):Week(1):lastdayofweek(maximum(dates))
      validdateseow = lastdayofweek.(validdates)
      mismatches = dateseow .≠ validdateseow
      sum(mismatches)==0 || throw(
        "Failed to validate weekly dates with $(sum(mismatches)) mismatches.
        List of mismatches (given, minimum(dates):Week(1):maximum(dates)):
          $(collect(zip(dateseow, validdateseow))[mismatches])")
    end
  elseif Symbol(lowercase(string(frequency))) == :month
    datesym = (d->year(d) .+ month(d) ./ 100).(dates)
    allunique(datesym) || throw("months and years of monthly data are not unique")
    if validateexactfrequency #we don't insist on equality of dates, just year-month
      validdates= firstdayofmonth(minimum(dates)):Month(1):lastdayofmonth(maximum(dates))
      validdatesym = (d->year(d) .+ month(d) ./ 100).(validdates)
      try
        mismatches = datesym .≠ validdatesym
        sum(mismatches)==0 || throw(error("
          Failed to validate monthly dates with $(sum(mismatches)) mismatches.
          List of mismatches (given, minimum(dates):Month(1):maximum(dates)):
            $(collect(zip(datesym, validdatesym))[mismatches])"))
      catch err
        debugdf = vcat(DataFrame(dates=dates, datesym=datesym),
          DataFrame(validdates=validdates, validdatesym=validdatesym), cols=:union)
        debugdf |> CSV.write("$(PARAM[:testpath])\\validdatesym_dump.csv")
        throw(err)
      end
    end
  else
    throw("unrecognized date frequency $frequency")
  end
  return nothing
end

#used to clean up names in a dataframe
function cleanname(dirty::T)::T where T
  local str::String
  local clean::T #holds the resulting answer

  str = string(dirty)

  str = replace(str, "-"=>"")
  str = replace(str, " "=>"")
  str = replace(str, "_"=>"")
  str = lowercase(str)

  T==String && return str
  return T(clean)
end

#same as above but always returns a symbol
cleanpropertyname(dirty) = Symbol(cleanname(dirty))

parseormissing(::Missing, s::String) = missing

#parses strings to a type or makes missing
function parseormissing(::Type{T}, s::String) where T<:Any
  val::Union{T,Nothing} = tryparse(T, s)

  isnothing(val) && (return missing)
  return val
end

#candidate for Finometrics
function winsorizewithin!(df::DataFrame;
    Fsource::Symbol=error("Fsource is required"),
    prop::Float64=error("prop is required"),
    Fgroup::Symbol=error("Fgroup is required"),
    twosided=true,
    Fnew::Symbol = Symbol("WX$Ftarget"))

  df[!, Fnew] = similar(df[!, Fsource])
  sdfs = groupby(df, Fgroup)
  Threads.@threads for k ∈ keys(sdfs)
    sdf = sdfs[k]
    (sum((!ismissing).(sdf[!,Fsource])) == 0) && continue
    #we can check the winsorization was performed as expected by recording the percentiles
    oldptile = (quantile(skipmissing(sdf[!, Fsource]),prop),
      quantile(skipmissing(sdf[!, Fsource]), 1-prop))

    sdf[:, Fnew] .= winsorizequantile(sdf[!, Fsource], prop; twosided)

    #perform the check
    newminmax = (minimum(skipmissing(sdf[!, Fnew])), maximum(skipmissing(sdf[!, Fnew])))
    @assert all(oldptile .≈ newminmax) || all(oldptile .≈ reverse(newminmax))
  end
  return nothing
end

#cleans up the names after a call to aggregate
@inline nofunction(s::Symbol, todelete::String) = Symbol(replace(string(s), todelete=>""))
@inline nofunction(s::String, todelete::String) = replace(s, todelete=>"")

finiteormissing(::Missing) = missing
function finiteormissing(f::Real)
  (isfinite(f)) && return f
  return missing
end



#converts gzip to lz4
function gz2zstd(file::String; filetype::String = "csv", location::Union{String,Nothing} = nothing,
  gzextension::String = "gz", zstdextension::String = "zstd",
  deleteold::Bool = false, level=1)

  local prefix::String

  prefix = isnothing(location) ?  "" : "$location\\"
  infile::String = "$prefix$file.$filetype.$gzextension"
  outfile::String = "$prefix$file.$filetype.$zstdextension"

  isfile(infile) || error("file not found. Full path: $prefix$file.$filetype.$gzextension")


  #Below doesnt work
  indata = open(infile)
  outdata = open(outfile, "w")
  s = GzipDecompressorStream(ZstdCompressorStream(outdata, level=level))

  try
    write(s, indata)
  finally
    close(s)
  end
end

#use this to clear a dictionary
function clearcollection(c::AbstractDict)
  for i ∈ 1:length(c)
    pop!(c)
  end
end

#check if all columns in a df or subdf are equal
#extra logic is to handle missing values
@inline function eachrowequal(df::AbstractDataFrame;
  cols::Vector{Symbol}=propertynames(df))::Bool
  for c ∈ eachcol(view(df, :, cols))
    completec = collect(skipmissing(c))
    (length(completec) == 0) && continue #all missing case

    #check for equality
    if  (length(completec) == length(c)) && (all(c[1] .== c))
     continue
    end

    return false
  end
  return true
end

@inline cleandate(d::Date) = replace("$d", "-"=>"")

#a vectorized exponential function
#handy for cuda work
vecexp(v::CuVector) = (CUDA.exp).(v)
vecexp(v::AbstractVector) = (exp).(v)

function pressakey()
  print("press any key to continue")
  read(stdin,1)
  return nothing
end


#this may allow for gpu toggling
maybegpu(obj::AbstractVector, ::Type{T}, ::Val{true}) where T = CuVector{T}(obj)
maybegpu(obj::AbstractMatrix, ::Type{T},  ::Val{true}) where T = CuMatrix{T}(obj)
maybegpu(obj, ::Type{T}, ::Val{false}) where T = obj
fluxgpu(obj, ::Type{T}=PARAM[:fluxtype]) where T<:Real = maybegpu(obj, T, Val(PARAM[:fluxgpu]))


#compute the jacobian from the pullback
function jacobian(f, x, Π)
  _, back = Zygote.pullback(()->f(x),Π)

  return back(ones(length(x))) #a vector of 1s evaluates df/dx|x
end

#maps the indices of a shape to its corresponding linear indices in Π
function map2linear(start::Int, dims::Tuple)
  N::Int = prod(dims)
  lins::Vector{Int} = collect(1:N) .+ (start-1)
  linmap = reshape(lins, dims)

  #this should be true, as vec should make this an identity
  @assert vec(linmap) == lins
  return start+N, linmap
end


#this utility re-interprets the object(!)
#the object MUST remain referenced, or bad things will happen
#tries to do some basic safety checks, but its still possible to screw things up
function less_unsafe_wrap(::Type{Tdest};
  source::Vector = error("src vector is required"),
  start::Int = error("start is required"),
  dims::Tuple)::Tdest where Tdest<:AbstractArray

  Ndat = prod(dims)

  T = eltype(Tdest)
  (eltype(source) == T) || error("eltype of src must be of same type as T")

  #this has the effect of making sure we are getting space from somewhere that exists
  protowrap::SubArray = view(source, start:(start + Ndat -1))
  @assert length(protowrap) == Ndat

  p = pointer(source, start)
  wrap::Tdest = unsafe_wrap(Tdest, p, dims)

  #final checks
  @assert sum(vec(protowrap) .== vec(wrap)) == Ndat

  linmap = map2linear(start, dims)
  #return the object
  return wrap
end

#CUDA version of the above
#WARNING- I think something is wrong- consider not using this
function less_unsafe_wrap(::Type{Tdest};

  source::AbstractVector = error("src vector is required"),
  start::Int = error("start is required"),
  dims::Tuple)::Tdest where Tdest<:CuArray
  @warn("Possible race error in less_unsafe_wrap- function is potentially broken.
    Consider not using.")

  Ndat = prod(dims)

  T = eltype(Tdest)
  (eltype(source) == T) || error("eltype of src must be of same type as T")

  #this has the effect of making sure we are getting space from somewhere that exists
  protowrap = view(source, start:(start + Ndat -1))
  @assert length(protowrap) == Ndat

  #a final check
  p = CUDA.pointer(source, start)
  #@info typeof(p)
  wrap::Tdest = CUDA.unsafe_wrap(CuArray, p, dims)
  @assert sum(vec(protowrap|>Array) .== vec(wrap|>Array)) == Ndat

  return wrap
end

#nice printing of matrices
printm(m) = show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
function printmln(m)
  show(IOContext(stdout, :limit=>false), MIME"text/plain"(), m)
  println()
end

#if createa a labelling string if the dates are constricted
function boundrangestr(;
    lowerdatebound::NDate=PARAM[:lowerdatebound],
    upperdatebound::NDate=PARAM[:upperdatebound])
  lowerdatestr = lowerdatebound === nothing ?  "" : "_$(Dates.format(lowerdatebound,"yymd"))"
  upperdatestr = upperdatebound === nothing ?  "" : "_$(Dates.format(upperdatebound,"yymd"))"
  return "$lowerdatestr$upperdatestr"
end

#adds zeros at the end of a matrix
#used in hyperreg
function pad0(M::TM, trows::Int, ::Type{T} = eltype(TM)) where {TM<:AbstractMatrix, T}
  Mrows = size(M,1)
  (Mrows == trows) && return (view(M,:,:), deepcopy(M))

  padded = vcat(M, TM(zeros(T, trows-size(M,1), size(M,2))))
  sM = view(padded, 1:size(M,1), 1:size(M,2))
  return (sM, padded)
end
function pad0(V::TV, trows::Int, ::Type{T} = eltype(TV)) where {TV<:AbstractVector, T}
  NV::Int = length(V)
  (NV == trows) && return (view(V,:), deepcopy(V))

  padded = vcat(V, TV(zeros(T, trows-NV)))
  sV = view(padded, 1:NV)
  return (sV, padded)
end

Base.ones(::Type{T}, n::Int, ::Type{<:AbstractVector}) where {T} = ones(T,n)
Base.ones(::Type{T}, n::Int, m::Int, ::Type{<:AbstractMatrix}) where {T} = ones(T,n,m)
Base.ones(::Type{T}, dims::Tuple, ::Type{<:AbstractArray}) where {T} = ones(T,dims)
Base.ones(::Type{T}, n::Int, ::Type{<:CuVector}) where {T} = CUDA.ones(T,n)
Base.ones(::Type{T}, n::Int, m::Int, ::Type{<:CuMatrix}) where {T} = CUDA.ones(T,n,m)
Base.ones(::Type{T}, dims::Tuple, ::Type{<:CuMatrix}) where {T} = CUDA.ones(T,dims)

Base.show(io, M::CUDA.CUSPARSE.CuSparseMatrixCSR{Float32}) = Base.show(io,Matrix(M))
Base.show(io, V::CUDA.CUSPARSE.CuSparseVector{Float32}) = Base.show(io,Matrix(V))

#these are basicalyl convert methods that only create a copy if necessary
Matrix!(m::Matrix) = m
Matrix!(m::TM, ::Type{T}=eltype(TM)) where {TM<:AbstractMatrix, T} = m |> Matrix{T}

Vector!(v::Vector) = v
Vector!(v::TV, ::Type{T}=eltype(TV)) where {TV<:AbstractVector, T} = v |> Vector{T}
T!(obj::T,::Type{T}) where T = obj
T!(obj::U, ::Type{T}) where {U,T} = obj |> T


#expensive versions of cumprod(_,dims=2) to work around Zygote's bugs and quirks
#these is terrible and should go away once cumprod starts working with Zygote
function _cumprodbyrow(M::CuMatrix)
  cols = (i->CUDA.prod(M[:,1:i], dims=2)).(1:size(M,2))
  #CUDA.reduce(hcat, cols)
  CUDA.hcat(cols...)
end

_cumprodbyrow(M::Matrix) = hcat((i->prod(M[:,1:i], dims=2)).(1:size(M,2))...)

function cumprodbyrow(M::TM, ::Type{T}=eltype(M);
  limitA::Bool = PARAM[:iterlimitA],
  growthtype::Symbol = PARAM[:itergrowthtype],
  largelimit = T <: Float64 ? LARGE_VAL : LARGE_VAL32,
  smalllimit = T <: Float64 ? SMALL_VAL : SMALL_VAL32
  ) where {TM, T}

  if limitA
    if growthtype === :log
      return limitcumprodbyrow(M, largelimit, smalllimit)
    else
      return limitcumprodbyrow(M, largelimit)
    end
  end

  return _cumprodbyrow(M)
end


#bounds cumprod
function limitcumprodbyrow(M::TM, limit::T, ::Type{T}=eltype(TM);
   init=ones(T,size(M,1))) where {TM<:Matrix,T}

  function checklimit(x, y)
    xy = x * y
    ifelse(abs(xy) < limit, xy, limit*sign(xy))
  end

  cumprodcol(c1,c2) = map(checklimit, c1, c2)

  cols = eachcol(M) |> collect
  return reduce(hcat, Iterators.accumulate(cumprodcol,cols, init=init))
end

#bounds cumprod if we are working on a log scale
function limitcumprodbyrow(M::TM, largelimit::T, smalllimit::T, ::Type{T}=eltype(TM);
   init=ones(T,size(M,1))) where {TM<:Matrix,T}

  checklimit(x, y)= max(min(x * y, largelimit), smalllimit)
  cumprodcol(c1,c2) = map(checklimit, c1, c2)

  cols = eachcol(M) |> collect
  return reduce(hcat, Iterators.accumulate(cumprodcol,cols, init=init))
end

limitcumprodbyrow(M::TM, limit::T, ::Type{T}=eltype(TM)
  ) where {TM<:CuMatrix,TV, T} = (
    limitcumprodbyrow(M |> Matrix{T}, limit) |> TM)
limitcumprodbyrow(M::TM, largelimit::T, smalllimit::T, ::Type{T}=eltype(TM)
  ) where {TM<:CuMatrix,TV, T} = (
    limitcumprodbyrow(M |> Matrix{T}, largelimit, smalllimit) |> TM)
#testing for cumprodbyrow
function testcumprodbyrow()
  M = rand(10,100)
  Mcu = M |> CuMatrix{Float64}

  print("Testing cumprod....")
  @assert cumprodbyrow(M, limitA=false) ≈ cumprod(M, dims=2)
  @assert cumprodbyrow(Mcu, limitA=false) ≈ cumprod(Mcu, dims=2)

  M = [1.0 10.0 -10.0 10.0; 1.0 10.0 -0.01 10.0]

  largelimit = 1.0
  growthtype=:identity
  Mver = hcat(Iterators.accumulate((x,y)->
    ((xᵢ,yᵢ)-> abs(xᵢ * yᵢ) < largelimit ? xᵢ * yᵢ : largelimit * sign(xᵢ * yᵢ)).(x,y),
    eachcol(M), init=[1.0,1.0])...)
  (cumprodbyrow(M, limitA=true; largelimit, growthtype) ≈ Mver) || error("cumprodbyrow ≈/ Mver
    cumprodbyrow = $(cumprodbyrow(M, limitA=true; largelimit, growthtype))\nM=$M\nMver=$Mver")

  smalllimit = -1.0
  growthtype=:log
  Mver = hcat(Iterators.accumulate((x,y)->
    ((xᵢ,yᵢ)-> max(min(xᵢ*yᵢ,largelimit), smalllimit)).(x,y),
    eachcol(M), init=[1.0,1.0])...)
  (cumprodbyrow(M, limitA=true; largelimit, smalllimit, growthtype) ≈ Mver) || error("cumprodbyrow ≈/ Mver
    cumprodbyrow = $(cumprodbyrow(M, limitA=true; largelimit, smalllimit, growthtype))\nM=$M\nMver=$Mver")

  @info "Passed tests in testcumprodbyrow"
end

#this is the sign command without a 0
signsans0(x::T) where T = (x≥zero(T))*T(2.0)-T(1.0)

#this is a function for converting a CuArray to an Array in a Zygote-safe manner
zygotedevice2host(M::TM, ::Type{T}=eltype(TM)
  ) where {TM, T<:Real} = copyto!(Zygote.Buffer(Array{T}(undef, size(M))), M) |> copy

zygotecuones(::Type{T}, N) where T<:Real = ones(N) |> CuVector{T}

#a recursive parsing routine for converting Named Tuples to dicts
function Base.Dict(nt::NamedTuple)
  parsevalue(x::Any) = x
  parsevalue(x::Union{NamedTuple,Dict}) = Dict(x)
  return Dict(k=>parsevalue(v) for (k,v) ∈ zip(keys(nt), values(nt)))
end


#a simple function to get the lead of a field within a group
function simpleleadwithin!(df::DataFrame,
  target::Symbol=error("target is required"),
  group::Symbol=error("group is required");
  Ntarget::Symbol = Symbol(:N, target),
  Ngroup::Symbol = Symbol(:N, group))

  issorted(df, group) || error("df must be sorted by $group, date")
  #shift up target and group
  df[!, Ntarget] = [df[!, target][2:end]; missing]
  df[!, Ngroup] = [df[!, group][2:end]; missing]

  df[df[!, Ngroup] .!== df[!, group] , Ntarget] .= missing
  select!(df, Not([Ngroup]))

  #check type stability
  @assert eltype(df[!, Ntarget]) <: Union{eltype(df[!, target]), Missing}

  #NOTE: can delete the below at some point
  #only for testing purposes
  lagwithin2!(df, [Ntarget], group)
  LNtarget = Symbol(:L, Ntarget)
  @assert all((df[!, target] .=== df[!, LNtarget]) .| (df[!, LNtarget] .=== missing))
  select!(df, Not([LNtarget]))

  return nothing
end

#simple function to add dates and date periods with missing values
incrementdt(d1,d2) = d1 + d2
incrementdt(::Missing, ::Any) = missing
incrementdt(::Any, ::Missing) = missing

#formats a dict for logging purposes
formatdict(d;
  title="",
  header="\nkeys\tvalues\n",
  ) = title * header * join(string.(keys(d)) .* ":\t" .* string.(values(d)) .* "\n")

function StatsBase.partialcor(X::AbstractArray{T,2}, Z::AbstractArray{T,2};
    zincludesintercept = false) where T
  N,K = size(X)
  @assert size(X,1) == size(Z,1)

  local Xe::Matrix
  if zincludesintercept
    Xₑ = X - Z * (cholesky!(Symmetric(Z'*Z))\(Z'*X))
  else
    Z1 = hcat(ones(T,N), Z)
    Xₑ = X - Z1 * (cholesky!(Symmetric(Z1'*Z1))\(Z1'*X))
  end
  return cor(Xₑ)
end

function testpartialcor(;N=252, K=500, R = 4)
  X = rand(N, K)
  Z = rand(N, R)
  #slow, manual version
  function partialcorver(X::AbstractArray{T,2}, Z::AbstractArray{T,2}) where T
    N,K = size(X)
    @assert size(X,1) == size(Z,1)

    ρ = ones(T, K,K)
    for c ∈ 1:K
      for r ∈ 1:(c-1)
        ρ[r,c] = partialcor(X[:,c], X[:,r],Z)
      end
    end

    ρ .= Symmetric(ρ)
    return ρ
  end

  #create an orthogonalized version for testing purposes
  Z1 = hcat(ones(N), Z)
  Xr = X .- Z1*(cholesky!(Z1'*Z1)\(Z1'*X))

  @info "testing partial correlations"
  @assert (partialcor(Xr, Z) .≈ cor(Xr)) |> all
  @assert (partialcor(X,Z) .≈ partialcor(X,Z1,zincludesintercept=true)) |> all
  #@btime partialcor($X, $Z1,zincludesintercept=true)
  @btime partialcor($X, $Z)

  if K ≤ 1000
    @assert (partialcorver(X,Z) .≈ cor(Xr)) |> all
    @assert (partialcor(X,Z) .≈ partialcorver(X,Z)) |> all
  end

end
