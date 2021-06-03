
struct ZSpec{Tscale<:Dict}
  originalms::MeasureSpec
  Zms::MeasureSpec
  scale::Tscale
end

#standardizes the data to make it easier to process
#WARNING- some field specific code here
function condition!(panel::AbstractDataFrame, ms::MeasureSpec)



  #scale the fields by the standard deviation
  #IMPORTANT: We need to use the same std dev on the lags in order to maintain integrity


  #assemble the new names and replace characteristics with standardized versions
  standardnames::Dict{NSymbol,NSymbol} = Dict{NSymbol,NSymbol}()
  standardnames[nothing]=nothing

  #create a dict to store the scale values
  scale::Dict{Symbol, Float64} = Dict()

  #start by standardizing the characteristics
  #NOTE: each characteristic weight and the lags, and controls are standardized by a single value
  #to avoid distorting changes in weights (so each characteristic uses a single standardizer)
  #QUESTION- could I use a seperate stnadardization value for the controls?
  Zξs::Vector{Characteristic} = Vector{Characteristic}()
  for ξ ∈ ms.ξs

    #IMPORTANT: use same scaling for lagged variables as unlagged equivelents
    standardnames[ξ.Fξ] = Symbol(:Z, ξ.Fξ)
    standardnames[ξ.X[:Fw]] = Symbol(:Z, ξ.X[:Fw])
    σ::Float64 = std(skipmissing(panel[!, ξ.X[:Fw]])) #use the clean panel to compute sds

    #now apply the same scaling to the lagged equivelents
    scale[standardnames[ξ.X[:Fw]]] = σ #do it this way to make sure we got the correct standard name
    for s ∈ [values(ξ.X)...; values(ms.FWbykt[ξ.Fξ])...;]
      standardnames[s] = Symbol(:Z, s)
      panel[!,standardnames[s]] = panel[!, s] ./ σ
      scale[standardnames[s]] = σ #record the scale
    end

    #construct the standardized characteristic
    #println("ZX: $(Dict(s=>standardnames[ξ.X[s]] for s ∈ keys(ξ.X)))")
    push!(Zξs, Characteristic(ξ, standardnames[ξ.Fξ],
      Dict(s=>standardnames[ξ.X[s]] for s ∈ keys(ξ.X))))
  end

  additionalfields::Vector{Symbol} = [ms.Fvol; ms.FW; ms.FWbyt]
  if !(ms.FDvol===nothing)
    push!(additionalfields, ms.FDvol)
  end

  #now standardize all other fields
  for F ∈ additionalfields
    standardnames[F] = Symbol(:Z, F)
    σ::Float64 = std(skipmissing(panel[!, F]))
    if σ + 0.01 ≈ 0.01 #check if σ is approaximately 0
      σ = 1.0#if so, leave the scale alone
    end
    scale[standardnames[F]] = σ
    panel[!, standardnames[F]] = panel[!, F] ./ σ
  end

  #create new names for the standardized controls
  ZFWbykt::Dict{Symbol, Vector{Symbol}} = Dict(standardnames[Fξ]=>(
    FW::Symbol->standardnames[FW]).(ms.FWbykt[Fξ])
    for Fξ ∈ keys(ms.FWbykt))
  Zcontrolnames = (;
    FW = (F->Symbol(:Z, F)).(ms.FW),
    FWbyt = (F->Symbol(:Z, F)).(ms.FWbyt),
    FWbykt = ZFWbykt,
  )

  #construct the standardized measurespec
  Zms = MeasureSpec(Zξs...,
    Fvol=standardnames[ms.Fvol],
    FDvol=standardnames[ms.FDvol],
    controlnames=Zcontrolnames)

  #check on the name integrity
  panelnames = propertynames(panel)
  for f ∈ [weightfields.(ms.ξs)...;(Fξ->ms.FWbykt[Fξ]).(ms.Fξs)...; additionalfields]
    (f ∈ panelnames) || error("$f not found in panel names")
    (Symbol(:Z, f) ∈ panelnames) || error("$(Symbol(:Z, f)) not found in panel names")
  end

  zs = ZSpec(ms, Zms, scale)


  #record the scaling values
  b = IOBuffer()
  #write(b, "Scaling info for $scale:")
  for p ∈ collect(scale)
    write(b, "$(p[1]) σ:\t $(p[2])\n")
  end
  standardizationinfo = String(take!(b))
  open("$(PARAM[:testpath])\\scaling\\$(iterresultsname()).txt", "w+") do f
    write(f, standardizationinfo)
  end


  return zs
end

function destandardize!(results::DataFrame, zs::ZSpec; intercept::Bool = true)
  #coefficients = intercept ? Fcoefficients : Fcoefficients0

  #this assumes the other controls are projected as opposed to calculated directly
  Fcoefficients(ms::MeasureSpec) = [ms.Fξs;#= ms.FW;=#]

  Fcoefs::Vector{Symbol} = Fcoefficients(zs.originalms)
  FZcoefs::Vector{Symbol} = Fcoefficients(zs.Zms)
  Fresults::Vector{Symbol} = propertynames(results)

  #make sure we have the expected coefficients
  if length(Fcoefs) ≠ length(FZcoefs) || length(setdiff(FZcoefs, Fresults)) ≠ 0
    error("Inconsistentcies found between Zspec and measure spec, and results:
      Fcoefs: $Fcoefs\nFZcoefs: $FZcoefs\nFresults: $Fresults")
  end

  #make sure the names are what we expect
  originalnames::Dict = Dict()
  for (Fcoef, FZcoef) ∈ zip(Fcoefs, FZcoefs)
    if FZcoef == Symbol(:Z, Fcoef)
      originalnames[FZcoef] = Fcoef
    else
      error("Unexpected field names found in Zms: \n Fcoefs: $Fcoefs \nFZcoefs: $FZcoefs")
    end
  end
  Symbol(:Z, zs.originalms.Fvol) == zs.Zms.Fvol || "Fvol and ZFvol have inconsistent names"
  originalnames[zs.Zms.Fvol] = zs.originalms.Fvol

  #now do the scaling
  @info "Beginning linear descaling of result fields $FZcoefs"
  #=NOTE: since V_nt/σᵥ =  sum_k 1/σ_kw * |ZA_tk * w_ntk - LZA_t * Lw_t*r_nt| + ε
    we have A_tk = σᵥ/σ_kw * ZA_tk=#

  σᵥ::Float64 = zs.scale[zs.Zms.Fvol]
  for (Fcoef, FZcoef) ∈ zip(Fcoefs, FZcoefs)
    resultcol::Vector = results[!, FZcoef]

    FZw::Symbol = zs.Zms[FZcoef].X[:Fw]
    σ_kw::Float64 = zs.scale[FZw]
    results[!, Fcoef] = resultcol .* σᵥ ./ σ_kw
  end

  standardizenames!(results, zs)
  return nothing

end

function standardizenames!(results::DataFrame, zs::ZSpec)
  Fcoefs::Vector{Symbol} = zs.originalms.Fξs
  FZcoefs::Vector{Symbol} = zs.Zms.Fξs
  FGcoefs::Vector{Symbol} = (s->Symbol(:G_,s)).(Fcoefs)

  @assert setdiff([:date; Fcoefs; FZcoefs; FGcoefs], propertynames(results)) |> isempty
  #setdiff(FGcoefs, propertynames(results)) |> isempty #"

#    Invalid setdiff: FGcoefs: $FGcoefs propertynames(results): $(propertynames(results))"
  #@assert setdiff(propertynames(results), [:date; Fcoefs; FZcoefs; FGcoefs]) |> isempty

  rename!(results, Dict(s=>Symbol(:A_,s) for s ∈ Fcoefs))
  rename!(results, Dict(s=>Symbol(:Z_,s) for s ∈ FZcoefs))



  return results
end

#some re-usable code regarding a variety of items needed to run the optimizatinos
function modelparts(panel::AbstractDataFrame,
    ::Type{T}=throw("3rd arg Type{T} is required!"),
    ::Type{TV}=throw("4th arg Type{TV} is required!"),
    ::Type{TM}=throw("5th arg Type{TM} is required!");
    zs::Union{ZSpec, Nothing} = nothing,
    ms::MeasureSpec = zs ≡ nothing ? throw("Need either zs::ZSpec or ms::MeasureSpec") : zs.Zms,
    weightdata::NamedTuple = error("weight data are required"),
    Fvol::Symbol = throw("Fvol is required"),
    FW::Vector{Symbol} = Vector{Symbol}(),
    FWgroup::Vector{Symbol} = Vector{Symbol}(),
    placebovol::Bool = PARAM[:placebovol],) where {
      T<:Real, TV<:AbstractVector{T}, TM<:AbstractMatrix{T}}




  #get some convenience quantities
  dims = getdims(panel, ms)
  tidx::Dict = gettidx(panel, ms)

  #get the time index for each date
  essentialfields::Vector{Symbol} = [[collect(values(weightdata))...;]...;
    FW...; FWgroup...; Fvol]
  cleanpanel::SubDataFrame = view(panel, completecases(panel, essentialfields), :)

  dat::Dict = Dict{Symbol, Any}()
  for k ∈ keys(weightdata)
    dat[k] = TM(reduce(hcat,(f->(T).(cleanpanel[!, f])).(weightdata[k])))
  end

  #test routine
  if placebovol
    dat[:Fvol] = (T).(cleanpanel[!, Fvol]) |> Random.shuffle |> TV
  else
    dat[:Fvol] = (T).(cleanpanel[!, Fvol]) |> TV
  end

  #acquire controls data. Creates an empty matrix if no controls
  emptymat = Matrix{T}(undef, size(cleanpanel,1), 0)
  dat[:FW] = reduce(hcat, (f->(T).(cleanpanel[!, f])).(FW), init=emptymat) |>TM
  dat[:FWgroup] = reduce(hcat, (f->(T).(cleanpanel[!, f])).(FWgroup), init=emptymat) |>TM
  dat[:permno] = cleanpanel[!, :permno] |> deepcopy
  dat[:date] = cleanpanel.date |> deepcopy

  ts::Vector{Int} = (d::Date->tidx[d]).(cleanpanel[!,:date])
  @assert issorted(ts)
  #V::TV = (T).(cleanpanel[!,ms.Fvol]) |> fluxgpu

  #now get the scaling values
  #Ascale = Vector{T}(undef, dims.K)
  #@assert length(Ascale) == size(dat[:Fw],2) #integrity check
  #FZcoefs::Vector{Symbol} = Fcoefficients0(zs.Zms)
  #@assert length(FZcoefs)  == length(Ascale)

  #=
  σᵥ::Float64 = zs.scale[zs.Zms.Fvol]
  for (i,FZcoef) ∈ enumerate(FZcoefs) #loop through each coefficient and identify the scaling value
    FZw::Symbol = zs.Zms[FZcoef].X[:Fw] #pull the standard weighting name (this is keyed to the scale)
    σ_kw::Float64 = zs.scale[FZw]
    Ascale[i] = σᵥ / σ_kw
  end

  @info "Scaling Values: (σᵥ, $σᵥ), $((zip((FZcoef->zs.Zms[FZcoef].X[:Fw]).(FZcoefs), Ascale)...))"
  =#
  return (dims=dims, tidx=tidx, cleanpanel=cleanpanel,
    dat=dat, ts=ts#=, Ascale=Ascale=#)
end
