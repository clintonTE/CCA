module CTModCopy

using  DataFrames, Distributions, StatsBase, GLM#, Base.BLAS, Base.LinAlg, RData

const DEBUG_CTMOD = false
const DEFAULT_STAR_LEGEND = "*p < 0.1, **p < 0.05, ***p < 0.01 s.t. p=Pr(>|T|)"

export CTLM, CT2SLS, project!, getResid!, getHomoskedΣ!, getModWhiteΣ!,
  pullModelMatrix, getCoeff!, getModWhiteΣSlow, getHomoskedΣSlow,
  dropNullsFromDF, CTExpr, CTSym, texTable, writeTables2File, CTQR,
  get1stStage, getTerm, getR

#this alias is useful when parsing inputs to formula objects
const CTExpr = Union{Symbol,Expr,Void}
const CTSym = Union{Symbol,Void}

abstract type CTModel end

##################Utility Functions ###############################

#helper function  to get the model matrix from a dataframe given a formula
#IN: A dataframe and formula object
#OUT: the model matrix
function getModelMatrix(df::T, f::Formula)::Matrix{Float64} where
  T <: AbstractDataFrame
  return ModelMatrix(ModelFrame(f, df)).m
end

#Same as above but allows for an expression
function getModelMatrix(df::T, exp::V)::Matrix{Float64} where
  {T <: AbstractDataFrame, V <: CTExpr}

  #special case which crashes Formula
  if exp == Symbol("")
    return ones(Float64,size(df,1),1)
  end

  return getModelMatrix(df, get1SidedFormula(exp))
end

#helper function  to get the list of symbols in an expression
#IN: A dataframe and a string object, typically representing a formula
#OUT: a vector of symbols in the string object
function getSymbolsFromExpr(exp::V)::Vector{Symbol} where {V <: CTExpr}
  symStrings::Vector{String} = split(string(exp),
    ['|','=','~',' ','+','*','&',')','(','-'], keep=false)

  filter!(s::String->typeof(parse(s))≠Int,symStrings)

  #convert the strings to symbols
  return Symbol.(symStrings)
end


#helper function to create a one-sided formula given an expression
#IN: an expression and dataframe
#OUT: A one-sided formula object
function get1SidedFormula(RHS::T)::Formula where
  {T <: CTExpr}

  return Formula(nothing, RHS)
end

#drops nulls  from dataframe given a list of symbols
#IN: a data frame and symbols
#OUT: dataframe less the null entires as dictated by the symbol list
dropNullsFromDF(df::DataFrame, syms::Vector{Symbol})::DataFrame =
  df[completecases(df[:,syms]),:]

#drops nulls  from dataframe given a list of symbols
#IN: a sub-dataframe and symbols, #OUT: subdataframe less the null entries
dropNullsFromDF(dfSub::SubDataFrame, syms::Vector{Symbol})::SubDataFrame =
    view(dfSub,completecases(dfSub[:,syms]))

#helper function for the above, takes in a single symbol
#IN: a subdataframe and symbol #OUT: subdataframe less the null entries
function dropNullsFromDF(df::T, s::Symbol)::T  where {T<:AbstractDataFrame}

  return dropNullsFromDF(df, collect([s]))
end

#helper function for the above, de-nulls entire dataframe
#IN: a subdataframe #OUT: subdataframe less the null entries
function dropNullsFromDF(df::T)::T where {T<:AbstractDataFrame}
  return dropNullsFromDF(df, names(df))::T
end

#####################CTQR Type###########################
#=This is a storage object specific to the QR optimizations
meant to hold the QR factorization and a handle for X
Values: X matrix that is decomposed, Q matrix, and R matrix
The X matrix has dimensions of NxK. R is the inverse, stored
for computational efficiency.=#

struct CTQR
  X::Matrix{Float64}
  Q::Matrix{Float64}
  R::Matrix{Float64}

  N::Int
  K::Int
  RInv::Matrix{Float64}
end

#helper method which takes the inverse of the R matrix (which is used in a lot)
#IN: X matrix and its Q R decomposition in matrix form
#OUT: CTQR Object
CTQR(X::Matrix{Float64},Q::Matrix{Float64},R::Matrix{Float64};keepX=true)::CTQR =
  CTQR(keepX?X:Matrix{Float64}(),Q,R,size(X,1),size(X,2),(R)\eye(size(X,2)))

#main constructor for the CTQR object
#IN: X Matrix
#OUT: CTQR object
function CTQR(X::Matrix{Float64}; keepX = true)::CTQR
  Q::Matrix{Float64},R::Matrix{Float64} = qr(X)
  return CTQR(X,Q,R)
end

#default constructor to clear memory
CTQR() = CTQR(Matrix{Float64}(),Matrix{Float64}(),Matrix{Float64}(),
    0,0,Matrix{Float64}())::CTQR


#####################CTLM Model Type#########################
#This is designed to hold the minimal info for a linear model
#Plus the associated CTQR object
 mutable struct CTLM <: CTModel
  xqr::CTQR #handle to the xqr object
  X::Matrix{Float64} #the data and x matrix (handle)
  Y::Vector{Float64} #Y data

  N::Int #Number of data points (rows)
  K::Int #Number of data columns
  β::Vector{Float64} #beta coefficient
  ε::Vector{Float64} #residuals

  XNames::Vector{Symbol} #names of the X variables
  YName::Symbol #names of the Y variables
  #rowLabels::Vector{Date}


end

#this function calculates the Pearson correlation coefficient of the predicted values
#relative to the realized values
getR(lin::CTLM)::Float64 = cor(lin.Y, lin.Y .- lin.ε)


#convenience function, calls CTQR constructor
#Standard constructor to be called from outside
#IN: X matrix of independent vars, Y vector for dependent var
#OUT: CTLM Object
function CTLM(xqr::Matrix{Float64}, Y::Vector{Float64};
    XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
    YName::Symbol = :Y)::CTLM

    #showall(X[1:100,1:2])
    xqr::CTQR = CTQR(X, keepX=false)
    β::Vector{Float64} = Vector{Float64}(xqr.K)
    getCoeff!(xqr,Y,β)

    ε::Vector{Float64} = Vector{Float64}(xqr.N)
    getResid!(X,Y,β,ε)

    return CTLM(xqr, X, Y, xqr.N, xqr.K, β, ε, XNames, YName)
end

#convenience function, calls CTQR constructor
#Standard constructor to be called from outside
#IN: X matrix of independent vars, Y vector for dependent var
#OUT: CTLM Object
function CTLM(X::Matrix{Float64}, Y::Vector{Float64};
    XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
    YName::Symbol = :Y)::CTLM

    #showall(X[1:100,1:2])
    xqr::CTQR = CTQR(X, keepX=false)
    β::Vector{Float64} = Vector{Float64}(xqr.K)
    getCoeff!(xqr,Y,β)

    ε::Vector{Float64} = Vector{Float64}(xqr.N)
    getResid!(X,Y,β,ε)

    return CTLM(xqr, X, Y, xqr.N, xqr.K, β, ε, XNames, YName)
end

#=Constructor and helper method which gets required info from DataFrame
IN: The source dataframe, a formula expression for the RHS,
    the dependent variable as a symbol, optionally names for X columns,
    Y and a switch to elliminate null values
OUT: A CTLM object
NOTE: RHS expression must be ordered with factors last, otherwise
    the factor names will be incorrect =#
function CTLM(df::DataFrame,  XExpr::T, YSym::Symbol;
    XNames::Vector{Symbol}=[Symbol(:Intercept)], YName::Symbol=:Y,
    eliminateNulls=true, fixedEffectsSym::V = nothing)::CTLM  where {T <: CTExpr, V<:Union{Symbol,Void}}

    #get the list of symbols from the expression for X
    XSym::Vector{Symbol} = getSymbolsFromExpr(XExpr)

    #make a view so if we drop nulls it doesn't affect the original
    if fixedEffectsSym == nothing || fixedEffectsSym ∈ [XSym; YSym]
        dfSub::SubDataFrame = view(df[:,[XSym; YSym]],1:size(df,1))
    else
        dfSub = view(df[:,[XSym; YSym; fixedEffectsSym]],1:size(df,1))
    end

    if eliminateNulls     # if we are going to drop nulls
        dfSub = dropNullsFromDF(dfSub)
    end

    #get the factor expanded model matrix (Adds the dummies)
    XModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, XExpr)

    #this is a check to make sure the factor expansion follows other columns
    if length(XSym) > 1
        for i ∈ 2:length(XSym) #all types
            if typeof(dfSub[:,XSym[i-1]])<:PooledDataVector &&
                !(typeof(dfSub[:,XSym[i]])<:PooledDataVector)
                warn("Type of $(typeof(dfSub[:,i])) is preceeded by type of
                    factor. Column names are likely MISALIGNED. Put factors
                    after numerical variables in RHS expressions.")
            end
        end
    end

    #Now assign the names as appropriate
    XNamesFull::Vector{Symbol} = Vector{Symbol}(size(XModelMatrix,2))
    for i ∈ 1:length(XNamesFull)
        XNamesFull[i] = i≤length(XNames)?XNames[i]:Symbol(:X,i)
    end

    #get the appropriate matrices and construct the CTLM object

    if fixedEffectsSym != nothing
            #println(typeof(Vector(dfSub[:,fixedEffectsSym])))
        return CTLM(XModelMatrix, Vector{Float64}(dfSub[:,YSym]), Vector(dfSub[:,fixedEffectsSym]),
            XNames=XNamesFull, YName=YName)
    else
        return CTLM(XModelMatrix,
            Vector{Float64}(dfSub[:,YSym]),
            XNames=XNamesFull, YName=YName)
    end
end

#One-way clustered SEs
function CTLM(X::Matrix{Float64}, Y::Vector{Float64}, fixedEffects::Vector{V};
    XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
    YName::Symbol = :Y)::CTLM where V<:Union{Float64,Int,Date,Symbol}

    #we get the unique values and make a two way table
    #tVars::Vector{T} = unique(tFI)
    iVars::Vector{V} = unique(fixedEffects)
    iN::Int = length(iVars)
    iTable::Dict = Dict(iVars[i] => i for i::Int ∈ 1:iN)

    K::Int = size(X,2)
    N::Int = size(X,1)

    #now create the codified versions
    #Xi::Matrix{Float64} = Matrix{Float64}(length(iVars),size(X,2))
    fixedEffectsCode::Vector{Int} = ((x::V)->iTable[x]).(fixedEffects)

    meanXY::Matrix{Float64} = zeros(iN,K+1)
    NXY::Vector{Int} = zeros(iN)

    #println("$N,$K,$iN,$(size(NXY))")
    @fastmath @inbounds for r ∈ 1:N
        meanXY[fixedEffectsCode[r],K+1] += Y[r]
        NXY[fixedEffectsCode[r]] += 1.0
    end

    @fastmath @inbounds for c ∈2:K, r ∈1:N
            meanXY[fixedEffectsCode[r],c] += X[r,c]
    end

    meanXY ./= NXY

    @fastmath @inbounds for c ∈2:K, r ∈1:N
            X[r,c] -= meanXY[fixedEffectsCode[r],c]
    end

    @fastmath @inbounds for r ∈ 1:N
         Y[r] -= meanXY[fixedEffectsCode[r],K+1]
    end


    #=@fastmath Threads.@threads for i ∈ 1:iN
        cluster::Vector{Int} = find((x::Int)->x==i,fixedEffectsCode)

        meanY::Float64 = mean(Y[cluster])
        #print(i<10?"$meanY\n":"")
        Y[cluster] .-= meanY

        #We get a dramatic benefit from the parallelization
        for j::Int ∈ 1:K
            if XNames[j] != :intercept
                meanX::Float64 = mean(X[cluster,j])
                X[cluster, j] .-= meanX
            end
        end
    end=#

    return CTLM(X, Y, XNames=XNames, YName=YName)
end

#blank default constructor for memory management
CTLM() = CTLM(CTQR(), Matrix{Float64}(), Vector{Float64}(), 0, 0,
        Vector{Float64}(), Vector{Float64}(), Vector{Symbol}(), Symbol())::CTLM

#fucntion for getting regression coefficients
getTerm(lm::CTLM, s::Symbol) = (lm.β[findfirst(lm.XNames, s)])::Float64

######2SLS constructor
#Holds the modelling componenets for a 2SLS model (except errors)
#Variable descriptions are in the struct
# Runs standard 2SLS using QR optimization
mutable struct CT2SLS <: CTModel
    zaqr::CTQR #QR decomposition of Za (Z|W)
    xaqr::CTQR #QR decomposition of fitted 1st stage values (X̃|W)

    X::Matrix{Float64} #Endogeneous covariates
    W::Matrix{Float64} #Exogeneous covariates
    Y::Vector{Float64} #Dependent var
    Z::Matrix{Float64} #IV data

    N::Int #number of data points
    K::Int # of endogeneous covariates
    L::Int # of instrumental variables
    KW::Int # of exogeneous covariates

    Π1::Matrix{Float64} #first stage results
    Ξ1::Matrix{Float64} #first stage residuals

    δ2::Vector{Float64} #2SLS coefficients
    ξ2::Vector{Float64} #2SLS residuals

    XNames::Vector{Symbol} #names of the X variables
    WNames::Vector{Symbol} #names of the W variables
    YName::Symbol #names of the Y variables
    ZNames::Vector{Symbol} #names of the Z variables

    #CT2SLS()
end


#= Convenience constructor for CT2SLS that allocates a 0 element array for the
case of no exogeneous covariates. NOTE: no exog covariates implies no intercept,
so there probably should be  some exogeneous covariates
#IN: X Matrix, Y Vector, and the Z matrix of Instrumental variables
#OUT: CT2SLS object=#
CT2SLS(X::Matrix{Float64},Y::Vector{Float64},Z::Matrix{Float64};
    XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
    YName::Symbol = :Y,
    ZNames::Vector{Symbol} = Symbol.(:Z, 1:size(Z,2)))::CT2SLS =
        CT2SLS(X,Matrix{Float64}(Size(X,1),0),Y,Z)

#=Main constructor for CT2SLS object
IN: Endogeneous matrix X, Exogeneous matrix W, Dependent vector Y,
instrumental zariable matrix Z
OUT: A 2SLS object
=#
function CT2SLS(X::Matrix{Float64}, W::Matrix{Float64},
  Y::Vector{Float64},Z::Matrix{Float64};
  XNames::Vector{Symbol} = Symbol.(:X, 1:size(X,2)),
  WNames::Vector{Symbol} = Symbol.(:W, 1:size(W,2)),
  YName::Symbol = :Y,
  ZNames::Vector{Symbol} = Symbol.(:Z, 1:size(Z,2)))::CT2SLS

  #memory pre-allocaiton
  N::Int = size(X,1) #number of data points
  K::Int = size(X,2) # of endogeneous variables
  L::Int = size(Z,2) # of instrumental variables
  KW::Int = size(W,2)# of exogeneous covariates

  Π1::Matrix{Float64} = Matrix{Float64}(L+KW,K)
  X̂::Matrix{Float64} = Matrix{Float64}(N,K) #holds fitted first stage values
  Ξ1::Matrix{Float64} = similar(X̂) #holds first stage residuals
  δ2::Vector{Float64} = Vector{Float64}(K+KW)
  ξ2::Vector{Float64} = similar(Y)
  #Xa::Matrix{Float64} = Matrix{Float64}(N,K+KW) #fitted 1st stage plus exog covariates

  #start with setting up and running the first stage
  ZW::Matrix{Float64} = [Z W]
  zaqr::CTQR = CTQR(ZW, keepX=false) #augment the matrices to ZW and get the QR
  getCoeff!(zaqr, X, Π1) #get the first stage coefficients
  BLAS.gemm!('N','N',1.0,ZW,Π1,0.0,X̂) #get the fitted first stage values
  getResid!(ZW, X, Π1, Ξ1) #get the residuals

  #now run the second stage
  xaqr::CTQR = CTQR([X̂ W], keepX = false) #augment X̂W and get the QR
  getCoeff!(xaqr, Y, δ2) #get the 2SLS coefficients
  getResid!([X W], Y, δ2, ξ2) #get the 2SLS residuals

  #Finally, return the constructed model
  return CT2SLS(zaqr, xaqr, X, W, Y, Z, N, K, L, KW, Π1, Ξ1, δ2, ξ2,
    XNames, WNames, YName, ZNames)
end




#Constructor and helper method which gets required info from DataFrame
#IN: The source dataframe, a formula expression/symbol for X, W, Y, and Z
# the dependent variable as a symbol
#OUT: A CTLM object
function CT2SLS(df::DataFrame, XExpr::T, WExpr::U, YSym::Symbol, ZExpr::V;
    XNames::Vector{Symbol} = Vector{Symbol}(),
    WNames::Vector{Symbol} = Vector{Symbol}([:Intercept]),
    YName::Symbol = :Y,
    ZNames::Vector{Symbol} = Vector{Symbol}(),
    eliminateNulls=true)::CT2SLS where {T <: CTExpr, U <: CTExpr, V <: CTExpr}

    #get the list of symbols for X, W and Z
    XSym::Vector{Symbol} = getSymbolsFromExpr(XExpr)
    WSym::Vector{Symbol} = getSymbolsFromExpr(WExpr)
    ZSym::Vector{Symbol} = getSymbolsFromExpr(ZExpr)

    #make a view so if we drop nulls it doesn't affect the original
    dfSub::SubDataFrame = view(df[:, [XSym; WSym; YSym; ZSym]],1:size(df,1))

    # if we are going to drop nulls
    if eliminateNulls
        dfSub = dropNullsFromDF(dfSub)
    end

    #get the model matrices, which will expand all factors
    XModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, XExpr)
    WModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, WExpr)
    ZModelMatrix::Matrix{Float64} = getModelMatrix(dfSub, ZExpr)

    #this is a check to make sure the factor expansion follows other columns
    #We must repeat the check on each variable matrix
    for symVec::Vector{Symbol} ∈ [XSym, WSym, ZSym]
        if length(symVec) > 1
            for i ∈ 2:length(symVec) #all types
                if typeof(dfSub[:,symVec[i-1]])<:PooledDataVector &&
                    !(typeof(dfSub[:,symVec[i]])<:PooledDataVector)
                    warn("Type of $(typeof(dfSub[:,i])) is preceeded by type of
                        factor. Column names are likely MISALIGNED. Put factors
                        after numerical variables in RHS expressions.")
                end
            end
        end
    end

    #now create the namesfor the expanded matrix
    XNamesFull::Vector{Symbol} = Vector{Symbol}(size(XModelMatrix,2))
    for i ∈ 1:length(XNamesFull)
        XNamesFull[i] = i≤length(XNames)?XNames[i]:Symbol(:X,i)
    end

    WNamesFull::Vector{Symbol} = Vector{Symbol}(size(WModelMatrix,2))
    for i ∈ 1:length(WNamesFull)
        WNamesFull[i] = i≤length(WNames)?WNames[i]:Symbol(:W,i)
    end

    ZNamesFull::Vector{Symbol} = Vector{Symbol}(size(ZModelMatrix,2))
    for i ∈ 1:length(ZNamesFull)
        ZNamesFull[i] = i≤length(ZNames)?ZNames[i]:Symbol(:Z,i)
    end

    #get the appropriate matrices and construct the CT2SLS object
    return CT2SLS(XModelMatrix, WModelMatrix, Vector{Float64}(dfSub[:,YSym]),
        ZModelMatrix, XNames=XNamesFull, WNames=WNamesFull,
        YName=YName, ZNames=ZNamesFull)

end

#default constructor to help with memory management
CT2SLS() = CT2SLS(CTQR(), CTQR(), Matrix{Float64}(), Matrix{Float64}(),
    Vector{Float64}(), Matrix{Float64}(), 0, 0, 0, 0, Matrix{Float64}(),
    Matrix{Float64}(), Vector{Float64}(), Vector{Float64}(),
    Vector{Symbol}(), Vector{Symbol}(), Symbol(), Vector{Symbol}())::CT2SLS

#=function get1stStage
#The purpose is to get the first stage of the IV as a linear model
#Unfortunately, we need to copy some of the data, so this can be expensive,
#although Z and W should only be copied once
#IN: an IV object
#OUT: a vector of CTLM objects, one for each focal (X) variable in the
#initial regression=#
function get1stStage(iv::CT2SLS)::Vector{CTLM}
    ZW::Matrix{Float64} =  [iv.Z iv.W]
    iv1st::Vector = Vector{CTLM}()

    for i ∈ 1:length(iv.XNames)
        push!(iv1st, CTLM(iv.zaqr, [iv.Z iv.W], iv.X[:,i],
            iv.zaqr.N, iv.zaqr.K, iv.Π1[:,i], iv.Ξ1[:,i],
            [iv.ZNames; iv.WNames], iv.XNames[i]))
    end

    return iv1st
end


#########################Project!##############################
#=gets the projection matrix

IN: THe CTQR decomposition and the memory for  the projection matrix
OUT: Writes and returns the projection matrix
NOTE: This method will allocate the projection matrix if no second
argument is supplied. =#
function project!(xqr::CTQR, P::Matrix{Float64}=Matrix{Float64}(xqr.N, xqr.N)
  )::Matrix{Float64}

  #use the BLAS library for fast matrix algebra
  return BLAS.gemm!('N','T',1.0,xqr.Q,xqr.Q,0.0,P)

end

#=This version gets only the diagonal of the projection matrix
IN: THe CTQR decomposition and the memory for  the projection diagonal
OUT: Writes and returns the projection matrix=#
function project!(xqr::CTQR, P::Vector{Float64})::Vector{Float64}
  #get the diagonal quickly
  #equivlenet to diag(Q*Q')
  P .= 0.0
  @fastmath for j::Int ∈ 1:xqr.K, i::Int ∈ 1:xqr.N
    P[i] += xqr.Q[i,j] * xqr.Q[i,j]
  end

  return P

end

#=This version runs either of the above methods after generating
the xqr object given an independent variable matrix
IN: Independent variable X matrix and the memory for  the projection matrix
or vector for the diagonal
OUT: Writes and returns the projection matrix=#
function project!(X::Matrix{Float64},
  P::Array{Float64}=Matrix{Float64}(size(X,1),size(X,1)))::Array{Float64}

  return project!(CTQR(X),P)
end

#=simple test function to verify the projection optimizations
output should be identical to the matrix function above
IN: Independent variable X matrix, memory for the projection matrix
OUT: Writes and returns teh projection matrix=#
function projectSlow!(X::Matrix{Float64},P::Matrix{Float64})
  P .= (X * ((X' * X)\eye(size(X,2))) * X')
  return P
end

#########################getResid!#########################
#=gets the residuals in the 2SLS framework

Note this is a wrapper function which can accept multiple independent
and multiple dependent variables. NOTE the change in the naming convention
relative to other verisons of this method (due to the 2SLS context)
IN: QR decomposition of the RHS variables, LHS dependent variable matrix,
a coefficient matrix, optionally the memory for the residual matrix
OUT: Writes and returns the residual matrix=#
function getResid!(ZW::Matrix{Float64}, X::Matrix{Float64}, Π1::Matrix{Float64},
  Ξ1::Matrix{Float64}=similar(X))::Matrix{Float64}

  #get the residuals by reference
  for i ∈ 1:size(X,2) #get the vector of reisduals multiple times
    getResid!(ZW, view(X,:,i), view(Π1,:,i), view(Ξ1,:,i))
  end
  return Ξ1
end

#=Main method to get residuals
IN: X matrix the RHS variables, LHS dependent variable,
a coefficient vector, optionally the memory for the residual vector
OUT: Writes and returns the residual vector=#
function getResid!(X::Matrix{Float64}, Y::T,
  β::T, ε::T=similar(Y))::T where T<:StridedVector{Float64}

  return BLAS.axpy!(1.0, Y, BLAS.gemv!('N',-1.0, X, β,0.0,ε)) #ε=Y-Xβ
end

#########################getCoef!###########################
# Gets the coefficients in a regression

#=This is a wrapper function which gets a matrix of coefficients for use in
# the first stage of 2SLS. NOTE the change in the naming convention
relative to other verisons of this method (due to the 2SLS context)
IN: QR decomposition of the RHS variables, LHS dependent variable matrix,
optionally the memory for the coefficient matrix
OUT: Writes and returns the coefficient matrix=#
function getCoeff!(zaqr::CTQR, X::Matrix{Float64},
  Π1::Matrix{Float64} = Matrix{Float64}(zaqr.K, size(X,2)))::Matrix{Float64}

  for i ∈ 1:size(X,2) # for each dependent variable vector
    getCoeff!(zaqr, view(X,:,i), view(Π1,:,i)) #run the regression as needed
  end

  return Π1
end

#=Main method to get the regression coefficients
IN: QR decomposition of X matrix the RHS variables, LHS dependent variable,
optionally the memory for the coefficient vector
OUT: Writes and returns the coefficient vector=#
function getCoeff!(xqr::CTQR, Y::T,
  β::T = Vector{Float64}(xqr.K))::T where  T<:StridedVector{Float64}

  return  BLAS.gemv!('N',1.0,BLAS.gemm('N', 'T', xqr.RInv, xqr.Q),Y,0.0,β)
end

#=Convenience method to get the regression coefficients
IN: X matrix the RHS variables, LHS dependent variable,
optionally the memory for the coefficient vector
OUT: Writes and returns the coefficeint vector=#
getCoeff!(X::Vector{Float64}, Y::Vector{Float64}) =
  getCoeff!(CTQR(X), Y, Vector{Float64}(size(X,2)))::Vector{Float64}


###################getHomosked!#######################
#Gets the covariance matrix under Homoskedastic assumptions σ^2[X'X]^-1

#main method to get the covariance matrix
#IN: The QR decomposition from the regression, a vector of residuals,
#optionally memory for the covariance matrix
#OUT: Writes and returns the covariance matrix
function getHomoskedΣ!(xqr::CTQR, ε::Vector{Float64},
  Σ::Matrix{Float64} = Matrix{Float64}(lin.K, lin.K))::Matrix{Float64}

  return BLAS.gemm!('N','T',(ε⋅ε)/(xqr.N-xqr.K),xqr.RInv,xqr.RInv,0.0,Σ) #[X'X]^-1*σ2
end

#=helper method for the above function which extracts the QR decomposition
from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: WRites and returns the covariance matrix=#
function getHomoskedΣ!(lin::CTLM,
  Σ::Matrix{Float64}=Matrix{Float64}(lin.K, lin.K))::Matrix{Float64}

  return getHomoskedΣ!(lin.xqr, lin.ε, Σ)
end

#=helper method for the above function which extracts the QR decomposition
from a CT2SLS object
IN: A linear model, optionally memory for the covariance matrix
OUT: WRites and returns the covariance matrix=#
function getHomoskedΣ!(iv::CT2SLS,
  Σ::Matrix{Float64}=Matrix{Float64}(iv.xaqr.K, iv.xaqr.K))::Matrix{Float64}

  return getHomoskedΣ!(iv.xaqr, iv.ξ2, Σ)
end

#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
getHomoskedΣSlow(X::Matrix{Float64}, ε::Vector{Float64})::Matrix{Float64} =
  (X.'*X)\eye(size(X,2))*(ε⋅ε)/(size(X,1)-size(X,2))

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
getHomoskedΣSlow(lin::CTLM)::Matrix{Float64} =
  getHomoskedΣSlow(lin.X, lin.Y.-lin.X * lin.β)

#Convenience method for above: IN: A CT2SLS object #OUT: A covariance matrix
getHomoskedΣSlow(iv::CT2SLS)::Matrix{Float64} =
  getHomoskedΣSlow([[iv.Z iv.W]*iv.Π1 iv.W], iv.Y.-[iv.X iv.W] * iv.δ2)

###################getMWErrors!##############################
#Main method to get the modified white SEs Does #[X'X]^-1X'ΛX[X'X]^-1

#=Main method to calculate the modified white SEs
IN: A CTQR decomposition object, a residuals vector,
optional memory for a covariance matrix
OUT: Writes and returns the covariance matrix=#
function getModWhiteΣ!(xqr::CTQR, ε::Vector{Float64},
  Σ::Matrix{Float64} = Matrix{Float64}(lin.K, lin.K))::Matrix{Float64}

  Λ::Vector{Float64} = similar(ε)

  project!(xqr,Λ) #get the projection diagonal
  Λ .= (ε./(1.0.-Λ)).^2.0

  QRtInv::Matrix{Float64} = BLAS.gemm('N','T',xqr.Q,xqr.RInv) #Q(R')^-1
  RInvQΛ::Matrix{Float64} = QRtInv.'

  #loop to scale the matrix by Λ (modified part of modified white)
  @fastmath for j ∈ 1:xqr.N, i∈1:xqr.K
    RInvQΛ[i,j] *= Λ[j]
  end

  #final multiplicaiton and assignment
  return BLAS.gemm!('N','N',1.0,RInvQΛ, QRtInv,0.0,Σ) #R^-1Q'ΛQ(R')^-1
end

#=helper method which extracts the required components from a linear model
IN: A linear model, optionally memory for the covariance matrix
OUT: Writes and allocates the covariance matrix=#
function getModWhiteΣ!(lin::CTLM,
  Σ::Matrix{Float64} = Matrix{Float64}(lin.K, lin.K))::Matrix{Float64}

  return getModWhiteΣ!(lin.xqr, lin.ε, Σ)
end

#=helper method which extracts the required components from a 2SLS model
IN: A 2SLS model, optionally memory for the covariance matrix
OUT: Writes and allocates the covariance matrix=#
function getModWhiteΣ!(iv::CT2SLS,
  Σ::Matrix{Float64}=Matrix{Float64}(iv.xaqr.K, iv.xaqr.K))::Matrix{Float64}

  return getModWhiteΣ!(iv.xaqr, iv.ξ2, Σ)
end

#=Simple test function to verify the QR-related optimizations
Output should be identical to the standard methods
IN: Independent variable X matrix
OUT: Returns the covariance matrix=#
function getModWhiteΣSlow(X::Matrix{Float64}, ε::Vector{Float64})
  XXInv::Matrix{Float64} = (X.' * X)\eye(size(X,2))

  P::Matrix{Float64} = X * (XXInv) * X.'
  return XXInv * X.' * diagm(ε  ./ (1.0 .- diag(P))) .^ 2.0 * X * XXInv

end

#Convenience method for above: IN: A linear model #OUT: A covariance matrix
getModWhiteΣSlow(lin::CTLM) = getModWhiteΣSlow(lin.X, lin.Y.-lin.X * lin.β)

#Convenience method for above: IN: A 2SLS object #OUT: A covariance matrix
getModWhiteΣSlow(iv::CT2SLS)::Matrix{Float64} =
  getModWhiteΣSlow([[iv.Z iv.W]*iv.Π1 iv.W], iv.Y.-[iv.X iv.W] * iv.δ2)



######################IO################


#texTable
#=    Creates a well-formed Latex table from string components
  The content matrix consists of an array of matrices, which are overlayed
  via a one line offset. For example, the first matrix might be the values,
  the second matrix the errors. Observational rows can be put in
  descContent, which is a vector of rows. Similarly, colNames is a vector
  of rows. So that column headers can span multiple columns, the widthColNames
  and widthDescContent provides the option to specify dimentions
  IN: caption (title), colNames, rowNames, content, descContent (rows of
      descriptive data), notes, optionally widthColNames and widthDescConent
      (which are the width of each entry in columns)
  OUT: The tex table as a string=#

function texTable(titleCaption::String,
  caption::String,
  colNames::Vector{Vector{String}},
  contentRowNames::Vector{String},
  content::Vector{Matrix{String}},
  descRowNames::Vector{String},
  descContent::Vector{Vector{String}},
  notes::Vector{String};
  arrayStretch::Float64 = 1.5,
  lineSpacer::String = "\\\\",
  summaryMathMode::Bool = true,
  widthColNames::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
  alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
    broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
  widthDescContent::Vector{Vector{Int}} = #contains the number of columns for each entry
    broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
  columnSepPt::Int = 0,
  colHeaderName::Vector{String} = ["" for i ∈ 1:length(colNames)])

  #size parameters for content
  numContentRows::Int = length(contentRowNames)

  numContentCols::Int = sum(widthColNames[1])

  numContentSubRows::Int = length(content)

      #intiate the stream
  b::IOBuffer = IOBuffer()
  write(b, """
      %This table was programatically generated

      \\begin{table} \\caption{$titleCaption} \\label{} \\centering
      $(length(caption)>0?"\\textit{$caption\\\\}":"")
      \\renewcommand{\\arraystretch}{$arrayStretch}
      \\begin{tabular}{l | """)

    for i ∈ 1:(numContentCols) #set the dimensions in tex
      write(b,"r")
    end
    #write(b, "{\\textwidth}{Xccccccc}")
    write(b,"}\n \\toprule")#filler tex

  for r::Int ∈ 1:length(colNames) #for each row of column headings
      #write(b,"\n\\\\[-1.8ex]") #use this if the below doesn't work
    write(b, colHeaderName[r])
    if length(colNames[r]) ≠ numContentCols
      for c::Int ∈ 1:length(colNames[r]) #for each column heading
          write(b,"\t&\t\\multicolumn{$(widthColNames[r][c])}{$(alignmentColNames[r][c])}{$(colNames[r][c])}")
      end
    else
      for c::Int ∈ 1:length(colNames[r]) #for each column heading
        write(b,"\t&\t$(colNames[r][c])")
      end
    end
      write(b,"\n\\\\")
  end
  write(b, " \\midrule\n ")

  #now write out the table content
    if numContentRows > 0
        write(b," \t \t ")
        for r::Int ∈ 1:numContentRows
          write(b, "$(contentRowNames[r])") #the row label
          for s::Int ∈ 1:numContentSubRows #print the sub-rows for each row
            for c::Int ∈ 1:numContentCols #print the columns for each sub-row
              #write(b,"\t&\t\\multicolumn{1}{r}{\$$((content[s])[r, c])\$}")
              write(b,"\t&\t$((content[s])[r, c])")
            end
            write(b,"\n $(lineSpacer) \t\t") #line-break stylistic formatting
          end
        end
        write(b, "\n \\midrule\n")
    end
    sumMathFlag::String = summaryMathMode?"\$":""

    if length(descContent)>0
        for r::Int ∈ 1:length(descContent) #for each description row
          write(b, "$(descRowNames[r])") #the row label
          if length(widthDescContent[r]) ≠ numContentCols
            for c::Int ∈ 1:length(descContent[r]) #for each descriptive row
              if length(descContent[r][c]) ≥ 1
                write(b,"\t&\t\\multicolumn{$(widthDescContent[r][c])}{r}{$(sumMathFlag)$(descContent[r][c])$(sumMathFlag)}")
              else
                write(b,"\t&\t\\multicolumn{$(widthDescContent[r][c])}{r}{}")
              end
            end
          else
            for c::Int ∈ 1:length(descContent[r]) #for each descriptive row
              if length(descContent[r][c]) ≥ 1
                write(b,"\t&\t$(sumMathFlag)$(descContent[r][c])$(sumMathFlag)")
              else
                write(b,"\t&\t")
              end
            end
          end
            write(b,"\n $(lineSpacer) \t\t") #line-break stylistic formatting
        end
    end

      write(b,""" \\bottomrule""")
  #if we have footnotes
  if length(notes) > 0
      write(b,"""\n \\\\[-1.0ex] \\textit{Notes:} \t """)
      for r ∈ length(notes)
          write(b, "\t&\t \\multicolumn{$numContentCols}{l}{$(notes[r])}\n \\\\")
      end
  end
  write(b, "\t \\end{tabular}\n \\\\ \\end{table}")
  return String(b)
end


#=this version formats a table from a set of linear models
While customizable, the idea here is to have sensible default arguments
where possible.
IN: A collection of models, the rows to capture (returns an empty string in
the table if a row is not present, optionally column names, their width,
names of the rows, names of the descriptive rows, descriptive content,
a boolean switch for including signficance stars, signficance levels for the stars,
a scaling factor, the number of digits in values, and customizable
notes
OUT: A latex table string
=#
function texTable(models::Vector{CTLM}, getΣ::Function, rows::Vector{Symbol};
    titleCaption::String = "",
    caption::String = "",
    colNames::Vector{Vector{String}} = [["C$i" for i ∈ 1:length(models)]],
    contentRowNames::Vector{String} = String.(rows),
    descRowNames::Vector{String}=Vector{String}(),
    descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    arrayStretch::Float64 = 1.5,
    lineSpacer::String = "\\\\",
    summaryMathMode::Bool = true,
    widthColNames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
    alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
    widthDescContent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
    stars::Bool=true,
    starLvls::Vector{Float64} = [.9, .95, .99],
    starLegend::String = stars?DEFAULT_STAR_LEGEND:"",
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimalDigits::Int = 2,
    columnSepPt::Int = 25-length(models)*5,
    colHeaderName::Vector{String} = ["" for i::Int ∈ 1:length(colNames)],
    clearMem=false)

  numCols = length(models)
  numContentRows = length(rows)

  #initialize a vector of content matrices
  content::Vector{Matrix{String}} = [fill("",numContentRows,numCols) for i::Int ∈ 1:2]

  #Pre-allocate the vector of SE errors
  modelsσ::Vector{Vector{Float64}} =
      [Vector{Float64}(models[i].K) for i∈1:numCols]

  #pull out the β coefficients and N
  modelsβ::Vector{Vector{Float64}} = [models[i].β for i::Int ∈ 1:numCols]
  modelsN::Vector{Int} = [models[i].N for i ∈ 1:numCols]
  modelsXNames::Vector{Vector{Symbol}} = [models[i].XNames for i::Int ∈ 1:numCols]

  #get the standard errors
  for c ∈ 1:numCols
      modelsσ[c] .= sqrt.(diag(getΣ(models[c])))
      if clearMem #clears the memory of a model no longer used. Hurts performance.
          models[c] .= CTLM()
          gc()
      end
  end

  if length(starLegend) > 0
      notes = [starLegend; notes]
  end

  #get the content matrices
  getContentMatrices!( modelsβ, modelsσ, modelsXNames, modelsN, rows, content,
      stars=stars, starLvls=starLvls, scaling=scaling,
      decimalDigits=decimalDigits, starStrings=starStrings)

  return texTable(titleCaption,caption, colNames, contentRowNames, content, descRowNames,
      descContent, notes, arrayStretch=arrayStretch, lineSpacer=lineSpacer, summaryMathMode=summaryMathMode, widthColNames=widthColNames,
      widthDescContent=widthDescContent, columnSepPt=columnSepPt, colHeaderName=colHeaderName,
      alignmentColNames=alignmentColNames)
end

function texTable(models::Vector{CT2SLS}, getΣ::Function, rows::Vector{Symbol};
    titleCaption::String = "",
    caption::String = "",
    colNames::Vector{Vector{String}} = [["C$i" for i ∈ 1:length(models)]],
    contentRowNames::Vector{String} = String.(rows),
    descRowNames::Vector{String}=Vector{String}(),
    descContent::Vector{Vector{String}}=Vector{Vector{String}}(),
    notes::Vector{String} = Vector{String}(),
    arrayStretch::Float64 = 1.5,
    lineSpacer::String = "\\\\",
    summaryMathMode::Bool = true,
    widthColNames::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(colNames[i])),1:length(colNames)),
    alignmentColNames::Vector{Vector{String}} = #contains the number of columns for each entry
      broadcast((i::Int)->["r" for i ∈ 1:length(colNames[i])],1:length(colNames)),
    widthDescContent::Vector{Vector{Int}} =
        broadcast((i::Int)->ones(Int,length(descContent[i])),1:length(descContent)),
    stars::Bool=true,
    starLvls::Vector{Float64} = [.9, .95, .99],
    starLegend::String = stars?DEFAULT_STAR_LEGEND:"",
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)),
    decimalDigits::Int = 2,
    columnSepPt::Int = 25-length(models)*5,
    colHeaderName::Vector{String} = ["" for i ∈ 1:length(colNames)],
    clearMem = false)

    numCols = length(models)
    numContentRows = length(rows)

    #initialize a vector of content matrices
    content::Vector{Matrix{String}} = [fill("",numContentRows,numCols) for i ∈ 1:2]

    #Pre-allocate the vector of SE errors
    modelsσ::Vector{Vector{Float64}} =
    [Vector{Float64}(models[i].K+models[i].KW) for i∈1:numCols]


    #pull out the β coefficients and N
    modelsβ::Vector{Vector{Float64}} = [models[i].δ2 for i ∈ 1:numCols]
    modelsN::Vector{Int} = [models[i].N for i ∈ 1:numCols]
    modelsXWNames::Vector{Vector{Symbol}} = [[models[i].XNames; models[i].WNames]  for i ∈ 1:numCols]

    #get the standard errors
    for c ∈ 1:numCols
        modelsσ[c] .= sqrt.(diag(getΣ(models[c])))
        if clearMem #clears the memory of a model no longer used. Hurts performance.
            models[c] .= CT2SLS()
            gc()
        end
    end

    if length(starLegend) > 0
        notes = [starLegend; notes]
    end

    #get the content matrices
    getContentMatrices!( modelsβ, modelsσ, modelsXWNames, modelsN, rows, content,
        stars=stars, starLvls=starLvls, scaling=scaling,
        decimalDigits=decimalDigits, starStrings=starStrings)

    return texTable(titleCaption, caption, colNames, contentRowNames, content, descRowNames,
        descContent, notes, arrayStretch=arrayStretch, lineSpacer=lineSpacer, summaryMathMode=summaryMathMode, widthColNames=widthColNames,
        widthDescContent=widthDescContent, columnSepPt=columnSepPt, alignmentColNames=alignmentColNames)
end

#=getContentMatrices!
this is a bit of utiltiy code for making tex tables. It is not model specific
hence why it was extracted into its own function
#IN: a vector of βs, σs, names, rows (coefficeints selected). Optional parameters
include a pre-allcoation of the string matrix, a switch for the inclusion of
stars, a scaling factor, and the number of digits (rounding level)=#
#OUT: Writes to and returns the content matrix
function getContentMatrices!(modelsβ::Vector{Vector{Float64}},
    modelsσ::Vector{Vector{Float64}},
    modelsXNames::Vector{Vector{Symbol}},
    modelsN::Vector{Int},
    rows::Vector{Symbol},
    content::Vector{Matrix{String}} = #will hold the coefficients and errors
        [fill("",length(rows),length(modelsCoef)) for i ∈ 1:2];
    stars::Bool=true, #whether to display signficance stars
    starLvls::Vector{Float64} = [.9, .95, .99],  #cutoffs for signficance (must be sorted)
    starStrings::Vector{String} =
      ["\\ensuremath{^\\text{*}}","\\ensuremath{^\\text{**}}","\\ensuremath{^\\text{***}}"],
    scaling::Vector{Float64}=ones(length(rows)), # an optional scaling factor
    decimalDigits::Int = 2) #number of decimal digits



    #iterate through the models
    for c ∈ 1:length(modelsβ)

        #build a dictionary of the names
        XNameTbl::Dict = Dict(modelsXNames[c][i] => i for i::Int ∈ 1:length(modelsβ[c]))

        for r ∈ 1:length(rows)
            if haskey(XNameTbl, rows[r]) #need to check if it exists
                ind::Int = XNameTbl[rows[r]]
                p::Float64 =
                    cdf(TDist(modelsN[c]), modelsβ[c][ind]/modelsσ[c][ind]) #get CDF from T distribution
                p = p > .5 ? 1-(1.0 - p)*2.0 : 1.0 - p*2.0 #calc 2-tailed p value
                sigLevel::Int = sum(p.>starLvls)
                if sigLevel > 0 && stars
                    starString::String = starStrings[sum(p.>starLvls)]
                else
                    starString = ""
                end

                #scale, round and write the β coefficeint and σ into the string matrices
                content[1][r,c] =
                  #"$(round(scaling[r]*modelsβ[c][ind],decimalDigits))^{$starString}"
                  "\$$(round(scaling[r]*modelsβ[c][ind],decimalDigits))\$$starString"
                content[2][r,c] =
                  #"($(round(scaling[r]*modelsσ[c][ind],decimalDigits)))"
                  "(\$$(round(scaling[r]*modelsσ[c][ind],decimalDigits))\$)"
            end
        end
    end

  return content
end

#=Writes the fully formed string to the file
  IN: A string with the output, optionally a path to the file, and the
  output name
  OUT: Writes the string to a file=#
function writeTables2File(tablesString::String;
    path::String=pwd(),
    outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

    oStream::IOStream = open("$path\\$outName","w+")
    write(oStream, tablesString)
    close(oStream)
end

#=Writes the table strings to a file
  Convenience method which takes an array of table strings and optionally
  header and footer strings and writes the table to a tex file
  IN: A string array of table strings, optionally a header string, a footer
  string, a path to the file, the output name
  OUT: Writes the tables to a file=#
function writeTables2File(tables::Vector{String};
    header::String = "\n", footer::String = "",
    path::String=pwd(),
    outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

    b::IOBuffer = IOBuffer() #make an IO buffer
    write(b,"\n$header\n\n")
    for t ∈ tables
        write(b,"$t\n\n\n") #write each table
    end
    write(b,footer) #write the footer

    #write the output via another instance
    writeTables2File(String(b), path=path, outName = outName)
end

#=Writes the table strings to a file
Convenience method which takes an array of table strings and the locations of
header and footer files.
IN: A string array of table strings, the name of a header file, the name of a
a footer file, optionally a path to the file, the output name
OUT: Writes the tables to a file=#
function writeTables2File(tables::Vector{String}, headerName::String,
  footerName::String;
  path::String=pwd(),
  outName::String = "$path\\results$(Dates.format(now(),"yymmdd_THMS")).tex")

  iStream::IOStream = open("$path\\$headerName")
  header::String = readstring(iStream)
  close(iStream)

  iStream = open("$path\\$footerName")
  footer::String = readstring(iStream)
  close(iStream)

  #write the output via another instance
  writeTables2File(tables, path=path, outName = outName, header=header,
      footer=footer)
end

#array2String
  #=takes content of an array, formats each cell into a string, appends a
  prefix and suffix, and sends back an array of strings with the same
  dimension as the input
  IN: A general array, optionally a prefix, suffix, scaling factor, and
   set number of decimal digits
  OUT: An array of strings with the same dimensions as the input=#
function array2String(m::Array{Float64}; prefix::String = "",
  suffix::String = "", scaling::Float64 = 1.0,
  decimalDigits::Int = 2)::Array{String}

  #broadcast the prefix, value in m, and suffix and concatenate respectively
  stringMat::Array{String} =
      ((x::Float64)->(prefix*"$(round(x,decimalDigits))"*suffix)).(m)

  return stringMat
end



function IOTest()
  #test parameters
  nrowsContent = 8
  ncols = 4
  nsecs = ncols ÷ 2
  nrowsDesc = 2
  nsubRows = 2

  path = pwd() * "\\results"
  outName = "test.tex"
  footerName = "footer.tex"
  headerName = "header.tex"

  #test table
  caption = "test table"
  colNames::Vector{Vector{String}} =
      [broadcast(i->"secs$i",1:nsecs),broadcast(i->"ncols$i",1:ncols)]
  contentRowNames::Vector{String} = ["nrows$i" for i ∈ 1:nrowsContent]
  content::Vector{Matrix{String}} =
      array2String.(broadcast((i::Int)->i .* ones(nrowsContent,ncols), 1:nsubRows))
  descRowNames::Vector{String} = ["Desc. row $i" for i ∈ 1:nrowsDesc]
  descContent::Vector{Vector{String}} =
      [broadcast(i->"desc-r1c$i",1:ncols),broadcast(i->"desc-r1c$i",1:ncols)]
  notes::Vector{String} = ["a note"]
  widthColNames::Vector{Vector{Int}} =
      [broadcast(i->(ncols ÷ nsecs), 1:nsecs),broadcast(i->1, 1:ncols)]

  s=texTable(caption, colNames,  contentRowNames, content, descRowNames,
      descContent, notes, widthColNames=widthColNames)

  writeTables2File([s,s], headerName, footerName, path=path, outName = outName)

  #println(array2String(fill(.45/π,3,3,4),decimalDigits=4))

print(s)
end

#################testing methods for this module
function LMtest()

  #Number of Observations
  N::Int = 200
  K::Int = 2

  #Allocate
  X::Matrix{Float64} = Matrix{Float64}(N,K)
  Y::Vector{Float64} = Vector{Float64}(N)
  ε::Vector{Float64} = similar(Y)

  #parameters for the simulation
  σ2::Vector{Float64} = [2.0^2.0 for i::Int ∈ 1:N]
  σ2[1:ceil(Int, N/2)] .= 0.5^2.0
  β::Vector{Float64} = [(1.0 * i) for i::Int ∈ 1:K]
  #println("β: ", β)
  #println("σ2: $σ2")

  #this will hold the sampled beta
  e::Vector{Float64} = similar(ε)

  #run the simulation
  X[:,1] .= 1.0 #intercept
  for i ∈ 2:K
    X[:,i] .= rand(Uniform(),N)
  end

  ε .= map((s2::Float64)->rand(Normal(0.0,s2^0.5)),σ2)
  Y .=  X*β .+ ε

  #get the linear model
  lin::CTLM = CTLM(X, Y)

  #get the homoskedastic SEs
  ΣHomosked::Matrix{Float64} = getHomoskedΣ!(lin)
  ΣHomoskedSlow::Matrix{Float64} = getHomoskedΣSlow(lin)

  #print the coefficients
  println("Coefficients: ",lin.β)
  println("Homoskedastic Errors: ", diag(ΣHomosked).^.5)
  println("Check: ", diag(ΣHomoskedSlow).^.5)

  #get the modified white SEs
  ΣWhite::Matrix{Float64} = similar(ΣHomosked)
  getModWhiteΣ!(lin, ΣWhite)
  ΣWhiteSlow::Matrix{Float64} = getModWhiteΣSlow(lin)

  #print the coefficients
  println("Modified White Errors: ",diag(ΣWhite).^.5)
  println("Check: ", diag(ΣWhiteSlow).^.5)

    #=Π1::Matrix{Float64} = Matrix{Float64}(K,2)
  getCoeff!(lin.xqr, [Y 2.*Y], Π1)
  println("1st stage Z: $(lin.X)")
  println("1st Stage X: $([Y 2.*Y])" )
  println("1st Stage Test Coef: $Π1" )
  Ξ::Matrix{Float64} = getResid!(lin.xqr, [Y 2.*Y], Π1)
  println("1st Stage Resid: $Ξ" )=#


  #test the project routines
  Pa::Vector{Float64} = Vector{Float64}(N)
  PM::Matrix{Float64} = Matrix{Float64}(N,N)
  PS::Matrix{Float64} = Matrix{Float64}(N,N)

  project!(X,PM)
  project!(X,Pa)

  projectSlow!(X,PS)
  println("P: ", Pa[1:5])
  println("P (from full matrix): ", PM[1:3,1:3])
  println("P Slow: ", diag(PS)[1:10])
end

function IVTest()

  srand(11)

  N::Int = 200 #Number of Observations
  K::Int = 3 # of endogneous covariates
  L::Int = 4 # of instruments
  KW::Int = 5# of exog covariates

  #Allocate
  X::Matrix{Float64} = Matrix{Float64}(N,K)
  Z::Matrix{Float64} = Matrix{Float64}(N,L)
  W::Matrix{Float64} = Matrix{Float64}(N,KW)
  Y::Vector{Float64} = Vector{Float64}(N)
  A::Vector{Float64} = Vector{Float64}(N)

  Π1::Matrix{Float64} = Matrix{Float64}(L+KW,K)
  π2::Vector{Float64} = Vector{Float64}(L+KW)
  Ξ1::Matrix{Float64} = similar(X)
  e2::Vector{Float64} = similar(Y)

  #parameters for the simulation
  σ2::Vector{Float64} = [2.0^2 for i::Int ∈ 1:N]
  σ2[1:ceil(Int, N/2)] .= 0.5^2.0

  #generate π2 coefficients (pattern is 1 times the index, ie 1 2 3 ...)
  π2 .= [(10*i) for i::Int ∈ 1:(L+KW)]
  π2[L+1] .= 1.0 #intercept term
  #println("π2: ",π2)

  #generate Π1, pattern is δ2 times 1/10 the index, ie ie .1δ2 .2δ2 .3δ2 ...)
  for i ∈ 1:(L+KW), j∈1:K
    Π1[i,j] .= i + j^2
  end
  Π1[L+1,:] .= 1.0 #set the intercept term
  #display(Π1)

  #generate the exogeneous variables
  for i ∈ 2:KW
    W[:,i] .= rand(Uniform(-10,10),N)
  end
  W[:,1] .= 1.0 #intercept

  #generate the instrumental variables
  for i ∈ 1:L
    Z[:,i] .= rand(Uniform(-10,10),N)
  end

  #map((s2::Float64)->rand(Normal(0.0,s2^0.5)),σ2)
  #generate the errors
  for i ∈ 1:N
    for j ∈ 1:K
      Ξ1[i,j] = rand(Normal(0.0,σ2[i]^0.5))
    end
    e2[i] = rand(Normal(0.0,σ2[i]^0.5))
  end

  #println(Ξ1)
  #populate Y and X
  X .= [Z W]*Π1 .+ Ξ1
  Y .= [Z W]*π2 .+ e2

  #add some endogeneity
  for i ∈ 1:K
    A .= rand(Uniform(-i,i),N)
    X[:,i] .+= A
    Y .+= A .* 10
  end


  #probably is a more graceful way to do this
  XNames::Vector{Symbol} = [Symbol("X$i") for i = 1:K]
  WNames::Vector{Symbol} = [Symbol("W$i") for i = 1:KW]
  ZNames::Vector{Symbol} = [Symbol("Z$i") for i = 1:L]


  iv::CT2SLS = CT2SLS(X,W,Y,Z, XNames=XNames, YName=:Y, WNames=WNames, ZNames=ZNames)
  println("coefficients: \n", [XNames; WNames] ,"\n", iv.δ2)

  #get the homoskedastic SEs
  ΣHomosked::Matrix{Float64} = getHomoskedΣ!(iv)
  ΣHomoskedSlow::Matrix{Float64} = getHomoskedΣSlow(iv)

  #print the coefficients
  println("Homoskedastic Errors: ", diag(ΣHomosked).^.5)
  println("Check: ", diag(ΣHomoskedSlow).^.5)

  #get the modified white SEs
  ΣWhite::Matrix{Float64} = getModWhiteΣ!(iv)
  ΣWhiteSlow::Matrix{Float64} = getModWhiteΣSlow(iv)

  #print the coefficients
  println("Modified White Errors: ",diag(ΣWhite).^.5)
  println("Check: ", diag(ΣWhiteSlow).^.5)

  iv1st::Vector{CTLM} = get1stStage(iv)
  println("1st Stage results: \n $(iv1st[1].XNames)")
  for iv1 ∈ iv1st
      println("coefficients $(iv1.YName):  $(iv1.β )")
      println("Error (Homosked $(iv1.YName)): $(diag(getHomoskedΣ!(iv1)).^.5)")
  end



  oStream::IOStream = open(
    "C:\\Users\\Clinton\\Dropbox\\Projects\\Summer17RTester\\Summer17RTester\\IVTestOut.csv","w+")


  write(oStream, [string.(iv.XNames); string.(iv.WNames)]::Vector{String})
  writecsv(oStream, [X W Y Z])
  close(oStream)
end



if DEBUG_CTMOD
  #LMtest()
  IVTest()
  #IOTest()
  #print(perfTestStr(10^6))
end



#end module
end
