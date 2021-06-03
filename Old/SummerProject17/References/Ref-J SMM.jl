module CTProj

using DataFrames, JuMP,  NLopt, Distributions,
  StatsBase, Base.LinAlg.BLAS, Gadfly
#pyplot()

#=markdown path
weave(Pkg.dir("C:\\Users\\Clinton\\Dropbox\\Coursework\\Theory\\PS3\\SMM.jl"),
  informat="script",
  out_path = "C:\\Users\\Clinton\\Dropbox\\Coursework\\Theory\\PS3\\outSMM.html",
  doctype = "md2html")
=#

#start the program
#in/out: nothing
function p5SMM()
  flush(STDOUT)
  #ProcessCompData() #uncomment for data analysis


  const iNSims=1250
  const iBurnIn = 250
  vState=rand(Float64,iNSims,1)
  vStateTemp=rand(Float64,iNSims,1)
  const iSimsForSE = 100

  vMom = ones(4)

  println("Running at θ=0.6:")
  converge!(0.6, vMom = vMom, iNSims=iNSims,iBurnIn=iBurnIn, vState=vState)

  println("\nMean I/K: ", vMom[1],
    "\nσ(I/K): ", vMom[2],
    "\nMean Q: ",vMom[3],
    "\nσ(Q): ",vMom[4])

  mW = eye(4)

  #first stage
  println("Running First Stage")
  getθ(iNSims,iBurnIn,vState,mW,verbose=true)

  #get the standard errors
  mB = zeros(iSimsForSE,4)
  for i in 1:iSimsForSE
    mB[i,:] .= (getθ(iNSims,iBurnIn,rand(Float64,iNSims,1),mW,verbose=false))[:]
  end

  mB .= mB .- ones(iSimsForSE)*
    [mean(mB[:,1]) mean(mB[:,2]) mean(mB[:,3]) mean(mB[:,4])]

  mS = cov(mB)
  dnormConst = mean(abs.(mS))
  mS .= mS/dnormConst
  mW = mS\eye(4)
  println("Normalized Weighting Matrix:\n\n", mW)

  #second stage
  println("Running Second Stage")
  getθ(iNSims,iBurnIn,vState,mW,verbose=true)




end

function getθ(iNSims,iBurnIn,vState::Array{Float64,2},mW::Array{Float64,2};
  verbose=false)
  #the moments from the data

  const dIKMean = 0.26215608182706984
  const dIKσ = 0.18066965524433168
  const dQMean = 1.7303211980649078
  const dQσ = 0.8678322338204875
  const vMomData = [dIKMean, dIKσ, dQMean, dQσ]
  const iM = 4

  const steps = 100
  const dtol = 10.0^-7

  θLow = 0.01
  θHigh = .99
  θMid = (θHigh+θLow)/2.0

  JPrev = 0.0 #for measurement
  JMid = 0.0
  JLow = 0.0
  JHigh = 0.0

  iter = 0

  vMomSim = zeros(iM)
  vMomΔ= zeros(iM)

  #objective function, could be faster if devectorized
  #for notational convenience
  @inline function measureConverge(x::Float64)
    converge!(x, vMom = vMomSim, iNSims=iNSims,
      iBurnIn=iBurnIn, vState=vState)
    vMomΔ .= vMomData.-vMomSim
    return (vMomΔ'*mW*vMomΔ)[1]
  end

  #calculate values for initial J
  JMid = measureConverge(θMid)
  JLow = measureConverge(θLow)
  JHigh = measureConverge(θHigh)
  #println(JMid)

  #main loop
  while ((abs(JLow-JHigh) >= dtol) & (iter < 1000))

    #test points

    JHigh= measureConverge((θMid+θHigh)/2.0)
    JLow = measureConverge((θLow+θMid)/2.0)


    if JHigh<JMid #did we do better? if so, use the upper window
      θLow = θMid
      θMid = (θLow+θHigh)/2.0
      JMid = JHigh
    else
      if JLow<JMid #did we do better? if so, use the lower window
        θHigh = θMid
        θMid = (θLow+θHigh)/2.0
        JMid = JLow
      else #otherwise, just shrink the current window
        θHigh = (θMid + θHigh)/2.0
        θLow= (θMid + θLow)/2.0
      end
    end

    iter+=1

  end



  if verbose
    println("\nIter:", iter,
      "\nθ: ", θMid,
      "\nJ: ", JMid
    )

    println("\nMean I/K: ", vMomSim[1],
      "\nσ(I/K): ", vMomSim[2],
      "\nMean Q: ",vMomSim[3],
      "\nσ(Q): ",vMomSim[4])

    θ=0.01:.01:0.99
    J=map(x->min(measureConverge(x),10.0),θ)

    plOut = plot(melt(DataFrame(θ=θ,
        J=J),:θ),
      x=:θ,y=:value,color=:variable,Geom.line,
      Guide.title("J vs θ"),
      Guide.xlabel("θ"),Guide.ylabel("J"))
    draw(PNG("SMMPlotSecondStage.png", 7inch, 5inch), plOut)
    display(plOut)

  end

  return vMomSim
end

function ProcessCompData()

  #pre-processing data
  dfSMM= readtable("CompData.csv")
  display(summary(dfSMM))

  #drop unneeded columns
  dfSMM = dfSMM[:,[:gvkey,:fyear,:at,:capx,:ceq,:csho,:ppent,:txdb,:prcc_c]]

  #filter out NA values
  dfSMM = dfSMM[(~isna(dfSMM[:,:gvkey])&
    ~isna(dfSMM[:,:ppent])),:]

  dfSMM[:,:txdb] .= map!(((x)->isna(x)?0.0:x),dfSMM[:,:txdb])

  #calculate the market value
  dfSMM[:,:MV] = dfSMM[:,:prcc_c].*dfSMM[:,:csho]
  sort!(dfSMM,cols = (:gvkey,:fyear))
  dfSMM[:,:lppent] = zeros(size(dfSMM,1))

  #make sure it is ok do to do the lagging procedure
  dfSMM[:,:ycnt] = 0
  dcnt = 0

  by(dfSMM, :gvkey) do dfSelect
    dcnt = size(dfSelect,1)
    if maximum(dfSelect[:,:fyear])-minimum(dfSelect[1,:fyear])!=dcnt-1
      dcnt = 0
    end
    dfSelect[:,:ycnt] = dcnt
  end

  dfSMM = dfSMM[dfSMM[:,:ycnt].>1,:]

  by(dfSMM, :gvkey) do dfSelect
    dfSelect[:,:lppent] = [0.0; dfSelect[1:end-1,:ppent]]
  end

  #filter out more NA values
  dfSMM = dfSMM[(~isna(dfSMM[:,:csho])&
    ~isna(dfSMM[:,:prcc_c])&
    ~isna(dfSMM[:,:at])&
    ~isna(dfSMM[:,:capx])&
    ~isna(dfSMM[:,:ceq])),:]

  #filter out non-sensical inputs
  dfSMM = dfSMM[dfSMM[:,:MV].>=100,:]
  dfSMM = dfSMM[dfSMM[:,:at].>.001,:]
  dfSMM = dfSMM[dfSMM[:,:ceq].>=0.001,:]
  dfSMM = dfSMM[dfSMM[:,:lppent].>0.001,:]

  #Calculate a Q
  dfSMM[:,:Q] = (dfSMM[:,:at].+dfSMM[:,:MV].-dfSMM[:,:ceq].-
    dfSMM[:,:txdb])./dfSMM[:,:at]
  dfSMM = dfSMM[dfSMM[:,:Q].>=0,:]

  dfSMM[:,:IOverK] =  dfSMM[:,:capx]./dfSMM[:,:lppent]
  dfSMM = dfSMM[:,[:gvkey,:fyear,:at,:capx,:ceq,:MV,:lppent,:txdb,:Q,:IOverK]]

  #windsorization
  dIKLower = percentile(dfSMM[:,:IOverK],5)
  dIKUpper = percentile(dfSMM[:,:IOverK],95)
  dQLower = percentile(dfSMM[:,:Q],5)
  dQUpper = percentile(dfSMM[:,:Q],95)

  dfSMM = dfSMM[dfSMM[:,:IOverK].>dIKLower,:]
  dfSMM = dfSMM[dfSMM[:,:IOverK].<dIKUpper,:]
  dfSMM = dfSMM[dfSMM[:,:Q].>dQLower,:]
  dfSMM = dfSMM[dfSMM[:,:Q].<dQUpper,:]

  #debug stuff
  println("IOverK mean: ", mean(dfSMM[:,:IOverK]))
  println("Q mean: ", mean(dfSMM[:,:Q]))
  println("IOverK sd: ", var(dfSMM[:,:IOverK])^.5)
  println("Q sd: ", var(dfSMM[:,:Q])^.5)

  display(summary(dfSMM))

end

function converge!(θ::Float64; vMom::Vector{Float64}=[1.0,1.0,1.0,1.0],
  iNSims=1250, iBurnIn = 250,vState::Array{Float64,2}=rand(Float64,iNSims,1),
  iSteps=50)

  θ = (θ==1.0?1.01:θ)
  #println("θ: ", θ)

  const r=0.04
  const β=1/(1+r)
  #const θ=.6
  const δ=0.15

  const vy=[exp(-0.2);exp(0.2)]
  const mpiy=[.6 .4;.4 .6]
  const mpiyβ=β.*mpiy
  const iM = length(vy)

  #the approximate steady state is:
  const dkss = ((r+δ)/(θ*0.5*(exp(-0.2)+exp(0.2))))^(1/(θ-1))
  #println("kss=",dkss)
  #finite dimensional state space
  const dStep=dkss/(iSteps/2.0)
  const dbmax=2.0*dkss

  vb=dStep:dStep:dbmax    #grid steps for bonds are .1*avg. income

  ################begin: My Code
  dTol=10.0^-8.0
  iN = length(vb)

  #first we create our operator
  m3R = ones(iN,iN,iM)
  mI = ones(Int,iN,iM)
  vOneiN = ones(iN)
  mZeroiN = zeros(iN,iN)
  mTPlusV = Array{Float64}(iN,iN)

  for j = 1:iM
    #display(ones(iN,1)*((vb).^θ.*vy[j]+vb*(1-δ))'-vb*ones(1,iN))
    m3R[:,:,j].=ones(iN,1)*(vb.^θ.*vy[j]+vb.*(1-δ))'-vb*ones(1,iN)
  end

  map!((x)->x<0?-10^10:x,m3R) #switch negative consumption to a penalty
  #display(m3R)

  dErr = 1.0
  iIter = 1
  mV = ones(iN,iM)
  mVm1 = ones(iN,iM)

  #BLAS.set_num_threads(12)

  #tic()
  while dErr> dTol && iIter < 10000
      #to get policy and V
      for j = 1:iM
          mTPlusV.=m3R[:,:,j]
          for k = 1:iM
            ger!(mpiyβ[j,k],mVm1[:,k],vOneiN,mTPlusV) #BLAS optimization αxy'+M
          end
          mV[:,j],mI[:,j] = findmax(mTPlusV,1) #functions like matlab max
      end

      dErr = maximum(abs(mVm1.-mV))
      mVm1.=mV

      iIter = iIter + 1;
  end
  #toc()

  mI .= ((mI-1) .% iN)+1  #prior indices were linear, this fixes that
  #display(iIter)

  #Calc our stats by MC

  iState=Int8(2)
  ibCur=Int16(round(iN/2))
  ibNext=Int16(1)
  vStateHist = ones(Int8,iNSims,1)
  vbHist = ones(iNSims,1)

  vInvestment = ones(iNSims,1)
  vIOverK=ones(iNSims,1)
  vQ = ones(iNSims,1)

  #assign the state
  fNextState(x::Float64,s::Int8) =  Int8(x<mpiy[s,1]?1:2)

  for i=1:iNSims
    #update the state

    iState = fNextState(vState[i],iState)

    vStateHist[i] = iState
    ibNext = mI[ibCur,iState]

    vIOverK[i] = (vb[ibNext]-(1-δ)*vb[ibCur])
    vQ[i] = mV[ibCur,iState]/vb[ibCur]

    vbHist[i] = vb[ibNext]
    ibCur = ibNext
  end
  vIOverK[2:end].=vIOverK[2:end]./vbHist[1:(end-1)]

  vbHist = vbHist[iBurnIn:end]
  vStateHist = vStateHist[iBurnIn:end]
  vIOverK = vIOverK[iBurnIn:end]
  vQ = vQ[iBurnIn:end]

  vMom[1] = mean(vIOverK)
  vMom[2] = var(vIOverK)^.5
  vMom[3] = mean(vQ)
  vMom[4] = var(vQ)^.5

end

p5SMM()

end
