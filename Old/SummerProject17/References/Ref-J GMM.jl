module CTProj3

using DataFrames, JuMP,  NLopt, Distributions, Gadfly#Plots
#pyplot()
flush(STDOUT)
Gadfly.push_theme(:dark)

#=markdown path
weave(Pkg.dir("C:\\Users\\Clinton\\Dropbox\\Coursework\\Theory\\PS3\\GMM.jl"),
  informat="script",
  out_path = "C:\\Users\\Clinton\\Dropbox\\Coursework\\Theory\\PS3\\outGMM.html",
  doctype = "md2html")
=#

#start the program
#in/out: nothing
function p5GMM()
  dfGMM= readtable("psgmmdata.csv")

  #get matrix for the dependent variable
  vDV = convert(Array,dfGMM[:ct_ct_1])
  vRM = Array(dfGMM[:,[:rm_rf]])
  vHML = Array(dfGMM[:,[:hml]])
  vSMB = Array(dfGMM[:,[:smb]])
  dERM = mean(dfGMM[:rm_rf])

  println("Means:")
  display(mean(Array(dfGMM[:,[:ct_ct_1,:rm_rf,:hml,:smb]]),1))

  println("Correlation Matrix:")
  display(cor(Array(dfGMM[:,[:ct_ct_1,:rm_rf,:hml,:smb]])))

  #this is a vectorized function for getting the g and d vectors
  ff(b) = vDV.^(-b)
  fdf(b) = -vDV.^(-b).*log(vDV)
  ffm(b,vM) = vM'*vDV.^(-b)
  vp = zeros(length(vDV))

  #lets do some graphs!
  γ=0:.1:100

  plOut = plot(melt(DataFrame(γ=γ,Rm=((x->(vRM'*ff(x))[1]).(γ)),
      HML=((x->(vHML'*ff(x))[1]).(γ)),
      SMB=((x->(vSMB'*ff(x))[1]).(γ))),:γ),
    x=:γ,y=:value,color=:variable,Geom.line,
    Guide.title("E(Ct+1/Ct*Re) vs γ"),
    Guide.xlabel("γ"),Guide.ylabel("E(Ct+1/Ct*Re)"),
    Guide.colorkey("Legend"))
  draw(PNG("InitialPlots.png", 7inch, 5inch), plOut)

  display(plOut)

  testSingleGMMEst!(ff,fdf,vp,vRM,
    eye(1),dERM , test="Rm Only")
  testSingleGMMEst!(ff,fdf,vp,vHML,
    eye(1),dERM , test="Hml Only")
  testSingleGMMEst!(ff,fdf,vp,vSMB,
    eye(1),dERM , test="SMB Only")
  testSingleGMMEst!(ff,fdf,vp,Array(dfGMM[:,[:rm_rf,:hml]]),
    eye(2),dERM , test="Rm_HML")
  testSingleGMMEst!(ff,fdf,vp,Array(dfGMM[:,[:rm_rf,:smb]]),
    eye(2),dERM , test="Rm_SMB")



end

#a wrapper for each of the moment sets
#in: dependent variable function, its derivtive, moment vectors, W matrix,
#the name of the test
#out: value 1 for success. Also, W is changed by reference
function testSingleGMMEst!(ff,fdf,vp,mMom, mW, dERM; test="")

  #some intermediate variables for future calculations
  γ=0:.1:100

  println("Running first stage for: ",test)
  db1,dVar1,valpha,valphaVar,dchi2,dchi2p,vgWg1,va =
     singleGMMEst!(ff,fdf,vp,mMom,mW,dERM=dERM,pltRng=γ)
  print("\nFirst stage ",test, " results\n",
    "\nγ: ", db1,
    "\nσ: ", dVar1^0.5,
    "\nχ2 Stat: ", size(mMom,2)>1 ? dchi2 : "NA",
    "\nχ2 p Val: ", size(mMom,2)>1 ? dchi2p : "NA",
    "\nα (E[mRe]-E[Re]): ",valpha,
    "\nα Var: ",valphaVar,
    "\nW (normalized by 1/mean(|S|)): ", mW,
    "\na: ", va,
    "\n")

  println("Running second stage for: ",test)
  db2,dVar,valpha,valphaVar,dchi2,dchi2p,vgWg2,va =
     singleGMMEst!(ff,fdf,vp,mMom,mW,dERM=dERM,pltRng=γ)

  print("\nSecond stage ",test, " results",
    "\nγ: ", db2,
    "\nσ: ", dVar^0.5,
    "\nχ2 Stat: ", size(mMom,2)>1 ? dchi2 : "NA",
    "\nχ2 p Val: ", size(mMom,2)>1 ? dchi2p : "NA",
    "\nα (E[mRe]-E[Re]): ",valpha,
    "\nα Var: ",valphaVar,
    "\nW (normalized by 1/mean(|S|)): ", mW,
    "\na: ", va,
    "\n")

  if(db1!=db2)
    plOut = plot(melt(DataFrame(γ=γ,FirstStage=vgWg1,SecondStage=vgWg2),:γ)
      , x=:γ,y=:value,color=:variable,Geom.line,
      Guide.title(test),
      Guide.xlabel("γ"),Guide.ylabel("gWg"),
      Guide.colorkey("Legend"))
    draw(PNG(test*".png", 7inch, 5inch), plOut)
    display(plOut)
  else
    plOut = plot(melt(DataFrame(γ=γ,gWg=vgWg1),:γ)
      , x=:γ,y=:value,Geom.line,
      Guide.title(test),
      Guide.xlabel("γ"),Guide.ylabel("gWg"))
    draw(PNG(test*".png", 5inch, 5inch), plOut)
    display(plOut)
  end

  return (1)
end

#gets the GMM estiamte fo a single parameter
#in:  dependent variable function, its derivtive, moment vectors, W matrix
#out: parameter value, its variance, alpha, its variance,
#     chisq stat (2> moments only), its p value,
#     updated weighting W (by ref)
function singleGMMEst!(ff,fdf,vp,mMom,mW; dERM = 0, pltRng = 0:1:100)

  #Different optimization models work: GN_ISRES, LN_COBYLA, LN_NELDERMEAD,
  #LN_BOBYQA, # GN_CRS2_LM, GN_DIRECT, LD_CCSAQ, LD_SLSQP,
  # LD_TNEWTON_PRECOND_RESTART, LD_LBFGS. Many others do not work.
  # I use LD_LBFGS due to prevelence and relative simplicity
  # I think LN_COBYLA is also a good choice
  modOp=Model(solver=NLoptSolver(algorithm=:LD_LBFGS))
  #modOp=Model(solver=NLoptSolver(algorithm=:LN_COBYLA))
  iN = Int(size(mMom,1))
  iM= Int(size(mMom,2))

  vg = Vector{Float64}(iM)
  vEp = sum(vp)/iN
  vf = Vector{Float64}(iN)

  #note this could be made faster through devectorization or BLAS routines
  #vectorized version
  function fgWg(b)
    vg = (mMom'*ff(b))./iN.-vEp
    return (vg'*mW*vg)[1]
  end

  #start the optimization macro
  @variable(modOp,-100<=jb<=100, start=1)
  JuMP.register(modOp,:ff,1,ff,autodiff=true)
  JuMP.register(modOp,:fgWg,1,fgWg,autodiff=true)
  @NLobjective(modOp, Min, fgWg(jb))
  output = solve(modOp)

  #get outputs and an average of f
  dbVal=getvalue(jb)
  vfBar=(mMom'*ff(dbVal)./iN)
  vf = ff(dbVal)
  vg = vfBar-vEp

  #fully devecotrized calculation for s
  mS=zeros(iM,iM)
  for i=1:iN
    for j=1:iM #this loop adds the de-meaned  covar matrix at the ith data point
      for k = 1:iM
        mS[k,j]+=(mMom[i,k]*vf[i]-vfBar[k])*(vf[i]*mMom[i,j]-vfBar[j])
      end
    end
  end
  mS=mS/iN

  #calculate the standard errors
  vd = (mMom'*fdf(dbVal))./iN
  ddWdInv = 1/((vd'*mW*vd)[1])

  #calculate the alpha stats
  valpha = vg #- mMom'*ones(iN)./iN #Not sure about this part
  vWd_dwdInv_d = mW*vd*ddWdInv*vd'
  valphaVar = (eye(iM)-vWd_dwdInv_d')*mS*(eye(iM)-vWd_dwdInv_d)/iN

  #Get the stats for the identification test
  if iM>1
    dchi2 = (vg'*(pinv(valphaVar))*vg*(iM-1))[1]
    dchi2p = 1-cdf(Chisq(iM-1),dchi2)
  else
    dchi2 = 0.0
    dchi2p = 0.0
  end

  #recalc the weighting matrix
  dNormConst = (1/abs(mean(mS))) #improves numerical stability/scaling
  vgWg = fgWg.(pltRng)
  mW.=(mS*dNormConst)\eye(iM)

  ddWdInv = 1/((vd'*mW*vd)[1])
  dse2 = ddWdInv^2*(vd'*mW*mS*mW*vd)[1]/iN

  va = mW*vd.*dNormConst

  return dbVal, dse2, valpha,valphaVar,dchi2,dchi2p, vgWg, va
end

#call the main wrapper function
p5GMM()

end
