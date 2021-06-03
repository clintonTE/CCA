rm(list = ls());

require(devtools);
require(DataAnalytics);
require(plyr);
require(ggplot2);
require(Matrix);
require(AER);
require(sandwich);

set.seed(11);
#allows for reproducability

#adapted from Rossi tsls funciton
#In: focal variables mX, exogeneous covariates mW, intstrumental variables mZ, y variables
#out: coeficients of the regression and a matrix of modified white standard errors
twoSLS = function(X, W, Z, Y) {
    #assumes both Z and W are exogeneous
    Za = cbind(Z, W);
    RInv = backsolve(chol(crossprod(Za)), diag(ncol(Za)));
    PZa = tcrossprod(Za %*% RInv);
    XTild = cbind(PZa %*% X, W);
    XspXsinv = chol2inv(chol(crossprod(XTild)));

    XW = cbind(X, W);
    B = XspXsinv %*% crossprod(XTild, Y);

    rownames(B) = c(colnames(X), colnames(W));
    resid = Y - XW %*% B;
    WhiteErrors = GetMWhiteErrors(XTild, resid);
    colnames(WhiteErrors) = rownames(B)

    return(list(
    Coef = B,
    MWSE = WhiteErrors,
    HSSE = rep(var(resid), ncol(XW)) * XspXsinv));
}

#in: X variables and residuals
#out: Modified White stnadard errors
GetMWhiteErrors = function(X, vResid) {
    #extract the temporary variables
    qrX = qr(X);
    mX = qr.X(qrX);
    mR = qr.R(qrX);
    mQ = qr.Q(qrX);

    mQRTInv = mQ %*% t(backsolve(mR, diag(ncol(mR))));
    return(t(mQRTInv * rep((vResid / (1 - rowSums(mQ * mQ))) ^ 2, ncol(mX))) %*% mQRTInv);
}


#in: a regression object
#out: a matrix object for the coefficients
GetCoefAsMatrix = function(reg) matrix(coef(reg),
    nrow = length(coef(reg)),
    dimnames = list(names(coef(reg))));

Rtsls = function(X, W, Z, y) {
    # model y=Xbeta + Wgamma + epsilon
    # Z is a matrix of instruments
    Za = cbind(Z, W) # note that the set of "instruments" is all exog var
    Rinv = backsolve(chol(crossprod(Za)), diag(ncol(Za)))
    Pza = tcrossprod(Za %*% Rinv) # project matrix X on Za, Za(Za'Za)-1Za'
    Xstar = cbind(Pza %*% X, W)
    XspXsinv = chol2inv(chol(crossprod(Xstar)))
    coef = XspXsinv %*% crossprod(Xstar, y)
    s = sqrt(sum((y - cbind(X, W) %*% coef) ^ 2) / length(y))
    stderr = s * sqrt(diag(XspXsinv))
    list(s = s, coef = coef, stderr = stderr)
}

testDat = read.csv("RHFData.csv", header = TRUE)
head(testDat)

iNumPts = length(testDat) #Number of Observations
#K = 3 # of endogneous covariates
#L = 4 # of instruments
#KW = 5 # of exog covariates

#X = data.matrix(testDat[1:K])
#W = data.matrix(testDat[(K+1):(K+KW)])
#Y = testDat$Y
#Z = data.matrix(testDat[(K + 2 + KW):(K + 1 + KW+L)])

#formula is X vars plus intercept
WSpec = paste(paste(paste("performanceLag", 1:12, sep = ""), collapse = " + "),
    "+ fDate + main_strategy")

specLM = as.formula(paste("performance ~ negLFlows + ", WSpec))

OLS = lm(specLM, data = testDat)
OLSCoef = GetCoefAsMatrix(OLS)
OLSResid = resid(OLS)

print("OLS Coeficients")
print(t(OLSCoef))

print("OLS Homoskedastic SE")
print(sqrt(diag(vcov(OLS))))


#print("OLS Modified White Standard Errors")
#OLSMWSE = sqrt(diag(GetMWhiteErrors(cbind(X, W), OLSResid)))
#print(OLSMWSE)

#Prints in order of X, W, intercept
#lSLSOut = twoSLS(X = X, W = W, Z = Z, Y = Y)

##Coeficients:
#print("2SLS Coeficients")
#print(t(lSLSOut$Coef))

#print("My MWSE:")
#print(sqrt(diag(lSLSOut$MWSE)))

print("2SLS")
spec2SLS = as.formula(paste("performance ~ negLFlows + ", WSpec, "| noticeLag1 + ", WSpec))



#spec2SlS = as.formula("Y ~ X1 + X2 + X3 + W1 + W2 + W3 + W4 + W5")
iv = ivreg(formula = spec2SLS, data = testDat, singular.ok = TRUE)
#iv = tsls(formula = Y ~ W1 + W2 + W3 + W4 + W5 + X1 + X2 + X3 - 1, instruments = ~W1 + W2 + W3 + W4 + W5 + Z1 + Z2 + Z3 + Z4 - 1, data = testDat)
print(summary(iv))

print("2SLS Modified White Standard Error")
print(sqrt(diag(vcovHC(iv, type = "HC3"))))
#print("Rossi SE")
#Rtsls(X, W, Z, Y)$stderr

