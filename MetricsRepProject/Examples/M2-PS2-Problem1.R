set.seed(11);
rm(list = ls());

require(devtools);
require(DataAnalytics);
require(plyr);
require(ggplot2);

################################Additional helper functions###################

#in: a regression object
#out: a matrix object for the coefficients
GetCoefAsMatrix = function(reg) matrix(coef(reg),
    nrow = length(coef(reg)),
    dimnames = list(names(coef(reg))));

#Obtains the modified white standard errors given X and residuals
#In: A QR decomposition object for X, a vector of residuals
#Out: a vector of white standard errors
GetMWhiteErrors = function(qrX, vResid) {
    #extract the temporary variables
    mX = qr.X(qrX);
    mR = qr.R(qrX);
    mQ = qr.Q(qrX);

    mQRTInv = mQ %*% t(backsolve(mR,diag(ncol(mR))));
    return(t(mQRTInv * rep((vResid / (1 - rowSums(mQ * mQ))) ^ 2, ncol(mX))) %*% mQRTInv);
}


#in: a matrix of X values and a vector of residuals
#out: a a matrix of white standard errors 
IneffGetMWhiteErrors = function(mX, vResid) {
    #extract the temporary variables

    return(solve(t(mX) %*% mX) %*% t(mX) %*%
        diag(as.vector((vResid / (1 - diag(mX %*% solve(t(mX) %*% mX) %*% t(mX)))) ^ 2)) %*%
             mX %*% solve(t(mX) %*% mX));

}

####################################Script Entry Point##########################
Problem1Script = function() {
    #Create a test data set
    iNumPoints = 100;
    iVars = 3;

    #simulate our data and add some noise
    mX = matrix(rnorm(iNumPoints*iVars),iNumPoints);
    vY = mX %*% matrix(1:iVars, iVars) + rnorm(iNumPoints) * 10;
    m1X = cbind(matrix(1, iNumPoints), mX);
    qrX = qr(m1X);

    reg = lm(vY ~ mX);
    vResid = vY - m1X %*% (GetCoefAsMatrix(reg = reg));

    ptm = proc.time();

    print("Efficient White Errors");
    print(GetMWhiteErrors(qrX = qrX, vResid = vResid));
    ptm2 = proc.time();
    print(ptm2 - ptm);

    print("Inefficient White Errors");
    print(IneffGetMWhiteErrors(mX = m1X, vResid = vResid));
    print(proc.time() - ptm2);

    return(1);
}

Problem1Script();







