#Allfunctions, pseudoclasses, and lower - level code for
    #    the replication




#in: X variables and residuals
#out: Modified White stnadard errors
GetMWhiteErrors = function(mRawX, vResid, SEOnly = TRUE, AddConst = TRUE, Efficient = TRUE) {
    #extract the temporary variables
    if (AddConst)
        mRawX = cbind(matrix(1, nrow = nrow(mRawX)), mRawX)
    if (Efficient) {
        qrX = qr(mRawX);
        mX = qr.X(qrX);
        mR = qr.R(qrX);
        mQ = qr.Q(qrX);

        mQRTInv = mQ %*% t(backsolve(mR, diag(ncol(mR))));
        out = t(mQRTInv * rep((vResid / (1 - rowSums(mQ * mQ))) ^ 2, ncol(mX))) %*% mQRTInv;
    } else {
        #inefficient mode for verification. Only use for testing.
        mRawX = as.matrix(mRawX);
        out = solve(t(mRawX) %*% mRawX) %*% t(mRawX) %*%
        diag(as.vector((vResid / (1 - diag(mRawX %*% solve(t(mRawX) %*% mRawX) %*% t(mRawX)))) ^ 2)) %*%
             mRawX %*% solve(t(mRawX) %*% mRawX);
    }
    
    if (SEOnly)
        out = sqrt(diag(out));
    return(out);
}

IneffGetMWhiteErrors = function(mX, vResid) {
    #extract the temporary variables

    return(solve(t(mX) %*% mX) %*% t(mX) %*%
        diag(as.vector((vResid / (1 - diag(mX %*% solve(t(mX) %*% mX) %*% t(mX)))) ^ 2)) %*%
             mX %*% solve(t(mX) %*% mX));

}

NumTypes = function(x) length(unique(x));



#in: the contemporary art data
#out: a cleaned, transformed, and sorted data table that we can use to produce the summary statistics



###some stuff that is not used
#lot numbering (Not used???)
    #dfCRet[, "newlot"] = ifelse(dfCRet[, "date"] == c(as.Date(0), dfCRet[1:(iCurLength - 1), "date"]), 0, 1)
    #for (i in 1:iCurLength) {
    #    if (dfCRet[i, "newlot"] != 1) { dfCRet[i, "newlot"] = dfCRet[i - 1, "newlot"] + 1 }
    #    }