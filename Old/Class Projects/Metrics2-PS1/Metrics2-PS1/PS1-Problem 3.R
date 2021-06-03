require(devtools);
require(DataAnalytics);

################################helper functions###################
#in: a regression object
#out: a matrix object for the coefficients
getCoefAsMatrix = function(reg) matrix(coef(reg),
    nrow = length(coef(reg)),
    dimnames = list(names(coef(reg))));

#this function takes in a data frame and outputs the beta coeficients for all pairwise combinations
#in: a data frame
#out: a matrix of beta coefficients for each pair of columns
getBetaMatrix = function(dat) {
    numCombo = ncol(dat);
    indices = expand.grid(1:numCombo, 1:numCombo);

    #the helper function below runs the pairwise regression
    return(matrix(mapply(function(x,y)
        getCoefAsMatrix(lm(dat[, y] ~ dat[, x], data = dat))[2, 1],
        indices[, 1], indices[, 2]), nrow = numCombo));
}

####################################Script Entry Point##########################
load("pl_share_demo_2005.rda");

PLData = pl_share_demo_2005[, c("pl_share", "lnIncome", "lnAge", "num_children")];
#get just the data we are analyzing

#first run the long regression
longPL = lm(pl_share ~ lnIncome + lnAge + num_children, data=PLData);

#extract the coefficients into a matrix
longCoef = getCoefAsMatrix(longPL);
print("Coeficients for part A")
print(longCoef);

betaMat = getBetaMatrix(PLData);
gMat = betaMat[2:4, 2:4]; # we get the matrix of beta coeficients where we regress column on row
bVecShort = betaMat[2:4, 1];
#this is a vector of our short regression coeficients

longCoefMethod2 = solve(gMat) %*% bVecShort;
print("Coeficients for part B");
print(longCoefMethod2);
print("Max Difference Between Methods:");
print(abs(max(longCoef[2:4] - longCoefMethod2)));

E13 = resid(lm(lnIncome ~ num_children, data = PLData))
E23 = resid(lm(lnAge ~ num_children, data = PLData))
E13d23 = resid(lm(E13 ~ E23))
print(lm(pl_share ~ E13d23, data = PLData))









