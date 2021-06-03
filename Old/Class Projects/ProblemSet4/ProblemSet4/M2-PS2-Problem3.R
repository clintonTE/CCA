rm(list = ls());
set.seed(11);

source("M2-PS2-Problem1.R");



################################helper functions###################
#in: a regression object
#out: a matrix object for the coefficients
GetCoefAsMatrix = function(reg) matrix(coef(reg),
    nrow = length(coef(reg)),
    dimnames = list(names(coef(reg))));


####################################Script Entry Point##########################
load("pl_share_demo_2005.rda");

#we drop irrelevant columns and rows with missing data
dfPLData = pl_share_demo_2005[complete.cases(pl_share_demo_2005[, c("pl_share", "lnIncome", "DMA2005")]), 
    c("pl_share", "lnIncome", "DMA2005")];

iN = nrow(dfPLData);

#generate the x matrix

m1X = cbind(rep(1, iN), dfPLData$lnIncome);

regOLS = lm(dfPLData$pl_share ~ dfPLData$lnIncome);
mOLSBeta = GetCoefAsMatrix(regOLS);
vResid = dfPLData$pl_share - m1X %*% mOLSBeta;
print("Results of OLS regression:")
print(regOLS);
print("OLS homoskedastic standard errors:");
print(sqrt(diag(as.numeric(var(vResid)) * solve(t(m1X) %*% m1X))));
print("OLS White corrected standard errors:")
print(sqrt(diag(GetMWhiteErrors(qr(m1X), vResid = vResid))));

#estimate the variances and apply the weights
mDMAVars = aggregate(x = vResid, by = list(dfPLData$DMA2005), FUN = var);
mDMAVars = rename(cbind(mDMAVars, 1 / mDMAVars$V1), c("1/mDMAVars$V1" = "weights"));

dfPLData = merge(dfPLData, mDMAVars, by.x = "DMA2005", by.y = "Group.1", all.x = TRUE);
m1X = cbind(rep(1, iN), dfPLData$lnIncome); #reassign due to different sorting order

mGLSBeta = matrix(solve(t(m1X * dfPLData$weights) %*% m1X) %*% t(m1X * dfPLData$weights) %*%
    dfPLData$pl_share, nrow = 2, ncol = 1);

print("The GLS intercept and beta coeficients are respectively:");
print(mGLSBeta);
print("The standard error of the two GLS coefficients is:")
print(sqrt(diag(solve(t(m1X * rep(dfPLData$weights, 2)) %*% m1X))));



