set.seed(11);
rm(list = ls());

source("M2-PS2-Problem1.R");
#Problem1Script();

#########################Script Entry Point
iN = 200;
#number of points
iNumSim = 2500;

mX = matrix(runif(iN), nrow = iN, ncol = 1);
mPopResid = rbind(matrix(rnorm(iN / 2 * iNumSim, mean = 0, sd = 0.5), nrow = iN / 2),
    matrix(rnorm(iN / 2 * iNumSim, mean = 0, sd = 2), nrow = iN / 2));
mY = matrix(rep(mX, iNumSim), nrow = iN, ncol = iNumSim) * 2 + mPopResid;
m1X = cbind(rep(1, iN), mX);

#first do the OLS calcualtions
mPreMult = solve(t(m1X) %*% m1X);
mBetaOLS = matrix(apply(mY, 2, function(v) mPreMult %*% t(m1X) %*% mY), nrow = 2, ncol = iNumSim);
mResid = matrix(rep(mX, iNumSim), nrow = iN, ncol = iNumSim) * rep(t(mBetaOLS[2,]), iN) - mY;

print("We calculate our beta coefficient using GLS:");

#estimate the standard deviations for the model among the groups
dVar1 = matrix(apply(mResid, 2, function(v) var(v[1:iN / 2])), nrow = 1, ncol = iNumSim);
dVar2 = matrix(apply(mResid, 2, function(v) var(v[(iN / 2):iN])), nrow = 1, ncol = iNumSim);

#generate the weights as a vector
mWeights = rbind(1 / matrix(rep(dVar1, iN / 2), nrow = iN / 2, ncol = iNumSim, byrow = TRUE),
    1 / matrix(rep(dVar2, iN / 2), nrow = iN / 2, ncol = iNumSim, byrow = TRUE));

#now print the GLS simulation stats
mBetaGLS = matrix(sapply(1:iNumSim, function(x) solve(t(m1X * mWeights[, x])
    %*% m1X) %*% t(m1X * mWeights[, x]) %*% mY[, x]), nrow = 2, ncol = iNumSim);
print("The average GLS beta is:")
print(rowMeans(mBetaGLS));

print("The simulation standard error of the GLS beta estimates is:")
print(sd(mBetaGLS[2,]));
print("The simulation variance of the GLS beta estimates is:")
print(var(mBetaGLS[2,]));

#now get the standard error for the GLS regressions
mBetaGLSSE = matrix(sapply(1:iNumSim, function(x) diag(solve(t(m1X * mWeights[, x])
    %*% m1X))), nrow = 2, ncol = iNumSim);
print("The RMS (sqrt of mean variance) of the GLS SE of the beta coefficient is:")
print(sqrt(mean(mBetaGLSSE[2,])));
print("The mean of the calculated GLS variance is:")
print(mean(mBetaGLSSE[2,]));

#now provide the OLS results
print("The average OLS beta is:")
print(rowMeans(mBetaOLS));
print("The simulation standard error of the OLS beta estimates is:")
print(sd(mBetaOLS[2,]));
print("The simulation variance of the OLS beta estimates is:")
print(var(mBetaOLS[2,]));

print("The RMS of the homoskedastic calculated standard error of OLS estimates is:")
mBetaOLSSE = matrix(apply(mResid, 2, var), nrow = 1, ncol = iNumSim) * diag(mPreMult)[2];
print(sqrt(mean(mBetaOLSSE)));
print("The mean of the calculated homoskedastic variances is:")
print(mean(mBetaOLSSE));


print("The true OLS standard error of beta is: ");
trueOLSVar = mPreMult %*% t(m1X) %*% diag(c(rep(.5 ^ 2, iN / 2), rep(2 ^ 2, iN / 2))) %*% m1X %*% t(mPreMult);
print(sqrt(trueOLSVar[2,2]));
print("The true OLS variance is: ")
print(trueOLSVar[2, 2]);

#Now provide results for White Standard Errors
qrX = qr(m1X);
mBetaWhiteSE = matrix(apply(mResid, 2, function(v) diag(GetMWhiteErrors(qrX=qrX, v))),ncol=iNumSim);
print("The RMS  of the calculated Modified White SE of OLS estimates is:")
print(sqrt(mean(mBetaWhiteSE[2,])));
print("The mean  of the calculated Modified White SE of OLS estimates is:")
print(mean(mBetaWhiteSE[2,]));

#Graphical Summary of results
ggplot(stack(data.frame(as.vector(mBetaGLS[2,]), as.vector(mBetaOLS[2,]))), aes(x = values)) +
    geom_density(aes(group = ind, color = ind)) + theme_bw() + ggtitle("Distribution of GLS and OLS Betas");

ggplot(stack(data.frame(as.vector(sqrt(mBetaGLSSE[2,])), as.vector(sqrt(mBetaOLSSE)), as.vector(sqrt(mBetaWhiteSE[2,])))), aes(x = values)) +
    geom_density(aes(group = ind, color = ind)) + theme_bw() +
    ggtitle("Distribution of GLS, OLS Homoskedastic, and OLS Modified White Standard Errors");







