#Clinton Tepper
#Problem 3-c
# This program creates QQ plots for the transformation X1_b and X2_b where X1_b and X2_b 
#are the sample averages of a bivariate normal. The draws are plotted against the asymptotic approximation.


rm(list = ls());

#Returns a matrix of draws from the multivariate normal distribution
#In: numSamples, mu (vector), Ut (transpose of cholesky decomp)
#Out: a matrix of numSamples vectors of multivariate normal draws
bivNorm = function(numSamples, mu, Ut) {
    l = length(mu)
    return(matrix(sapply(1:numSamples, function(x) mu + Ut %*% rnorm(l)),
         numSamples,
         l,
         byrow = TRUE));
}

#Returns an R^2 -> R transformation of the sample average of bivariate normals
#In: transform function g, number of draws, sequence number, mu (vector), cov (matrix)
#Out: a vector of samples of the transformation
getBivTransDraws = function(g, numSamples, n, mu, cov) {
tCholMat = t(chol(cov));
bivDraws = matrix(sapply(rep(n, numSamples), function(x) colSums(bivNorm(x, muVec, tCholMat))),
    numDraws, length(muVec), byrow = TRUE);
return (g(bivDraws[, 1], bivDraws[, 2]));
}

########################Script entry point############################
numDraws = 10 ^ 4;
covMat = matrix(c(1, 1 / 2, 1 / 2, 1 / 3), 2, 2);
# Note- Hilbert matrices are positive definite. 
muVec = c(10, 1/3);
tranVec = matrix(c(1 / muVec[2], - muVec[1] / (muVec[2] ^ 2)), 2, 1);
sigmaAsympt = sqrt(t(tranVec) %*% covMat %*% tranVec);

seqN = 10;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function (x,y) x/y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=10, sigma22=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");

seqN = 100;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function(x, y) x / y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=100, sigma22=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");

seqN = 1000;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function(x, y) x / y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=1000, sigma22=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");


#####Now try different mu_2 and sigma_22
muVec = c(10, 10);
tranVec = matrix(c(1 / muVec[2], - muVec[1] / (muVec[2] ^ 2)), 2, 1);
sigmaAsympt = sqrt(t(tranVec) %*% covMat %*% tranVec);

seqN = 10;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function(x, y) x / y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=10, sigma22*30=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");

seqN = 100;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function(x, y) x / y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=100, sigma22*30=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");

seqN = 1000;
asymptDraws = rnorm(numDraws) * sigmaAsympt / sqrt(seqN) + muVec[1] / muVec[2];
actDraws = getBivTransDraws(function(x, y) x / y, numDraws, seqN, muVec, covMat);
qqplot(actDraws,
asymptDraws,
main = "QQPlot for N=1000, sigma22*30=mu",
xlab = "Actual X1_b/X2_b Distribution",
ylab = "Sample Distribution");

