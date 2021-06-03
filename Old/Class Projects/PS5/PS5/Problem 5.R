require(grDevices)

#The likelihood function for a two column data set where data is censored below C
#In: C, sigma, beta, xVec, yVec
#Out: likelihood
tobitLike = function(sigma, beta, cutoff, xVec, yVec) {
    singleLike = function(x, y) 
        if (y > cutoff) pnorm((y - beta * x) / sigma)
        else dnorm((cutoff - beta * x) / sigma);
    return(prod(mapply(singleLike, xVec, yVec)));
}

########################Problem 5 Script entry point############################
numVals = 10 ^ 3; #Define the parameters
edgeSize = 10 ^ 1;
sigmaMin = 0.01;
sigmaMax = .6;
betaMin = -3;
betaMax = 3;

#sigBetaGrid = matrix(expand.grid(seq(sigmaMin, sigmaMax, (sigmaMax - sigmaMin) / edgeSize),
#    seq(betaMin, betaMax, (betaMax - betaMin) / edgeSize)),nrow=numVals, ncol=2); #establish the parameter grid
xSeq = seq(sigmaMin, sigmaMax, (sigmaMax - sigmaMin) / edgeSize);
ySeq = seq(betaMin, betaMax, (betaMax - betaMin) / edgeSize);
sigBetaGrid = expand.grid(xSeq, ySeq); #establish the parameter grid
#print(sigBetaGrid);

xVector = runif(numVals, -1, 1);
err = rnorm(numVals);
likelihoods = mapply(function(s, b) tobitLike(s, b, 0, xVector, max(b * xVector + s * err, 0)),
    sigBetaGrid[, 1], sigBetaGrid[, 2]);
likeLevels = quantile(likelihoods, seq(0, 1, .1));
#print(matrix(likelihoods, length(xSeq), length(ySeq)));
#print(likeLevels);
filled.contour(xSeq, ySeq, matrix(likelihoods, length(xSeq), length(ySeq)), levels = likeLevels, color.palette = heat.colors,lty = "solid");




