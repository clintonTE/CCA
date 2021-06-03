#require(grDevices)

#The likelihood function for a two column data set where data is censored below C
#In: C, sigma, beta, xVec, yVec
#Out: likelihood
tobitLike = function(beta, sigma, cutoff, xVec, yVec) {
    minLike = -10^8; #minimum likelihood to prevent -inf results
    singleLike = function(x, y) 
        log(if (y > cutoff) dnorm(y-beta*x, 0, sigma) else pnorm(cutoff-beta*x, 0, sigma));
        #print(mapply(singleLike, xVec, yVec));
    return(max(sum(mapply(singleLike, xVec, yVec)),minLike));
}

########################Problem 5 Script entry point############################
numVals = 4000; #Define the parameters
edgeSize = 50;
sigmaMin = 0.1;
sigmaMax = 1;
betaMin = -3;
betaMax = 3;

# the true population paramaters
betaPop = 1; 
sigmaPop = 0.6;

xSeq = seq(betaMin, betaMax, (betaMax - betaMin) / edgeSize); #establish the parameter grid
ySeq = seq(sigmaMin, sigmaMax, (sigmaMax - sigmaMin) / edgeSize);
sigBetaGrid = data.matrix(expand.grid(xSeq, ySeq)); #establish the parameter grid

xVector = runif(numVals, -1, 1);
yVector = mapply(function(x,e) max(betaPop*x + sigmaPop*e,0),xVector,rnorm(numVals));

likelihoods = mapply(function(b, s) tobitLike(b, s, 0, xVector, yVector), #pass the appropriate parameters to tobitlike
    sigBetaGrid[, 1], sigBetaGrid[, 2]);
likeLevels = quantile(likelihoods, seq(0, 1, .04));
#print(likelihoods);
contour(xSeq, ySeq, matrix(likelihoods, length(xSeq), length(ySeq)), levels = likeLevels,
main = "Problem 5 Contour Map", xlab="beta", ylab="sigma");





