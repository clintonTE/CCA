
#The likelihood function for a two column data set where data is censored below C
#In: C, sigma, beta, xVec, yVec
#Out: likelihood
tobitLogLike = function(beta, sigma, cutoff, xVec, yVec) {
    minLike = -10 ^ 8;
    #minimum likelihood to prevent -inf results
    singleLike = function(x, y)
        log(if (y > cutoff) dnorm(y - beta * x, 0, sigma) else pnorm(cutoff - beta * x, 0, sigma));
    return(max(sum(mapply(singleLike, xVec, yVec)), minLike));
}

#The gradiant of the tobit likelihood function for a two column data set where data is censored below C
#In: C, sigma, beta, xVec, yVec
#Out: gradient of the tobit likelihood
tobitGradLogLike = function(beta, sigma, cutoff, xVec, yVec) {
    minLike = -10 ^ 8;
    #minimum likelihood to prevent -inf results
    singleGrad = function(x, y) matrix(
            if (y > cutoff) c(x * (y - beta * x) / sigma ^ 2, (y - beta * x) ^ 2 / sigma ^ 3 - 1 / sigma)
            else c(0, - x / sigma ^ 2 * dnorm( - x / sigma) / pnorm( - x / sigma)), ncol = 2);
            #calculates one value for the likelihood

    return(apply(mapply(singleGrad, xVec, yVec),1, function(x) max(sum(x), minLike)));
}


#Simulates y|x provided that y|x is normally distributed about beta*x
#In: alpha, beta, sigma, a vector of x values, (optional) a cutoff for censored values
#Out: a vector of simulated y values
simLinYData = function(alpha, beta, sigma, xVec, cutOff) {
    if (missing(cutOff)) 
        return(mapply(function(x, e) alpha + beta * x + sigma * e, xVec, rnorm(length(xVec))));
    return(mapply(function(x, e) max(alpha + beta * x + sigma * e,cutOff), xVec, rnorm(length(xVec))));
}

##########################Problem 6 Entry Point##########################
#define the simulation parameters
numVals = 1000;
numSampleSets = 50;

#boundaries for optimization
sigmaMin = 0.1;
sigmaMax = 2;
betaMin = -3;
betaMax = 3;

# the true population paramaters
betaPop = 1;
sigmaPop = 0.6;
alphaPop = 0;

#create a matrix of x values
xVector = runif(numVals, -1, 1); 

#simulate the y values for each sample set
yVectors = apply(matrix(rep(xVector, numSampleSets), nrow = numVals, ncol = numSampleSets),
2, function(x) simLinYData(alphaPop, betaPop, sigmaPop, x, 0));

#simulate the maximum likelihood for each data set

optimFunc = function(y) optim(c(.5, .5), function(betaSigma) -tobitLogLike(betaSigma[1], betaSigma[2], 0,xVector, y),
    function(betaSigma) - tobitGradLogLike(betaSigma[1], betaSigma[2], 0, xVector, y), 
    method = "L-BFGS-B", lower = c(betaMin, sigmaMin), upper = c(betaMax, sigmaMax));

mleResults = apply(yVectors, 2, optimFunc);

mleParResults = matrix(NA, numSampleSets, 2)
for (j in 1:numSampleSets) {
    mleParResults[j,] = mleResults[[j]]$par;
}

#boxplots
boxplot(mleParResults[, 1], main = "Problem 6 Box Plot of Beta", ylab = "beta");
boxplot(mleParResults[, 2], main = "Problem 6 Box Plot of Sigma", ylab = "sigma");






