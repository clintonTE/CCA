set.seed(11); #allows for reproducibility

#The likelihood function for a multinomial distribution 
#In: A set of parameters and a data set
#Out: The sum of the log likelihoods for each data point
logitLogLike = function(V, cMat) {
    minLike = -10 ^ 9;
    maxLike = 10 ^ 9;
    denom = matrix(apply(V, 2, function(x) log(sum(exp(x)))), nrow = 1);
    #calculate the denominator for each price
    #calculate the log probability for each value of v
    logPVec = t(apply(V, 1, function(x) x - denom));

    #The below command looks up the logged probability for each entry in cMat and sums the total
    return(sum(sapply(1:ncol(V), function(x) sum(logPVec[cMat[, x], x]))));
}

#A wrapper function for logitLike amenable to the setup in problem 6
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, and the  data set
#Out: The sum of the log likelihoods for each data point
multinomChoiceModelLike = function(alpha, beta, lPrice, cMat) {

    #form the matrix of Vs
    return(logitLogLike(apply(lPrice, 2, function(x) alpha + beta * x), cMat));
}

#A gradient function specific to problem 1
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, price, and the  data set. Alpha[0] and price[0] are assumed to be zero.
#Out: The gradient of the likelihood
multinomChoiceModelGrad = function(alpha, beta, lPrice, cMat) {
    L = nrow(cMat);
    W = ncol(cMat);
    J = length(alpha);
    eta = matrix(sapply(1:J,
        function(x) apply(cMat, 2, function(y) sum(y == x))), nrow = J, ncol = W, byrow = TRUE);
        #calculate the denominator for each price series    
    denom = apply(lPrice, 2, function(x) sum(exp(alpha + x * beta)));
    #the below return calculates and sums the gradients accross prices

    return(rbind(matrix(sapply(2:J,
        function(x) sum(eta[x,] - L * exp(alpha[x] + beta * lPrice[x,]) / denom)), ncol = 1),
        sum(sapply(1:J,
            function(x) sum(lPrice[x,] * (eta[x,] - L * exp(alpha[x] + beta * lPrice[x,]) / denom))))));
}

##########################Problem 1 Entry Point##########################

#Define the Parameters
numSamples = 50;
numPriceVectrs = 20;
alpha = matrix(NA, 2);
alpha[1] = 0;
alpha[2] = 0;
beta = 500;
lowerBound = 0.1;
upperBound = 5000;
maxIter = 100;

numPrices = length(alpha);
#get the prices
logPrices = rbind(matrix(rep(0, numPriceVectrs), nrow = 1),
    matrix(runif(numPriceVectrs * (numPrices - 1), -1, 1), nrow = numPrices - 1));

#Simulate the data
unifPts = 0
iter=0
while ((unifPts != numSamples * numPriceVectrs) && (iter<maxIter)) {
    popPVec = exp(rep(alpha, numPriceVectrs) + logPrices * beta);
    popPVec = apply(popPVec, 2, function(x) x / sum(x));
    sampMat = matrix(sapply(1:numPriceVectrs,
        function(x) sample.int(numPrices, numSamples, replace = TRUE, popPVec[, x])),
        nrow = numSamples, ncol = numPriceVectrs);
    unifPts = sum(sapply(1:numPriceVectrs, function(x) sum(sampMat[1, x] == sampMat[, x])));
    iter = iter + 1;
}
#print(sampMat);
if (iter == maxIter) print("Warning: Conformed data set not generated");
#Get the maximum likelihood of the parameters
likeFunc = function(b)
    - multinomChoiceModelLike(alpha, b, logPrices, sampMat);

gradFunc = function(b)
    -multinomChoiceModelGrad(alpha, b, logPrices, sampMat)[numPrices];

gMin = 100;
gMax = 1000;
xVec = gMin:gMax;
yVec = sapply(gMin:gMax, likeFunc)


plot(xVec, yVec, main = "Problem 1 Plot", xlab = "beta", ylab = "likelyhood");

out = optim(700, likeFunc, gradFunc,
    method = "L-BFGS-B",
    lower = c(rep(lowerBound, numPrices)),
    upper = c(rep(upperBound, numPrices)));

print(out);

#multinomChoiceModelLike(out$par[1:numPrices], out$par[numPrices + 1], logPrices, sampVec)
