#The likelihood function for a multinomial distribution 
#In: A set of parameters and a data set
#Out: The sum of the log likelihoods for each data point
logitLogLike = function(V, cMat) {
    minLike = -10 ^ 8;
    denom = matrix(apply(V, 2, function(x) log(sum(exp(x)))),nrow=1);
    #calculate the denominator for each price
    #calculate the log probability for each value of v
    logPVec = t(apply(V, 1, function(x) x - denom));

    #The below command looks up the logged probability for each entry in cMat and sums the total
    return(max(sum(sapply(1:ncol(V), function(x) sum(logPVec[cMat[, x], x])), minLike)));
}

#A wrapper function for logitLike amenable to the setup in problem 6
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, and the  data set
#Out: The sum of the log likelihoods for each data point
multinomChoiceModelLike = function(alpha, beta, lPrice, cMat) {

    #form the matrix of Vs
    return(logitLogLike(apply(lPrice, 2, function(x) alpha + beta * x), cMat));
}

#A gradient function specific to problem 6
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
        function(x) sum(eta[x,] - L * exp(alpha[x] + beta * lPrice[x,]) / denom)),ncol=1),
        sum(sapply(1:J,
            function(x) sum(lPrice[x,] * (eta[x,] - L * exp(alpha[x] + beta * lPrice[x,]) / denom))))));
}

##########################Problem 5 Entry Point##########################

#Define the Parameters
numSamples = 100;
numPriceVectrs = 200;
alpha = matrix(NA, 3);
alpha[1] = 0;
alpha[2] = 0.5;
alpha[3] = 1.5;
beta = 2;
lowerBound = 0.1;
upperBound = 5;

numPrices = length(alpha);
#get the prices
logPrices = rbind(matrix(rep(0, numPriceVectrs), nrow = 1),
    matrix(runif(numPriceVectrs * (numPrices - 1), -1, 1), nrow = numPrices - 1));

#Simulate the data
popPVec = exp(rep(alpha, numPriceVectrs) + logPrices * beta);
popPVec = apply(popPVec, 2, function(x) x/sum(x));
sampMat = matrix(sapply(1:numPriceVectrs,
    function(x) sample.int(numPrices, numSamples, replace = TRUE, popPVec[, x])),
    nrow = numSamples, ncol = numPriceVectrs);

#Get the maximum likelihood of the parameters
likeFunc = function(alphaBeta)
    - multinomChoiceModelLike(rbind(0, matrix(alphaBeta[1:(numPrices - 1)], nrow = numPrices-1)),
        alphaBeta[numPrices], logPrices, sampMat);

gradFunc = function(alphaBeta)
    - multinomChoiceModelGrad(rbind(0, matrix(alphaBeta[1:(numPrices - 1)], nrow = numPrices-1)),
        alphaBeta[numPrices], logPrices, sampMat);

out = optim(c(rep(1, numPrices)), likeFunc, gradFunc,
    method = "L-BFGS-B",
    lower = c(rep(lowerBound, numPrices)),
    upper = c(rep(upperBound, numPrices)));

print(out);
#multinomChoiceModelLike(out$par[1:numPrices], out$par[numPrices + 1], logPrices, sampVec)
