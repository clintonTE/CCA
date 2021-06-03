


#The likelihood function for a multinomial distribution 
#In: A set of parameters and a data set
#Out: The sum of the log likelihoods for each data point
logitLogLike = function(V, cVec) {
    minLike = -10 ^ 8;
    denom = log(sum(exp(V)));
    logPVec = sapply(V, function(x) x - denom);
    #print(sum(sapply(cVec, function(x) logPVec[x])));
    return(max(sum(sapply(cVec, function(x) logPVec[x])),minLike));
}

#A wrapper function for logitLike amenable to the setup in problem 6
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, and the  data set
#Out: The sum of the log likelihoods for each data point
multinomChoiceModelLike = function(alpha, beta, lPrice, cVec) {
    #print(logitLogLike(rbind(beta + min(lPrice), alpha + beta * lPrice), cVec));
    return(logitLogLike(rbind(beta*min(lPrice), alpha + beta * lPrice),cVec));
}

#A gradient function specific to problem 6
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, price, and the  data set
#Out: The gradient of the likelihood
multinomChoiceModelGrad = function(alpha, beta, lPrice, cVec) {
    N = length(cVec);
    J = length(alpha) + 1;
    alpha = rbind(0, matrix(alpha, nrow = J - 1));
    lPrice = rbind(min(lPrice), matrix(lPrice, nrow = J - 1));
    eta = matrix(sapply(1:J, function(x) sum(cVec == x)), nrow = J);
    denom = sum(sapply(1:J, function(x) exp(alpha[x] + beta * lPrice[x])));
    return(rbind(matrix(sapply(2:J, function(x) eta[x] - N * exp(alpha[x] + beta * lPrice[x]) / denom), nrow = J-1),
        sum(sapply(1:J, function(x) lPrice[x] * eta[x] - N * lPrice[x] * exp(alpha[x] + beta * lPrice[x]) / denom))));
}


##########################Problem 5 Entry Point##########################

#Define the Parameters
numSamples=200;
alpha = matrix(NA, 3);
alpha[1] = 0.5;
alpha[2] = 1.5;
alpha[3] = 2;
beta = 2;
lowerBound = 0.1;
upperBound = 5;

numPrices = length(alpha);
logPrices = matrix(runif(numPrices, -1, 1),nrow=numPrices); #get the prices

#Simulate the data
popPVec = rbind(0, alpha) + rbind(min(logPrices), logPrices) * beta;
popPVec = exp(popPVec) / sum(exp(popPVec));
sampVec = sample.int(numPrices + 1, numSamples, replace = TRUE, popPVec);

print(multinomChoiceModelLike(alpha, beta, logPrices, sampVec));

#Get the maximum likelihood of the parameters
out=optim(c(rep(1, numPrices + 1)), 
    function(alphaBeta) - multinomChoiceModelLike(alphaBeta[1:numPrices], alphaBeta[numPrices + 1], logPrices, sampVec),
    function(alphaBeta) - multinomChoiceModelGrad(alphaBeta[1:numPrices], alphaBeta[numPrices + 1], logPrices, sampVec),
    method = "L-BFGS-B",
    lower = c(rep(lowerBound, numPrices + 1)),
    upper = c(rep(upperBound, numPrices + 1)));

print(out$par);
#multinomChoiceModelLike(out$par[1:numPrices], out$par[numPrices + 1], logPrices, sampVec)







