#The likelihood function for a multinomial distribution 
#In: A set of parameters and a data set
#Out: The sum of the log likelihoods for each data point
logitLogLike = function(V, cVec) {
    minLike = -10 ^ 8;
    denom = log(sum(exp(V)));
    logPVec = sapply(V, function(x) x - denom);
    #print(sum(sapply(cVec, function(x) logPVec[x])));
    return(max(sum(sapply(cVec, function(x) logPVec[x])), minLike));
}

#A wrapper function for logitLike amenable to the setup in problem 6
# Note that a null alternative is added. The null alternative is assumed to be J=1.
#In: Alpha, Beta, and the  data set
#Out: The sum of the log likelihoods for each data point
multinomChoiceModelLike = function(alpha, lPrice, cVec) {
    #print(logitLogLike(rbind(beta + min(lPrice), alpha + beta * lPrice), cVec));
    #print(logitLogLike(alpha + lPrice, cVec));
    return(logitLogLike(alpha + lPrice, cVec));
}


##########################Problem 5 Entry Point##########################

#Define the Parameters
numSamples = 20000;
alpha = matrix(NA, 3);
alpha[1] = 0.25;
alpha[2] = 0.4;
alpha[3] = 0.6;
lowerBound = 0.1;
upperBound = 3;

numPrices = length(alpha);
logPrices = matrix(runif(numPrices, -1, 1), nrow = numPrices);
#get the prices

#Simulate the data
popPVec = exp(alpha + logPrices);
popPVec = popPVec / sum(popPVec);
sampVec = sample.int(numPrices, numSamples, replace = TRUE, popPVec);


#Get the maximum likelihood of the parameters
out = optim(c(rep(0.5, numPrices)),
    function(param) - multinomChoiceModelLike(param, logPrices, sampVec),
    method = "L-BFGS-B", lower = c(rep(lowerBound, numPrices)),
    upper = c(rep(upperBound, numPrices)));

print(out);
