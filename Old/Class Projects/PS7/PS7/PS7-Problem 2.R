set.seed(11); #allows for reproducability

#The likelihood function for a multinomial distribution 
#In: A set of parameters and a data set
#Out: The sum of the binary likelihoods for each data point
binLogLike = function(beta, xMat, yVec) {
    expon = xMat %*% beta;
    #calculate the exponent

    return(t(yVec) %*% expon - sum(log(1 + exp(expon))));
}

#A gradient function specific to problem 7
#In: beta coeficients, matrix of xvalues, y values. The first beta value is assumed to be the intercept
#Out: The gradient of the likelihood
binLogGrad = function(beta, xMat, yVec) {
    DoF = length(beta);
    expon = xMat %*% beta;

    return(rbind(sum(yVec) - sum(exp(expon) / (1 + exp(expon))),
        matrix(sapply(2:DoF,
            function(k) t(xMat[, k]) %*% (yVec - (exp(expon) / (1 + exp(expon))))),
        nrow = dof - 1)));
}

##########################Problem 2 Entry Point##########################

#Define the Parameters
numSamples = 10^5;
beta = matrix(seq(0.5,2.5,.5), nrow = 5);
#beta = matrix(c(.38, 0, 0, 0, 0), nrow = 5);
initBeta = matrix(rep(rep(1, 5), nrow = 5));
dof = length(beta);
lowerBound = -1;
upperBound = 5;

#Simulate the x values
xMat = cbind(rep(1, numSamples),
    matrix(runif((dof - 1) * numSamples, -1, 1), nrow = numSamples, ncol = dof - 1));

#simulate the y values
yVec = runif(numSamples);
expon = xMat %*% beta;
yVec = (yVec<exp(expon) / (1 + exp(expon)));

#Get the maximum likelihood of the parameters
likeFunc = function(b) - binLogLike(b, xMat, yVec);
gradFunc = function(b) - binLogGrad(b, xMat, yVec);

#get the mle data
out = optim(initBeta, likeFunc, gradFunc,
    method = "L-BFGS-B",
    lower = c(rep(lowerBound, dof)),
    upper = c(rep(upperBound, dof)), 
    hessian = TRUE);

print(out);

mlePar = matrix(out$par, nrow = dof);
mleVal = out$value;
mleHess = out$hessian[1:dof, 1:dof];
nullH = rbind(sum(yVec) / (sum(yVec) + numSamples), matrix(rep(0, dof - 1), nrow = dof - 1));
lambda = 2*(-mleVal + likeFunc(nullH));
print("Likelihood Ratio Test");
print(1-pchisq(lambda, dof-1));

rVec = mlePar[2:dof] - nullH[2:dof];
W = t(rVec) %*% mleHess[2:dof,2:dof] %*% rVec;
print("Wald Test");
print(1-pchisq(W, dof-1));


