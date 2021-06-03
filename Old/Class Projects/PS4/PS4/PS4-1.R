#Clinton Tepper
#Problem 1
#Creates QQ plots for distributions of averages of ensembles of arbitrary distributions
#against rs normal distribution

rm(list = ls());

#Samples the average of an ensemble of numInstances arbitrary identical distributions
#In: function for sampling random variables (n times), 
#number of distributions, number of samples
#Out: a vector of samples, each of which is the average of numInstances distributions
sampleOfAvgs = function(rndFunc, numInstances, samples) {

    return(rep(1, samples) * sapply(rep(1, samples) * numInstances,
    function(n) sum(rndFunc(n))/numInstances));
}

########################Script entry point############################

#begin with distribution of 0.5
numSamples = 10 ^ 5;
p = .5;

rBern = function(n) rbinom(n, 1, p);
sigma = (p * (1 - p)) ^ .5

nInstances = 10;
qqplot(rnorm(numSamples)*sigma/nInstances^.5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=10",
xlab = "R norm Distribution",
ylab = "Sample Distribution");

nInstances = 100;
qqplot(rnorm(numSamples) * sigma / nInstances ^ .5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=100",
xlab = "R norm Distribution",
ylab = "Sample Distribution");

nInstances = 1000;
qqplot(rnorm(numSamples) * sigma / nInstances ^ .5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=1000",
xlab = "R norm Distribution",
ylab = "Sample Distribution");


#Now do distribution of 0.05
p = .05;
rBern = function(n) rbinom(n, 1, p);
sigma = (p*(1-p))^.5

nInstances = 10;
qqplot(rnorm(numSamples) * sigma / nInstances ^ .5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=10",
xlab = "R norm Distribution",
ylab = "Sample Distribution");

nInstances = 100;
qqplot(rnorm(numSamples) * sigma / nInstances ^ .5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=100",
xlab = "R norm Distribution",
ylab = "Sample Distribution");

nInstances = 1000;
qqplot(rnorm(numSamples) * sigma / nInstances ^ .5+p,
sampleOfAvgs(rBern, nInstances, numSamples),
main = "QQPlot for N=1000",
xlab = "R norm Distribution",
ylab = "Sample Distribution");










