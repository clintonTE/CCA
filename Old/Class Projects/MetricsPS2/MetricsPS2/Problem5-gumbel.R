#Clinton Tepper
#Problem 5-Gumbel
#Includes a method for taking the maximum of multiple CDFs and calculating the pdf

require(ggplot2);
require(reshape2);

#Calculates the density function for the max of several sub distributions
#In: x, f, F, n
#Out: the probabilty distribution of the maximum of the n input variables
getMaxProbDist=function(x,n,f, fCum) {
    if (n == 1) {
        return(f(x));
    }

    return(n * f(x) * fCum(x) ^ (n - 1));
}

##############################Script entry point##################


betaParam = 2;
xVal = seq(-10, 10, .001);
f_in = function(x) 1/betaParam*exp(-x/betaParam)*exp(-exp(-x/betaParam));
fCumIn = function(x) exp(-exp( - x / betaParam));
dat = xVal;

for (j in 1:5) {
    dat = cbind(dat, sapply(xVal, getMaxProbDist, 2 ^ j, f_in, fCumIn));
}
colnames(dat) = c("x", "n=2", "n=4", "n=8", "n=16", "n=32");

datFrame = data.frame(dat);
longList = melt(datFrame, id = "x");
ggplot(longList, aes(x, value, color = variable)) + geom_line();

#Conclude that increasing n shifts the distribution to the right 
# without signficantly affecting the shape















