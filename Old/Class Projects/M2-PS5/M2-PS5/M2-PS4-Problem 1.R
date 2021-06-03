rm(list = ls());

require(devtools);
require(DataAnalytics);
require(plyr);
require(ggplot2);

set.seed(11);
#allows for reproducability

#A function to get the truncated expectation given a distribution and break points 
#In: The distribution and break points
#Out: The expectation
TruncMean = function(Dist, dLo, dHi) {
    if (missing(dLo)) dLo = -Inf;
    if (missing(dHi)) dHi = Inf;
    return(integrate(function(x) x * Dist(x), lower = dLo, upper = dHi)$val);
}

#A function to get the weights using the weighted derivitive interpretation of a regression
#In: The distribution, cumulative distribution, and point for the weight
#Out: The expectation
regWeight = function(Dist, CumDist, dT, dLowerBound, dUpperBound) {
    if (missing(dLowerBound)) dLowerBound = -Inf;
    if (missing(dUpperBound)) dUpperBound = Inf;
    dProb = 1-CumDist(dT);
    return((TruncMean(Dist, dLo = dT, dHi = dUpperBound)/dProb -
        TruncMean(Dist, dLo = dLowerBound, dHi = dT)/(1-dProb)) * dProb * (1 - dProb));
}

#########################################Script Entry Point###########################
Dist = function(x) dchisq(x, 4);
CumDist = function(x) pchisq(x, 4);
dLo = 0.1;
dHi = 25;
dInterval = .1;

dXVal = seq(from = dLo, to = dHi, by = dInterval);
dYVal = sapply(dXVal, function(x) regWeight(Dist = Dist, CumDist = CumDist, dT = x, dLowerBound = 0));

ggplot(data.frame(dXVal, dYVal), aes(x = dXVal, y = dYVal)) +
    geom_point() + theme_bw() +
    ggtitle("Chi-Sq (4 dof) Weight vs X Values");




