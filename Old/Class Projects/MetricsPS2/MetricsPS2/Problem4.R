#Problem 4
#tests the memory of three distributions

#Computes p(s<X<s+ds)|X>t
#In: the cumulative distribution function, s, ds
#Out: The corresponding X value
intervalGivenS = function(s, ds, cumDist) {
    return((cumDist(s + ds) - cumDist(s)) / (1 - cumDist(s)));
}

#Computes the maximum memory error
#In: the cumulative distribution function, the test minimum, test maximum, and interval width
#Out: The corresponding X value
testMemory = function(ds, sMin, sMax, cumDist) {
    pts = seq(sMin, sMax, ds);
    err = sapply(pts, intervalGivenS, ds, cumDist); #apply the interval test to each test point
    err = abs(err - cumDist(ds) + cumDist(0) )/ (cumDist(ds) - cumDist(0)) * 100; #generate the percent error
    return(max(err));
}


#####test with the exponential function
#set parameters
#lambda = 10;
#cat("The max percent difference between the interval at time t and at zero\n for the exponential distribution is:\n");
#print(testMemory(.01, 0, 2, function(x) pexp(x, lambda)));

#####check the normal
#set parameters
lambda = 10;
cat("The max percent difference between the interval at time t and at zero\n for the standard normal distribution is:\n");
print(testMemory(.01, -3, 3, pnorm));

#####now check the chi-squared
#set parameters
dof = 3;
cat("The max percent difference between the interval at time t and at zero\n for the chi-squared distribution is:\n");
print(testMemory(.01, 0, 5, function(x) pchisq(x, dof)));


