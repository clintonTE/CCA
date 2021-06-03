#Clinton Tepper
#Problem 2
#Computes nth order tailor series expansion for e^x

#method: ePowerSeries
#Inputs: x, the number of series terms
#Out: approximation of e^x at 0 for n terms
ePowerSeries = function(n,x) {
    return (sum(sapply(0:n,function (pow) x^pow/factorial(pow))))
}

################################################R Script- Entry Point #############################
#script to get the error value at 10 for the expansion around zero
approxTerms = 1:100;
xToTest = 10;
mAns = matrix(sapply(approxTerms, ePowerSeries, xToTest), nrow = length(approxTerms)); #get the approximations
mAns = cbind(approxTerms, exp(1)^xToTest, mAns, abs(exp(1)^xToTest - mAns));
colnames(mAns) = c("nTerms", "act", "approx", "abs % error");
cat("Approximation is:\n")
print(mAns);


