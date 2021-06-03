#Clinton Tepper
#Problem 3
#creates a qq plot for the exponential distribution using a determiend lambda
#inverse CDF: X=ln(1-F)/-Lambda


#Computes the inverse of the CDF
#In: Cumualtive probability
#Out: The corresponding X value
expInverse = function(pCum,lambda) {
    return(log(1 - pCum) / ( - lambda));
}



################################################R Script- Entry Point #############################
#performs a qqplot of expInverse against R's exponential distribution
lambdaIn = 1;
qqplot(rexp(10 ^ 5, rate = lambdaIn), 
    expInverse(runif(10 ^ 5), lambdaIn), 
    main = "QQPlot", 
    xlab = "R Exp Distribution", 
    ylab = "Sample Distribution");

#Conclude:
#1) The "custom" inverse function is likely the same distrbution as the R distribution
#2) The test loses accuracy at distribution values greater then ~8.5 


