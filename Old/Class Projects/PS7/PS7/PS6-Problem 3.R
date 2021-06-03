#Creates a non-parametric bootstrap distribution from a population distribution given an
#arbitrary sample statistic (1 arg)
#In: population data, a stat function (1 arg), number of samples per draw, number of draws
#Out: a vector of bootstraps
boot.NPBootStrapper = function(datSet, statFunc, samplesPerDraw, numDraws) {

    return(matrix(sapply(1:numDraws,
    function(x) statFunc(sample(datSet, samplesPerDraw, replace = TRUE))),
    numDraws));
}


#Creates a parametric bootstrap distribution from a population distribution given an
#arbitrary sample statistic (1 arg)
#In: a function with parameters estimated (1arg, a stat function (1 arg), 
# number of samples per draw, number of draws
#Out: a vector of bootstraps
boot.PBootStrapper = function(estFunc, statFunc, samplesPerDraw, numDraws) {

    return(matrix(sapply(1:numDraws, function(x) statFunc(estFunc(samplesPerDraw))), numDraws));
}

#Creates a sample distribution from a population function given an
#arbitrary sample statistic (1 arg)
#In: population fucntion (1 arg which takes in number of samples per draw), a stat function (1 arg), 
# number of samples per draw, number of draws
#Out: a vector of sample calculations
boot.sampleDist = function(popFunc, statFunc, samplesPerDraw, numDraws) {
    return(matrix(sapply(rep(samplesPerDraw, numDraws), function(x) statFunc(popFunc(x))), numDraws));
}

#######################Problem 3 Script Entry Point###############
draws = 10000;
datSize = 50;
popGen = function(x) runif(x);
#define the statistic and population generator function
statGen = function(x) max(x);

dataSet = popGen(datSize);
#acquire the data set
datStat = statGen(dataSet);
problemDist = boot.NPBootStrapper(dataSet, statGen, datSize, draws);

cat("\n\nThe problem 3 non-parametric standard error is:\n");
cat(var(problemDist) ^ .5);
cat("\n\nThe pivotal non-parametric confidence interval is:\n");
CIPivotal = c(2 * datStat - quantile(problemDist, .025), 2 * datStat - quantile(problemDist, .975));
cat(CIPivotal);
hist(problemDist, 50, main = "Problem 3 Non-parametric Bootstrap Distribution");


problemDist = boot.PBootStrapper(function(x) runif(x, 0, datStat), statGen, datSize, draws);

cat("\n\nThe problem 3 parametric standard error is:\n");
cat(var(problemDist) ^ .5);
cat("\n\nThe pivotal parametric confidence interval is:\n");
CIPivotal = c(2 * datStat - quantile(problemDist, .025), 2 * datStat - quantile(problemDist, .975));
cat(CIPivotal);
hist(problemDist, 50, main = "Problem 3 Parametric Bootstrap Distribution");


sampleDist = boot.sampleDist(popGen, statGen, datSize, draws)
hist(sampleDist, 50, main = "Problem 3 Sample Distribution");