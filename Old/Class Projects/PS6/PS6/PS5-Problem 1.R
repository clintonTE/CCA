

#Creates a bootstrap distribution from a population distribution given an
#arbitrary sample statistic (1 arg)
#In: population data, a stat function (1 arg), number of samples per draw, number of draws
#Out: a vector of bootstraps
boot.bootStrapper = function(popSet, statFunc, samplesPerDraw, numDraws) {

    return(matrix(sapply(1:numDraws,
    function(x) statFunc(sample(popSet, samplesPerDraw, replace = TRUE))),
    numDraws));
}

#Creates a sample distribution from a population function given an
#arbitrary sample statistic (1 arg)
#In: population fucntion (1 arg which takes in number of samples per draw), a stat function (1 arg), 
# number of samples per draw, number of draws
#Out: a vector of sample calculations
boot.sampleDist = function(popFunc, statFunc, samplesPerDraw, numDraws) {
    return(matrix(sapply(rep(samplesPerDraw, numDraws), function(x) statFunc(popFunc(x))),numDraws));
}


########################Problem 1 Script entry point############################
population = 100;
draws = 10000;
popGen = function(x) rnorm(x, 5, 1); #define the statistic and population generator function
statGen = function(x) exp(mean(x));

dataSet = popGen(population); #acquire the data set
popStat = statGen(dataSet);
problemDist = boot.bootStrapper(dataSet, statGen, population, draws);
hist(problemDist, 50, main = "Problem 1 Bootstrap Distribution");

cat("\nThe problem 1 standard error is:\n");
cat(var(problemDist) ^ .5);
cat("\n\nThe pivotal confidence interval is:\n")
CIPivotal = c(2*popStat - quantile(problemDist, .025), 2*popStat- quantile(problemDist, .975));
cat(CIPivotal);

problemDist = boot.sampleDist(popGen, statGen, population, draws)
hist(problemDist, 50, main = "Problem 1 Sample Distribution");

#######################Problem 2 Script Entry Point###############
population = 50;
popGen = function(x) runif(x);
#define the statistic and population generator function
statGen = function(x) max(x);

dataSet = popGen(population);
#acquire the data set
popStat = statGen(dataSet);
problemDist = boot.bootStrapper(dataSet, statGen, population, draws);
CIPivotal = c(2 * popStat - quantile(problemDist, .025), 2 * popStat - quantile(problemDist, .975));

cat("\n\nThe problem 2 standard error is:\n");
cat(var(problemDist) ^ .5);
cat("\n\nThe pivotal confidence interval is:\n")
cat(CIPivotal);

hist(problemDist, 50, main = "Problem 2 Bootstrap Distribution");
problemDist = boot.sampleDist(popGen, statGen, population, draws)
hist(problemDist, 50, main = "Problem 2 Sample Distribution");


