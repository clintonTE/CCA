require(ggplot2) #for graphs
require(sandwich) #standard error
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data


set.seed(11) #A seed for me



#holds constants and program parameters
CONST = list(
    NUM_ROWS = 1000,
    NUM_SAMPLES = 10000,
    EPSILON = .Machine$double.eps, #machine precision
    NUM_WORKERS = max(round(detectCores() * .5), 2), #just a heuristic for multi-threading
    NBSIZE = 1 / 3

    )

#generates a test sample from the asymtotic distribution
binomSample = function(R = CONST$NUM_ROWS, nbsize = CONST$NBSIZE) {
    #pre-allocate
    X = rnorm(R)

    #create the Y vector
    Y = exp(X / 100)
    Y = rnbinom(R, size = CONST$NBSIZE, mu = Y)

    return(list(Y = Y, X = X))
}

glmPoisson = function(Y, X) {
    pois = glm(Y ~ X, family = poisson())
    beta = pois$coefficients

    #Get the standard SE
    AInv = vcov(pois)
    SE = diag(AInv) ^ 0.5
    names(SE) = names(beta)


    #Get the robust SE
    score = estfun(pois)
    B = t(score) %*% score
    AInvBAInv = AInv %*% B %*% AInv
    SERobust = diag(AInvBAInv) ^ 0.5
    names(SERobust) = names(beta)

    return(list(beta = beta, SE = SE, SERobust = SERobust))
}

compareModelsOnce = function() {
    S = binomSample() #get the main sample
    pois = glmPoisson(S$Y, S$X) #get the model output
    print(pois) #print it
}

compareModelsOnce()

simulateModel = function(K = CONST$NUM_SAMPLES) {
    #prepare to parallelize
    cl = makeCluster(CONST$NUM_WORKERS)
    clusterExport(cl = cl, varlist = c("glmPoisson", "binomSample", "CONST"))
    clusterEvalQ(cl, require(sandwich))

    #first get the samples
    samples = parLapply(cl, 1:K, function(x) binomSample())
    simModels = parLapply(cl, samples, function(x) glmPoisson(x$Y, x$X))
    betas = data.table(method = "simulation", b1 = sapply(simModels, function(x) x$beta[1]),
        b2 = sapply(simModels, function(x) x$beta[2]))

    #plot
    p1 = ggplot(betas, aes(x = b1)) +
        geom_density(aes(group = method, color = method)) + theme_bw() +
        ggtitle("Distribution of simulated beta-1")

    p2 = ggplot(betas, aes(x = b2)) +
        geom_density(aes(group = method, color = method)) + theme_bw() +
        ggtitle("Distribution of simulated beta-2")

    print(p1)
    print(p2)

    cat("Cross-sectional standard deviation of b1: ", sd(betas[, b1]), "\n")
    cat("Cross-sectional standard deviation of b2: ", sd(betas[, b2]), "\n")

    #cleannup
    stopCluster(cl)
}


system.time(simulateModel())