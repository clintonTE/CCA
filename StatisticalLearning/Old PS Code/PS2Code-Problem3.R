require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
set.seed(11) #A seed for me



#holds constants and program parameters
CONST = list(
    NUM_ROWS = 200,
    NUM_COLS = 3,
    NUM_SAMPLES = 4000,
    EPSILON = .Machine$double.eps, #machine precision
    NUM_WORKERS = max(round(detectCores()*.5),2) #just a heuristic for multi-threading
)

#This is the probit likelihood function
llikelihoodProbit = function(b, Y, X) {
    epsilon = CONST$EPSILON
    argvec = X %*% b

    #avoid numerical issues with logs of small numbers
    pnorms = pmin(pmax(pnorm(argvec), epsilon), 1.0 - epsilon) 

    #Use vectorized ifelse
    likes = ifelse(Y, log(pnorms), log(1-pnorms))

    return (sum(likes))
}

#this is the gradient of the previous
llikelihoodProbitGrad = function(b, Y, X) {
    epsilon = CONST$EPSILON
    argvec = X %*% b
    pnorms = pnorms = pmin(pmax(pnorm(argvec), epsilon), 1.0 - epsilon)
    dnorms = dnorm(argvec)

    #Use vectorized ifelse
    premults = ifelse(Y, (1 / pnorms), - (1 / (1 - pnorms))) * dnorms

    #R's equivelent to broadcast
    grads = apply(X,MARGIN = 2, function(x) x * premults)
    return(colSums(grads))
}


probitModel = function(Y, X, suppressIntercept = FALSE) {
    #make the intercept as needed
    if (!suppressIntercept) {
        if (min(X[, ncol(X)]) != 1 || max(X[, ncol(X)]) != 1) X = cbind(X, rep(1, nrow(X)))
    }

    # Get some convenience constants 
    R = nrow(Y)
    C = ncol(X)

    #initial value of b
    b = rep(1, C)

    #make single argument versions for optim
    ll = function(x) -1.0*llikelihoodProbit(x, Y, X)
    llgrad = function(x) -1.0*llikelihoodProbitGrad(x, Y, X)

    #call the optimizer
    opt = optim(b, ll, gr = llgrad, method = "BFGS", hessian = TRUE)
    if(opt$convergence != 0) print("WARNING! Optimizer did not converge")

    #Efficient matrix inversion
    U = chol(opt$hessian)
    UInv = solve(chol(opt$hessian))
    Sigma = t(UInv) %*% UInv

    #form the info we want into a named list
    prob = list(B = opt$par, llikelihood = opt$value, varB = diag(Sigma), seB = diag(Sigma)^0.5)
    return(prob)
}

#generates a test sample from the asymtotic distribution
testSample = function(R = CONST$NUM_ROWS, C = CONST$NUM_COLS, beta = 1 / (1:C)) {
    #pre-allocate
    X = matrix(rnorm(R * C), nrow = R, ncol = C)

    #create the Y vector
    Y = apply(X, 1, function(x) pnorm(x %*% beta))
    Y = rbinom(R, 1, Y)

    return (list(Y=Y,X=X))
}

#tests the model a single time and prints the results
testProbitModelOnce = function() {

    S = testSample()
    prob = probitModel(S$Y, S$X)
    print(prob)
}

testProbitModelOnce()

#this generates a multi-variate bootstrap sample
bootSample = function(Y, X) {
    R = nrow(X)

    #first pick the rows we will sample
    sampledRows = sample(1:R, R, replace = TRUE)

    #sample the rows
    Y = sapply(sampledRows, function(r) Y[r])
    X = matrix(sapply(sampledRows, function(r) X[r,]), nrow = R, byrow = TRUE)

    return(list(Y = Y, X = X))
}

examineProbitDistributions = function(N=CONST$NUM_SAMPLES) {
    #maybe this will take a while, so lets multi-thread (process)
    cl = makeCluster(CONST$NUM_WORKERS)
    clusterExport(cl = cl,
        varlist = c("llikelihoodProbit", "llikelihoodProbitGrad", "probitModel",
        "testSample", "CONST", "bootSample"))

    #get the primary sample and model
    S = testSample()
    prob = probitModel(S$Y, S$X)
    betasAsymp = data.table(method = "asymp", b1 = rnorm(N, mean = prob$B[1], sd = (prob$varB[1] ^ 0.5)),
        b2 = rnorm(N, mean = prob$B[2], sd = (prob$varB[2] ^ 0.5)),
        b3 = rnorm(N, mean = prob$B[3], sd = (prob$varB[3] ^ 0.5))
    )

    #Get the bootstrap samples and solve for the MLE
    bootSamples = parLapply(cl, 1:N, function(x) bootSample(S$Y, S$X))
    bootModels = parLapply(cl, bootSamples, function(s) probitModel(s$Y, s$X))
    betasBoot = data.table(method = "boot", b1 = sapply(bootModels, function(x) x$B[1]),
        b2 = sapply(bootModels, function(x) x$B[2]),
        b3 = sapply(bootModels, function(x) x$B[3]))
        
    #get the true samples
    trueSamples = parLapply(cl, 1:N, function(x) testSample())
    trueModels = parLapply(cl, trueSamples, function(s) probitModel(s$Y, s$X))
    betasTrue = data.table(method = "true", b1 = sapply(trueModels, function(x) x$B[1]),
        b2 = sapply(trueModels, function(x) x$B[2]),
        b3 = sapply(trueModels, function(x) x$B[3]))

    #combine into a ggplot2 friendly structure
    betas = rbind(betasAsymp, betasBoot, betasTrue)

    #plot the densities of the estimates
    p1 = ggplot(betas, aes(x = b1)) +
        geom_density(aes(group = method, color = method)) + theme_bw() +
        ggtitle("Distribution of asymptotic, bootstrap, and simulated true beta-1")

    p2 = ggplot(betas, aes(x = b2)) +
        geom_density(aes(group = method, color = method)) + theme_bw() +
        ggtitle("b2 distribution of asymptotic, bootstrap, and simulated true beta-2")
    
    p3 = ggplot(betas, aes(x = b3)) +
        geom_density(aes(group = method, color = method)) + theme_bw() +
        ggtitle("Distribution of asymptotic, bootstrap, and simulated true beta-3")
    print(p1)
    print(p2)
    print(p3)

    #cleanup
    stopCluster(cl)

}

#system.time(testProbitModelOnce())
system.time(examineProbitDistributions())