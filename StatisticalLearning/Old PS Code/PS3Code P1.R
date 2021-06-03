require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
require(xtable)
set.seed(10) #A seed for me



#holds constants and program parameters
CONST = list(
    EPSILON = .Machine$double.eps, #machine precision
    MIN_DOUBLE = .Machine$double.xmin,
    NUM_WORKERS = max(round(detectCores() * .5), 2), #just a heuristic for multi-threading
    X_NAMES = c("l.physint", "l.polity2", "l.lneconaid", "l.lnrgdp", "l.lnpop", "colony"),
    Y_NAME = "lneconaid",
    CI95 = 1.96,
    FILE_NAME = "nielsenaid.csv"
)

#This is the tobit loglikelihood function
llikelihoodTobit = function(bs, Y, X) {
    epsilon = CONST$EPSILON
    k = length(bs)
    b = bs[1:(k - 1)]
    s = max(abs(bs[k]), epsilon)
    #s = b[k]

    likes = ifelse(Y > 0,
        dnorm((Y - (X %*% b)) / s) / s, 1 - pnorm((X %*% b) / s))
    likes = log(likes)
    likes = ifelse(is.finite(likes), likes, -10^10)

    return(sum(likes))
}




tobitModel = function(Y, X, yname=NULL, xnames=NULL, suppressintercept = FALSE) {

    pnames = xnames

    #make the intercept as needed
    if (!suppressintercept) {
        if (min(X[, ncol(X)]) != 1 || max(X[, ncol(X)]) != 1) {
            X = cbind(X, rep(1, nrow(X))) #bind the intercept column
            if (!is.null(pnames)) {
                pnames = c(pnames, "intercept") #adjust the names
            }
        }
    }
    pnames = c(pnames, "sigma")

    # Get some convenience constants 
    R = nrow(Y)
    C = ncol(X)

    #initial value of b
    bs = rep(.001, C+1)

    #make single argument versions for optim
    ll = function(x) - 1.0 * llikelihoodTobit(x, Y, X)

    #call the optimizer
    opt = optim(bs, fn=ll, 
        method = "BFGS", hessian = TRUE)
    if (opt$convergence != 0) print("WARNING! Optimizer did not converge")

    print(opt$par)
    #Efficient matrix inversion
    U = chol(opt$hessian)
    UInv = solve(chol(opt$hessian))
    Sigma = t(UInv) %*% UInv

    #form the info we want into a named list
    tob = list(B = opt$par, llikelihood = opt$value, varB = diag(Sigma), seB = diag(Sigma) ^ 0.5, k = length(opt$par))

    if (!is.null(pnames)) {
        names(tob$B) = pnames #name the coefficients
        names(tob$seB) = pnames
    }

    return(tob)
}

prepsample = function() {
    d = fread(file = CONST$FILE_NAME)

    X = as.matrix(d[, CONST$X_NAMES, with = FALSE])
    Y = as.matrix(d[, CONST$Y_NAME, with=FALSE])

    S = list(Y = Y, X = X, yname = CONST$Y_NAME, xnames = CONST$X_NAMES)

    return (S)
}

#tests the model a single time and prints the results
runTobitModelOnce = function() {

    S=prepsample()
    tob = tobitModel(Y = S$Y, X = S$X, xnames = S$xnames, yname = S$yname)
    #print(tob)

    testvals = as.matrix(sapply(1:ncol(S$X), function(x) median(S$X[, x]))) #set the test vector as the median of the X values
    testvals[length(testvals)] = 0 #mode for colony
    testvals[length(testvals) + 1] = 1

    predictTobit = function(xb, s) pnorm(xb/s)*xb + dnorm(xb/s) * s

    b = as.matrix(tob$B[1:(tob$k - 1)])
    s = tob$B["sigma"]

    print(b)
    print(s)
    print(testvals)
    print(t(b) %*% testvals)

    results = list(
    predicted = predictTobit(t(b) %*% testvals, s),
    CI = c(predictTobit(t(b) %*% testvals - s * CONST$CI95, s), predictTobit(t(b) %*% testvals + s * CONST$CI95, s)))
    print(results)

    outdf = data.table(coef = c(S$xnames, "intercept", "sigma"), value = tob$B, SE = tob$seB)
    xtable(t(outdf))
}

runTobitModelOnce()