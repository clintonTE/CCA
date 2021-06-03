require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
require(xtable)
set.seed(10) #A seed for me



#holds constants and program parameters
CONST = list(
    EPSILON = .Machine$double.eps, #machine precision
    MIN_DOUBLE = .Machine$double.xmin,
    NUM_WORKERS = 7#max(round(detectCores() * .5), 2) #just a heuristic for multi-threading
)

predict0 = function(m, dt) {
    options(warn = -1) #disable the rank defficient warning
    p = predict(m, dt)
    options(warn = 1)

    return(p)
}

#univariate k-fold polynomial cross-validation
polygnome = function(dt, degree, k, yname = "Y", xname = "X") {
    print("got here")
    X = dt[[xname]]
    n = length(X)
    
    #create the formula
    rhs = "intercept"
    coefs = rhs
    if (degree > 0) {
        rhs = paste0(c(rhs, " + ", xname))
        coefs = append(coefs, rhs)
    }


    if (degree > 1) {
        for (i in 1:degree) {

            coefname = paste0(xname, i, collapse = "")
            dt[, eval(coefname)] = X ^ i
            rhs = paste0(c(rhs, " + ", coefname))
            coefs = append(coefs, rhs)
        }
    }

    #handle the intercept ourselves
    rhs = paste0(c(rhs, " + 0"))

    #create the formula
    dt$intercept = rep(1,n)
    spec = formula(paste0(c(yname, " ~ ", rhs), collapse = ""))


    #create the folds
    foldwidth = trunc(n / k)
    dt$fold = rep(k,n)
    r = 1
    for (i in 1:(k - 1)) {
        for (j in 1:foldwidth) {
            dt$fold[r] = i
            r = r + 1
        }
    }

    #cross-validate
    folderrors = rep(0.0, k)
    for (i in 1:k) {
        m = lm(spec, dt[fold != i])
        p = predict0(m, dt[fold==i])
        folderrors[i] = sum((dt[fold == i,Y] - p)^2)
    }
    secross = sum(folderrors)/n

    #compute the in-sample variance
    m = lm(spec, dt)
    p = predict0(m, dt)
    sein = mean((p-dt[,Y])^2)

    return(list(model=m, secross=secross, sein = sein))
}

#form the data
formdgps = function(n) {
    X = runif(n, min = -4, max = 4)
    epsilon = rnorm(n)

    Y1 = -2 * (X < -3) + 2.55 * (X > -2) - 2 * (X > 0) + 4 * (X > 2) - 1 * (X > 3) + epsilon
    Y2 = 6 + .4 * X - .36 * X ^ 2 + .005 * X ^ 3 + epsilon
    Y3 = 2.83 * sin(pi / 2 * X) + epsilon
    Y4 = 4 * sin(3 * pi * X) * (X > 0) + epsilon

    return(list(dgp1 = data.table(Y = Y1, X = X),
      dgp2 = data.table(Y = Y2, X = X),
      dgp3 = data.table(Y = Y3, X = X),
      dgp4 = data.table(Y = Y4, X = X)))
}

testpolignome = function(n, showplots = TRUE) {
    cl = makeCluster(CONST$NUM_WORKERS)
    clusterExport(cl = cl,
        varlist = c("testpolignome", "polygnome", "CONST", "predict0", "formdgps", "data.table"))

    dgpn = formdgps(n)

    processpar = function(dgp) {
        return(parLapply(cl, 0:10, function(d) polygnome(dgp, degree = d, k = 10)))
    }

    #clusterExport(cl=cl, varlist = c("processpar"))


    #run polygnome for each polynomial size and data generating process
    results010list = lapply(dgpn, processpar )

    #recast the results in a usable form
    results010 = lapply(results010list,
        function(results) melt(data.table(d = 0:10,
            secross = sapply(0:10, function(x) results[[x + 1]]$secross),
            sein = sapply(0:10, function(x) results[[x + 1]]$sein)), id.vars = "d"))

    #make the plots
    d = rep(-1,4)
    for (i in 1:4) {
        p = ggplot(results010[[i]], aes(x = d, y = value, color = variable)) + geom_line() + theme_bw() +
            ggtitle(paste0(c("dgp-", i, " ", "(n=", n, ")"), collapse = ""))

        if (showplots) print(p)
        d[i] = which.min(results010[[i]][variable == "secross", value])
        print(paste0("dsg-",i, " min mse: degree=", d[i]-1, "  for n = ", n))
    }

    #extract the best models
    models = mapply(function(resultlist, degree) resultlist[[degree]]$model, results010list, d, SIMPLIFY = FALSE)

    #cleanup
    stopCluster(cl)

    return(list(dgpn = dgpn, models = models, d = d))
}

####problem 3.2
p32results = testpolignome(n = 100)

plotpredictions = function(dgpn, models, degrees) {

    n = nrow(dgpn[[1]])

    #helper function to make the dt containing the predictions
    predictdt = function(dgp, model, degree) {
        X = dgp[,X]
        if (degree > 1) {
            for (i in 1:degree) {
                coefname = paste0("X", i, collapse = "")
                dgp[, eval(coefname)] = X ^ i
            }
        }

        dgp$intercept = rep(1, n)

        p = predict0(model, dgp)
        dt = data.table(X = dgp[, X], actual = dgp[, Y], predictions = p)

        return (melt(dt, id.vars="X"))
    }


    predtables = mapply(predictdt, dgpn, models, degrees, SIMPLIFY = FALSE)

    for (i in 1:4) {
        p = ggplot(predtables[[i]], aes(x = X, y = value, color = variable)) + geom_line() + theme_bw() +
            ggtitle(paste0(c("dgp-", i, " ", "(n=", n, ") Predicted and Actual Values, degree=", degrees[i]-1), collapse = ""))
        print(p)
    }

    return(NULL)
}

###problem 3.3
plotpredictions(dgpn = p32results$dgpn, models = p32results$models, degrees=p32results$d)

###problem 3.4
system.time(testpolignome(n = 10000, showplots = FALSE))

