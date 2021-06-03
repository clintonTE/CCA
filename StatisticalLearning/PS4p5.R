require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
require(xtable)
require(randomForest)
require(glmnet)
require(reshape2)





#holds constants and program parameters
CONST = list(
    EPSILON = .Machine$double.eps, #machine precision
    MIN_DOUBLE = .Machine$double.xmin,
    NUM_WORKERS = max(round(detectCores() * .5), 2), #just a heuristic for multi-threading
    FILE_NAME = "Heart.csv",
    Y_NAME = "AHD",
    X_NAMES = ".",
    FORMULA = formula("AHD ~ ."),
    TRAIN = 150,
    MTRY = 3,
    NTREE = 1000
)

#reads in the sample and makes any data adjustmetns needed
prepsample = function(ntrain = CONST$TRAIN) {
    d = fread(file = CONST$FILE_NAME)

    if (CONST$X_NAMES == ".")
        xnames = setdiff(names(d), CONST$Y_NAME)
    else
        xnames = CONST$X_NAMES

    #break into test and training samples
    d = na.omit(d)

    train = sample(nrow(d), ntrain)
    test = setdiff(1:nrow(d), train)
    d$AHD = as.factor(d$AHD)
    d$Thal = as.factor(d$Thal)
    d$ChestPain = as.factor(d$ChestPain)

    S = list(d = d, dtrain = d[train,], dtest = d[test,], train = train, test = test, spec = CONST$FORMULA, yname = CONST$Y_NAME, xnames = xnames)

    return(S)
}

#primary funciton for estimating a forest
estimateForest = function(S, mtry = CONST$MTRY, ntree = CONST$NTREE) {

    rf = randomForest(formula = CONST$FORMULA, data = S$dtrain, mtry = mtry, ntree = ntree, importance = TRUE)
    rfpredict = predict(rf, newdata = S$d[!S$train])

    #get the appropriate error metrics
    ooberr = tail(rf$err.rate[, "OOB"], 1)
    ooserr = sum(rfpredict != S$dtest[, AHD]) / nrow(S$dtest)

    return(list(mod = rf, ooberr = ooberr, ooserr = ooserr))

}

p3c = function(S, verbose = true) {
    set.seed(11) #A seed for me
    S = prepsample()
    results = estimateForest(S)

    #this gets the results of problem
    print(results)

    return(S)
}

p3c()

p3d = function() {
    set.seed(11) #A seed for me
    S = prepsample()

    allmtrys = 1:(ncol(S$d) - 1)
    estimates = lapply(allmtrys, function(mtry) estimateForest(S, mtry = mtry, ntree = 500))

    ooberrs = sapply(estimates, function(r) r$ooberr)
    ooserrs = sapply(estimates, function(r) r$ooserr)

    resultswide = data.table(mtry = allmtrys, ooberrs = ooberrs, ooserrs = ooserrs)
    results = melt(resultswide, id.vars = "mtry")

    p = ggplot(results, aes(x = mtry, y = value, color = variable)) + geom_line() + theme_bw() +
            ggtitle("OOB and OOS errors vs mtry")

    print(p)
}

p3d()

estimateNet = function(S, alpha=1.0) {

    #setup training and test
    xmattrain = model.matrix(~. , S$dtrain[, S$xnames, with = FALSE])
    yvectrain = as.factor(S$dtrain[[S$yname]])
    xmattest = model.matrix(~., S$dtest[, S$xnames, with = FALSE])
    yvectest = as.factor(S$dtest[[S$yname]])


    #alpha=1 corresponds to full lasso
    net = cv.glmnet(xmattrain, yvectrain, alpha = alpha, family = "binomial", type.measure = "mse")
    netpredict = predict(net, newx = xmattest, type="class")

    #get the appropriate error metrics
    ooserr = sum(netpredict != yvectest) / nrow(S$dtest)

    return(list(mod = net, ooserr = ooserr))

}

p3e = function() {
    set.seed(11) #A seed for me
    S = prepsample()
    rfresults = estimateForest(S)
    #print(rfresults)

    lassoresults = estimateNet(S)
    cat("Lasso OOS Error: ", lassoresults$ooserr)


}

p3e()
