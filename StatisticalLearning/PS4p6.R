require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
require(xtable)
require(randomForest)
require(glmnet)
require(e1071)
require(reshape2)
require(rpart)

setwd("C:\\Users\\Clinton\\Dropbox\\Projects\\StatisticalLearning\\StatisticalLearning")




#holds constants and program parameters
CONST = list(
    NUM_WORKERS = max(round(detectCores() * .5), 2), #just a heuristic for multi-threading
    DATA_NAME = "pset4words.RData",
    TABLE_NAME = "words.rbin",
    Y_NAME = "vote",
    X_NAMES = ".",
    FORMULA = formula("vote ~ ."),
    NTREE = 1000, #set to 1000 in production
    PER_TRAINING = .75, #set to .75 in production
    SVM_COST = 1,
    RF_W = 3,
    LASSO_W = 2,
    SVM_W = 2,
    TREE_W = 1

)

#no need to call this every time
prepDataSet = function() {
    load(CONST$DATA_NAME)

    vote = pmax(vote, 0)
    d = cbind(data.table(vote = vote), data.table(wordMatrix))
    d$vote = pmax(d$vote, 0) #gets rid of the annoying negative values

    #this was cruel to debug
    setnames(d, "else", "elsenew")
    setnames(d, "break", "breaknew")
    setnames(d, "function", "functionnew")
    setnames(d, "repeat", "repeatnew")

    #names(mydf)[names(mydf) == "MyName.1"] = "MyNewName"

    #setDT(data)
    for (j in names(d)) {
        set(d, i = NULL, j = j, value = as.factor(d[[j]])) #convert all columns to factors
        if (length(unique(d[[j]])) == 1) #delete columns with the same value
            d[[j]] == NULL
    }

    save(d, file = CONST$TABLE_NAME)



}

#reads in the sample and makes any data adjustmetns needed
prepsample = function(ntrain = CONST$TRAIN, prep = FALSE) {
    if (prep) prepDataSet()

    load(CONST$TABLE_NAME, verbose = FALSE)



    #break into test and training samples
    d = na.omit(d)

    train = sample(nrow(d), round(nrow(d) * CONST$PER_TRAINING))
    test = setdiff(1:nrow(d), train)

    for (j in names(d)) {
        if (length(unique(d[[j]][train] )) == 1) #delete columns with the same value
        {
            d[, c(j) := NULL]
        }
            
    }

    #get the names of x for when we need them
    if (CONST$X_NAMES == ".")
        xnames = setdiff(names(d), CONST$Y_NAME)
    else
        xnames = CONST$X_NAMES

    S = list(d = d, dtrain = d[train,], dtest = d[test,], train = train, test = test, spec = CONST$FORMULA, yname = CONST$Y_NAME, xnames = xnames)

    return(S)
}

#primary funciton for estimating a forest
estimateForest = function(S, mtry = round(sqrt(ncol(S$d))), ntree = CONST$NTREE) {
    set.seed(11) #A seed for me

    rf = randomForest(formula = CONST$FORMULA, data = S$dtrain, mtry = mtry, ntree = ntree)
    rfpredict = predict(rf, newdata = S$d[!S$train])

    #get the appropriate error metrics
    ooberr = tail(rf$err.rate[, "OOB"], 1)
    ooserr = sum(rfpredict != S$dtest[[S$yname]]) / nrow(S$dtest)

    cat("rfresult OOS Error: ", ooserr, "\n")
    return(list(modtype = "randomforest", pred = as.numeric(as.character(rfpredict)), mod = rf, ooberr = ooberr, ooserr = ooserr))

}

estimateNet = function(S, alpha = 1.0) {
    set.seed(11) #A seed for me
    #setup training and test
    xmattrain = model.matrix(~., S$dtrain[, S$xnames, with = FALSE])
    yvectrain = as.factor(S$dtrain[[S$yname]])
    xmattest = model.matrix(~., S$dtest[, S$xnames, with = FALSE])
    yvectest = as.factor(S$dtest[[S$yname]])


    #alpha=1 corresponds to full lasso)
    net = cv.glmnet(xmattrain, yvectrain, alpha = alpha, family = "binomial", type.measure = "mse")
    netpredict = predict(net, newx = xmattest, type = "class")

    #get the appropriate error metrics
    ooserr = sum(netpredict != yvectest) / nrow(S$dtest)

    cat("Lasso OOS Error: ", ooserr, "\n")
    return(list(modtype = "lasso", mod = net, pred = as.numeric(as.character(netpredict)), ooserr = ooserr))

}

estimateSVM = function(S, tunesvm = FALSE, kern = "linear") {
    set.seed(11) #A seed for me
    if (tunesvm) {
        tuning = tune(svm, S$spec, data = S$dtrain, kernel = kern)
        mod = tuning$best.model
    } else {
        mod = svm(S$spec, data = S$dtrain, kernel = kern, cost = CONST$SVM_COST)
    }
    #setup training and test
    svmpredict = predict(mod, S$dtest)

    ooserr = sum(svmpredict != S$dtest[[S$yname]]) / nrow(S$dtest)

    cat("SVM OOS Error: ", ooserr, "\n")
    return(list(modtype = "svm", mod = mod, pred = as.numeric(as.character(svmpredict)), ooserr = ooserr))
}

#estimates the tree
estimatetree = function(S) {
    set.seed(11) #A seed for me
    inittree = rpart(S$spec, data = S$dtrain)
    idxmin = which.min(inittree$cptable[, "xerror"])
    besttree = prune(inittree, cp = inittree$cptable[idxmin, "CP"])
    treepredict = predict(besttree, newdata=S$dtest, type="class")

    #get the appropriate error metrics
    ooserr = sum(treepredict != S$dtest[[S$yname]]) / nrow(S$dtest)

    cat("Tree OOS Error: ", ooserr, "\n")
    return(list(modtype = "tree", mod = besttree, pred = as.numeric(as.character(treepredict)), ooserr = ooserr))
}

#runs all the models
runmodels = function(S) {
    set.seed(11) #A seed for me

    rfresults = estimateForest(S)
    lassoresults = estimateNet(S)
    svmresults = estimateSVM(S)
    treeresults = estimatetree(S)

    ##now build the ensemble
    maxw = CONST$RF_W + CONST$LASSO_W + CONST$SVM_W + CONST$TREE_W
    threshold = maxw / 2
    fallbackthreshold = (CONST$RF_W + CONST$LASSO_W + CONST$SVM_W)/2

    rawensemble = (rfresults$pred * CONST$RF_W +
        lassoresults$pred * CONST$LASSO_W +
        svmresults$pred * CONST$SVM_W +
        treeresults$pred * CONST$TREE_W)

    tiebreakerfallback = (rfresults$pred * CONST$RF_W +
        lassoresults$pred * CONST$LASSO_W +
        svmresults$pred * CONST$SVM_W)

    ensemble = ifelse(rawensemble == threshold, tiebreakerfallback > fallbackthreshold, rawensemble > threshold)
    ensemblepred = ifelse(ensemble, 1, 0)
    ooserr = sum(ensemblepred != S$dtest[[S$yname]]) / nrow(S$dtest)

    ensembleresults = list(modtype = "ensemble", pred = ensemblepred, ooserr = ooserr)
    cat("Ensemble OOS Error: ", ooserr, "\n")

    modeltypes = c(rfresults$modtype, lassoresults$modtype, svmresults$modtype, treeresults$modtype, ensembleresults$modtype)
    ensembleWeights = c(CONST$RF_W, CONST$LASSO_W, CONST$SVM_W, CONST$TREE_W, NA)
    ooserrs = c(rfresults$ooserr, lassoresults$ooserr, svmresults$ooserr, treeresults$ooserr, ensembleresults$ooserr)

    outtab = data.table(type = modeltypes, ooserr = ooserrs, ensembleW = ensembleWeights)
    xtable(outtab, type = "html")



}

p6 = function() {
    set.seed(3) #A seed for me
    S = prepsample()

    results = runmodels(S)

}

p6()



