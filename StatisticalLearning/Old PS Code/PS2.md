---
title: "Problem Set 2"
output: html_document
author: Clinton Tepper
---

Note: On the previous assignment you asked why I used the notation
$\int\left(\cdot\right)p\left(s\right)ds$. Generally I find using
a new variable as a dummy inside the integrand enhances clarity. If
it makes the results less clear, let me know and I will stop using
dummy variables.


## Problem 1

### 1. 

* Write down the total likelihood:
$$
\begin{aligned}
L\left(\hat{\theta}|X\right)= & \prod_{i\in1:n}\iota\left(\hat{\theta}\ge x_{i}\right)\frac{1}{\hat{\theta}}\\
l\left(\hat{\theta}|X\right)= & \begin{cases}
\sum_{i\in1:n}\frac{1}{\hat{\theta}} & \hat{\theta}\ge\max\left\{ x_{i}\right\} _{1}^{N}\\
-\infty & otherwise
\end{cases}
\end{aligned}
$$
* Since $l$ is strictly decreasing for all $\hat{\theta}\ge\max\left\{ x_{i}\right\} _{1}^{N}$,
we have $\hat{\theta}_{MLE}=\max\left\{ x_{i}\right\} _{1}^{N}$

### 2.

* Denote $\hat{\theta}\equiv\max\left\{ x_{i}\right\} _{1}^{N}$. Note
also that $n>1$. Note we can re-write the problem in the following
form:
$$
\begin{aligned}
0.95= & \gamma\int_{0}^{\kappa}L\left(\tau|X\right)d\tau\\
= & \gamma\int_{0}^{\kappa}\prod_{i\in1:n}\iota\left(\tau\ge x_{i}\right)\frac{1}{\tau}d\tau=0.95\\
= & \gamma\int_{\hat{\theta}}^{\kappa}\frac{1}{\tau^{n}}d\tau\\
= & \frac{\gamma}{\left(n-1\right)}\left(\frac{1}{\hat{\theta}^{n-1}}-\frac{1}{\kappa^{n-1}}\right)\\
\kappa^{n-1}= & \frac{1}{\frac{1}{\hat{\theta}^{n-1}}-\frac{0.95\left(n-1\right)}{\gamma}}
\end{aligned}
$$
* We need to solve for $\gamma$. Impose the properties of a CDF:
$$
\begin{aligned}
1 & =\lim_{\kappa\to\infty}\frac{\gamma}{\left(n-1\right)}\left(\frac{1}{\hat{\theta}^{n-1}}-\frac{1}{\kappa^{n-1}}\right)\\
\hat{\theta}^{n-1} & =\frac{\gamma}{\left(n-1\right)}\\
\gamma & =\hat{\theta}^{n-1}\left(n-1\right)\\
\kappa^{n-1}= & \frac{\hat{\theta}^{n-1}}{1-0.95}
\end{aligned}
$$
* Plugging in, $\hat{\theta}=2.85$ and $\gamma=263.9$, so $\kappa=6.03$

### 3.

* We have:
$$
\begin{aligned}
\frac{L\left(\theta|X\right)}{L\left(\hat{\theta}|X\right)}= & \frac{\hat{\theta}^{n}}{\theta^{n}}\forall\theta>\hat{\theta}\\
p\left(\frac{L\left(\theta|X\right)}{L\left(\hat{\theta}|X\right)}<c\right)= & 1-p\left(\frac{L\left(\theta|X\right)}{L\left(\hat{\theta}|X\right)}>c\right)\\
= & \gamma_{r}\int_{\hat{\theta}}^{\kappa}\frac{\hat{\theta}^{n}}{\tau^{n}}d\tau\\
= & \frac{\gamma_{r}\hat{\theta}^{n}}{n-1}\left(\frac{1}{\hat{\theta}^{n-1}}-\frac{1}{\kappa^{n-1}}\right)\\
= & \frac{\gamma_{r}\hat{\theta}^{n}}{n-1}\left(\frac{1}{\hat{\theta}^{n-1}}-\frac{1-c}{\hat{\theta}^{n-1}}\right)\\
= & \frac{\gamma_{r}\hat{\theta}c}{n-1}
\end{aligned}
$$
* Solve for $\gamma_{r}$ as before:
$$
\begin{aligned}
1= & \gamma_{r}\lim_{\kappa\to\infty}\int_{\hat{\theta}}^{\kappa}\frac{\hat{\theta}^{n}}{\tau^{n}}d\tau\\
1= & \frac{\gamma_{r}\hat{\theta}}{n-1}\\
\gamma_{r}= & \frac{n-1}{\hat{\theta}}
\end{aligned}
$$
* So
$$
\begin{aligned}
p\left(\frac{L\left(\theta|X\right)}{L\left(\hat{\theta}|X\right)}>c\right)= & 1-c\checkmark
\end{aligned}
$$


## Problem 2

### 1.

* This is an odds ratio of odds ratios. 
* The numerator is the ratio of the probability of someone voting to
the probability of someone not voting given that they had a high level
of education and covariates $w$. 
* The denominator is the ratio of the probability of someone voting
to the probability of someone not voting given that they had low education
and covariates $w$. 
* The overall expression is the ratio of the odds ratio of voting given
high education to the odds ratio given low education, all given covariates
$w$.

### 2.

* This is simply plugging into the provided assumption $p\left(Y_{i}=1|X_{i}\right)=\frac{\exp\left(X_{i}'\beta\right)}{1+\exp\left(X_{i}'\beta\right)}$:
$$
\begin{aligned}
p\left(Y_{i}=1|T_{i}=1,W_{i}=w\right)= & \frac{\exp\left(\alpha+\gamma+\delta'w\right)}{1+\exp\left(\alpha+\gamma+\delta'w\right)}\\
p\left(Y_{i}=0|T_{i}=1,W_{i}=w\right)= & \frac{1}{1+\exp\left(\alpha+\gamma+\delta'w\right)}\\
p\left(Y_{i}=1|T_{i}=0,W_{i}=w\right)= & \frac{\exp\left(\alpha+\delta'w\right)}{1+\exp\left(\alpha+\delta'w\right)}\\
p\left(Y_{i}=1|T_{i}=0,W_{i}=w\right)= & \frac{1}{1+\exp\left(\alpha+\delta'w\right)}
\end{aligned}
$$
* Plugging in:
$$
\begin{aligned}
OR\left(w\right)= & \frac{p\left(Y_{i}=1|T_{i}=1,W_{i}=w\right)/p\left(Y_{i}=0|T_{i}=1,W_{i}=w\right)}{p\left(Y_{i}=1|T_{i}=0,W_{i}=w\right)/p\left(Y_{i}=0|T_{i}=0,W_{i}=w\right)}\\
= & \frac{\exp\left(\alpha+\gamma+\delta'w\right)}{\exp\left(\alpha+\delta'w\right)}\\
= & \exp\left(\gamma\right)
\end{aligned}
$$
* This allows us to interpret the estimated coefficient $\hat{\gamma}$
as an estimate of the log of the odds ratio of interest.

### 3.

* By the continuous mapping theorem, $\hat{\sigma}_{n}\stackrel{p}{\to}\sigma$
implies $\hat{\sigma}_{n}^{2}\stackrel{p}{\to}\sigma^{2}$
* As the standard error exists and $\hat{\gamma}_{n}\stackrel{p}{\to}\gamma$
and $\hat{\sigma}_{n}^{2}\stackrel{p}{\to}\sigma^{2}$, we can apply
the central limit theorem:
$$
\begin{aligned}
\sqrt{n}\left(\hat{\gamma}_{n}-\gamma\right)\stackrel{d}{\rightarrow} & N\left(0,\;\sigma^{2}\right)
\end{aligned}
$$
* Hence the delta method provides the asymptotic distribution:
$$
\begin{aligned}
\sqrt{n}\left(g\left(\hat{\gamma}_{n}\right)-g\left(\gamma\right)\right)\stackrel{d}{\rightarrow} & N\left(0,\;\sigma^{2}g'\left(\gamma\right)^{2}\right)\\
\sqrt{n}\left(e^{\hat{\gamma}_{n}}-e^{\gamma}\right)\stackrel{d}{\rightarrow} & N\left(0,\;\sigma^{2}e^{2\gamma}\right)
\end{aligned}
$$


\subsection*{Problem 2}

### 4.

* We have $X_{i}'\beta=\pi_{i}$, so the link function is given by $g\left(\mu\right)=\mu$.

### 5.

* We have
$$
\begin{aligned}
E\left[\varepsilon_{i}^{2}|X_{i}\right]= & E\left[\left(Y_{i}-X_{i}'\beta\right)^{2}|X_{i}\right]\\
= & E\left[Y_{i}^{2}|X_{i}\right]-E\left[Y_{i}X_{i}\beta|X_{i}\right]+\left(X_{i}'\beta\right)^{2}\\
= & E\left[Y_{i}^{2}|X_{i}\right]-\left(X_{i}'\beta\right)^{2}
\end{aligned}
$$
* Since $p\left(Y_{i}|X_{i}\right)\sim B\left(\pi_{i}\right)$, we have
$$
\begin{aligned}
E\left[\varepsilon_{i}^{2}|X_{i}\right]= & V\left[Y_{i}|X_{i}\right]+\pi_{i}^{2}-\pi_{i}^{2}\\
= & \pi_{i}\left(1-\pi_{i}\right)
\end{aligned}
$$
* Thus outside of some degenerate cases (e.g. $\pi_{i}$ is a constant),
the error terms are conditionally heteroskedastic.
* Per White 1980, if exogeneity holds such that $E\left[\varepsilon|X_{i}\right]=0$,
corrected estimators of the standard error are asymptotically consistent. 

### 6.

* The likelihood is given by:
$$
\begin{aligned}
L\left(\theta|X\right)= & \prod_{i\in1:P}\pi_{i}^{Y_{i}}\left(1-\pi_{i}\right)^{1-Y_{i}}\\
= & \prod_{i\in1:P}\left(X_{i}\beta\right)^{Y_{i}}\left(1-X_{i}'\beta\right)^{1-Y_{i}}
\end{aligned}
$$

### 7.

* Let $n_{Y}=\sum Y_{i}$. Then the log likelihood is:
$$
\begin{aligned}
l\left(\theta|X\right)= & \sum_{i\in1:P}\left[Y_{i}ln\left(X_{i}'\beta\right)+\left(1-Y_{i}\right)ln\left(1-X_{i}'\beta\right)\right]
\end{aligned}
$$
* The score is thus
$$
\begin{aligned}
S\left(\beta|X\right)=\nabla l\left(\theta|X\right)= & \sum_{i\in1:P}\left[\frac{Y_{i}}{X_{i}'\beta}X_{i}-\frac{1-Y_{i}}{1-X_{i}'\beta}X_{i}\right]\\
= & \sum_{i\in1:P}\left[\frac{Y_{i}\left(1-X_{i}'\beta\right)-X_{i}'\beta\left(1-Y_{i}\right)}{X_{i}'\beta\left(1-X_{i}'\beta\right)}X_{i}\right]\\
= & \sum_{i\in1:P}\left[\frac{Y_{i}-X_{i}'\beta}{X_{i}'\beta\left(1-X_{i}'\beta\right)}X_{i}\right]\checkmark
\end{aligned}
$$

### 8.

* Under correct specification, 
$$
\begin{aligned}
\beta_{MLE}\sim & N\left(\beta,\;-E\left[H\right]^{-1}\right)
\end{aligned}
$$
* But under correct specificaiton, 
$$
\begin{aligned}
-E\left[H\right]^{-1}= & I_{N}\left(\beta|X\right)^{-1}\\
= & E\left[S\left(\beta|X_{i}\right)S\left(\beta|X_{i}\right)'\right]^{-1}
\end{aligned}
$$
* Then
$$
\begin{aligned}
S\left(\beta|X_{i}\right)S\left(\beta|X_{i}\right)'= & \left[\frac{Y_{i}-X_{i}'\beta}{X_{i}'\beta\left(1-X_{i}'\beta\right)}\right]^{2}X_{i}X_{i}'\\
E\left[S\left(\beta|X_{i}\right)S\left(\beta|X_{i}\right)'\right]= & E\left[\left[\frac{Y_{i}-X_{i}'\beta}{X_{i}'\beta\left(1-X_{i}'\beta\right)}\right]^{2}X_{i}X_{i}'|X_{i}\right]\\
= & E\left[\left[\frac{Y_{i}-X_{i}'\beta}{X_{i}'\beta\left(1-X_{i}'\beta\right)}\right]^{2}|X_{i}\right]X_{i}X_{i}'\\
= & \frac{V\left(\varepsilon|X_{i}\right)}{\left[X_{i}'\beta\left(1-X_{i}'\beta\right)\right]^{2}}X_{i}X_{i}'\\
= & \frac{X_{i}'\beta\left(1-X_{i}'\beta\right)}{\left[X_{i}'\beta\left(1-X_{i}'\beta\right)\right]^{2}}X_{i}X_{i}'\text{ (from Q5)}\\
= & \frac{X_{i}X_{i}'}{X_{i}'\beta\left(1-X_{i}'\beta\right)}
\end{aligned}
$$
* Plugging in, we thus have
$$
\begin{aligned}
\beta_{MLE}\sim & N\left(\beta,\;\left[\frac{X_{i}X_{i}'}{X_{i}'\beta\left(1-X_{i}'\beta\right)}\right]^{-1}\right)
\end{aligned}
$$


\subsection*{Problem 3}

### 9. 

* The likelihood and log-likelihood are given by:
$$
\begin{aligned}
L\left(\theta|X\right)= & \prod_{i\in1:P}\Phi\left(X_{i}'\beta\right)^{Y_{i}}\left(1-\Phi\left(X_{i}'\beta\right)\right)^{1-Y_{i}}\\
l\left(\theta|X\right)= & \sum_{i\in1:P}\left[Y_{i}\ln\left[\Phi\left(X_{i}'\beta\right)\right]+\left(1-Y_{i}\right)\ln\left(1-\Phi\left(X_{i}'\beta\right)\right)\right]
\end{aligned}
$$

### 10. and 11.

(Sorry for the wacky R code. Julia is my main language.)

* Need gradient for reliable optimization:
$$
\begin{aligned}
\nabla l\left(\theta|X\right)= & \left[\frac{Y_{i}}{\Phi\left(X_{i}'\beta\right)}-\frac{\left(1-Y_{i}\right)}{1-\Phi\left(X_{i}'\beta\right)}\right]\phi\left(X_{i}'\beta\right)\beta
\end{aligned}
$$


```r
require(ggplot2) #for graphs
require(parallel) #good for bootstrapping
require(data.table) #this and the below package are needed to work with data
require(knitr)
set.seed(11) #A seed for me



#holds constants and program parameters
CONST = list(
    NUM_ROWS = 200,
    NUM_COLS = 3,
    NUM_SAMPLES = 2000,
    EPSILON = .Machine$double.eps, #machine precision
    NUM_WORKERS = max(round(detectCores() * .5), 2) #just a heuristic for multi-threading
)

#This is the probit likelihood function
llikelihoodProbit = function(b, Y, X) {
    epsilon = CONST$EPSILON
    argvec = X %*% b

    #avoid numerical issues with logs of small numbers
    pnorms = pmin(pmax(pnorm(argvec), epsilon), 1.0 - epsilon)

    #Use vectorized ifelse
    likes = ifelse(Y, log(pnorms), log(1 - pnorms))

    return(sum(likes))
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
    grads = apply(X, MARGIN = 2, function(x) x * premults)
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
    ll = function(x) - 1.0 * llikelihoodProbit(x, Y, X)
    llgrad = function(x) - 1.0 * llikelihoodProbitGrad(x, Y, X)

    #call the optimizer
    opt = optim(b, ll, gr = llgrad, method = "BFGS", hessian = TRUE)
    if (opt$convergence != 0) print("WARNING! Optimizer did not converge")

    #Efficient matrix inversion
    U = chol(opt$hessian)
    UInv = solve(chol(opt$hessian))
    Sigma = t(UInv) %*% UInv

    #form the info we want into a named list
    prob = list(B = opt$par, llikelihood = opt$value, varB = diag(Sigma), seB = diag(Sigma) ^ 0.5)
    return(prob)
}

#generates a test sample from the asymtotic distribution
testSample = function(R = CONST$NUM_ROWS, C = CONST$NUM_COLS, beta = 1 / (1:C)) {
    #pre-allocate
    X = matrix(rnorm(R * C), nrow = R, ncol = C)

    #create the Y vector
    Y = apply(X, 1, function(x) pnorm(x %*% beta))
    Y = rbinom(R, 1, Y)

    return(list(Y = Y, X = X))
}

#tests the model a single time and prints the results
testProbitModelOnce = function() {

    S = testSample()
    prob = probitModel(S$Y, S$X)
    print(prob)
}

testProbitModelOnce()
```

```
## $B
## [1]  1.294512138  0.429558179  0.302285828 -0.008555375
## 
## $llikelihood
## [1] 85.56511
## 
## $varB
## [1] 0.02820121 0.01650585 0.01472355 0.01217555
## 
## $seB
## [1] 0.1679322 0.1284751 0.1213406 0.1103429
```

* Note that the true betas are 1.0, 0.5 and 0.33
* The results seem reasonably close to the true betas given the small
sample size and binary nature of the dependent variable. 
* From a frequentest standpoint, we cannot reject any of the true betas
using the estimates.

### 12.


```r
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

examineProbitDistributions = function(N = CONST$NUM_SAMPLES) {
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

system.time(examineProbitDistributions())
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-2.png)![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-3.png)

```
##    user  system elapsed 
##    1.27    2.13   12.33
```

* The distribution of the bootstrap and asymtotic error seem reasonably
close.
* Qualitatively, the median of the asymptotic and bootstrap distribution
at least occurs within a reasonable part of the true beta distribution.
* Because the true distribution seems leptokurtic, I generally trust
the bootstrap more in this situation. 


## Problem 4

### 13.

* Because the conditional mean of both specifications is the same, the
MLE estimates of beta are consistent.

### 14.

* As discussed in the slides, the sandwich estimator does not simplify.
    + The estimator is thus distributed $\hat{\theta}\sim N\left(\theta,\,E\left[H^{-1}\right]E\left[S\left(\theta\right)S\left(\theta\right)'\right]E\left[H^{-1}\right]\right)$

* Only under correct specification does $I\left(\theta|X_{i}\right)=E\left[S\left(\theta|X_{i}\right)S\left(\theta|X_{i}\right)'\right]$
* PROOF (Univariate case, borrowing from MLE2\_handout.pdf slides 8
and 9):

    + First write down the square of the score, but using the true probability
distribution to compute the expectation:
$$
\begin{aligned}
E\left[S\left(\theta|Y_{i}\right)^{2}\right]= & \int S\left(\theta|Y_{i}\right)^{2}q\left(Y_{i}|\theta\right)dY_{i}\\
= & \int\left[\frac{\partial lnp\left(Y_{i}|\theta\right)}{\partial\theta}\right]^{2}q\left(Y_{i}|\theta\right)dY_{i}\\
= & \int\frac{1}{p\left(Y_{i}|\theta\right)^{2}}\left[\frac{\partial p\left(Y_{i}|\theta\right)}{\partial\theta}\right]^{2}q\left(Y_{i}|\theta\right)dY_{i}
\end{aligned}
$$
    + Do the same for the Hessian (doesn't quite match up due to typo in
bottom of slide 8):
$$
\begin{aligned}
-E\left[H\left(\theta|Y_{i}\right)\right]= & -\int\frac{\partial^{2}lnp\left(Y_{i}|\theta\right)}{\partial\theta^{2}}q\left(Y_{i}|\theta\right)dY_{i}\\
\frac{\partial^{2}lnp\left(Y_{i}|\theta\right)}{\partial\theta^{2}}= & \frac{\partial}{\partial\theta}\left[\frac{1}{p\left(Y_{i}|\theta\right)}\frac{\partial p\left(Y_{i}|\theta\right)}{\partial\theta}\right]\\
= & \frac{-1}{p\left(Y_{i}|\theta\right)^{2}}\left(\frac{\partial p\left(Y_{i}|\theta\right)}{\partial\theta}\right)^{2}+\frac{1}{p\left(Y_{i}|\theta\right)}\frac{\partial^{2}p\left(Y_{i}|\theta\right)}{\partial^{2}\theta}\\
-E\left[H\left(\theta|Y_{i}\right)\right]= & -\int\left[\frac{-q\left(Y_{i}|\theta\right)}{p\left(Y_{i}|\theta\right)^{2}}\left(\frac{\partial p\left(Y_{i}|\theta\right)}{\partial\theta}\right)^{2}+\frac{q\left(Y_{i}|\theta\right)}{p\left(Y_{i}|\theta\right)}\frac{\partial^{2}p\left(Y_{i}|\theta\right)}{\partial^{2}\theta}\right]dY_{i}\\
= & E\left[S\left(\theta|Y_{i}\right)^{2}\right]-\int\frac{q\left(Y_{i}|\theta\right)}{p\left(Y_{i}|\theta\right)}\frac{\partial^{2}p\left(Y_{i}|\theta\right)}{\partial^{2}\theta}dY_{i}\\
= & E\left[S\left(\theta|Y_{i}\right)^{2}\right]-\int\frac{q\left(Y_{i}|\theta\right)}{p\left(Y_{i}|\theta\right)}\frac{\partial^{2}p\left(Y_{i}|\theta\right)}{\partial\theta^{2}}dY_{i}\checkmark
\end{aligned}
$$
    + Note if q=p we achieve the desired simplification.


### 15.

* {[}CODE HERE{]}

### 16.

* As expected, the standard errors are higher using the more robust
technique.

### 17.

* The true standard errors seem reasonably close to the standard errors
from the robust estimation technique. 
* They are substantially more than the standard errors computed assuming
correct specification.
