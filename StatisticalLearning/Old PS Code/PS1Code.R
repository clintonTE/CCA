require(data.table)
set.seed(11) #A seed for me


## ---- p1q1
p1q1 = function(N = 100) {
    
    X = matrix(c(rep(1, N), rnorm(N), rnorm(N)), nrow = N, ncol = 3)
    ep = rnorm(N)
    B = c(1, 2, 3)
    Y = X %*% B + ep

    Xqr = qr(X)
    Q = qr.Q(Xqr)
    R = qr.R(Xqr)
    Best = solve(R) %*% t(Q) %*% Y
    cat("Our estimate: ", Best, "\n")
    cat("lm estimate:\n")
    print(coef(lm(Y ~ X + 0)))
}

## ---- Execute
p1q1(100)

