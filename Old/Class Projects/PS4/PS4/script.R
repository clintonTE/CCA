lambda = 2;
mu = 1;

unxRV = runif(10 ^ 5);
unyRV = runif(10 ^ 5);
xRV = sapply(unxRV, function(x) qexp(x, lambda));
yRV = mapply(function(x, y) qnorm(x, mu, y), unyRV, xRV);
#plot(unyRV, yRV);
print(mean(yRV));
print(var(yRV));
