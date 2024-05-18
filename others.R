rm(list = ls())

a1 <- function(u) {
    .138 + (.316 + .982 * u) * exp(-3.89 * u^2)    
}

a2 <- function(u) {
    -.437 - (.659 + 1.26 * u) * exp(-3.89 * u^2)    
}

sim.expar <- function(n, a1, a2, burn = 1000) {
    y <- rnorm(n + burn)
    e <- rnorm(n + burn, sd = .2)
    for (i in 3:(n + burn)) {
        y[i] <- a1(y[i-1]) * y[i-1] + a2(y[i-1]) * y[i-2] + e[i]
    }
    y[-(1:burn)]
}

set.seed(31415921)
y <- sim.expar(400, a1, a2)
X <- matrix(0, ncol = 2, nrow = length(y) - 2)
for (i in 1:nrow(X)) {
    X[i, ] <- y[(i+1):i]
}
u <- X[, 1]
y <- matrix(y[-(1:2)], ncol = 1)

fit.list <- fcar.fit(y, X, u, epanechnikov, .41)

par(mfrow = c(1, 1))
ts.plot(y)
lines(fit.list[[1]], col = "red")


par(mfrow = c(1, 2))
a1.orig <- a1(fit.list[[3]])
a1.hat <- fit.list[[2]][1, ]

plot(fit.list[[3]], a1.orig, ylim = c(min(a1.orig, a1.hat), max(a1.orig, a1.hat)))
lines(fit.list[[3]], a1.hat, col = "red")

a2.orig <- a2(fit.list[[3]])
a2.hat <- fit.list[[2]][2, ]

plot(fit.list[[3]], a2.orig, ylim = c(min(a2.orig, a2.hat), max(a2.orig, a2.hat)))
lines(fit.list[[3]], a2.hat, col = "red")
