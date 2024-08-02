generate.X <- function(Y, p) {
    Y1 <- Y[, 1]; Y2 <- Y[, 2]
    n <- length(Y1)
    X <- matrix(0, nrow = n - p, ncol = 2 * p)
    even <- seq(2, 2 * p, by = 2)
    odd <- seq(1, 2 * p, by = 2)
    for (i in 1:nrow(X)) {
        X[i, odd] <- Y1[(i + p - 1):i]
        X[i, even] <- Y2[(i + p - 1):i]
    }
    X
}

simulate.far <- function(n, d, u, p, p.u, mixing.matrices, sigma.noise = 1, burn = 1000) {
    p.star <- max(p, p.u)
    Y <- matrix(0, ncol = burn + n, nrow = d)
    ref <- numeric(burn + n - p.star)
    Y[, 1:p.star] <- rnorm(p.star * d)
    for (i in (p.star + 1):(burn + n)) {
        for (j in 1:length(mixing.matrices)) {
            ref[i - p.star] <- ref.point <- ifelse(u == 0, rnorm(1), Y[u, i - p.u]) 
            Y[, i] <- Y[, i] + mixing.matrices[[j]](ref.point) %*% Y[, i - j] + matrix(rnorm(d, sd = sigma.noise), ncol = 1)
        }
    }
    list(ts = t(Y[, -(1:burn)]), ref = ref[-(1:burn)], p.star = p.star)
} 

m1 <- function(u) {
    matrix(c(
        -.3, 0,
        0, -0.4*exp(-0.45 * u^2)
    ), ncol = 2, byrow = T)
}

m2 <- function(u) {
    matrix(c(
       0.8 * ((exp(5 * u)) / (1 + exp(5 * u))) - 0.3, 0.2,
       -0.9 * ((exp(5 * u)) / (1 + exp(5* u))) + 0.5, 0.3
    ), ncol = 2, byrow = T)
}

ex <- simulate.far(256, 2, 0, 2, 1, list(m1, m2))
dim(ex$ts)
length(ex$ref)
ex$p.star

X.ex <- generate.X(ex$ts, ex$p.star)
u.ex <- ex$ref[1:(256 - 2)]
h.ex <- (max(u.ex) - min(u.ex)) * .2
y1.ex <- ex$ts[-(1:2), 1]
y2.ex <- ex$ts[-(1:2), 2]
f1 <- fcar.fit(y1.ex, X.ex, u.ex, gaussian, h.ex, ex$p.star)
f2 <- fcar.fit(y2.ex, X.ex, u.ex, gaussian, h.ex, ex$p.star)

ts.plot(y1.ex)
lines(f1$yhat, col = "red")

layout(matrix(1:4, ncol = 2, byrow = T))
plot(f1$coeffs[1, ], type = "l", ylim = c(-1, 1))
plot(f1$coeffs[2, ], type = "l", ylim = c(-1, 1))
plot(f1$coeffs[3, ], type = "l", ylim = c(-1, 1))
plot(f1$coeffs[4, ], type = "l", ylim = c(-1, 1))

ts.plot(y2.ex)
lines(f2$yhat, col = "red")

layout(matrix(1:4, ncol = 2, byrow = T))
plot(f2$coeffs[1, ], type = "l", ylim = c(-1, 1))
plot(f2$coeffs[2, ], type = "l", ylim = c(-1, 1))
plot(f2$coeffs[3, ], type = "l", ylim = c(-1, 1))
plot(f2$coeffs[4, ], type = "l", ylim = c(-1, 1))


mse(y1.ex, f1$yhat)
mse(y2.ex, f2$yhat)

gc.ex <- permutation.test(y1.ex, y2.ex, u.ex, X.ex, gaussian, h.ex, 10, P = 1000)
