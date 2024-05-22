u <- rnorm(100)
Mybounds <- as.numeric(quantile(u, probs = c(0.05, 0.95)))
fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = 50) #getting a sequence of evenly spaced numbers between the quantile 5 and the quantile 95
fhat_points <- stats::filter(fhat_int, c(0.5, 0.5))[-length(fhat_int)] #getting the midpoint on every subinterval of the sequence defined above
comdiff <- fhat_points[2] - fhat_points[1] #distance between all points (they are evenly spaced)
#adding one point at the beginning and one at the end (still evenly distributed)
fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  


fit1 <- fcar.fit(Y1, X, u, epanechnikov, h)

plot(fit1$functional_points, f11(fit1$functional_points), 
    ylim = c(min(f11(fit1$functional_points), fit1$coeffs[1,]), max(f11(fit1$functional_points), fit1$coeffs[1,])))
lines(fit1$functional_points, fit1$coeffs[1, ])

plot(fit1$functional_points, f12(fit1$functional_points), 
    ylim = c(min(f12(fit1$functional_points), fit1$coeffs[1,]), max(f12(fit1$functional_points), fit1$coeffs[1,])))
lines(fit1$functional_points, fit1$coeffs[2, ])

fit2 <- fcar.fit(Y2, X, u, epanechnikov, h)

plot(fit2$functional_points, f21(fit2$functional_points),  ylim = c(-.5, .5))
lines(fit2$functional_points, fit2$coeffs[1, ])

plot(fit2$functional_points, f22(fit2$functional_points), ylim = c(.3 - .5, .3 + .5))
lines(fit2$functional_points, fit2$coeffs[2, ])



lol <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 100)
min(lol$ref.distribution1)
abline(v = lol$null.stat1)

Y <- mydata$y
X.mul <- Y[-nrow(Y), ]
Y <- Y[-1, ]
U <- mydata$ref[-1]
uo <- U[1]
X.tilde <- cbind(X.mul, X.mul * (U - uo))
h.mul <- (max(u) - min(u)) * .2
W <- diag(epanechnikov(uo - U, h))
theta1 <- solve(t(X.tilde) %*% W %*% X.tilde) %*% t(X.tilde) %*% W %*% Y[, 1]
theta1

solve(t(X.tilde) %*% W %*% X.tilde) %*% t(X.tilde) %*% W %*% Y

fit.tmp <- fvcar.fit(Y, X.mul, U, epanechnikov, h.mul)

plot(fit.tmp$functional_points, f12(fit.tmp$functional_points))
lines(fit.tmp$functional_points, fit.tmp$coeffs[1,2,])

ts.plot(Y[, 1])
lines(fit.tmp$yhat[,1], col = "red")
mse(Y[, 1], fit.tmp$yhat[, 1])

k.tmp <- function(u, h) {
    exp(-.5 * (u/h)^2) / (h * sqrt(2 * pi))
}

sizas <- fvcar.fit(mydata[[1]], mydata[[2]], k.tmp, 1, 2)
plot(sizas$functional_points, f21(sizas$functional_points), ylim = c(-.5, .5))
plot(sizas$functional_points, sizas$coeffs[2, 1,])
