rm(list = ls())
source("auxiliary_functions.R")

##### Functional coefficients - EXPAR (only Y1 GC Y2)
f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){rep(0, length(x))} #0.6*exp(-(0.30 + theta)*x^2)
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}


mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)

Y1 <- mydata[[1]][, 1]
Y2 <- mydata[[1]][, 2]

u <- mydata[[2]][1:(length(Y1) - 2)]
h <- (max(u) - min(u)) * .2
X <- matrix(0, ncol = 2, nrow = length(Y1) - 2)
for (i in 1:nrow(X)) {
    X[i, ] <- c(Y1[i+1], Y2[i+1])
}

Y1 <- matrix(Y1[-(1:2)], ncol = 1)
Y2 <- matrix(Y2[-(1:2)], ncol = 1)

set.seed(314159)
permutation.test(Y1, Y2, u, X, epanechnikov, h, c(2, NA), P = 1000)
