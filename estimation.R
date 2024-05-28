rm(list = ls())
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)
source("auxiliary_functions.R")

niters <- 100
numcores <- detectCores()
##### Functional coefficients - EXPAR 

#Both GC

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){0.6*exp(-(0.30 + theta)*x^2)} 
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}

set.seed(1234)
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)
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

Y1 <- Y1[-(1:2)]
Y2 <- 
fit <- fcar.fit(Y1, X, u, epanechnikov, h)
ts.plot(Y1)
lines(fit$yhat, col = "red")

plot(fit$functional_points, f22(fit$functional_points))
lines(fit$functional_points, fit$coeffs)

