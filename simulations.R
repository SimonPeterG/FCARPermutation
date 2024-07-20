rm(list = ls())
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)
source("auxiliary_functions.R")

niters <- 2
numcores <- detectCores()
##### Functional coefficients - EXPAR 

# only Y2 GC Y1

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){rep(0, length(x))} 
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}


set.seed(1234)
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


registerDoParallel(numcores)

set.seed(i)
tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)

expar1 <- foreach(i = 1:1) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}

saveRDS(expar1, "results/expar1.Rds")

#There is no GC

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){rep(0, length(x))} 
f21 <- function(x,theta=0){rep(0, length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}

set.seed(7)
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)

Y1 <- mydata[[1]][, 1]
Y2 <- mydata[[1]][, 2]

u <- mydata[[2]][1:(length(Y1) - 2)]
h <- (max(u) - min(u)) * .1
X <- matrix(0, ncol = 2, nrow = length(Y1) - 2)
for (i in 1:nrow(X)) {
    X[i, ] <- c(Y1[i+1], Y2[i+1])
}

Y1 <- matrix(Y1[-(1:2)], ncol = 1)
Y2 <- matrix(Y2[-(1:2)], ncol = 1)

expar2 <- foreach(i = 1:niters) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}

saveRDS(expar2, "results/expar2.Rds")

#Both GC

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){0.6*exp(-(0.30 + theta)*x^2)} 
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}

set.seed(1234)
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

expar3 <- foreach(i = 1:niters) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}

saveRDS(expar3, "results/expar3.Rds")

##### Functional coefficients - Logistic

# Only Y1 GC Y2

f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0.2,length(x)) + theta}
f21 <- function(x,theta=0){rep(0, length(x)) + theta}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}

set.seed(4321)
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 2,
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

logis1 <- foreach(i = 1:niters) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}
saveRDS(logis1, "results/logis1.Rds")

# There is no GC

f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0.2,length(x)) + theta}
f21 <- function(x,theta=0){-0.9*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))+0.5}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}


set.seed(951413)
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 2,
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

logis2 <- foreach(i = 1:niters) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}

saveRDS(logis2, "results/logis2.Rds")


# Both GC

f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0.2,length(x)) + theta}
f21 <- function(x,theta=0){-0.9*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))+0.5}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}

set.seed(1212)
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 2,
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

logis3 <- foreach(i = 1:niters) %dopar% {
    set.seed(i)
    tmp <- permutation.test(Y1, Y2, u, X, epanechnikov, h, P = 500)
    print(paste0("Iteration ", i, " completed"))
    tmp
}

saveRDS(logis3 , "results/logis3.Rds")