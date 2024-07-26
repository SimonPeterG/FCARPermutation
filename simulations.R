rm(list = ls())
library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)
source("auxiliary_functions.R")

niters <- 100
numcores <- detectCores()
block.sizes <- c(10, 50, 100, 256 - 2)
##### Functional coefficients - EXPAR 

# only Y2 GC Y1

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){rep(0, length(x))} 
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}


set.seed(1)
mydata <- gendata(Tlength = 256,d = 2,Y_d = 0,
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


ex1 <- permutation.test(Y1, Y2, u, X, gaussian, h, 10, P = 500)

ex1$null.stat1
range(ex1$ref.distribution1)
hist(ex1$ref.distribution1)
abline(v = ex1$null.stat1, col = "red")

ex1$null.stat2
range(ex1$ref.distribution2)
hist(ex1$ref.distribution2, xlim = c(-2, 55))
abline(v = ex1$null.stat2, col = "red")


#simulation of first scenario

expar1b1 <- gc.test(block.sizes[1], niters, numcores, 256, 1000)
saveRDS(expar1b1, paste0("results/expar1_", block.sizes[1], ".Rds"))

expar1b2 <- gc.test(block.sizes[2], niters, numcores, 256, 1000)
saveRDS(expar1b1, paste0("results/expar1_", block.sizes[2], ".Rds"))

expar1b3 <- gc.test(block.sizes[3], niters, numcores, 256, 1000)
saveRDS(expar1b1, paste0("results/expar1_", block.sizes[3], ".Rds"))

expar1b4 <- gc.test(block.sizes[4] - 2, niters, numcores, 256, 1000)
saveRDS(expar1b1, paste0("results/expar1_", block.sizes[4], ".Rds"))

#There is no GC

f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){rep(0, length(x))} 
f21 <- function(x,theta=0){rep(0, length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}

set.seed(1)
mydata <- gendata(Tlength = 256,d = 2,Y_d = 0,
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

set.seed(10)
ex2 <- permutation.test(Y1, Y2, u, X, epanechnikov, h, block.sizes[1], P = 1000)

ex2$null.stat1
range(ex2$ref.distribution1)
hist(ex2$ref.distribution1)
abline(v = ex2$null.stat1, col = "red")

ex2$null.stat2
range(ex2$ref.distribution2)
hist(ex2$ref.distribution2)
abline(v = ex2$null.stat2, col = "red")

#simulation of the second scenario
expar2b1 <- gc.test(block.sizes[1], niters, numcores, 256, 1000)
saveRDS(expar2b1, paste0("results/expar2_", block.sizes[1], ".Rds"))

expar2b2 <- gc.test(block.sizes[2], niters, numcores, 256, 1000)
saveRDS(expar2b2, paste0("results/expar2_", block.sizes[2], ".Rds"))

expar2b3 <- gc.test(block.sizes[3], niters, numcores, 256, 1000)
saveRDS(expar2b3, paste0("results/expar2_", block.sizes[3], ".Rds"))

expar2b4 <- gc.test(block.sizes[4] - 2, niters, numcores, 256, 1000)
saveRDS(expar2b4, paste0("results/expar2_", block.sizes[4], ".Rds"))

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
mydata <- gendata(Tlength = 256,d = 2,Y_d = 2,
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

set.seed(10)
ex4 <- permutation.test(Y1, Y2, u, X, gaussian, h, block.sizes[1], P = 1000)

ex4$null.stat1
range(ex4$ref.distribution1)
hist(ex4$ref.distribution1)
abline(v = ex4$null.stat1, col = "red")

ex4$null.stat2
range(ex4$ref.distribution2)
hist(ex4$ref.distribution2)
abline(v = ex4$null.stat2, col = "red")

#simulation of first scenario

logis1b1 <- gc.test(block.sizes[1], niters, numcores, 256, 1000)
saveRDS(logis1b1, paste0("results/logis1_", block.sizes[1], ".Rds"))

logis1b2 <- gc.test(block.sizes[2], niters, numcores, 256, 1000)
saveRDS(logis1b1, paste0("results/logis1_", block.sizes[2], ".Rds"))

logis1b3 <- gc.test(block.sizes[3], niters, numcores, 256, 1000)
saveRDS(logis1b1, paste0("results/logis1_", block.sizes[3], ".Rds"))

logis1b4 <- gc.test(block.sizes[4] - 2, niters, numcores, 256, 1000)
saveRDS(logis1b1, paste0("results/logis1_", block.sizes[4], ".Rds"))

# There is no GC

f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0,length(x)) + theta}
f21 <- function(x,theta=0){rep(0,length(x)) + theta}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}


set.seed(951413)
mydata <- gendata(Tlength = 256,d = 2,Y_d = 2,
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

set.seed(1771)
ex5 <- permutation.test(Y1, Y2, u, X, gaussian, h, block.sizes[1], P = 1000)

ex5$null.stat1
range(ex5$ref.distribution1)
hist(ex5$ref.distribution1)
abline(v = ex5$null.stat1, col = "red")

ex5$null.stat2
range(ex5$ref.distribution2)
hist(ex5$ref.distribution2)
abline(v = ex5$null.stat2, col = "red")

#simulation of second scenario

logis2b1 <- gc.test(block.sizes[1], niters, numcores, 256, 1000)
saveRDS(logis2b1, paste0("results/logis2_", block.sizes[1], ".Rds"))

logis2b2 <- gc.test(block.sizes[2], niters, numcores, 256, 1000)
saveRDS(logis2b1, paste0("results/logis2_", block.sizes[2], ".Rds"))

logis2b3 <- gc.test(block.sizes[3], niters, numcores, 256, 1000)
saveRDS(logis2b1, paste0("results/logis2_", block.sizes[3], ".Rds"))

logis2b4 <- gc.test(block.sizes[4] - 2, niters, numcores, 256, 1000)
saveRDS(logis2b1, paste0("results/logis2_", block.sizes[4], ".Rds"))


# Both GC

f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0.2,length(x)) + theta}
f21 <- function(x,theta=0){-0.9*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))+0.5}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}

set.seed(1212)
mydata <- gendata(Tlength = 256,d = 2,Y_d = 2,
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

set.seed(7117)
ex6 <- permutation.test(Y1, Y2, u, X, gaussian, h, block.sizes[1], P = 1000)

ex6$null.stat1
range(ex6$ref.distribution1)
hist(ex6$ref.distribution1)
abline(v = ex6$null.stat1, col = "red")

ex6$null.stat2
range(ex6$ref.distribution2)
hist(ex6$ref.distribution2)
abline(v = ex6$null.stat2, col = "red")

#simulation of third scenario

logis3b1 <- gc.test(block.sizes[1], niters, numcores, 256, 1000)
saveRDS(logis3b1, paste0("results/logis3_", block.sizes[1], ".Rds"))

logis3b2 <- gc.test(block.sizes[2], niters, numcores, 256, 1000)
saveRDS(logis3b1, paste0("results/logis3_", block.sizes[2], ".Rds"))

logis3b3 <- gc.test(block.sizes[3], niters, numcores, 256, 1000)
saveRDS(logis3b1, paste0("results/logis3_", block.sizes[3], ".Rds"))

logis3b4 <- gc.test(block.sizes[4] - 2, niters, numcores, 256, 1000)
saveRDS(logis3b1, paste0("results/logis3_", block.sizes[4], ".Rds"))