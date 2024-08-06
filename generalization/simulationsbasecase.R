rm(list = ls())
setwd("./generalization")
source("./auxiliary_functions.R")

niters <- 100
numcores <- detectCores()
obs.per.block <- c(10, 50, 100)
##### Functional coefficients - EXPAR 

# only Y1 GC Y2
expar1 <- function(x) {
  matrix(c(
    -.3, 0,
    -.2, -.4 * exp(-.45 * x ^2)
  ), byrow = T, ncol = 2)
}

# set.seed(1234)
# ex1 <- simulate.far(256, 2, 0, 1, 1, list(expar1))

# per.test1 <- permutation.test(ex1$ts, ex1$ref, gaussian, .2, ex1$p, ex1$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test1$null.stat1
# range(per.test1$ref.distribution1)
# mean(per.test1$ref.distribution1 >= per.test1$null.stat1)
# mean(per.test1$ref.distribution1 >= per.test1$null.stat1) < 0.05

# per.test1$null.stat2
# range(per.test1$ref.distribution2)
# mean(per.test1$ref.distribution2 >= per.test1$null.stat2)
# mean(per.test1$ref.distribution2 >= per.test1$null.stat2) < 0.05

#simulation of first scenario

expar1b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar1)), P = 1000)

saveRDS(expar1b1, paste0("../results/basecasegen/expar1_", obs.per.block[1], ".Rds"))

expar1b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar1)), P = 1000)
saveRDS(expar1b2, paste0("../results/basecasegen/expar1_", obs.per.block[2], ".Rds"))

expar1b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar1)), P = 1000)
saveRDS(expar1b3, paste0("../results/basecasegen/expar1_", obs.per.block[3], ".Rds"))

#There is no GC

expar2 <- function(x) {
  matrix(c(
    -.3, 0,
    0, -.4 * exp(-.45 * x ^2)
  ), byrow = T, ncol = 2)
}

# set.seed(314)
# ex2 <- simulate.far(256, 2, 0, 1, 1, list(expar2))

# per.test2 <- permutation.test(ex2$ts, ex2$ref, gaussian, .2, ex2$p, ex2$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test2$null.stat1
# range(per.test2$ref.distribution1)
# mean(per.test2$ref.distribution1 >= per.test2$null.stat1)
# mean(per.test2$ref.distribution1 >= per.test2$null.stat1) < 0.05

# per.test2$null.stat2
# range(per.test2$ref.distribution2)
# mean(per.test2$ref.distribution2 >= per.test2$null.stat2)
# mean(per.test2$ref.distribution2 >= per.test2$null.stat2) < 0.05

#simulation of the second scenario
expar2b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar2)), P = 1000)

saveRDS(expar2b1, paste0("../results/basecasegen/expar2_", obs.per.block[1], ".Rds"))

expar2b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar2)), P = 1000)
saveRDS(expar2b2, paste0("../results/basecasegen/expar2_", obs.per.block[2], ".Rds"))

expar2b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar2)), P = 1000)
saveRDS(expar2b3, paste0("../results/basecasegen/expar2_", obs.per.block[3], ".Rds"))

#Both GC
expar3 <- function(x) {
  matrix(c(
    -.3, 0.6*exp(-0.30*x^2),
    -.2, -0.4*exp(-0.45*x^2)
  ), byrow = T, ncol = 2)
}


# set.seed(69)
# ex3 <- simulate.far(256, 2, 0, 1, 1, list(expar3))

# per.test3 <- permutation.test(ex3$ts, ex3$ref, gaussian, .2, ex3$p, ex3$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test3$null.stat1
# range(per.test3$ref.distribution1)
# mean(per.test3$ref.distribution1 >= per.test3$null.stat1)
# mean(per.test3$ref.distribution1 >= per.test3$null.stat1) < 0.05

# per.test3$null.stat2
# range(per.test3$ref.distribution2)
# mean(per.test3$ref.distribution2 >= per.test3$null.stat2)
# mean(per.test3$ref.distribution2 >= per.test3$null.stat2) < 0.05

#simulation of the third scenario
expar3b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar3)), P = 1000)

saveRDS(expar3b1, paste0("../results/basecasegen/expar3_", obs.per.block[1], ".Rds"))

expar3b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar3)), P = 1000)
saveRDS(expar3b2, paste0("../results/basecasegen/expar3_", obs.per.block[2], ".Rds"))

expar3b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar3)), P = 1000)
saveRDS(expar3b3, paste0("../results/basecasegen/expar3_", obs.per.block[3], ".Rds"))

##### Functional coefficients - Logistic

# only Y1 GC Y2
logis1 <- function(x) {
  matrix(c(
    0.8*(exp(5*x)/(1+exp(5*x)))-0.3, 0,
    -0.9*(exp(5*x)/(1+exp(5*x)))+0.5, .3
  ), byrow = T, ncol = 2)
}

# set.seed(7)
# ex4 <- simulate.far(256, 2, 0, 1, 1, list(logis1))

# per.test4 <- permutation.test(ex4$ts, ex4$ref, gaussian, .2, ex4$p, ex4$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test4$null.stat1
# range(per.test4$ref.distribution1)
# mean(per.test4$ref.distribution1 >= per.test4$null.stat1)
# mean(per.test4$ref.distribution1 >= per.test4$null.stat1) < 0.05

# per.test4$null.stat2
# range(per.test4$ref.distribution2)
# mean(per.test4$ref.distribution2 >= per.test4$null.stat2)
# mean(per.test4$ref.distribution2 >= per.test4$null.stat2) < 0.05

#simulation of first scenario

logis1b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis1)), P = 1000)

saveRDS(logis1b1, paste0("../results/basecasegen/logis1_", obs.per.block[1], ".Rds"))

logis1b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis1)), P = 1000)
saveRDS(logis1b2, paste0("../results/basecasegen/logis1_", obs.per.block[2], ".Rds"))

logis1b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis1)), P = 1000)
saveRDS(logis1b3, paste0("../results/basecasegen/logis1_", obs.per.block[3], ".Rds"))

#There is no GC

logis2 <- function(x) {
  matrix(c(
    0.8*(exp(5*x)/(1+exp(5*x)))-0.3, 0,
    0, .3
  ), byrow = T, ncol = 2)
}

# set.seed(1221)
# ex5 <- simulate.far(256, 2, 0, 1, 1, list(logis2))

# per.test5 <- permutation.test(ex5$ts, ex5$ref, gaussian, .2, ex5$p, ex5$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test5$null.stat1
# range(per.test5$ref.distribution1)
# mean(per.test5$ref.distribution1 >= per.test5$null.stat1)
# mean(per.test5$ref.distribution1 >= per.test5$null.stat1) < 0.05

# per.test5$null.stat2
# range(per.test5$ref.distribution2)
# mean(per.test5$ref.distribution2 >= per.test5$null.stat2)
# mean(per.test5$ref.distribution2 >= per.test5$null.stat2) < 0.05

#simulation of the second scenario
logis2b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis2)), P = 1000)

saveRDS(logis2b1, paste0("../results/basecasegen/logis2_", obs.per.block[1], ".Rds"))

logis2b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis2)), P = 1000)
saveRDS(logis2b2, paste0("../results/basecasegen/logis2_", obs.per.block[2], ".Rds"))

logis2b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis2)), P = 1000)
saveRDS(logis2b3, paste0("../results/basecasegen/logis2_", obs.per.block[3], ".Rds"))

#Both GC
logis3 <- function(x) {
  matrix(c(
    0.8*(exp(5*x)/(1+exp(5*x)))-0.3, 0.2,
    -0.9*(exp(5*x)/(1+exp(5*x)))+0.5, .3
  ), byrow = T, ncol = 2)
}


# set.seed(69)
# ex6 <- simulate.far(256, 2, 0, 1, 1, list(logis3))

# per.test6 <- permutation.test(ex6$ts, ex6$ref, gaussian, .2, ex6$p, ex6$p.u, 
#                               function(x) block_permute1(x, 10), P = 1000)

# per.test6$null.stat1
# range(per.test6$ref.distribution1)
# mean(per.test6$ref.distribution1 >= per.test6$null.stat1)
# mean(per.test6$ref.distribution1 >= per.test6$null.stat1) < 0.05

# per.test6$null.stat2
# range(per.test6$ref.distribution2)
# mean(per.test6$ref.distribution2 >= per.test6$null.stat2)
# mean(per.test6$ref.distribution2 >= per.test6$null.stat2) < 0.05

#simulation of the third scenario
logis3b1 <- gc.test(niters, numcores, 256, obs.per.block[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis3)), P = 1000)

saveRDS(logis3b1, paste0("../results/basecasegen/logis3_", obs.per.block[1], ".Rds"))

logis3b2 <- gc.test(niters, numcores, 256, obs.per.block[2], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis3)), P = 1000)
saveRDS(logis3b2, paste0("../results/basecasegen/logis3_", obs.per.block[2], ".Rds"))

logis3b3 <- gc.test(niters, numcores, 256, obs.per.block[3], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis3)), P = 1000)
saveRDS(logis3b3, paste0("../results/basecasegen/logis3_", obs.per.block[3], ".Rds"))

setwd("./..")
rm(list = ls())
