source("./generalization/auxiliary_functions.R")

simulate.far <- function(n, d, u, p, p.u, mixing.matrices, sigma.noise = 1, burn = 1000) {
    p.star <- max(p, p.u)
    Y <- matrix(0, ncol = burn + n, nrow = d)
    ref <- numeric(burn + n)
    Y[, 1:p.star] <- rnorm(p.star * d)
    for (i in (p.star + 1):(burn + n)) {
        for (j in 1:length(mixing.matrices)) {
            ref[i - p.star] <- ref.point <- ifelse(u == 0, rnorm(1), Y[u, i - p.u]) 
            Y[, i] <- Y[, i] + mixing.matrices[[j]](ref.point) %*% Y[, i - j] + matrix(rnorm(d, sd = sigma.noise), ncol = 1)
        }
    }
    
    list(ts = t(Y[, -(1:burn)]), ref = ref[-(1:burn)], p = p, p.u = p.u, p.star = p.star)
} 

m1 <- function(u) {
    matrix(c(
        -.3, 0,
        1, -0.4*exp(-0.45 * u^2)
    ), ncol = 2, byrow = T)
}

m2 <- function(u) {
    matrix(c(
       0.8 * ((exp(5 * u)) / (1 + exp(5 * u))) - 0.3, 0.2,
       -0.9 * ((exp(5 * u)) / (1 + exp(5* u))) + 0.5, 0.3
    ), ncol = 2, byrow = T)
}

m0 <- function(u) {
    matrix(0, ncol = 2, nrow = 2)
}

#modelo 1
ex <- simulate.far(256, 2, 1, 1, 1, list(m1))

per.test <- permutation.test(ex$ts, ex$ref, gaussian, .2, ex$p, ex$p.u, permute = function(x) block_permute1(x, 10))

per.test$null.stat1
range(per.test$ref.distribution1)
mean(per.test$ref.distribution1 >= per.test$null.stat1)

per.test$null.stat2
range(per.test$ref.distribution2)
mean(per.test$ref.distribution2 >= per.test$null.stat2)

gc.test1 <- gc.test(2, 2, 256, 10, function(n) simulate.far(n, 2, 1, 1, 1, list(m1)), P = 500)

gc.test1[[1]]$null.stat1
range(gc.test1[[1]]$ref.distribution1)
mean(gc.test1[[1]]$ref.distribution1 >= gc.test1[[1]]$null.stat1)

gc.test1[[1]]$null.stat2
range(gc.test1[[1]]$ref.distribution2)
mean(gc.test1[[1]]$ref.distribution2 >= gc.test1[[1]]$null.stat2)

gc.test1[[2]]$null.stat1
range(gc.test1[[2]]$ref.distribution1)
mean(gc.test1[[2]]$ref.distribution1 >= gc.test1[[2]]$null.stat1)

gc.test1[[2]]$null.stat2
range(gc.test1[[2]]$ref.distribution2)
mean(gc.test1[[2]]$ref.distribution2 >= gc.test1[[2]]$null.stat2)

#mmodel 2
ex <- simulate.far(256, 2, 1, 1, 1, list(m2))

per.test <- permutation.test(ex$ts, ex$ref, gaussian, .2, ex$p, ex$p.u, permute = function(x) block_permute1(x, 10), P = 1000)

per.test$null.stat1
range(per.test$ref.distribution1)
mean(per.test$ref.distribution1 >= per.test$null.stat1)

per.test$null.stat2
range(per.test$ref.distribution2)
mean(per.test$ref.distribution2 >= per.test$null.stat2)

gc.test2 <- gc.test(2, 1, 256, 10, function(n) simulate.far(n, 2, 1, 1, 1, list(m2)), P = 500)

gc.test2[[1]]$null.stat1
range(gc.test2[[1]]$ref.distribution1)
mean(gc.test2[[1]]$ref.distribution1 >= gc.test2[[1]]$null.stat1)

gc.test2[[1]]$null.stat2
range(gc.test2[[1]]$ref.distribution2)
mean(gc.test2[[1]]$ref.distribution2 >= gc.test2[[1]]$null.stat2)

gc.test2[[2]]$null.stat1
range(gc.test2[[2]]$ref.distribution1)
mean(gc.test2[[2]]$ref.distribution1 >= gc.test2[[2]]$null.stat1)

gc.test2[[2]]$null.stat2
range(gc.test2[[2]]$ref.distribution2)
mean(gc.test2[[2]]$ref.distribution2 >= gc.test2[[2]]$null.stat2)

#m3
ex <- simulate.far(256, 2, 0, 2, 3, list(m1, m2))

per.test <- permutation.test(ex$ts, ex$ref, gaussian, .2, ex$p, ex$p.u, 
                             permute = function(x) block_permute1(x, 10), P = 1000)

per.test$null.stat1
range(per.test$ref.distribution1)
mean(per.test$ref.distribution1 >= per.test$null.stat1)

per.test$null.stat2
range(per.test$ref.distribution2)
mean(per.test$ref.distribution2 >= per.test$null.stat2)

#expar 1

expar <- function(x) {
  matrix(c(
    -0.3, 0,
    -0.2, -0.4*exp(-(0.45)*x^2)
    
  ), byrow = T, ncol = 2)
}

set.seed(1234)
ex <- simulate.far(256, 2, 0, 1, 1, list(expar))
layout(matrix(1:2, nrow = 1))
ts.plot(ex$ts[, 1])
ts.plot(ex$ts[, 2])


per.test <- permutation.test(ex$ts, ex$ref, gaussian, .2, ex$p, ex$p.u, 
                             permute = function(x) block_permute1(x, 10), P = 1000)

per.test$null.stat1
range(per.test$ref.distribution1)
mean(per.test$ref.distribution1 >= per.test$null.stat1) < 0.05

per.test$null.stat2
range(per.test$ref.distribution2)
mean(per.test$ref.distribution2 >= per.test$null.stat2) < 0.05

