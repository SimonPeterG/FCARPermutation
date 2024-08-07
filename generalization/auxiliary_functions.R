library(parallel)
library(foreach)
library(doParallel)

simulate.far <- function(n, d, u, p, p.u, mixing.matrices, sigma.noise = 1, burn = 1000) {
  p.star <- max(p, p.u)
  Y <- matrix(0, ncol = burn + n, nrow = d)
  ref <- numeric(burn + n)
  Y[, 1:p.star] <- rnorm(p.star * d)
  for (i in (p.star + 1):(burn + n)) {
    for (j in 1:length(mixing.matrices)) {
      ref[i - p.star] <- ref.point <- ifelse(u == 0, rnorm(1), Y[u, i - p.u]) 
      Y[, i] <- Y[, i] + mixing.matrices[[j]](ref.point) %*% Y[, i - j]
    }
    Y[, i] <- Y[, i] + matrix(rnorm(d, sd = sigma.noise), ncol = 1)
  }
  
  list(ts = t(Y[, -(1:burn)]), ref = ref[-(1:burn)], p = p, p.u = p.u, p.star = p.star)
} 

generate.X <- function(Y, p, null.idx = NULL) {
    n <- nrow(Y)
    X <- matrix(0, nrow = n - p, ncol = 2 * p)
    for (i in 1:nrow(X)) {
        X[i, 1:p] <- Y[(i + p - 1):i, 1]
        X[i, -(1:p)] <- Y[(i + p - 1):i, 2]
    }
    if (!is.null(null.idx)) {
        X <- X[, 1:p + p*(null.idx == 2)]
    }
    if (is.null(dim(X))) {
      return(matrix(X, ncol = 1))
    }
    X
}

epanechnikov <- Vectorize(
    function(u, h) {
        .75 * max(0, 1 - (u / h)^2) / h
    }, vectorize.args = c("u")
)

gaussian <- Vectorize({
  function(u, h) {
    dnorm(u/h)/h
  }
})

get.coefficients <- function(uo, y, X, u, kernel.function, h) {
    W <- diag(kernel.function(u - uo, h))
    X.tilde <- cbind(X, X * (u - uo))
    tmp <- MASS::ginv(t(X.tilde) %*% W %*% X.tilde) %*% t(X.tilde) %*% W %*% y
    tmp[1:(nrow(tmp)/2), ]
}

block_permute1 <- function(x, nblocks){
    n <- length(x)
    block_size <- n %/% nblocks
    res <- numeric(n)
    curr_idx <- 1
    permuted_blocks <- sample(1:nblocks)
    # print("Order of the permuted blocks") 
    # print(permuted_blocks)
    # print(paste(rep("=", 20), collapse = ""))
    for (i in permuted_blocks) {
        start_block <- (i - 1) * block_size + 1
        end_block <- start_block + block_size - 1
        if (i == nblocks) end_block <- n
        block <- x[start_block:end_block]
        res[curr_idx:(curr_idx + length(block) - 1)] <- block
        curr_idx <- curr_idx + length(block)
    }
    res
}

block_permute2 <- function(x, block_size){
    n <- length(x)
    nblocks <- n %/% block_size
    res <- numeric(n)
    curr_idx <- 1
    permuted_blocks <- sample(1:nblocks)
    for (i in permuted_blocks) {
        start_block <- (i - 1) * block_size + 1
        end_block <- start_block + block_size - 1
        if (i == nblocks) end_block <- n
        block <- x[start_block:end_block]
        res[curr_idx:(curr_idx + length(block) - 1)] <- block
        curr_idx <- curr_idx + length(block)
    }
    res
}

block_permute <- function(x, block_size){
    n <- length(x)
    splits <- n %/% block_size
    remainder <- n %% block_size
    nblocks <- splits + (remainder > 0)
    res <- numeric(n)
    curr_idx <- 1
    permuted_blocks <- sample(1:nblocks)
    for (i in permuted_blocks) {
        start_block <- (i - 1) * block_size + 1
        end_block <- start_block + block_size - 1
        if (i == nblocks) end_block <- n
        block <- x[start_block:end_block]
        res[curr_idx:(curr_idx + length(block) - 1)] <- block
        curr_idx <- curr_idx + length(block)
        
    }
    res
}

gen_permutation <- function(X, permute = function(x) x, cols.to.permute = NULL) {
  
  if (is.null(dim(X))) {
    return(matrix(X, ncol = 1))
  }
  
  if (!is.null(cols.to.permute)) {
        for (i in cols.to.permute) {
        X[, i] <- permute(X[, i])
    }
  }
  
  X
}

fcar.fit <- function(Y, u, kernel.function, h, idx, p = 1, p.u = 1, null.idx = NULL, 
                     npoints = 75, permute = function(x) x, cols.to.permute = NULL) {
    
    Mybounds <- as.numeric(quantile(u, probs = c(.05, .95)))
    fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = npoints)
    fhat_points <- stats::filter(fhat_int, c(.5, .5))[-length(fhat_int)]
    comdiff <- fhat_points[2] - fhat_points[1]
    fhat_points <- c(fhat_points[1] - comdiff, fhat_points, fhat_points[length(fhat_points)] + comdiff)
    
    ajs <- matrix(0, nrow = ifelse(is.null(null.idx), 2 * p, p), ncol = length(fhat_points))
    
    p.star <- max(p, p.u)
    X <- generate.X(Y, p, null.idx)
    if (p.star > p) {
      X <- X[-(1:(p.star - p)), ]
    }
    X <- gen_permutation(X, permute, cols.to.permute)
    n <- nrow(Y)
    y <- Y[-(1:p.star), idx]
    u <- u[1:(n - p.star)]
    h <- (max(u) - min(u)) * h

    for (i in 1:length(fhat_points)) {
        ajs[, i] <- get.coefficients(fhat_points[i], y, X, u, kernel.function, h)
    }

    cat_u <- as.numeric(cut(u, breaks = c(-Inf, fhat_points, Inf)))
    cat_u[cat_u > length(fhat_points)] <- length(fhat_points)
    yhat <- numeric(length(y))
    resids <- numeric(length(y))
    
    for (i in 1:length(y)) {
        yhat[i] <- X[i, ] %*% ajs[, cat_u[i]]
        resids[i] <- y[i] - yhat[i]
    }
    list(yhat = yhat, coeffs = ajs, functional_points = fhat_points, residuals = resids)
}

mse <- function(x, y) {
    mean((x - y)^2)
}



permutation.test <- function(Y, u, kernel.function, h, p, p.u, permute, npoints = 75, P = 500) {
  
  fit1 <- fcar.fit(Y, u, kernel.function, h, idx = 1, p = p, p.u = p.u, null.idx = 1, npoints = npoints)
  fit2 <- fcar.fit(Y, u, kernel.function, h, idx = 2, p = p, p.u = p.u, null.idx = 2, npoints = npoints)
  fit.full1 <- fcar.fit(Y, u, kernel.function, h, 1, p, p.u)
  fit.full2 <- fcar.fit(Y, u, kernel.function, h, 2, p, p.u)
  null.resid1 <- sum((fit1$residuals)^2)
  null.resid2 <-  sum((fit2$residuals)^2)
  full.resid1 <- sum((fit.full1$residuals)^2)
  full.resid2 <-  sum((fit.full2$residuals)^2)
  Fobs1 <- (null.resid1 - full.resid1)/(2 * p * full.resid1/(nrow(Y) - 4 * p))
  Fobs2 <- (null.resid2 - full.resid2)/(2 * p * full.resid2/(nrow(Y) - 4 * p))
  F1 <- numeric(P)
  F2 <- numeric(P)
  
  for (i in 1:P) {
    if (i %% 100 == 0) print(paste0("Permutation: ", i))
    temp.fit1 <- fcar.fit(Y, u, kernel.function, h, idx = 1, p = p, p.u = p.u,
                          npoints = npoints, permute = permute, cols.to.permute = (p+1):(2*p))
    temp.fit2 <- fcar.fit(Y, u, kernel.function, h, idx = 2, p = p, p.u = p.u,
                          npoints = npoints, permute = permute, cols.to.permute = 1:p)
    per.resid1 <- sum(temp.fit1$residuals^2)
    per.resid2 <- sum(temp.fit2$residuals^2)
    F1[i] <- (null.resid1 - per.resid1)/(2 * p * per.resid1/(nrow(Y) - 4 * p))
    F2[i] <- (null.resid2 - per.resid2)/(2 * p * per.resid2/(nrow(Y) - 4 * p))
  }

  list(null.stat1 = Fobs1, null.stat2 = Fobs2,
  ref.distribution1 = F1,
  ref.distribution2 = F2)
}

gc.test <- function(niters, numcores, Tlength, obs.per.block, sim.function, P) {
  registerDoParallel(numcores)
  foreach(i = 1:niters) %dopar% {
      set.seed(i)
      gen.data <- sim.function(Tlength)
      
      per.test <- permutation.test(gen.data$ts, gen.data$ref, gaussian, .2, 
                                   gen.data$p, gen.data$p.u, 
                                   permute = function(x) block_permute(x, obs.per.block), P = P)
      print(paste0("Iteration ", i, " completed"))
      per.test
  }
}

#function was used on simulationslag1old.R
# gc.test.old <- function(niters, numcores, Tlength, nblock, sim.function, P) {
#   registerDoParallel(numcores)
#   foreach(i = 1:niters) %dopar% {
#       set.seed(i)
#       gen.data <- sim.function(Tlength)
      
#       per.test <- permutation.test(gen.data$ts, gen.data$ref, gaussian, .2, 
#                                    gen.data$p, gen.data$p.u, 
#                                    permute = function(x) block_permute1(x, nblock), P = P)
#       print(paste0("Iteration ", i, " completed"))
#       per.test
#   }
# }