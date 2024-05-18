gendata <- function(Tlength,d,Y_d,f11,f12,f21,f22,g11=0,g12=0,g21=0,g22=0){
  y <- matrix(NA,nrow = 500+Tlength,ncol = 2)
  y1 <- c(1,1)
  ylag <- rep(1,d)
  
  e <- MASS::mvrnorm(n = 500+Tlength,mu = rep(0,3),Sigma = diag(rep(1,3)))
  
  for(t in 1:(500+Tlength)){
    y[t,1] <- f11(ylag[d],g11)*y1[1] + f12(ylag[d],g12)*y1[2] + e[t,1]
    y[t,2] <- f21(ylag[d],g21)*y1[1] + f22(ylag[d],g22)*y1[2] + e[t,2]
    y1 <- c(y[t,1],y[t,2])
    for(i in d:2){
      ylag[i] <- ylag[i-1]
    }
    if(Y_d == 0){
      ylag[1] <- e[t,3]
    } else {
      ylag[1] <- y[t,Y_d] 
    }
  }
  y <- y[-c(1:500),]
  e <- e[-c(1:500),]
  
  if(Y_d == 0){
    output <- list(y = y, ref = e[,3])
  } else {
    output <- list(y = y, ref = y[,Y_d]) 
  }
  return(output)
}

epanechnikov <- Vectorize(
    function(u, h) {
        .75 * max(0, 1 - (u / h)^2) / h
    }, vectorize.args = c("u")
)

get.coefficients <- function(uo, y, X, u, kernel.function, h) {
    W <- diag(kernel.function(u - uo, h))
    X.tilde <- cbind(X, X * (u - uo))
    tmp <- cbind(diag(rep(1, 2)), diag(rep(0, 2)))
    tmp %*% solve(t(X.tilde) %*% W %*% X.tilde) %*% t(X.tilde) %*% W %*% y
}

fcar.fit <- function(y, X, u, kernel.function, h, npoints = 75) {
    
    Mybounds <- as.numeric(quantile(u, probs = c(.05, .95)))
    fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = npoints)
    fhat_points <- stats::filter(fhat_int, c(.5, .5))[-length(fhat_int)]
    comdiff <- fhat_points[2] - fhat_points[1]
    fhat_points <- c(fhat_points[1] - comdiff, fhat_points, fhat_points[length(fhat_points)] + comdiff)
    ajs <- matrix(0, nrow = 2, ncol = length(fhat_points))
    for (i in 1:length(fhat_points)) {
        ajs[, i] <- get.coefficients(fhat_points[i], y, X, u, epanechnikov, h)
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

plot.functions <- function(real, hats, pts) {
    par(mfrow = c(2, 1))
    for (i in 1:length(real)) {
        tmp.real <- real[[i]](pts)
        tmp.hat <- hats[i, ]
        plot(pts, tmp.real, ylim = c(min(tmp.real, tmp.hat), max(tmp.real, tmp.hat)))
        lines(pts, tmp.hat, col = "red")
    }
}

mse <- function(x, y) {
    mean((x - y)^2)
}

gen_permutation <- function(X, cols.to.permute = NA) {
  if (!is.na(cols.to.permute)) {
    idx <- 1:nrow(X)
    for (i in cols.to.permute) {
      X[, i] <- X[, i][sample(idx)]
    }
  }
  X
}

permutation.test <- function(Y1, Y2, u, X, kernel.function, h, cols.to.permute = c(NA, NA), npoints = 75, P = 500) {
  
  fit1 <- fcar.fit(Y1, X, u, kernel.function, h, npoints)
  fit2 <- fcar.fit(Y2, X, u, kernel.function, h, npoints)
  null.resid <- sum((fit1$residuals)^2) + sum((fit2$residuals)^2)
  resids <- numeric(P)
  
  for (i in 1:P) {
    if (i %% 10 == 0) print(paste0("Iteration: ", i))
    temp.fit1 <- fcar.fit(Y1, gen_permutation(X, cols.to.permute[1]), u, epanechnikov, h)
    temp.fit2 <- fcar.fit(Y2, gen_permutation(X, cols.to.permute[2]), u, epanechnikov, h)
    resids[i] <- sum(temp.fit1$residuals^2) + sum(temp.fit2$residuals^2)
  }

  list(null.resid = null.resid, ref.distribution = resids)
}
