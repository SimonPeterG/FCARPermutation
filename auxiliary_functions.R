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

fcar.fit <- function(y, X, u, kernel.function, h, p = 1, k = 2, npoints = 75) {
    
    Mybounds <- as.numeric(quantile(u, probs = c(.05, .95)))
    fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = npoints)
    fhat_points <- stats::filter(fhat_int, c(.5, .5))[-length(fhat_int)]
    comdiff <- fhat_points[2] - fhat_points[1]
    fhat_points <- c(fhat_points[1] - comdiff, fhat_points, fhat_points[length(fhat_points)] + comdiff)
    ajs <- matrix(0, nrow = p * k, ncol = length(fhat_points))
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

gen_permutation <- function(X, nblocks, cols.to.permute = NA) {
    for (i in cols.to.permute) {
      X[, i] <- block_permute1(X[, i], nblocks)
    }
  X
}

permutation.test <- function(Y1, Y2, u, X, kernel.function, h, nblocks, npoints = 75, P = 500) {
  
  fit1 <- fcar.fit(Y1, matrix(X[, 1], ncol = 1), u, kernel.function, h, npoints = npoints, k = 1)
  fit2 <- fcar.fit(Y2, matrix(X[, 2], ncol = 1), u, kernel.function, h, npoints = npoints, k = 1)
  fit.full1 <- fcar.fit(Y1, X, u, kernel.function, h)
  fit.full2 <- fcar.fit(Y2, X, u, kernel.function, h)
  null.resid1 <- sum((fit1$residuals)^2)
  null.resid2 <-  sum((fit2$residuals)^2)
  full.resid1 <- sum((fit.full1$residuals)^2)
  full.resid2 <-  sum((fit.full2$residuals)^2)
  Fobs1 <- (null.resid1 - full.resid1)/(full.resid1/(length(Y1) - 2))
  Fobs2 <- (null.resid2 - full.resid2)/(full.resid2/(length(Y1) - 2))
  F1 <- numeric(P)
  F2 <- numeric(P)
  
  for (i in 1:P) {
    if (i %% 10 == 0) print(paste0("Iteration: ", i))
    temp.fit1 <- fcar.fit(Y1, gen_permutation(X, nblocks, 2), u, kernel.function, h, npoints = npoints)
    temp.fit2 <- fcar.fit(Y2, gen_permutation(X, nblocks, 1), u, kernel.function, h, npoints = npoints)
    per.resid1 <- sum(temp.fit1$residuals^2)
    per.resid2 <- sum(temp.fit2$residuals^2)
    F1[i] <- (null.resid1 - per.resid1)/(per.resid1/(length(Y1) - 2))
    F2[i] <- (null.resid2 - per.resid2)/(per.resid2/(length(Y1) - 2))
  }

  print("after the foor loop")

  list(null.stat1 = Fobs1, null.stat2 = Fobs2,
  fit1 = fit1, fit2 = fit2,
  ref.distribution1 = F1,
  ref.distribution2 = F2)
}

individual.permutation <- function(Y1, u, X, kernel.function, h, cols, p = 1, k = 2, npoints = 75, P = 100) {
  
  fit1 <- fcar.fit(Y1, X, u, kernel.function, h, p, k, npoints)
  null.resid1 <- sum((fit1$residuals)^2)
  resids1 <- numeric(P)
  
  for (i in 1:P) {
    if (i %% 250 == 0) print(paste0("Iteration: ", i))
    temp.fit1 <- fcar.fit(Y1, gen_permutation(X, cols), u, epanechnikov, h, p, k, npoints)
    resids1[i] <- sum(temp.fit1$residuals^2)
  }

  as.numeric(null.resid1 < quantile(resids1, probs = c(0.05)))
}


##### Estimating multivariate FAR model locally at u0
## Inputs
#    y    : N x K matrix for the multivariate time series (N - timepoints,K - number of time series)
#    u    : vector of length N consisting the reference variable
#    p    : order of the FAR model
#    d    : lag of reference variable
#   bwp   : bandwidth for kernel (numeric in (0,1))
#    u0   : specific value of the reference variable for local linear approximation
## Output : K x 2KP matrix for the functional parameter estimates (local at u0) 

# cppFunction('arma::mat FAR_u0(arma::mat& y,arma::vec& u,double p,double d,double bwp,double u0){
    
#     arma::vec pd = {p,d};
#     double lagrm = max(pd);
#     int Tlength = y.n_rows;
#     int k = y.n_cols;
    
#     arma::mat mat_Y = y(arma::span(lagrm,Tlength-1),arma::span(0,k-1));
    
#     arma::mat mat_X(Tlength-lagrm,k*p);
#     for(int l = 0; l < p; l++){
#       for(int j = 0; j < k; j++){
#       mat_X(arma::span::all,k*l+j) = y(arma::span(lagrm-l-1,Tlength-l-2),j);
#       }
#     }

#     arma::mat mat_U = u(arma::span(lagrm-d,Tlength-d-1));
#     double bw = bwp * as_scalar(range(mat_U));
#     arma::mat W = ((1/arma::datum::sqrt2pi)*exp(-0.5*((mat_U-u0)/bw)%((mat_U-u0)/bw)))/bw;
#     arma::mat mat_W = diagmat(W);
    
#     arma::mat mat_UX(Tlength-lagrm,k*p);
#     for(int l = 0; l < k*p; l++){
#       mat_UX(arma::span::all,l) = mat_X(arma::span::all,l)%(mat_U-u0);
#     }
    
#     arma::mat mat_Xcurl = join_rows(mat_X,mat_UX);
#     arma::mat mat_XW = trans(mat_Xcurl) * mat_W;
#     arma::mat mat_Q = mat_XW * mat_Xcurl;
#     arma::mat mat_XWY = mat_XW * mat_Y;
#     arma::mat mat_Qinv = inv(mat_Q);
#     arma::mat fhat = mat_Qinv * mat_XWY;

#     return fhat;
#   }',depends="RcppArmadillo")

# APEq <- function(y,u,p,d, q, bwp = 0.1,intnum = 50){
  
#   lagrm <- max(p,d)
#   k <- ncol(y)
#   Tlength <- nrow(y)
#   r <- floor(0.1 * Tlength)

#   y.new <- y[1:(Tlength - r * q), ]
#   u.new <- u[1:(Tlength - r * q)]
  
#   Mybounds <- as.numeric(quantile(u.new,probs = c(0.05,0.95)))
#   fhat_int <- seq(Mybounds[1],Mybounds[2],length = intnum)
#   fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
#   comdiff <- fhat_points[2] - fhat_points[1]
#   fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  
#   fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)),dim = c(k,k*p,length(fhat_points)))

#   for(i in 1:length(fhat_points)){
#     try({
#       fhat_mean[,,i] <- t(FAR_u0(y.new,u.new,p,d,bwp,u0 = fhat_points[i]))[,1:(k*p)]
#     },silent = TRUE)
#     # cat(paste("Progress: ",i,"/",length(fhat_points))," \r")
#   }
  
#   mat_Y <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
#   mat_X <- matrix(NA,nrow = Tlength-lagrm,ncol = k*p)
#   mat_U <- matrix(NA,nrow = r, ncol = 1)
  
#   mat_Y[1:(Tlength-lagrm),] <- y[(lagrm+1):Tlength,]
#   mat_U[1:r,] <- u[(Tlength - r * q + 1 - d):(Tlength - r * q + r - d)]
  
#   for(l in 1:p){
#     for(j in 1:k){
#       mat_X[,k*(l-1)+j] <- y[(lagrm+1-l):(Tlength-l),j]
#     }
#   }  
  
#   cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
#   mat_Yhat <- matrix(NA,nrow = r ,ncol = k) #(Tlength - r * q + 1):(Tlength - r * q + r)
#   mat_Y <- y[(Tlength - r * q + 1):(Tlength - r * q + r), ]
  
#   for(t in 1:r){
#     fhat_t <- fhat_mean[,,cat_U[t]]
#     idx <- Tlength - r * q - lagrm + t
#     mat_Yhat[t,] <- t(fhat_t%*%matrix(mat_X[idx,],ncol=1))
#   }
  
#   sum(apply(mat_Y - mat_Yhat, 1, function(x) sum(x^2)))

# }

# Ape <- function(y,u,p,d) {
#     sum(sapply(1:4, function(x) APEq(y, u, p, d, x)))
# }

# FAR_est <- function(y,u,p,d,bwp = 0.1,intnum = 50){
  
#   lagrm <- max(p,d)
#   k <- ncol(y)
#   Tlength <- nrow(y)
  
#   Mybounds <- as.numeric(quantile(u,probs = c(0.05,0.95)))
#   fhat_int <- seq(Mybounds[1],Mybounds[2],length = intnum)
#   fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
#   comdiff <- fhat_points[2] - fhat_points[1]
#   fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  
#   fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)),dim = c(k,k*p,length(fhat_points)))
  
#   print(length(fhat_points))

#   for(i in 1:length(fhat_points)){
#     try({
#       fhat_mean[,,i] <- t(FAR_u0(y,u,p,d,bwp,u0 = fhat_points[i]))[,1:(k*p)]
#     },silent = TRUE)
#     # cat(paste("Progress: ",i,"/",length(fhat_points))," \r")
#   }
  
#   mat_Y <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
#   mat_X <- matrix(NA,nrow = Tlength-lagrm,ncol = k*p)
#   mat_U <- matrix(NA,nrow = Tlength-lagrm,ncol = 1)
  
#   mat_Y[1:(Tlength-lagrm),] <- y[(lagrm+1):Tlength,]
#   mat_U[1:(Tlength-lagrm),] <- u[(lagrm+1-d):(Tlength-d)]
#   for(l in 1:p){
#     for(j in 1:k){
#       mat_X[,k*(l-1)+j] <- y[(lagrm+1-l):(Tlength-l),j]
#     }
#   }  
#   cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
#   print(cat_U)
#   mat_Yhat <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
  
#   print(Tlength-lagrm)
#   for(t in 1:(Tlength-lagrm)){
#     fhat_t <- fhat_mean[,,cat_U[t]]
#     mat_Yhat[t,] <- t(fhat_t%*%matrix(mat_X[t,],ncol=1))
#   }
  
#   mat_resid <- mat_Y - mat_Yhat
  
#   output <- list(mat_Y = mat_Y,mat_Yhat = mat_Yhat,mat_U = mat_U,
#                  mat_resid = mat_resid, mat_X = mat_X, fhat_points = fhat_points,
#                  fhat_mean = fhat_mean)
  
#   return(output)
  
# }
