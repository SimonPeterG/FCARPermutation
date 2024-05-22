rm(list = ls())
##### Loading needed packages
library(Rcpp)
library(MASS)
library(vars)

###################################
##### Data Generating Process #####
###################################

##### Generating bivariate FAR(p=1)
## Inputs
# Tlength : number of time points
#    d    : lag of reference variable
#   Y_d   : index for the reference variable (0 - exogenous N(0,1), 1 - Y1, 2 - Y2)
# f11-f22 : functional coefficients
# g11-g22 : subject-specific parameter for f11-f22 (set to 0 for single subject)
## Output : list of length 2 ([[1]] - bivariate series, [[2]] - reference variable)

gendata <- function(Tlength,d,Y_d,f11,f12,f21,f22,g11=0,g12=0,g21=0,g22=0){
  y <- matrix(NA,nrow = 500+Tlength,ncol = 2)
  y1 <- c(1,1)
  ylag <- rep(1,d)
  
  e <- mvrnorm(n = 500+Tlength,mu = rep(0,3),Sigma = diag(rep(1,3)))
  
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

##### Estimating multivariate FAR model locally at u0
## Inputs
#    y    : N x K matrix for the multivariate time series (N - timepoints,K - number of time series)
#    u    : vector of length N consisting the reference variable
#    p    : order of the FAR model
#    d    : lag of reference variable
#   bwp   : bandwidth for kernel (numeric in (0,1))
#    u0   : specific value of the reference variable for local linear approximation
## Output : K x 2KP matrix for the functional parameter estimates (local at u0) 

cppFunction('arma::mat FAR_u0(arma::mat& y,arma::vec& u,double p,double d,double bwp,double u0){
    
    arma::vec pd = {p,d};
    double lagrm = max(pd);
    int Tlength = y.n_rows;
    int k = y.n_cols;
    
    arma::mat mat_Y = y(arma::span(lagrm,Tlength-1),arma::span(0,k-1));

    for(int l = 0; l < p; l++){
      for(int j = 0; j < k; j++){
      mat_X(arma::span::all,k*l+j) = y(arma::span(lagrm-l-1,Tlength-l-2),j);
      }
    }

    arma::mat mat_U = u(arma::span(lagrm-d,Tlength-d-1));
    double bw = bwp * as_scalar(range(mat_U));
    arma::mat W = ((1/arma::datum::sqrt2pi)*exp(-0.5*((mat_U-u0)/bw)%((mat_U-u0)/bw)))/bw;
    arma::mat mat_W = diagmat(W);
    
    arma::mat mat_UX(Tlength-lagrm,k*p);
    for(int l = 0; l < k*p; l++){
      mat_UX(arma::span::all,l) = mat_X(arma::span::all,l)%(mat_U-u0);
    }
    
    arma::mat mat_Xcurl = join_rows(mat_X,mat_UX);
    arma::mat mat_XW = trans(mat_Xcurl) * mat_W;
    arma::mat mat_Q = mat_XW * mat_Xcurl;
    arma::mat mat_XWY = mat_XW * mat_Y;
    arma::mat mat_Qinv = inv(mat_Q);
    arma::mat fhat = mat_Qinv * mat_XWY;

    return fhat;
  }',depends="RcppArmadillo")

##### Estimating multivariate FAR model for the entire support of the reference variable
## Inputs
#    y    : N x K matrix for the multivariate time series (N - timepoints,K - number of time series)
#    u    : vector of length N consisting the reference variable
#    p    : order of the FAR model
#    d    : lag of reference variable
#   bwp   : bandwidth for kernel (numeric in (0,1))
# intnum  : number of local bins to be considered in estimation
## Output : list of length 6 ([[1]] - multivariate time series input,
#                             [[2]] - predicted values,
#                             [[3]] - input reference signal,
#                             [[4]] - residuals from model fit
#                             [[5]] - set of local points for estimating the model
#                             [[6]] - estimated functional coefficients evaluated at [[5]])

FAR_est <- function(y,u,p,d,bwp = 0.1,intnum = 50){
  
  lagrm <- max(p,d)
  k <- ncol(y)
  Tlength <- nrow(y)
  
  Mybounds <- as.numeric(quantile(u,probs = c(0.05,0.95)))
  fhat_int <- seq(Mybounds[1],Mybounds[2],length = intnum)
  fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
  comdiff <- fhat_points[2] - fhat_points[1]
  fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  
  fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)),dim = c(k,k*p,length(fhat_points)))
  
  print(length(fhat_points))

  for(i in 1:length(fhat_points)){
    try({
      fhat_mean[,,i] <- t(FAR_u0(y,u,p,d,bwp,u0 = fhat_points[i]))[,1:(k*p)]
    },silent = TRUE)
    # cat(paste("Progress: ",i,"/",length(fhat_points))," \r")
  }
  
  mat_Y <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
  mat_X <- matrix(NA,nrow = Tlength-lagrm,ncol = k*p)
  mat_U <- matrix(NA,nrow = Tlength-lagrm,ncol = 1)
  
  mat_Y[1:(Tlength-lagrm),] <- y[(lagrm+1):Tlength,]
  mat_U[1:(Tlength-lagrm),] <- u[(lagrm+1-d):(Tlength-d)]
  for(l in 1:p){
    for(j in 1:k){
      mat_X[,k*(l-1)+j] <- y[(lagrm+1-l):(Tlength-l),j]
    }
  }  
  cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
  print(cat_U)
  mat_Yhat <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
  
  print(Tlength-lagrm)
  for(t in 1:(Tlength-lagrm)){
    fhat_t <- fhat_mean[,,cat_U[t]]
    mat_Yhat[t,] <- t(fhat_t%*%matrix(mat_X[t,],ncol=1))
  }
  
  mat_resid <- mat_Y - mat_Yhat
  
  output <- list(mat_Y = mat_Y,mat_Yhat = mat_Yhat,mat_U = mat_U,
                 mat_resid = mat_resid, fhat_points = fhat_points,
                 fhat_mean = fhat_mean)
  
  return(output)
  
}

####################
##### Examples #####
####################

##### Functional coefficients - EXPAR
f11 <- function(x,theta=0){rep(-0.3,length(x)) + theta}
f12 <- function(x,theta=0){0.6*exp(-(0.30 + theta)*x^2)}
f21 <- function(x,theta=0){rep(-0.2,length(x)) + theta}
f22 <- function(x,theta=0){-0.4*exp(-(0.45 + theta)*x^2)}

# Simulate data
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)

# Time plots
par(mfrow = c(2,1))
plot.ts(mydata[[1]][,1],ylab = "Y1")
plot.ts(mydata[[1]][,2],ylab = "Y2")

# Estimate
myest <- FAR_est(y = mydata[[1]],u = mydata[[2]],p = 1,d = 2,bwp = 0.1)
par(mfrow = c(2,2))
plot(myest$fhat_points, myest$fhat_mean[1,1,],ylab = "f11",xlab = "u",ylim = c(-0.7,0.1))
lines(myest$fhat_points,f11(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[1,2,],ylab = "f12",xlab = "u",ylim = c(0,0.7))
lines(myest$fhat_points,f12(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,1,],ylab = "f21",xlab = "u",ylim = c(-0.6,0.2))
lines(myest$fhat_points,f21(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,2,],ylab = "f22",xlab = "u",ylim = c(-0.6,0.2))
lines(myest$fhat_points,f22(myest$fhat_points),col = "red")

par(mfrow = c(2, 1))
plot(mydata[[1]][, 1], type = 'l', 
  ylim = c(min(mydata[[1]][, 1], myest$mat_Yhat[, 1]), max(mydata[[1]][, 1], myest$mat_Yhat[, 1])))
lines(myest$mat_Yhat[, 1], col = 'red')

plot(mydata[[1]][, 2], type = 'l', 
  ylim = c(min(mydata[[1]][, 2], myest$mat_Yhat[, 2]), max(mydata[[1]][, 2], myest$mat_Yhat[, 2])))
lines(myest$mat_Yhat[, 2], col = 'red')


##### Functional coefficients - Logistic
f11 <- function(x,theta=0){0.8*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))-0.3}
f12 <- function(x,theta=0){rep(0.2,length(x)) + theta}
f21 <- function(x,theta=0){-0.9*((exp((5 + theta)*x))/(1+exp((5 + theta)*x)))+0.5}
f22 <- function(x,theta=0){rep(0.3,length(x)) + theta}

# Simulate data
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)

# Time plots
par(mfrow = c(2,1))
plot.ts(mydata[[1]][,1],ylab = "Y1")
plot.ts(mydata[[1]][,2],ylab = "Y2")

# Estimate
myest <- FAR_est(y = mydata[[1]],u = mydata[[2]],p = 1,d = 2,bwp = 0.1)
par(mfrow = c(2,2))
plot(myest$fhat_points, myest$fhat_mean[1,1,],ylab = "f11",xlab = "u",ylim = c(-0.8,0.8))
lines(myest$fhat_points,f11(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[1,2,],ylab = "f12",xlab = "u",ylim = c(-0.2,0.5))
lines(myest$fhat_points,f12(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,1,],ylab = "f21",xlab = "u",ylim = c(-0.8,0.8))
lines(myest$fhat_points,f21(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,2,],ylab = "f22",xlab = "u",ylim = c(0,0.6))
lines(myest$fhat_points,f22(myest$fhat_points),col = "red")



##### Functional coefficients - Constant
f11 <- function(x,theta=0){rep(-0.65,length(x)) + theta}
f12 <- function(x,theta=0){rep(0.30,length(x)) + theta}
f21 <- function(x,theta=0){rep(-0.75,length(x)) + theta}
f22 <- function(x,theta=0){rep(0.50,length(x)) + theta}

# Simulate data
mydata <- gendata(Tlength = 1000,d = 2,Y_d = 0,
                  f11 = f11, f12 = f12, f21 = f21, f22 = f22)

# Time plots
par(mfrow = c(2,1))
plot.ts(mydata[[1]][,1],ylab = "Y1")
plot.ts(mydata[[1]][,2],ylab = "Y2")

# Estimate
myest <- FAR_est(y = mydata[[1]],u = mydata[[2]],p = 1,d = 2,bwp = 0.1)
par(mfrow = c(2,2))
plot(myest$fhat_points, myest$fhat_mean[1,1,],ylab = "f11",xlab = "u",ylim = c(-1.15,-0.15))
lines(myest$fhat_points,f11(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[1,2,],ylab = "f12",xlab = "u",ylim = c(-0.1,0.6))
lines(myest$fhat_points,f12(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,1,],ylab = "f21",xlab = "u",ylim = c(-1.25,-0.25))
lines(myest$fhat_points,f21(myest$fhat_points),col = "red")

plot(myest$fhat_points, myest$fhat_mean[2,2,],ylab = "f22",xlab = "u",ylim = c(0,1))
lines(myest$fhat_points,f22(myest$fhat_points),col = "red")


