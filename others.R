APEq <- function(y,u,p,d, q, bwp = 0.1,intnum = 50){
  
  lagrm <- max(p,d)
  k <- ncol(y)
  Tlength <- nrow(y)
  r <- floor(0.1 * Tlength)

  y.new <- y[1:(Tlength - r * q), ]
  u.new <- u[1:(Tlength - r * q)]
  
  Mybounds <- as.numeric(quantile(u.new,probs = c(0.05,0.95)))
  fhat_int <- seq(Mybounds[1],Mybounds[2],length = intnum)
  fhat_points <- as.numeric(stats::filter(fhat_int, c(0.5,0.5)))[-length(fhat_int)]
  comdiff <- fhat_points[2] - fhat_points[1]
  fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  
  fhat_mean <- array(rep(NA,(k^2)*p*length(fhat_points)),dim = c(k,k*p,length(fhat_points)))

  for(i in 1:length(fhat_points)){
    try({
      fhat_mean[,,i] <- t(FAR_u0(y.new,u.new,p,d,bwp,u0 = fhat_points[i]))[,1:(k*p)]
    },silent = TRUE)
    # cat(paste("Progress: ",i,"/",length(fhat_points))," \r")
  }
  
  mat_Y <- matrix(NA,nrow = Tlength-lagrm,ncol = k)
  mat_X <- matrix(NA,nrow = Tlength-lagrm,ncol = k*p)
  mat_U <- matrix(NA,nrow = r, ncol = 1)
  
  mat_Y[1:(Tlength-lagrm),] <- y[(lagrm+1):Tlength,]
  mat_U[1:r,] <- u[(Tlength - r * q + 1 - d):(Tlength - r * q + r - d)]
  
  for(l in 1:p){
    for(j in 1:k){
      mat_X[,k*(l-1)+j] <- y[(lagrm+1-l):(Tlength-l),j]
    }
  }  
  
  cat_U <- as.numeric(cut(mat_U,breaks = c(-Inf,fhat_int,Inf)))
  mat_Yhat <- matrix(NA,nrow = r ,ncol = k) #(Tlength - r * q + 1):(Tlength - r * q + r)
  mat_Y <- y[(Tlength - r * q + 1):(Tlength - r * q + r), ]
  
  for(t in 1:r){
    fhat_t <- fhat_mean[,,cat_U[t]]
    idx <- Tlength - r * q - lagrm + t
    mat_Yhat[t,] <- t(fhat_t%*%matrix(mat_X[idx,],ncol=1))
  }
  
  sum(apply(mat_Y - mat_Yhat, 1, function(x) sum(x^2)))

}

Ape <- function(y,u,p,d) {
    sum(sapply(1:4, function(x) APEq(y, u, p, d, x)))
}

Ape(mydata[[1]], mydata[[2]], 1, 2)

ananai <- apply(cbind(expand.grid(p = 1:4, d = 1:10), 1:40), 1,
    function(x) Ape(mydata[[1]], mydata[[2]], x[1], x[2]))

expand.grid(p = 1:4, d = 1:10)[which.min(ananai), ]

APEq(mydata[[1]], mydata[[2]], 1, 2, 1)
APEq(mydata[[1]], mydata[[2]], 1, 2, 2)
APEq(mydata[[1]], mydata[[2]], 1, 2, 3)
APEq(mydata[[1]], mydata[[2]], 1, 2, 4)

APEq(mydata[[1]], mydata[[2]], 2, 2, 1)
APEq(mydata[[1]], mydata[[2]], 2, 2, 2)
APEq(mydata[[1]], mydata[[2]], 2, 2, 3)
APEq(mydata[[1]], mydata[[2]], 2, 2, 4)

lol <- FAR_est(mydata[[1]], mydata[[2]], 2, 2)

mean(apply(lol$mat_resid, 1, function(x) sum(x^2)))
mean(apply(myest$mat_resid, 1, function(x) sum(x^2)))

plot(lol$fhat_points, f11(lol$fhat_points))
plot(lol$fhat_points, lol$fhat_mean[1, 1, ])

plot(lol$fhat_points, lol$fhat_mean[1, 2, ])

ts.plot(mydata[[1]][, 1])
lines(lol$mat_Yhat[, 1], col = "red")
myest
sum(sapply(1:4, function(Q) APEq(mydata[[1]], mydata[[2]], 1, 2, Q)))
1000 - 100*1 
