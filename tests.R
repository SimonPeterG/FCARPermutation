u <- rnorm(100)
Mybounds <- as.numeric(quantile(u, probs = c(0.05, 0.95)))
fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = 50) #getting a sequence of evenly spaced numbers between the quantile 5 and the quantile 95
fhat_points <- stats::filter(fhat_int, c(0.5, 0.5))[-length(fhat_int)] #getting the midpoint on every subinterval of the sequence defined above
comdiff <- fhat_points[2] - fhat_points[1] #distance between all points (they are evenly spaced)
#adding one point at the beginning and one at the end (still evenly distributed)
fhat_points <- c(fhat_points[1]-comdiff,fhat_points,fhat_points[length(fhat_points)]+comdiff)  

y
u
