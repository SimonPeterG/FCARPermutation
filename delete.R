get.vcoefficients <- function(uo, y, X, u, kernel.function, h) {
    W <- diag(kernel.function(u - uo, h))
    X.tilde <- cbind(X, X * (u - uo))
    XW <- t(X.tilde) %*% W 
    Q <- solve(XW %*% X.tilde)
    XWy <- XW %*% y
    tmp <- Q %*% XWy
    tmp[1:(nrow(tmp)/2), ]
}

fvcar.fit <- function(y, u, kernel.function, p, d, h = NULL, npoints = 75) {
    
    p.star <- max(p, d)
    init.x <- p.star + 1 - p; end.x <- nrow(y) - p.star + init.x - 1
    init.u <- p.star + 1 - d; end.u <- nrow(y) - p.star + init.u - 1
    
    X <- y[init.x:end.x, ]
    y <- y[-(1:p.star), ]
    u <- u[init.u:end.u]

    if (is.null(h)) h = (max(u) - min(u)) * .2

    Mybounds <- as.numeric(quantile(u, probs = c(.05, .95)))
    fhat_int <- seq(Mybounds[1], Mybounds[2], length.out = npoints)
    fhat_points <- stats::filter(fhat_int, c(.5, .5))[-length(fhat_int)]
    comdiff <- fhat_points[2] - fhat_points[1]
    fhat_points <- c(fhat_points[1] - comdiff, fhat_points, fhat_points[length(fhat_points)] + comdiff)

    ajs <- array(0, dim = c(ncol(y), ncol(y), length(fhat_points)))
    for (i in 1:length(fhat_points)) {
        ajs[,,i] <- get.vcoefficients(fhat_points[i], y, X, u, kernel.function, h)
    }
    cat_u <- as.numeric(cut(u, breaks = c(-Inf, fhat_points, Inf)))
    cat_u[cat_u > length(fhat_points)] <- length(fhat_points)
    yhat <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    resids <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for (i in 1:nrow(y)) {
        yhat[i, ] <- X[i, ] %*% ajs[,,cat_u[i]]
        resids[i, ] <- y[i, ] - yhat[i, ]
    }
    list(yhat = yhat, coeffs = ajs, functional_points = fhat_points, residuals = resids)
}
