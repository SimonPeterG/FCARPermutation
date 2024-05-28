rm(list = ls())
library(Rcpp)
source("auxiliary_functions.R")

df <- read.csv("subject_data/Subject_01.csv")
cols <- c("Fp1", "Fp2", "TP7", "TP8", "Cz")

get.epoch.state <- function(x, state = "Alert", epoch = 1) {
    mask <- x$state == state & x$epoch == 1
    as.matrix(x[mask, cols])
}


Y.alert <- get.epoch.state(df)
Y.alert <- as.matrix(scale(Y.alert))
U.alert <- Y.alert[, ncol(Y.alert)]
Y.alert <- Y.alert[, -ncol(Y.alert)]

param.grid <- expand.grid(p = 1:4, d = 1:10)

apes <- apply(param.grid, 1, function(x) Ape(Y.alert, U.alert, x[1], x[2]))

params <- as.numeric(param.grid[which.min(apes), ])

alert.fit <- FAR_est(Y.alert, U.alert, params[1], params[2])

X.alert <- alert.fit$mat_X

p.star.alert <- max(params)

length(Y.alert[-(1:p.star.alert), 1])


lol <- fcar.fit(matrix(Y.alert[-(1:p.star.alert), 1], ncol = 1), X.alert, U.alert[p.star.alert:(length(U.alert)-1)],
     epanechnikov, .5, 4, 4)

y <- matrix(Y.alert[-(1:p.star.alert), 1], ncol = 1)
u <- U.alert[p.star.alert:(length(U.alert)-1)]
X <- X.alert
h <- (max(u) - min(u)) * .1


adj.alert <- matrix(0, nrow = params[1], ncol = ncol(Y.alert))

for (i in 1:4) {
    for (j in 1:4) {
        if (i != j) {
            tmp.y <- matrix(Y.alert[-(1:p.star.alert), i], ncol = 1)
            cols <- seq(j, ncol(X), by = params[1])
            adj.alert[i, j] <- individual.permutation(tmp.y, u, X, epanechnikov, .1, cols, params[1], ncol(Y.alert))
        }
    }
}

adj.alert

permutation.test()
