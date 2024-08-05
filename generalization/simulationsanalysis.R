tmp <- readRDS("./results/univariate/expar2_10.Rds")

p.value.decision <- function(ref, ref.dist, alpha) {
    p.value <- mean(ref.dist >= ref)
    rejection <- p.value < alpha
    list(p.value, rejection)
}


# f <- lapply(1:100, function(x) p.value.decision(tmp[[x]]$null.stat1, tmp[[x]]$ref.distribution1, 0.05))
# mean(sapply(f, function(x) x[[2]]))


# f <- lapply(1:100, function(x) p.value.decision(tmp[[x]]$null.stat2, tmp[[x]]$ref.distribution2, 0.05))
# mean(sapply(f, function(x) x[[2]]))

#analysis on univariate case 
root <- "./results/univariate/"
cases <- 1:3
block.sizes <- c(10, 50, 100, 254)
expar.paths <- expand.grid(block.sizes, cases, "expar")
logis.paths <- expand.grid(block.sizes, cases, "logis")


percentage.of.interest <- function(tmp, x, case, model) {
    
    #returning the percentage of rejections
    f <- lapply(1:100, function(x) p.value.decision(tmp[[x]]$null.stat1, tmp[[x]]$ref.distribution1, 0.05))
    r1 <- mean(sapply(f, function(x) x[[2]]))
    f <- lapply(1:100, function(x) p.value.decision(tmp[[x]]$null.stat2, tmp[[x]]$ref.distribution2, 0.05))
    r2 <- mean(sapply(f, function(x) x[[2]]))
    if (model == "logis" & case == 1) {
        c(r2, r1)
    }
    c(r1, r2)
}

lapply(1:nrow(expar.paths), function(i) {
    x <- expar.paths[i, ]
    path <- gsub("\\s", "", paste0(root, x[[3]], x[[2]], "_", x[[1]], ".Rds"))
    print(path)
    tmp <- readRDS(path)
    percentage.of.interest(tmp, i, x[[2]], x[[3]])
    }
)

lapply(1:nrow(logis.paths), function(i) {
    x <- logis.paths[i, ]
    path <- gsub("\\s", "", paste0(root, x[[3]], x[[2]], "_", x[[1]], ".Rds"))
    print(path)
    tmp <- readRDS(path)
    percentage.of.interest(tmp, i, x[[2]], x[[3]])
    }
)

