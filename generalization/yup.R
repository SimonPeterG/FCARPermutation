foo <- function(i, x) {
  c(mean(x[[i]]$null.stat1 <= x[[i]]$ref.distribution1),
    mean(x[[i]]$null.stat2 <= x[[i]]$ref.distribution2)) < 0.05
}

expar1b1 <- gc.test(10, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar1)), P = 1000)

aygono <- gc.test(10, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar1)), P = 1000)

expar2b1 <- gc.test(5 , numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar2)), P = 1000)

expar3b1 <- gc.test(10, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(expar3)), P = 1000)

logis1b1 <- gc.test(5, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis1)), P = 1000)

logis2b1 <- gc.test(5, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis2)), P = 1000)

logis3b1 <- gc.test(5, numcores, 256, block.sizes[1], 
                    function(n) simulate.far(n, 2, 0, 1, 1, list(logis3)), P = 1000)

lapply(1:10, function(i) foo(i, expar1b1))
lapply(1:10, function(i) foo(i, aygono))
lapply(1:5, function(i) foo(i, expar2b1))
lapply(1:10, function(i) foo(i, expar3b1))
lapply(1:5, function(i) foo(i, logis1b1))
lapply(1:5, function(i) foo(i, logis2b1))
lapply(1:5, function(i) foo(i, logis3b1))

identical(expar1b1, aygono)
