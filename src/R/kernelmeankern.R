library(kernlab)

kernelmeankernel <- function(x1,x2, sigma=1.0) {
    mmd <- kmmd(x1, x2, ntimes=0, replace=FALSE)@mmdstats[1]
    return(exp(-mmd / (2*sigma)))
}

exp.mean <- function(x, sigma=1.0) {
    n <- length(x)
    kernel.mat <- matrix(0, n, n)
    for(i in 1:n) {
        for(j in i:n) {
            mmd <- kernelmeankernel(x[[1]], x[[2]])
            kernel.mat[i,j] <- mmd
            kernel.mat[j,i] <- mmd
        }
    }
    return(exp(-kernel.mat / (2*sigma)))
}

fit.gp <- function(X, y, kernel=rbfdot(), add.kernels=list(), var=1.0) {
    total.features <- ncol(X) + length(add.kernels)
    K <- (ncol(X) / total.features) * kernelMatrix(kernel, X)
    if(length(add.kernels) > 0) {   
        # uniform combination for now, we can get fancier later
        for(i in 1:length(add.kernels)) {
            K <- K + (1 / total.features) * add.kernels[[i]]
        }
    }
    K <- (K / 1 + length(add.kernels))
    my.alpha <- solve(K + diag(rep(var, length = ncol(K)))) %*% y
    return(list(alpha=my.alpha, K=K, X=X, add.kernels=add.kernels))
}

predict.gp <- function(gp, new.X, kernel=rbfdot()) {
    K <- kernelMatrix(kernel, gp$X, new.X)
    if(length(gp$add.kernels) > 0) {
        for(i in 1:length(gp$add.kernels[[i]])) {
            K <- K + gp$add.kernels[[i]]
        }
    }
    return(K %*% gp$alpha)
}

fit.relational.gp <- function(X, y, kernel=rbfdot(), relational.features=list(), var=1.0) {
    # learn the relational kernels first
    relational.kernels <- list()
    if(length(relational.features) > 0) {
        for(i in 1:length(relational.features)) {
            relational.kernels[[i]] <- exp.mean(relational.features[[i]])
        }
    }
    gp <- fit.gp(X, y, kernel=kernel, add.kernels=relational.kernels, var=var)
    return(gp) 
}
