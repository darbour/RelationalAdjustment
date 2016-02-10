library(kernlab)

kernelmeankernel <- function(x1,x2, sigma=1.0) {
    kpar <- list(sigma=sigest(rbind(x1,x2),scaled=FALSE, frac=1)[2])
    mmd <- kmmd(x1, x2, ntimes=0, replace=FALSE, kpar=kpar)@mmdstats[1]
    return(exp(-mmd / (2*sigma)))
}

exp.mean <- function(x, sigma=1.0) {
    require(hash)
    n <- length(x)
    kernel.mat <- matrix(0, n, n)
    # has because this thing takes too long
    cache <- hash()
    for(i in 1:n) {
        for(j in i:n) {
            hash.name <- paste(length(x[[i]]), sum(x[[i]]), length(x[[j]]), sum(x[[j]]), sep=",")
            if(has.key(hash.name, cache)) {
                mmd <- cache[[hash.name]]
            } else {
                mmd <- kernelmeankernel(matrix(x[[i]]), matrix(x[[j]]))
                cache[[hash.name]] <- mmd
                # mmd is symmetric so store it both ways
                cache[[paste(length(x[[j]]), sum(x[[j]]), length(x[[i]]), sum(x[[i]]), sep=",")]] <- mmd
            }
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

# Set the values of friends to meet a certain proportion. 
# NB: does NOT return a valid global configuration
set.to.proportion <- function(binary.friend.vals, prop) {
    N <- length(binary.friend.vals)
    new.vals <- binary.friend.vals
    for(i in 1:N) {
        num.friends <- length(new.vals[[i]])
        num.treated <- round(prop * num.friends)
        new.vals[[i]] <- c(rep(1, num.treated), rep(0, num.friends - num.treated))
    }
    return(new.vals)
}

# make predictions using the kernel mean embedding with a certain proportion of friends treated
predict.prop.treatment.gp <- function(gp, new.X, relational.features, prop) {
    new.features <- set.to.proportion(relational.features, prop)
    new.kernel <- exp.mean(new.features)
    K <- kernelMatrix(kernel, gp$X, new.X) + new.kernel
    return(K %*% gp$alpha)
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

fit.relational.gp <- function(X, y, kernel=rbfdot(), relational.features=list(), var=1.0, num.relational.features=1) {
    # learn the relational kernels first
    relational.kernels <- list()
    if(num.relational.features > 0) {
        if(num.relational.features == 1) 
            relational.kernels[[1]] <- exp.mean(relational.features)
        else {
            for(i in 1:num.relational.features) {
                relational.kernels[[i]] <- exp.mean(relational.features[[i]])
            }
        }
    }
    gp <- fit.gp(X, y, kernel=kernel, add.kernels=relational.kernels, var=var)
    return(gp) 
}
