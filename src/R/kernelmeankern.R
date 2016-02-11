library(kernlab)
library(hash)
library(compiler)
fast.mmd <- kmmd

kernelmeankernel <- function(x1,x2, sigma=1.0) {
    if(all(x1 == x2)) {
        return(1)
    }
    kpar <- list(sigma=sigest(rbind(x1,x2),scaled=FALSE, frac=1)[2])
    mmd <- kmmd(x1, x2, ntimes=0, replace=FALSE, kpar=kpar)@mmdstats[1]
    return(mmd)
}

cache <- hash()
exp.mean <- function(x, y=NULL, sigma=NULL) {
    n <- length(x)
    kernel.mat <- matrix(0, n, n)
    # has because this thing takes too long
    compute.val <- function(p, q) {
        hash.name <- paste(length(p), sum(p), length(q), sum(q), sep=",")
        if(has.key(hash.name, cache)) {
            return(cache[[hash.name]])
        } else {
            mmd <- kernelmeankernel(matrix(p), matrix(q))
            cache[[hash.name]] <- mmd
            # mmd is symmetric so store it both ways
            cache[[paste(length(q), sum(q), length(p), sum(p), sep=",")]] <- mmd
            return(mmd)
        }
    }
    if(is.null(y)) {
        for(i in 1:n) {
            for(j in i:n) {
                mmd <- compute.val(x[[i]], x[[j]])
                kernel.mat[i,j] <- mmd
                kernel.mat[j,i] <- mmd
            }
        }
    } else {
        for(i in 1:n) {
           for(j in 1:n) {
                mmd <- compute.val(x[[i]], y[[j]])
                kernel.mat[i,j] <- mmd
           }
        }
    }
    if(is.null(sigma)) {
        sigma <- median(kernel.mat)
    } 
    
    return(exp(-kernel.mat/(2*sigma)))
}

fast.mat <- cmpfun(kernelMatrix)

fit.gp <- function(X, y, kernel=rbfdot(), add.kernels=list(), var=1.0) {
    total.features <- ncol(X) + length(add.kernels)
    K <- fast.mat(kernel, matrix(X))
    if(length(add.kernels) > 0) {   
        # uniform combination for now, we can get fancier later
        for(i in 1:length(add.kernels)) {
            K <- K + (1 / total.features) * add.kernels[[i]]
        }
    }
    #K <- (K / 1 + length(add.kernels))
    
    my.alpha <- solve(K + diag(rep(var, length = ncol(K)))) %*% y
    return(list(alpha=my.alpha, K=K, X=matrix(X), add.kernels=add.kernels))
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
predict.prop.treatment.gp <- function(gp, new.X, relational.features, prop, kernel=rbfdot()) {
    print("I am here")
    new.features <- set.to.proportion(relational.features, prop)
    # NB this is only going to work for experimental data right now
    new.kernel <- exp.mean(new.features, gp$relational.features[[1]])
    K <- fast.mat(rbfdot(), gp$X, matrix(new.X)) 
    K <- K + new.kernel
    return(K %*% gp$alpha)
}

predict.gp <- function(gp, new.X, kernel=rbfdot()) {
    K <- fast.mat(kernel, gp$X, matrix(new.X))
    if(length(gp$add.kernels) > 0) {
        K <- K + gp$add.kernels
    }
    return(K %*% gp$alpha)
}

fit.relational.gp <- function(X, y, kernel=rbfdot(), relational.features=list(), var=1.0, num.relational.features=1) {
    # learn the relational kernels first
    relational.kernels <- list()
    if(num.relational.features > 0) {
            for(i in 1:num.relational.features) {
                relational.kernels[[i]] <- exp.mean(relational.features[[i]])
            }
    }
    print(length(relational.kernels))
    gp <- fit.gp(X, y, kernel=kernel, add.kernels=relational.kernels, var=var)
    gp$relational.features <- relational.features
    return(gp) 
}
