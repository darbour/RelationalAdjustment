rks.estimator <- function(x,y,s=1/6,f=sin, k=20, alpha=0.001) {
    x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
    transform.mat <- matrix(rnorm(ncol(x)*k),ncol(x))
    x <- s/ncol(x)*x%*%transform.mat
    x <- f(x)
    print(dim(x))
    W <- solve(t(x) %*% x + alpha * diag(k)) %*% t(x) %*% y
    obj <- list(weights=W, transform.mat=transform.mat, s=s, f=f)
    class(obj) <- "rks"
    return(obj)
}

predict.rks <- function(obj, newdata=NULL) {
    x <- as.matrix(newdata)
    x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
    x <- obj$s / ncol(x) * x %*% obj$transform.mat
    x <- obj$f(x)
    return(x %*% obj$weights)
}

