# k is a fraction of treated neighbors
ugander.horvitz.thompson <- function(adj.mat, data, clusters, k, exposure=compute.ugander.exposure.prob, treatment.prob=0.5) {
    treat.probs <- exposure(adj.mat, clusters, treatment.prob, k)
    control.probs <- exposure(adj.mat, clusters, treatment.prob, k, control=TRUE)
    degrees <- rowSums(adj.mat)
    
    exposed <- (degrees * k) < (adj.mat %*% data$t)
    nonexposed <- (degrees * (1 - k)) > (adj.mat %*% data$t)
    
    outcome.estimator <- function(indicator, outcome, prob) {
      1 / nrow(data) * sum((outcome / prob)[indicator])
    }
    
    tmean <- outcome.estimator(exposed, data$o, treat.probs)
    cmean <- outcome.estimator(nonexposed, data$o, treat.probs)
    
    return(list(tmean=tmean, cmean=cmean))
}


compute.ugander.exposure.prob <- function(adj.mat, clusters, cluster.treatment.prob, k, control=FALSE) {
    cluster.counts <- matrix(0, nrow(adj.mat), max(clusters))
    for(i in 1:nrow(adj.mat)) {
        row <- adj.mat[i, ]
        vals <- (row * clusters)
        for(cluster in unique(vals)) {
            if(cluster == 0) {
              next
            }
            cluster.counts[i, cluster] <- sum(vals == cluster)
        }
    } 
    probs <- sapply(1:nrow(adj.mat), function(i)  {
        # exclude your own cluster
        intra.cluster.count <- cluster.counts[i, clusters[i]]
        ind.cluster.count <- cluster.counts[i,-clusters[i]]
        ind.cluster.count <- ind.cluster.count[which(ind.cluster.count > 0)]
        #cat("intra cluster count: ", intra.cluster.count, "\n")
        #print(ind.cluster.count)

        # compute the number you'd need
        degree <- sum(adj.mat[i,])
        min.subjects <- k * sum(adj.mat[i,])
        if(length(ind.cluster.count) == 0) {
          if(control) {
            return((1 - cluster.treatment.prob) * (min.subjects <= degree))
          } else {
            return(cluster.treatment.prob * (min.subjects <= degree))
          }
        }
        
        if(control) {
            dprob <- (1 - compute.prob(length(ind.cluster.count), 
                                       degree - min.subjects + 1, cluster.treatment.prob, ind.cluster.count))
            return((1 - cluster.treatment.prob) * dprob)
        } else {
            dprob <- compute.prob(length(ind.cluster.count), min.subjects - intra.cluster.count, cluster.treatment.prob, ind.cluster.count)
            return(cluster.treatment.prob * dprob)
        }
    })
    return(probs)
}

compute.prob <- function(s, T, p, w) {
    if(s == 0) {
        stop("shouldn't get here")
    } else if(s == 1) {
        return(p * (T < w[s]))
    } else {
        return(p * compute.prob(s-1, T - w[s], p, w) +  (1 - p) * compute.prob(s-1, T, p, w))
    }
}

lam.I <- function(adj.mat, data) {
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    # fraction of treated friends
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    reg.model <- lm(o ~ t + frac.treated, data=reg.df)

    hyp.dat <- data.frame(reg.df)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            hyp.dat$frac.treated <- friendt
        }
        return(mean(predict(reg.model, newdata=hyp.dat)))
    }
    return(mean.po)
}

lam.II <- function(adj.mat, data) {
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    # fraction of treated friends
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    
    regmodel <- lm(o ~ t:frac.treated, data=reg.df)

    hyp.dat <- data.frame(reg.df)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            hyp.dat$frac.treated <- friendt
        }
        return(mean(predict(regmodel, newdata=hyp.dat)))
    }
    return(mean.po)
}

gbm.estimate <- function(adj.mat, data) {
    require(gbm)
    degrees <- apply(adj.mat, 1 ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
    model <- gbm(o ~ t + frac.treated, data=reg.df, cv.folds=10, n.trees=2000)
    opt.iter <- gbm.perf(model, plot.it=FALSE)
    model <- gbm(o ~ t + frac.treated, data=reg.df, cv.folds=10, n.trees=opt.iter)

    hyp.dat <- data.frame(reg.df)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            hyp.dat$frac.treated <- friendt
        }
        return(mean(predict(model, newdata=hyp.dat, n.trees=opt.iter)))
    }
    
    return(mean.po)
}

gp.estimate <- function(adj.mat, data,var=1, sigma=NULL) {
    require(kernlab)
    source('kernelmeankern.R') 
    degrees <- apply(adj.mat, 1 ,sum)
    relational.vars <- lapply(1:nrow(adj.mat), function(id) data$t[which(adj.mat[id,] == 1)])
    gp <- fit.relational.gp(matrix(data$t), matrix(data$o), relational.features=relational.vars)
    #reg.df <- data
    #reg.df$frac.treated <- as.vector((adj.mat %*% data$t) / degrees)
    #reg.df$t[which(reg.df$t == 0)] <- -1
    #if(is.null(sigma)) {
    #    gp <- gausspr(o ~ t + frac.treated, data=reg.df,var=var, kernel="laplacedot")
    #} else {
    #    gp <- gausspr(o ~ t + frac.treated, data=reg.df,var=var, kpar=list(sigma=sigma))
    #}

    hyp.dat <- data.frame(data)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            relational.features <- lapply(1:nrow(adj.mat), function(row) data$t[which(adj.mat[row,] == 1)])
            return(predict.prop.treatment.gp(gp, hyp.dat,  relational.features, friendt))
        } else {
            return(mean(predict(gp, hyp.dat)))
        }
    }
    
    return(mean.po)
}
