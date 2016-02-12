library(plyr)


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

get.po.func <- function(model, data, ...) {
    hyp.dat <- data.frame(data)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            hyp.dat$frac.treated <- friendt
        }
        return(mean(predict(model, newdata=hyp.dat, ...)))
    }
    return(mean.po)
}

obs.linear.simple <- function(adj.mat, data) {
    degrees <- rowSums(adj.mat)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    reg.df$mean.fc1 <- as.numeric((adj.mat %*% data$t) / degrees)
    reg.df$mean.fc2 <- as.numeric((adj.mat %*% data$t) / degrees)

    if(nrow(adj.mat) > 1024) {
        require(biglm)
        reg.model <- lm(o ~ t + c1 + c2 + frac.treated + mean.fc1 + mean.fc2, data=reg.df)
    } else {
        reg.model <- lm(o ~ t + c1 + c2 + frac.treated + mean.fc1 + mean.fc2, data=reg.df)
    }
    return(get.po.func(reg.model, reg.df))
}

obs.linear.sufficient <- function(adj.mat, data) {
    degrees <- rowSums(adj.mat)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    reg.df$mean.fc1 <- as.numeric((adj.mat %*% data$c1) / degrees)
    reg.df$mean.fc2 <- as.numeric((adj.mat %*% data$c2) / degrees)
    reg.df$var.fc1 <- as.vector((adj.mat %*% data$c1^2) - reg.df$mean.fc1^2) / degrees
    reg.df$var.fc2 <- as.vector((adj.mat %*% data$c2^2) - reg.df$mean.fc2^2) / degrees


    reg.model <- lm(o ~ t + c1 + c2 + frac.treated + mean.fc1 + mean.fc2 + var.fc2 + var.fc2+ mean.fc1:var.fc1 + mean.fc2:var.fc2, data=reg.df)
    return(get.po.func(reg.model, reg.df))
}

#obs.rks.sufficient <- function(adj.mat, data) {
#    degrees <- rwoSums(adj.mat)
#    reg.df <- data
#    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
#    reg.df$mean.fc1 <- as.numeric((adj.mat %*% data$c1) / degrees)
#    reg.df$mean.fc2 <- as.numeric((adj.mat %*% data$c2) / degrees)
#    reg.df$var.fc1 <- as.vector((adj.mat %*% data$c1^2) - reg.df$mean.fc1^2) / degrees
#    reg.df$var.fc1 <- as.vector((adj.mat %*% data$c2^2) - reg.df$mean.fc2^2) / degrees
#
#    
#    x <- cbind(apply(as.matrix(reg.df[,-c('o')]),2,function(u)rank(u)/length(u)),1)
#    reg.model <- lm(o ~ t + c1 + c2 + frac.treated + mean.fc1 + mean.fc2 + var.fc1 + var.fc2+ mean.fc1:var.fc1 + mean.fc2:var.fc1, data=reg.df)
#    return(get.po.func(reg.model, reg.df))
#}

obs.gp.sufficient <- function(adj.mat, data) {
    require(kernlab)
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    reg.df$mean.fc1 <- as.numeric((adj.mat %*% data$c1) / degrees)
    reg.df$mean.fc2 <- as.numeric((adj.mat %*% data$c2) / degrees)
    reg.df$var.fc1 <- as.vector((adj.mat %*% data$c1^2) - reg.df$mean.fc1^2) / degrees
    reg.df$var.fc2 <- as.vector((adj.mat %*% data$c2^2) - reg.df$mean.fc2^2) / degrees
    
    reg.model <- gausspr(o ~ t + c1 + c2 + frac.treated + mean.fc1 + mean.fc2 + var.fc1 + var.fc2, data=reg.df)
    return(get.po.func(reg.model, reg.df))
}

lam.I <- function(adj.mat, data) {
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)

    reg.model <- lm(o ~ t + frac.treated, data=reg.df)
    return(get.po.func(reg.model, reg.df))
}

lam.II <- function(adj.mat, data) {
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)

    regmodel <- lm(o ~ factor(t):frac.treated, data=reg.df)
    return(get.po.func(regmodel, reg.df))
}

obs.gbm.sufficient <- function(adj.mat, data) {
    require(gbm)
    degrees <- rowSums(adj.mat)
    reg.df <- data
    reg.df$frac.treated <- as.numeric((adj.mat %*% data$t) / degrees)
    print(names(data))
    reg.df$mean.fc1 <- as.numeric((adj.mat %*% data$c1) / degrees)
    reg.df$mean.fc2 <- as.numeric((adj.mat %*% data$c2) / degrees)
    reg.df$var.fc1 <- as.vector((adj.mat %*% data$c1^2) - reg.df$mean.fc1^2) / degrees
    
    reg.df$var.fc2 <- as.vector((adj.mat %*% data$c2^2) - reg.df$mean.fc2^2) / degrees
    
    if(nrow(adj.mat) > 1024) {
        model <- gbm(o ~ t + frac.treated + mean.fc1 + mean.fc2 + var.fc1 + var.fc2, data=reg.df[sample(1:nrow(adj.mat), 2048),], cv.folds=10, n.trees=2000, distribution="gaussian", n.cores=15)
        opt.iter <- gbm.perf(model, plot.it=FALSE)

        model <- gbm(o ~ t + frac.treated + mean.fc1 + mean.fc2 + var.fc1 + var.fc2, data=reg.df[sample(1:nrow(adj.mat), 2048),], n.trees=opt.iter, distribution="gaussian")
    } else {
        model <- gbm(o ~ t + frac.treated + c1 + c2 + mean.fc1 + mean.fc2 + var.fc1 + var.fc1, data=reg.df, cv.folds=10, n.trees=2000, distribution="gaussian")
        opt.iter <- gbm.perf(model, plot.it=FALSE)

        model <- gbm(o ~ t + frac.treated + c1 + c2 + mean.fc1 + mean.fc2 + var.fc1 + var.fc1, data=reg.df, n.trees=opt.iter, distribution="gaussian")
    }
    return(get.po.func(model, reg.df, n.trees=opt.iter))
}

gbm.estimate <- function(adj.mat, data) {
    require(gbm)
    degrees <- apply(adj.mat, 1 ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees

    model <- gbm(o ~ t + frac.treated, data=reg.df, cv.folds=10, n.trees=2000, distribution="gaussian")
    opt.iter <- gbm.perf(model, plot.it=FALSE)
    model <- gbm(o ~ t + frac.treated, data=reg.df, n.trees=opt.iter, distribution="gaussian")
    return(get.po.func(model, reg.df, n.trees=opt.iter))
}

obs.gp.kme <- function(adj.mat, data) {
    # TODO: call KME code
    require(kernlab)
    degrees <- apply(adj.mat, 1 ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees

    gp <- gausspr(o ~ t + frac.treated, data=reg.df)

    return(get.po.func(gp, reg.df))
}

gp.estimate <- function(adj.mat, data) {
    require(kernlab)
    degrees <- apply(adj.mat, 1 ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees

    gp <- gausspr(o ~ t + frac.treated, data=reg.df)

    return(get.po.func(gp, reg.df))
}
