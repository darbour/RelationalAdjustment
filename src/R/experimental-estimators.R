ugander.sample.mean <- function(adj.mat, data, clusters, k, exposure=compute.ugander.exposure.prob, treatment.prob=0.5) {
    treat.probs <- exposure(adj.mat, clusters, treatment.prob, k)
    control.probs <- exposure(adj.mat, clusters, treatment.prob, k, control=TRUE)
    return( 1/nrow(data) * ((data$y*data$t / treat.probs) - (data$y*(1 - data$t) / control.probs)))
}


compute.ugander.exposure.prob <- function(adj.mat, clusters, cluster.treatment.prob, k, control=FALSE) {
    cluster.counts <- matrix(0, nrow(adj.mat), max(clusters))
    for(i in 1:nrow(adj.mat)) {
        row <- adj.mat[i, ]
        row[i] <- 1
        vals <- (row * clusters)
        for(cluster in unique(vals)[-1]) {
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
        min.treated <- ceiling(k * sum(adj.mat[i,]))
        if(control) {
            return((1 - cluster.treatment.prob) * (1 - compute.prob(length(ind.cluster.count), sum(adj.mat[i,]) - min.treated + 1, cluster.treatment.prob, ind.cluster.count)))
        } else {
            dprob <- compute.prob(length(ind.cluster.count), min.treated - intra.cluster.count, cluster.treatment.prob, ind.cluster.count)
            #cat("dprob:", dprob, ", ", cluster.treatment.prob, "\n")
            return(cluster.treatment.prob * dprob)
        }
        return(1)
    } )
    return(probs)
}

compute.prob <- function(s, T, p, w) {
    if(s == 0) {
        return(1)
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
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
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
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
    
    regmodel <- lm(o ~ frac.treated, data=reg.df)

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

gp.estimate <- function(adj.mat, data, sigma=0.5) {
    require(kernlab)
    degrees <- apply(adj.mat, 1 ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
    gp <- gausspr(o ~ t + frac.treated, data=reg.df, kpar=list(sigma=sigma))

    hyp.dat <- data.frame(reg.df)
    mean.po <- function(myt=NULL, friendt=NULL) {
        if(!is.null(myt)) {
            hyp.dat$t <- myt
        }
        if(!is.null(friendt)) {
            hyp.dat$frac.treated <- friendt
        }
        return(mean(predict(gp, newdata=hyp.dat)))
    }
    
    return(mean.po)
}
