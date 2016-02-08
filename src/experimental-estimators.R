ugander.sample.mean <- function(adj.mat, data, clusters, k, exposure=compute.ugander.exposure.prob, treatment.prob=0.5) {
    treat.probs <- exposure(adj.mat, clusters, k)
    control.probs <- exposure(adj.mat, clusters, k, control=TRUE)
    return((data$y*data$t / probs) / sum(data$t == 1) - (data$y*(1 - data$t) / probs) / sum(data$t == 0))
}


compute.ugander.exposure.prob <- function(adj.mat, clusters, cluster.treatment.prob, k, control=FALSE) {
    cluster.counts <- matrix(0, nrow(adj.mat), max(clusters))
    for(i in 1:nrow(adj.mat)) {
        row <- adj.mat[i, ]
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
        # compute the number you'd need
        min.treated <- k #ceiling(k * sum(adj.mat[i,]))
        if(control) {
            return((1 - cluster.treatment.prob) * (1 - compute.prob(length(ind.cluster.count), sum(adj.mat[i,]) - min.treated + 1, p, ind.cluster.count)))
        } else {
            print(paste("Min treated", min.treated, "count", intra.cluster.count))
            return(cluster.treatment.prob * compute.prob(length(ind.cluster.count), min.treated - intra.cluster.count, p, ind.cluster.count))
        }
        return(1)
    } )
    return(probs)
}

compute.prob <- function(s, T, p, w) {
    print(s)
    if(s == 0) {
        stop("what the heck")
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
    regression <- lm(o ~ t + frac.treated, data=reg.df)
    return(sum(regression$coefficients[-1]))
}

lam.II <- function(adj.mat, data) {
    degrees <- apply(adj.mat, 1, sum)
    reg.df <- data
    # fraction of treated friends
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
    regression.t <- lm(o ~ frac.treated, data=reg.df[which(reg.df$t == 1),])
    regression.c <- lm(o ~ frac.treated, data=reg.df[which(reg.df$t == 0),])
    return(regression.t$coefficients[1] + regression.t$coefficients[2] - regression.c$coefficients[1])
}

gp.estimate <- function(adj.mat, data) {
    require(kernlab)
    degrees <- apply(adj.mat, 1, ,sum)
    reg.df <- data
    reg.df$frac.treated <- (adj.mat %*% data$t) / degrees
    gp <- gausspr(o ~ ., data=reg.df)
    treatment.vals <- predict(gp, data.frame(t=1, frac.treated=1))
    control.vals <- predict(gp, data.frame(t=0, frac.treated=0))
    return(mean(treatment.vals) - mean(control.vals))
}
