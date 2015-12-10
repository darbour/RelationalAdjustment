library(plyr)
library(ggplot2)

source("exposure.R")
source("ht-estimator.R")

adjmat <- as.matrix(read.csv("network1000.csv", header=FALSE))
tsamples <- as.matrix(read.csv("samples1000.csv", header=FALSE))
data <- read.csv("network_attrs1000.csv")

# TODO: need exposure model that considers intrinsic treatment status. 
#       currently only looking at peers (marginal effect?)

# figure out what exposure category we actually observed for each subject
obs.t <- t(data[, "t"])
obs.exposure <- calculate.exposure.dists(obs.t, adjmat, prop.exposed)
obs.exposure <- cut(obs.exposure,
    breaks=seq(0, 1, length.out=11), include.lowest=TRUE)

# figure out what exposure category each sample corresponds to for each subject
exposure.dists <- calculate.exposure.dists(tsamples, adjmat, prop.exposed)

discrete.exposure <- cut(c(exposure.dists), 
    breaks=seq(0, 1, length.out=11), include.lowest=TRUE)
exposure.cats <- unique(discrete.exposure)
exposure.mat <- matrix(discrete.exposure, ncol=ncol(exposure.dists))

# compute marginal exposure probabilities
eprobs <- matrix(0, nrow(adjmat), length(exposure.cats))
colnames(eprobs) <- exposure.cats
for(i in 1:nrow(exposure.mat)) {
    probs <- prop.table(table(exposure.mat[, i]))
    eprobs[i, names(probs)] <- probs
}

# TODO: need to compute pairwise probs 
# (tensor of size subjects x subjects x num exposures)
eprob.fun <- function(i, j, e) {
  if(i == j) {
    return(eprobs[i, e])
  } else {
    return(NA)
  }
}
estimates <- ht(obs.exposure, data[, "o"], eprob.fun)
print(estimates)
d <- data.frame(exposure=factor(names(estimates$mean), levels=names(estimates$mean)), outcome=unlist(estimates$mean))

# TODO: need variance estimates for confidence intervals
g <- ggplot(d, aes(x=exposure, y=outcome)) + 
    geom_point() + theme_bw(base_size=16) + labs(x="Exposure", y="Mean Outcome")
print(g)
png("exposure-trend.png", width=800, height=400)
print(g)
dev.off()


