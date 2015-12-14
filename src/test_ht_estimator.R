library(plyr)
library(ggplot2)

source("exposure.R")
source("ht-estimator.R")

adjmat <- as.matrix(read.csv("network1000.csv", header=FALSE))
tsamples <- as.matrix(read.csv("samples1000.csv", header=FALSE))
data <- read.csv("network_attrs1000.csv")

map.exposures <- function(neighborhood.prop, intrinsic.treatments) {
    neighborhood.prop.bin <- cut(neighborhood.prop, 
                breaks=seq(0, 1, length.out=11), include.lowest=TRUE)
    return(paste(intrinsic.treatments, neighborhood.prop.bin, sep=","))
}

# figure out what exposure category we actually observed for each subject
obs.t <- t(data[, "t"])
obs.exposure <- calculate.exposure.dists(obs.t, adjmat, prop.exposed)
obs.exposure <- map.exposures(obs.exposure, obs.t)

# figure out what exposure category each sample corresponds to for each subject
exposure.dists <- calculate.exposure.dists(tsamples, adjmat, prop.exposed)

# matrix -> vector for efficiency of 'cut'
discrete.exposure <- map.exposures(c(exposure.dists), c(tsamples)) 
exposure.cats <- unique(discrete.exposure)
ordered.exposure.cats <- exposure.cats[order(str_replace_all(exposure.cats, "\\[|\\(|\\]|\\)", ""))]

exposure.mat <- matrix(discrete.exposure, ncol=ncol(exposure.dists))

# compute marginal exposure probabilities
eprobs <- matrix(0, nrow(adjmat), length(exposure.cats))
colnames(eprobs) <- exposure.cats
for(i in 1:nrow(exposure.mat)) {
    probs <- prop.table(table(exposure.mat[i, ]))
    eprobs[i, names(probs)] <- probs
}

# compute joint exposure probabilities, 
# store in a list indexed by exposure category

# this is too slow:
#subjects <- 1:nrow(exposure.mat)
#pairwise.probs <- llply(exposure.cats, function(exposure.cat) {
#    outer(subjects, subjects, function(i, j) {
#        sum(exposure.mat[i, ] == exposure.cat && exposure.mat[j, ] == exposure.cat) / ncol(exposure.mat)
#    })
#})


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
d <- data.frame(exposure=factor(names(estimates$mean), levels=ordered.exposure.cats), outcome=unlist(estimates$mean))

# TODO: need variance estimates for confidence intervals
g <- ggplot(d, aes(x=exposure, y=outcome)) + 
    geom_point() + theme_bw(base_size=16) + labs(x="Exposure", y="Mean Outcome")
print(g)
png("exposure-trend.png", width=800, height=400)
print(g)
dev.off()


