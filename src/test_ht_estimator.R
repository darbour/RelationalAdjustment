library(plyr)

source("exposure.R")
source("ht-estimator.R")

adjmat <- as.matrix(read.csv("network1000.csv", header=FALSE))
tsamples <- as.matrix(read.csv("samples1000.csv", header=FALSE))
data <- read.csv("network_attrs1000.csv")

obs.t <- t(data[, "t"])
obs.exposure <- calculate.exposure.dists(obs.t, adjmat, prop.exposed)
obs.exposure <- cut(obs.exposure,
    breaks=seq(0, 1, length.out=11), include.lowest=TRUE)
exposure.dists <- calculate.exposure.dists(tsamples, adjmat, prop.exposed)

discrete.exposure <- cut(c(exposure.dists), 
    breaks=seq(0, 1, length.out=11), include.lowest=TRUE)
exposure.cats <- unique(discrete.exposure)
exposure.mat <- matrix(discrete.exposure, ncol=ncol(exposure.dists))

eprobs <- matrix(0, nrow(adjmat), length(exposure.cats))
colnames(eprobs) <- exposure.cats
for(i in 1:nrow(exposure.mat)) {
    probs <- prop.table(table(exposure.mat[, i]))
    eprobs[i, names(probs)] <- probs
}

eprob.fun <- function(i, j, e) {
  if(i == j) {
    return(eprobs[i, e])
  } else {
    return(NA)
  }
}
print(ht(obs.exposure, data[, "o"], eprob.fun))


