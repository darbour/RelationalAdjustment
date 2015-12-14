library(plyr)
library(stringr)
library(ggplot2)
library(Matrix)

source("ht-estimator.R")

adjmat <- as.matrix(read.csv("network1000.csv", header=FALSE))
data <- read.csv("network_attrs1000e.csv")


exposures <- outer(c(0, 1), unique(cut(seq(0, 1, length.out=11), 
            breaks=seq(0, 1, length.out=11), include.lowest=TRUE)), function(a, b) paste(a, b, sep=","))
ordered.exposure.cats <- c(exposures[1, ], exposures[2, ])
obs.exposure <- ordered.exposure.cats[as.numeric(data[, "exposure"])]

# store pairwise exposure probs in a list indexed by exposure category
pairwise.exposures <- list()
exp.prob.files <- Sys.glob("network1000_exposure*.csv")
for(f in exp.prob.files) {
    # files are tagged by factor level, not name
    exposure.cat <- as.numeric(str_replace(f, "network.*_exposure([0-9]+).csv", "\\1"))
    exposure.cat.str <- ordered.exposure.cats[exposure.cat]
    dat <- read.csv(f)
    pairwise.exposures[[exposure.cat.str]] <- sparseMatrix(i=dat$row, j=dat$col, x=dat$prob, dims=c(1000, 1000))
}

eprob.fun <- function(i, j, e) {
  if(i == j) {
    return(pairwise.exposures[[e]][i, j])
  } else {
    return(pairwise.exposures[[e]][i, j])
  }
}
estimates <- ht(obs.exposure, data[, "o"], eprob.fun)
print(estimates)
d <- data.frame(exposure=factor(names(estimates$mean), levels=ordered.exposure.cats), 
                outcome=unlist(estimates$mean), variance=unlist(estimates$variance))

g <- ggplot(d, aes(x=exposure, y=outcome)) + geom_errorbar(aes(ymin=outcome-sqrt(variance), ymax=outcome+sqrt(variance))) +
    geom_point() + theme_bw(base_size=16) + labs(x="Exposure", y="Mean Outcome") + theme(axis.text.x=element_text(angle=45))
print(g)
png("exposure-trend.png", width=1200, height=400)
print(g)
dev.off()


