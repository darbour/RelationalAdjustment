#!/usr/bin/env Rscript

#commandArgs <- function(trailingOnly=TRUE) { return(c(10, "~/repos/RelationalICausalInference/experiments/all_configurations.csv", "temp.csv"))} 

source("3-net.R")
source("experimental-estimators.R")
source("generate_data.R")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) {
  stop("Usage: run_instance.R <config ID> <all config file> <output file>")
}

config_id <- as.numeric(args[1])
config_file <- args[2]
output_file <- args[3]

for(trial in 1:50) {
  cat(paste0("Running instance ", config_id, " trial ", trial, " from configuration file  ", config_file, "\n"))
  cat(paste0("Output file: ", output_file, "\n"))
  
  # generate data
  configs <- read.csv(config_file)
  config <- configs[config_id, ]
  random.seed <- config$random.seed * trial
  gendata <- generate.by.index(configs, config_id, random.seed=random.seed, noise.sd=1)
  
  # estimate effects using GP, linear models, and horvitz-thompson
  gp.model <- gp.estimate(gendata$adj.mat, gendata$data)
  gbm.model <- gbm.estimate(gendata$adj.mat, gendata$data)
  lm.model <- lam.I(gendata$adj.mat, gendata$data)
  lm.modelII <- lam.II(gendata$adj.mat, gendata$data)
  if(config$graph.cluster.randomization) {
    htmeans <- ugander.horvitz.thompson(gendata$adj.mat, gendata$data, gendata$clusters, 0.75)
  } else {
    htmeans <- list(tmean=NA, cmean=NA)
  }
  
  estimates <- data.frame(method=c("Actual", "GP", "GBM", "LM-IND", "LM-INT", "HT"), config=config_id, trial=trial)
  
  hyp.friends.values <- seq(0, 1, length.out=11)
  for(hyp.friends in hyp.friends.values) {
    
    gp.est.outcome <- gp.model(friendt=hyp.friends)
    gbm.est.outcome <- gbm.model(friendt=hyp.friends)
    lmI.est.outcome <- lm.model(friendt=hyp.friends)
    lmII.est.outcome <- lm.modelII(friendt=hyp.friends)
    
    act.outcome <- gendata$outcome.function(friendt=hyp.friends, sdnoise=0)
    
    # NA for Horvitz-Thompson which only estimates 'global' effects
    estimates[, paste0("mf_", hyp.friends)] <- c(act.outcome, gp.est.outcome, gbm.est.outcome, lmI.est.outcome, lmII.est.outcome, NA)
  }
  estimates$mt_0 <- c(gendata$outcome.function(myt=0, sdnoise=0), gp.model(myt=0), gbm.model(myt=0), lm.model(myt=0), lm.modelII(myt=0), NA)
  estimates$mt_1 <- c(gendata$outcome.function(myt=1, sdnoise=0), gp.model(myt=1), gbm.model(myt=1), lm.model(myt=1), lm.modelII(myt=1), NA)
  estimates$global_treatment <- c(gendata$outcome.function(myt=1, friendt=1, sdnoise=0), gp.model(myt=1, friendt=1), gbm.model(myt=1, friendt=1),
                                    lm.model(myt=1, friendt=1), lm.modelII(myt=1, friendt=1), htmeans$tmean)
  estimates$global_control <- c(gendata$outcome.function(myt=0, friendt=0, sdnoise=0), gp.model(myt=0, friendt=0), gbm.model(myt=1, friendt=1),
                                    lm.model(myt=0, friendt=0), lm.modelII(myt=0, friendt=0), htmeans$cmean)
  
  # synchronize access to the results file
  system(paste0("lockfile results.lock"))
  append <- FALSE
  if(file.exists(output_file) && file.info(output_file)$size > 0) {
    append <- TRUE
  }
  write.table(estimates, output_file, append=append, col.names=!append, row.names=FALSE, sep=",")
  file.remove("results.lock")
  
}
