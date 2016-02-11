#!/usr/bin/env Rscript

#commandArgs <- function(trailingOnly=TRUE) { return(c(10, "~/repos/RelationalICausalInference/experiments/all_configurations.csv", "temp.csv"))} 

source("3-net.R")
source("estimators.R")
source("generate_data.R")

main <- function() {
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

      estimates <- run.one(config_file, config_id, trial)

      # synchronize access to the results file
      system(paste0("lockfile results.lock"))
      append <- FALSE
      if(file.exists(output_file) && file.info(output_file)$size > 0) {
        append <- TRUE
      }
      write.table(estimates, output_file, append=append, col.names=!append, row.names=FALSE, sep=",")
      file.remove("results.lock")
    }
}

run.one <- function(config_file, config_id, trial) {
  # generate data
  configs <- read.csv(config_file)
  config <- configs[config_id, ]
  random.seed <- config$random.seed * trial
  gendata <- generate.by.index(configs, config_id, random.seed=random.seed, noise.sd=1)
  
  # estimate effects using GP, linear models, and horvitz-thompson

  if(config$graph.cluster.randomization) {
      methods <- list("Actual"=function(junk1, junk2) gendata$outcome.function, 
                      "Exp-GBM"=gbm.estimate, 
                      "Exp-LM-IND"=lam.I, 
                      "Exp-LM-INT"=lam.II, 
                      "EXP-GP"=gp.estimate, 
                      "Exp-HT"=function(adj.mat, data) {
                        htmeans <- ugander.horvitz.thompson(adj.mat, data, gendata$clusters, 0.75)
                        return(function(myt=NULL, friendt=NULL) {
                          if(is.null(myt) || is.null(friendt)) {
                            return(NA)
                          } else if(myt == 1 & friendt == 1) {
                            return(htmeans$tmean)
                          } else if(myt == 0 & friendt == 0) {
                            return(htmeans$cmean)
                          } else {
                            return(NA)
                          }
                        })
                      })
  } else {
      methods <- list("Actual"=function(junk1, junk2) gendata$outcome.function, 
                      "Obs-GBM-Sufficient"=obs.gbm.sufficient,
                      #"Obs-GP-KME"=obs.gp.kme,
                      "Obs-LM-Simple"=obs.linear.simple, 
                      "Obs-LM-Sufficient"=obs.linear.sufficient, 
                      "Obs-GP-Sufficient"=obs.gp.sufficient

                      # include these as examples of unadjusted analyses
                      "Exp-GBM"=gbm.estimate, 
                      "Exp-LM-IND"=lam.I, 
                      "Exp-LM-INT"=lam.II) 
  }
  
  estimates <- data.frame(method=names(methods), config=config_id, trial=trial)
  # fit each of the models
  outcome.funcs <- alply(names(methods), 1, function(name) {
      cat("Fitting ", name, "\n")
      return(methods[[name]](gendata$adj.mat, gendata$data))
  })
  names(outcome.funcs) <- names(methods)

  cat("Getting estimates\n")
  # sweep across friend configurations if the design calls for it
    hyp.friends.values <- seq(0, 1, length.out=11)
    for(hyp.friends in hyp.friends.values) {
      est.outcomes <- aaply(names(methods), 1, function(name) {
          return(outcome.funcs[[name]](friendt=hyp.friends))
      })
      estimates[, paste0("mf_", hyp.friends)] <- est.outcomes
    }
      
    estimates$mt_0 <- aaply(names(methods), 1, function(name) outcome.funcs[[name]](myt=0))
    estimates$mt_1 <- aaply(names(methods), 1, function(name) outcome.funcs[[name]](myt=1))

  estimates$global_treatment <- aaply(names(methods), 1, function(name) outcome.funcs[[name]](myt=1, friendt=1))
  estimates$global_control <- aaply(names(methods), 1, function(name) outcome.funcs[[name]](myt=0, friendt=0))
  
  return(estimates)
}

main()
