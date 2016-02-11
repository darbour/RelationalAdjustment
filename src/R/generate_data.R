library(igraph)
library(plyr)

source("3-net.R")

# This function creates a collection of run configurations in a specified directory
create.configurations <- function(base.dir) {
  #sizes <- c(100, 400, 900)
  sizes <- c(1024)
  graph.settings <- expand.grid(graph.type=c("small-world"), degree=5, p=c(0.0, 0.01, 0.10, 0.50, 1.00), power=NA, size=sizes)
  #graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("erdos-renyi"), p=c(0.05, 0.1, 0.2), degree=NA, power=NA, size=sizes))
  graph.settings <- rbind(graph.settings, expand.grid(graph.type=c("barabasi-albert"), power=c(0.1, 0.5, 1), degree=NA, p=NA, size=sizes))
  
  experimental.function.settings <- expand.grid(exposure.type=c("linear", "sigmoid", "exponential", "rbf-friends"), 
                                   of.beta=c(0, 0.5, 1, 1.5), ot.beta=c(0, 0.5, 1, 1.5), 
                                   confounding.coeff=c(0), treatment.autocorr.coeff=0, graph.cluster.randomization=TRUE)
  observational.function.settings <- expand.grid(exposure.type=c("linear", "sigmoid", "exponential", "rbf-friends"), 
                                                of.beta=c(1), ot.beta=c(1), 
                                                confounding.coeff=c(0, 1, 2, 3), treatment.autocorr.coeff=0, graph.cluster.randomization=FALSE)
  observational.function.settings <- rbind(observational.function.settings, 
                                           expand.grid(exposure.type=c("linear", "sigmoid", "exponential", "rbf-friends"), 
                                                 of.beta=c(1), ot.beta=c(1), 
                                                 confounding.coeff=c(0, 1), treatment.autocorr.coeff=c(1,2,5,10), graph.cluster.randomization=FALSE))
  function.settings <- rbind(experimental.function.settings, observational.function.settings)
  # exaggerate effects for some classes of models
  multipliers <- list("sigmoid"=10, "exponential"=10)
  for(type in names(multipliers)) {
    function.settings[function.settings$exposure.type == type, ]$of.beta <- function.settings[function.settings$exposure.type == type, ]$of.beta * multipliers[[type]]
    function.settings[function.settings$exposure.type == type, ]$ot.beta <- function.settings[function.settings$exposure.type == type, ]$ot.beta * multipliers[[type]]
  }
  
  all.settings <- merge(function.settings, graph.settings)
  all.settings$random.seed <- 1:nrow(all.settings)
  write.csv(all.settings, file.path(base.dir, "all_configurations.csv"))
}

# Generates data corresponding to a row in the configuration file
#   config.data -- A data frame like that output by 'create.configurations'
#   idx -- An integer row number corresponding to the instance of 'config.data' to run
generate.by.index <- function(config.data, idx, random.seed=NULL, verbose=FALSE, ...) {
  if(is.null(random.seed)) {
    random.seed <- idx
  }
  config <- config.data[idx, ]
  generate.data(nsubjects=config$size, random.seed=random.seed, 
                graph.type=config$graph.type, 
                graph.parameters=list(degree=config$degree, p=config$p, power=config$power), 
                exposure.type=config$exposure.type, 
                confounding.coeff=config$confounding.coeff,
                treatment.autocorr.coeff=config$treatment.autocorr.coeff,
                of.beta=config$of.beta,
                ot.beta=config$ot.beta,
                t.binary=TRUE,
                graph.cluster.randomization=config$graph.cluster.randomization,
                result.dir=".",
                basename=idx, ...)
}

# Generates data.
#   nsubjects -- Number of nodes in the network
#   random.seed -- Used for reproducible randomness
#   graph.type -- Graph generation algorithm, one of 'small-world', 'erdos-renyi', or 'barabasi-albert'. 
#                 Each network expects specific specific keys in 'graph.parameters'
#   graph.parameters -- Specific to the 'graph.type' specified. 
#                         For 'small-world', expects 'degree' and 'p'
#                         For 'erdos-renyi', expects 'p'
#                         For 'barabasi-albert', expects 'power'
#   exposure.type -- Defines the nature of the relationship between treatments and outcomes.
#                       linear -- Outcome is a linear combination of treatment, friend's proportion of treatment (fp), and confounders
#                       rbf-friends -- Outcome is linear in treatment and confounders, but is related to fp via a radial basis function
#                       sigmoid -- Outcome is related to a linear combination of treatment, fp, and confounders through a logistic link
#                       exponential -- Outcome is related to a linear combination of treatment, fp, and confounders through an exponential link
#   graph.cluster.randomization -- Indicates whether treatment will be assigned experimentally with graph cluster randomization, or observed passively.
generate.data <- function(nsubjects, random.seed, graph.type, graph.parameters, exposure.type, confounding.coeff, 
                          treatment.autocorr.coeff, of.beta, ot.beta, t.binary, graph.cluster.randomization, 
                          result.dir, basename, noise.sd=1, verbose=FALSE) {
    set.seed(random.seed)
    adjacency.file = FALSE
    if(graph.type == "small-world") {
      if(floor(sqrt(nsubjects))**2 != nsubjects) {
        stop(paste("Expected the number of subjects to be perfect square for small world network generation. Got ", nsubjects, ", floor(sqrt(nsubjects))**2: ", floor(sqrt(nsubjects))**2, "."))
      }
      net <- simplify(watts.strogatz.game(1, nsubjects, graph.parameters$degree, graph.parameters$p))
    } else if(graph.type == "erdos-renyi") {
      net <- erdos.renyi.game(nsubjects, graph.parameters$p)
    } else if(graph.type == "barabasi-albert") {
      net <- barabasi.game(nsubjects, power=graph.parameters$power, directed=FALSE)
    } else {
      require(Matrix)
      adj.df <- read.table(as.character(graph.type), sep='\t', skip = 4)
      adj.mat <- sparseMatrix(i=adj.df$V1, j=adj.df$V2, symmetric=TRUE, index1=FALSE, dims=c(max(adj.df) + 1, max(adj.df) + 1))
      nsubjects <- max(adj.df)+1
      adjacency.file = TRUE
    }
    
    if(!(adjacency.file)) { 
      adj.mat <- as.matrix(get.adjacency(net))
    }
    c1 <- rnorm(nsubjects)
    c2 <- rnorm(nsubjects)
    degrees <- rowSums(adj.mat)
    c1fmean <- (adj.mat %*% c1) / degrees + rnorm(nsubjects)
    c2fmean <- (adj.mat %*% c2) / degrees + rnorm(nsubjects)
    c1fvariance <- aaply(1:nsubjects, 1, function(i) var(c1[as.logical(adj.mat[1, ])]))
    c2fvariance <- aaply(1:nsubjects, 1, function(i) var(c2[as.logical(adj.mat[1, ])]))
    confounding.terms <- cbind(c1, c2, c1fmean, c2fmean, c1fvariance, c2fvariance, c1fmean * c1fvariance, c2fmean * c2fvariance)
    confounding.beta <- runif(ncol(confounding.terms))
    confounding.term <- scale(confounding.terms %*% confounding.beta)
    if(verbose) {
      plot(density(confounding.term))
    }
    
    if(graph.cluster.randomization) {
      clustering <- three.net(net)
      clusters <- clustering$clusters
      treatment <- clustering$treatment.assignments[clustering$clusters]
      t.friends <- adj.mat %*% treatment
      friend.prop <- t.friends / rowSums(adj.mat)      
    } else {
      clusters <- NULL
      
      # possibly generate autocorrelated treatments
      scale.friend.prop <- 0
      for(i in 1:3) { #three Gibbs steps
        if(t.binary) {
          t.prob <- 1 / (1 + exp(-(confounding.coeff * confounding.term  + treatment.autocorr.coeff * scale.friend.prop) + rnorm(nsubjects)))
          treatment <- rbinom(nsubjects, 1, prob=t.prob)      
        } else {
          treatment <- confounding.coeff * confounding.term + rnorm(nsubjects) +  treatment.autocorr.coeff * scale.friend.prop
        }      
        t.friends <- adj.mat %*% treatment
        friend.prop <- t.friends / rowSums(adj.mat)      
        scale.friend.prop <- scale(friend.prop)
      }      
    }
    
    if(exposure.type == "linear") {
      ofun <- function(treat, friend.pt, sdnoise=1) {
        ot.beta * treatment + of.beta * friend.pt + confounding.coeff * confounding.term + rnorm(nsubjects, sd=sdnoise)
      }
    } else if(exposure.type == "rbf-friends") {
      ofun <- function(treat, friend.pt, sdnoise=1) {
        ot.beta * treatment + of.beta * exp(-(friend.pt - 0.5)**2) + confounding.coeff * confounding.term + rnorm(nsubjects, sd=sdnoise)
      }
    } else if(exposure.type == "sigmoid") {
      ofun <- function(treat, friend.pt, sdnoise=1) {
        1 / (1 + exp(-(ot.beta * treat + of.beta * friend.pt + confounding.coeff * confounding.term + rnorm(nsubjects, sd=sdnoise))))
      }
    } else if(exposure.type == "exponential") {
      ofun <- function(treat, friend.pt, sdnoise=1) {
        exp(-(ot.beta * treat + of.beta * friend.pt + confounding.coeff * confounding.term + rnorm(nsubjects, sd=sdnoise)))
      }
    } else {
      stop(paste0("Unrecognized exposure model", exposure.type))
    }
    o <- ofun(treatment, friend.prop, sdnoise=noise.sd)
    
    # a function for mean potential outcome
    po.fun <- function(myt=NULL, friendt=NULL, sdnoise=0) {
      if(is.null(myt)) {
        myt <- treatment
      }
      if(is.null(friendt)) {
        friendt <- friend.prop
      }
      return(mean(ofun(myt, friendt, sdnoise=sdnoise)))
    }
    
    if(verbose) {
      plot(friend.prop, o)
      
      actual.dose.response <- data.frame()
      for(hypothetical.friend.prop in seq(0, 1, length.out=200)) {
        potential.o <- mean(po.fun(friendt=hypothetical.friend.prop, sdnoise=0))
        actual.dose.response <- rbind(actual.dose.response, data.frame(hypothetical.friend.prop, potential.o))
      }
      with(actual.dose.response, plot(hypothetical.friend.prop, potential.o, type="l"))
    }
    
    df <- data.frame(c1, c2, t=treatment, o)
    return(list(data=df, adj.mat=adj.mat, outcome.function=po.fun, clusters=clusters))
}

# This function creates a collection of run configurations in a specified directory
create.rw.configurations <- function(base.dir) {
  graph.settings <- data.frame(graph.type=c('../../data/com-lj.ungraph.txt', '../../data/email-Enron.txt', '../../data/roadNet-CA.txt', '../../data/web-Stanford.txt'))
  
  experimental.function.settings <- expand.grid(exposure.type=c("linear", "sigmoid","rbf-friends"), 
                                   of.beta=c(0, 2), ot.beta=c(0, 2), 
                                   confounding.coeff=c(0), treatment.autocorr.coeff=0, graph.cluster.randomization=TRUE)
  observational.function.settings <- expand.grid(exposure.type=c("linear", "sigmoid",  "rbf-friends"), 
                                                of.beta=c(1), ot.beta=c(1), 
                                                confounding.coeff=c(3), treatment.autocorr.coeff=0, graph.cluster.randomization=FALSE)
  observational.function.settings <- rbind(observational.function.settings, 
                                           expand.grid(exposure.type=c("linear", "sigmoid",  "rbf-friends"), 
                                                 of.beta=c(1), ot.beta=c(1), 
                                                 confounding.coeff=c(3), treatment.autocorr.coeff=c(1), graph.cluster.randomization=FALSE))
  function.settings <- rbind(experimental.function.settings, observational.function.settings)
  # exaggerate effects for some classes of models
  multipliers <- list("sigmoid"=10)
  for(type in names(multipliers)) {
    function.settings[function.settings$exposure.type == type, ]$of.beta <- function.settings[function.settings$exposure.type == type, ]$of.beta * multipliers[[type]]
    function.settings[function.settings$exposure.type == type, ]$ot.beta <- function.settings[function.settings$exposure.type == type, ]$ot.beta * multipliers[[type]]
  }
  
  all.settings <- merge(function.settings, graph.settings)
  all.settings$random.seed <- 1:nrow(all.settings)
  write.csv(all.settings, file.path(base.dir, "all_configurations_rw.csv"))
}

#create.configurations("~/repos/RelationalICausalInference/experiments/")
