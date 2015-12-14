
ht <- function(exposure, outcome, eprob.fun) {
  subject.id <- 1:length(exposure)
  splits <- split(subject.id, factor(exposure))
  
  exposure.outcome <- list()
  exposure.variance <- list()
  for(cur.exposure in names(splits)) {
    valid.subject <- aaply(1:length(exposure), 1, function(s) eprob.fun(s, s, cur.exposure) > 0)
    n.valid.subjects <- sum(valid.subject)
    subjects <- splits[[cur.exposure]]
    subjects <- subjects[subjects %in% (1:length(exposure))[valid.subject]]
    cur.outcomes <- outcome[subjects]
    
    subject.probs <- aaply(subjects, 1, function(s) eprob.fun(s, s, cur.exposure))
    exposure.outcome[[cur.exposure]] <- 1/n.valid.subjects * sum(cur.outcomes / subject.probs)
    
    pairwise.probs <- outer(subjects, subjects, Vectorize(function(i, j) eprob.fun(i, j, cur.exposure)))
    corr.comp <- ((pairwise.probs - subject.probs %*% t(subject.probs)) / pairwise.probs) * ((cur.outcomes / subject.probs) %*% t(cur.outcomes / subject.probs))
    diag(corr.comp) <- 0            
    exposure.variance[[cur.exposure]] <- 1/(n.valid.subjects**2) * 
        (sum((1 - subject.probs) * (cur.outcomes / subject.probs) ** 2) + 
           sum(corr.comp * is.finite(corr.comp), na.rm=TRUE))
    
    cat("variance", exposure.variance[[cur.exposure]], "\n")
    # need to correct the variance if there are zero-valued pairwise probabilities
    zero.vals <- which(pairwise.probs == 0, arr.ind=TRUE)
    if(nrow(zero.vals) > 0) {
      bias.correction <- aaply(zero.vals, 1, function(r) {
        row <- r["row"]
        col <- r["col"]
        return(((exposure[row] == cur.exposure) * outcome[row] ** 2) / (2 * subject.probs[[as.character(row)]]) + 
                 ((exposure[col] == cur.exposure) * outcome[col] ** 2) / (2 * subject.probs[[as.character(col)]]))
      })
      cat("Bias correction", sum(bias.correction), "\n")
      exposure.variance[[cur.exposure]] <- exposure.variance[[cur.exposure]] + sum(bias.correction)
    }
  }
  return(list("mean"=exposure.outcome, "variance"=exposure.variance))
}

