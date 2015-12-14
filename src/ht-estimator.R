
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
    corr.comp <- (1 - (subject.probs %*% t(subject.probs) / pairwise.probs)) * ((cur.outcomes / subject.probs) %*% t(cur.outcomes / subject.probs))
    diag(corr.comp) <- 0            
    exposure.variance[[cur.exposure]] <- 1/(n.valid.subjects**2) * 
        (sum((1 - subject.probs) * (cur.outcomes / subject.probs) ** 2) + sum(corr.comp))
  }
  return(list("mean"=exposure.outcome, "variance"=exposure.variance))
}

