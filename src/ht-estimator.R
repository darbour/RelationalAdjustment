
ht <- function(exposure, outcome, eprob.fun) {
  subject.id <- 1:length(exposure)
  splits <- split(subject.id, factor(exposure))
  
  exposure.outcome <- list()
  exposure.variance <- list()
  for(cur.exposure in names(splits)) {
    subjects <- splits[[cur.exposure]]
    cur.outcomes <- outcome[subjects]
    
    subject.probs <- aaply(subjects, 1, function(s) eprob.fun(s, s, cur.exposure))
    exposure.outcome[[cur.exposure]] <- 1/length(exposure) * sum(cur.outcomes / subject.probs)
    
    pairwise.probs <- outer(subjects, subjects, Vectorize(function(i, j) eprob.fun(i, j, cur.exposure)))
    corr.comp <- (1 - (subject.probs %*% t(subject.probs) / pairwise.probs)) * ((cur.outcomes / subject.probs) %*% t(cur.outcomes / subject.probs))
    diag(corr.comp) <- 0            
    exposure.variance[[cur.exposure]] <- 1/(length(exposure)**2) * 
        (sum((1 - subject.probs) * (cur.outcomes / subject.probs) ** 2) + sum(corr.comp))
  }
  return(list("mean"=exposure.outcome, "variance"=exposure.variance))
}

n <- 1000
exposure <- sample(c("a", "b"), size=n, replace=TRUE)

eprob.fun <- function(i, j, e) {
  if(i == j) {
    return(0.5)
  } else {
    return(0.25)
  }
}
outcome <- ifelse(exposure=="a", rnorm(n), rnorm(n, mean=1, sd=10))
print(ht(exposure, outcome, eprob.fun))
