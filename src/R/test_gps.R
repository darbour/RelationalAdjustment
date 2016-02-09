library(ggplot2)
library(ranger)
library(caret)
library(kernlab)

# Demonstration of generalized propensity score
# adapted from "The Propensity Score with Continuous Treatments", Hirano and Imbens
# http://scholar.harvard.edu/imbens/files/hir_07feb04.pdf


# treatment is a linear function of some covariates
nsubjects <- 1000
x <- matrix(runif(nsubjects*5), ncol=5)
coeffs <- runif(5)
beta0 <- runif(1)
treat <- rnorm(nsubjects, mean=x %*% coeffs + beta0)

# y is a linear function of covariates and treatment cubed
ycoeffs <- runif(5)
ybeta0 <- runif(1)
yfun <- function(ct) {
  as.numeric(10 * x %*% ycoeffs + ybeta0 + ct^3)
}
y <- yfun(treat)
plot(treat, y)

# fit a linear model for treatment, consistent with the generating process
model <- lm(treat~ x)
model.summ <- summary(model)
s2 <- sd(model.summ$residuals)

propensity <- function(tvals) {
    # estimate normal density
  return(1/sqrt(2 * pi * s2) * exp(-1/(2*s2) * (tvals - model$fitted)**2))
}

# fit conditional expectation of y
density <- propensity(t)
y.model <- gausspr(y ~ treat + density)
print(summary(y.model))

# estimate response values at hypothetical values of treatment (dose)
dose.response <- NULL
for(tvals in seq(min(t), max(t), length.out=30)) {
  prop <- propensity(tvals)
  resp <- mean(predict(y.model, data.frame(treat=tvals, density=prop))
  actual <- mean(yfun(tvals))
  dose.response <- rbind(dose.response, data.frame(treat=tvals, resp=resp, actual=actual))
}

# plot true dose-response in red, estimated dose-response in black
m <- max(dose.response$resp)
g <- ggplot(dose.response, aes(x=treat)) + geom_point(aes(y=resp)) + geom_line(aes(y=resp)) + 
  geom_point(aes(y=actual), color="red") + geom_line(aes(y=actual), color="red") + 
  geom_density(data=data.frame(x=treat), aes(x=x, y=..scaled..), alpha=0.1) + theme_bw()
print(g)
