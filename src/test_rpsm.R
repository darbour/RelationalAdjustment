library(class)
library(optmatch)
library(ranger)
library(plyr)
library(caret)
library(e1071)
library(kernlab)

dat <- read.csv("data5k/data.csv")
adj_mat <- read.csv("data5k/network.csv")
adj_mat <- matrix(as.numeric(adj_mat == "true"), nrow=nrow(adj_mat))

dat$degree <- rowSums(adj_mat)
dat$degree2 <- rowSums(adj_mat %*% adj_mat)

dat$t_friends <- adj_mat %*% dat$t / dat$degree
dat$num_t_friends <- adj_mat %*% dat$t

# family binomial is cheating, but it's fine for now
prop.model <- glm(t ~ t_friends + c1 + c2 + c3, data=dat, family="binomial")
print(summary(prop.model))

dist <- match_on(prop.model)
dat$groupid <- fullmatch(dist, data=dat)
print(summary(lm(o ~ t + groupid, data=dat)))

adj2 <- adj_mat %*% adj_mat
adj2 <- ifelse(adj2 > 1, 1, adj2)
dat$f_c1 <- adj_mat %*% dat$c1 / dat$degree
dat$f_c2 <- adj_mat %*% dat$c2 / dat$degree
dat$f_c3 <- adj_mat %*% dat$c3 / dat$degree

ff_features <- c("c1", "c2", "c3", "t")
friends.friends.aggs <- adply(1:nrow(dat), 1, function(i) {
  d <- adply(which(adj_mat[i, ] == 1), 1, function(j) {
    nodes <- which(adj_mat[j, ] == 1)
    ffdat <- list()
    for(feature in ff_features) {
      ffdat[[paste0("ff_", feature, "_mean")]] <- mean(dat[nodes, feature])
      ffdat[[paste0("ff_", feature, "_mean")]] <- sd(dat[nodes, feature])
    }
    ffdat[["p_treat"]] <- predict(prop.model, dat[j, ], type="response")
    return(as.data.frame(ffdat))
  })
  return(colMeans(d[, 2:ncol(d)]))
})
friends.friends.aggs <- subset(friends.friends.aggs, select=-X1)
dat <- cbind(dat, friends.friends.aggs)

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Conditional Density Estimation as in Hirano & Imbens' GPS paper\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

laplace.density <- function(scale, location, x) {
  if(min(scale) <= 0) {
    warning(paste0(sum(scale <= 0),  " Laplace scale values are less than zero. That's weird"))
  }
  scale <- ifelse(scale <= 0, 0.001, scale)
  return(1 / (2 * scale) * exp(-abs(x - location) / scale))
}
gaussian.density <- function(scale, location, x) {
  1 / (scale * sqrt(2 * pi)) * exp(-(x - location)^2/(2 * scale**2))
}
model <- "gp"
if(model =="svm") {
  friends.model <- svm(t_friends ~ f_c1 + f_c2 + f_c3 + ff_c1_mean + ff_c2_mean + ff_c3_mean + ff_t_mean + degree + degree2, data=dat, probability=TRUE, C = 0.25)
  friends.scale <- predict(friends.model, type="probability")
  densfun <- laplace.density
} else if(model == "gp") {
  friends.model <- gausspr(t_friends ~ f_c1 + f_c2 + f_c3 + ff_c1_mean + ff_c2_mean + ff_c3_mean + ff_t_mean + degree + degree2 + p_treat, data=dat, variance.model=TRUE)
  friends.scale <- predict(friends.model, dat, type="sdeviation")
  densfun <- gaussian.density
}
friends.density <- densfun(friends.scale, dat$t_friends, predict(friends.model, dat, type="response"))
gp.model <- gausspr(o ~ t_friends + friends.density, data=dat)

# estimate response values at hypothetical values of friend proportion (dose)
dose.response <- NULL
for(friend.prop in unique(dat$t_friends)) {
  density <- densfun(friends.scale, dat$t_friends, friend.prop)
  
  resp <- mean(predict(gp.model, data.frame(t_friends=friend.prop, friends.density=density)))
  actual <- mean(dat$t + friend.prop + dat$c1 + dat$c2 + dat$c3)
  dose.response <- rbind(dose.response, data.frame(prop=friend.prop, resp=resp, actual=actual))
}

# plot true dose-response in red, estimated dose-response in black
m <- max(dose.response$resp)
g <- ggplot(dose.response, aes(x=prop)) + geom_point(aes(y=resp)) + geom_line(aes(y=resp)) + 
  geom_density(data=data.frame(x=dat$t_friends), aes(x=x, y=..scaled..), alpha=0.1, fill="gray") +
  geom_point(aes(y=actual), color="red") + geom_line(aes(y=actual), color="red") + theme_bw()
print(g)

cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Conditional Density Weighting as in marginal structural models\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

marginal.density <- density(dat$t_friends)
friends.model <- gausspr(t_friends ~ f_c1 + f_c2 + f_c3 + ff_c1_mean + ff_c2_mean + ff_c3_mean + ff_t_mean + degree + degree2 + p_treat, data=dat, variance.model=TRUE)
friends.scale <- predict(friends.model, dat, type="sdeviation")
friends.density <- laplace.density(friends.scale, dat$t_friends, predict(friends.model, dat, type="response"))
weight <- gausspr(t_friends)/densfun(friends.scale, dat$t_friends, predict(friends.model, dat, type="response"))
gp.model <- gausspr(o ~ t_friends + friends.density, data=dat, kernel="laplacedot")

# estimate response values at hypothetical values of friend proportion (dose)
dose.response <- NULL
for(friend.prop in unique(dat$t_friends)) {
  density <- densfun(friends.scale, dat$t_friends, friend.prop)
  
  resp <- mean(predict(gp.model, data.frame(t_friends=friend.prop, friends.density=density)))
  actual <- mean(dat$t + friend.prop + dat$c1 + dat$c2 + dat$c3)
  dose.response <- rbind(dose.response, data.frame(prop=friend.prop, resp=resp, actual=actual))
}

# plot true dose-response in red, estimated dose-response in black
m <- max(dose.response$resp)
g <- ggplot(dose.response, aes(x=prop)) + geom_point(aes(y=resp)) + geom_line(aes(y=resp)) + 
  geom_density(data=data.frame(x=dat$t_friends), aes(x=x, y=..scaled..), alpha=0.1, fill="gray") +
  geom_point(aes(y=actual), color="red") + geom_line(aes(y=actual), color="red") + theme_bw()
print(g)


cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Network sampling\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

samples <- data.matrix(read.csv("prop_results/samples.csv"))
subject.prop.samples <- alply(1:nrow(dat), 1, function(i) {
  return((samples %*% (adj_mat[i, ] * dat$t)) / dat$degree[i])
})
dat$t_friends_density <- aaply(1:nrow(dat), 1, function(i) {
  props <- subject.prop.samples[[i]]
  sum(props == dat$t_friends[i]) / length(props)
})
gp.model <- gausspr(o ~ t_friends + t_friends_density, data=dat)

# estimate response values at hypothetical values of friend proportion (dose)
dose.response2 <- NULL
for(friend.prop in unique(dat$t_friends)) {
  density <- aaply(1:nrow(dat), 1, function(i) {
    props <- subject.prop.samples[[i]]
    sum(props == friend.prop) / length(props)
  })
  resp <- mean(predict(gp.model, data.frame(t_friends=friend.prop, t_friends_density=density)))
  actual <- mean(dat$t + friend.prop + dat$c1 + dat$c2 + dat$c3)
  dose.response2 <- rbind(dose.response2, data.frame(prop=friend.prop, resp=resp, actual=actual))
}

# plot true dose-response in red, estimated dose-response in black
m <- max(dose.response2$resp)
g <- ggplot(dose.response2, aes(x=prop)) + geom_point(aes(y=resp)) + geom_line(aes(y=resp)) + 
  geom_point(aes(y=actual), color="red") + geom_line(aes(y=actual), color="red") + 
  geom_density(data=data.frame(x=friend.prop), aes(x=x, y=..scaled..), alpha=0.1) + theme_bw()
print(g)
