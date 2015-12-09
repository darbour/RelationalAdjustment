library(igraph)
library(plyr)

nsubjects <- 1000
net <- sample_smallworld(1, nsubjects, 3, 0.1)
adj.mat <- as.matrix(get.adjacency(net))

c1 <- rnorm(nsubjects)
c2 <- (adj.mat %*% c1) / rowSums(adj.mat) + rnorm(nsubjects)
c3 <- (adj.mat %*% c2) / rowSums(adj.mat) + rnorm(nsubjects)
t.prob <- 1 / (1 + exp(-scale(c1 + c2 + c3 + rnorm(nsubjects))))
t <- rbinom(nsubjects, 1, prob=t.prob)
o <- c1 + c2 + c3 + t + rnorm(nsubjects)

write.csv(data.frame(c1, c2, c3, t, o), file="network_attrs.csv", row.names = FALSE)
write.table(adj.mat, file="network.csv", row.names = FALSE, col.names=FALSE, sep=",")
