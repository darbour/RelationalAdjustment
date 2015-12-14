library(igraph)
library(plyr)

generate.data <- function(nsubjects) {
    nsubjects <- 1000
    net <- sample_smallworld(1, nsubjects, 3, 0.1)
    adj.mat <- as.matrix(get.adjacency(net))

    c1 <- rnorm(nsubjects)
    c2 <- (adj.mat %*% c1) / rowSums(adj.mat) + rnorm(nsubjects)
    c3 <- (adj.mat %*% c2) / rowSums(adj.mat) + rnorm(nsubjects)
    t.prob <- 1 / (1 + exp(-scale(c1 + c2 + c3 + rnorm(nsubjects))))
    treatment <- rbinom(nsubjects, 1, prob=t.prob)

    burn.in <- 100
    iter <- 0
    while(iter <= burn.in) {
        iter <- iter + 1
        t.friends <- (adj.mat %*% treatment) / rowSums(adj.mat)
        t.prob <- 1 / (1 + exp(-(scale(c1 + c2 + c3 + rnorm(nsubjects)) + scale(t.friends)) ) )
        treatment <- rbinom(nsubjects, 1, prob=t.prob)
    }

    print(summary(glm(t ~ ., data=data.frame(c1, c2, c3, t=treatment, t.friends))))

    o <- c1 + c2 + c3 + treatment + t.friends + rnorm(nsubjects)

    #write.csv(data.frame(c1, c2, c3, t=treatment, o), file="network_attrs.csv", row.names = FALSE)
    #write.table(adj.mat, file="network.csv", row.names = FALSE, col.names=FALSE, sep=",")
    df <- data.frame(c1, c2, c3, t=treatment, o)
    return list(df, adj.mat)
}
