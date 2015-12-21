
burn.in <- 1000
iter <- 0
while(iter <= burn.in) {
    iter <- iter + 1
    t.friends <- adj.mat %*% treatment
    for(node in 1:nrow(adj.mat)) {
        t.friends <- adj.mat[node,] * treatment / sum(adj.mat[node,])
    }
}
