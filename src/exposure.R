prop.exposed <- function(adj.row, treatments) {
    sum(adj.row * treatments) / sum(adj.row)
}

calculate.exposure <- function(treatments, adj.mat, exposure.function) {
    apply(adj.mat, 1, exposure.function, treatments=treatments)
}

calculate.exposure.dists <- function(treatment.matrix, adj.mat, exposure.function) {
    apply(treatment.matrix, 1, calculate.exposure, adj.mat=adj.mat, exposure.function=exposure.function)
}
