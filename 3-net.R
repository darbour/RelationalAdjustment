#Initially all vertices are unmarked.
#While there are unmarked vertices, in step j find an arbitrary
#unmarked vertex v, selecting v to be vertex vj and marking all
#vertices in B2(vj ).
#
#Suppose k such vertices are defined, and let S = {v1, v2, ..., vk}.
#
#For every vertex w of G, assign w to the closest vertex vi 2 S,
#breaking ties consistently (e.g. in order of lowest index).
#
#For every vj , let Cj be the set of all vertices assigned to vj .
three.net <- function(graph, p=0.5) {
    # get the shortest path distances
    dists <- distances(graph)
    N <- nrow(dists)
    B.2 <- dists <= 2
    unmarked <- rep(TRUE, N)
    seeds <- c()
    while(!any(unmarked)) {
        # randomly choose from the unmarked nodes
        new.seed <- sample(which(unmarked), 1)
        unmarked[B.2[new.seed, ]] <- FALSE
        seeds <- c(seeds, new.seed)
    }
    
    # now assign to clusters
    clusters <- apply(dists, 1, function(row) { seeds[which.min(row[seeds])] }) 
    treatment.assignments <- rbinom(length(unique(clusters)), 1, p)
    return(list(clusters=clusters, treatment.assignmets=treatment.assignments))
}

