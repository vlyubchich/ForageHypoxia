# Get information from in-connections of a directed network
if (FALSE) {
    seed = 1
    n.wave = 3
}

AvgImpactOnOneNode <- function(AM, Alag, seed) {

}


influencers <- function(AM, Alag, seed, n.wave = 2) {
    if (n.wave < 1) {
        stop("Number of waves, n.wave, should be >= 1.")
    }
    w <- 1
    n.used <- 1
    waves <- as.list(c(seed, rep(NA, n.wave)))
    while (w <= n.wave & nrow(AM) >= n.used) {
        # Select influencers of the previous wave
        sapply(waves[[1]], function(x) AvgImpactOnOneNode(AM, Alag, x))


        waves[[w]] <- unlist(lapply(waves[[w - 1]], function(nd) which(Alag[, nd] > 0)))
        lags <- unlist(lapply(waves[[w]]$ids, function(nd) Alag[which(Alag[, nd] > 0), nd]))



    }





    effEdges <- net$edges
    nodes.waves <- as.list(c(seed, rep(NA, n.wave)))
    wave <- 1
    while (wave <= n.wave & nrow(effEdges) >= 0) {
        tmp <- is.element(effEdges, nodes.waves[[wave]])
        if (any(tmp)) {
            tmp <- which(matrix(tmp, dim(effEdges)[1], 2), arr.ind = TRUE)
            nodes.waves[[wave + 1]] <- sort(effEdges[cbind(tmp[,
                                                               1], sapply(tmp[, 2], FUN = switch, 2, 1))])
            effEdges <- effEdges[-tmp[, 1], ]
            if (is.vector(effEdges)) {
                effEdges <- t(effEdges)
            }
        }
        wave <- wave + 1
    }
    nodes.waves
}
