# Get information from in-connections of a directed network

influencers <- function(AM, Alag, seed, n.wave = 3) {
    if (n.wave < 1) {
        stop("Number of waves, n.wave, should be >= 1.")
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
