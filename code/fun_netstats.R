# Output multiple network stats by running a single function

if (FALSE) {
    G <- graphGranger
    fast = TRUE
}

require(igraph)

netstats <- function(G,
                     fast = TRUE)
{
    directed <- igraph::is_directed(G)

    # Computing centrality measures for each vertex
    indegree <- igraph::degree(G, mode = "in")
    if (directed) {
        outdegree <- igraph::degree(G, mode = "out")
    } else {
        outdegree <- FALSE
    }
    if (!fast) {
        closeness <- igraph::closeness(G, mode = "total")
        betweenness <- igraph::betweenness(G, normalized = TRUE)
        triads <- igraph::triad_census(G)
    } else {
        closeness <- FALSE
        betweenness <- FALSE
        triads <- FALSE
    }
    components <- igraph::components(G, mode = "strong")

    # Graph-level variables
    G_order <- igraph::vcount(G)
    G_size <- igraph::ecount(G)
    G_density <- igraph::edge_density(G)
    G_reciprocity <- igraph::reciprocity(G)
    G_centralization <- igraph::centr_betw(G)$centralization
    G_pathLen <- igraph::mean_distance(G)

    G_components <- igraph::count_components(G, mode = "strong")

    # Graph-level variables based on the vertex stats
    G_isolated <- sum(indegree == 0L)

    list(indegree = indegree,
         outdegree = outdegree,
         closeness = closeness,
         betweenness = betweenness,
         triads = triads,
         components = components,
         G_order = G_order,
         G_size = G_size,
         G_density = G_density,
         G_reciprocity = G_reciprocity,
         G_centralization = G_centralization,
         G_pathLen = G_pathLen,
         G_components = G_components,
         G_isolated = G_isolated)
}
