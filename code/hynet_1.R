# Hypoxia network for each 1 year


# 1. Pairwise tests ----

## Example 1986 ----
# This subsection is an example file "hynet_1pair_1986.R"
# to be run on cluster (total 30 files, for 1986--2015)

# Hypoxia network for 1 year

# Packages ---
rm(list = ls())
library(data.table)
library(dplyr)
library(BigVAR)
library(parallel)

year = 1986
lagmax = 30L


# Data ---
rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1) %>%
    mutate(CellID = paste0("x", CellID)) %>%
    mutate(Atlantic = (LAT < 37.22 & LON > -75.96) |
               (LAT < 37.5 & LON > -75.8) |
               (LAT < 38.0 & LON > -75.6) |
               (LAT < 36.96 & LON > -75.99)
    )
AtlanticCells <- rca_cells %>%
    filter(Atlantic) %>%
    pull(CellID)
if (FALSE) {
    library(ggplot2)
    library(plotly)
    p <- ggplot(rca_cells, aes(x = LON, y = LAT, color = Atlantic)) +
        geom_point()
    ggplotly(p)
}

D <- data.table::fread(paste0("./data_rca/rca_ts_", year, "_2.csv"),
                       select = c("DOAVEG_avg", "CellID", "Date")) %>%
    mutate(CellID = paste0("x", CellID)) %>%
    filter(!is.element(CellID, AtlanticCells))
# head(D)
Dwide <- data.table::dcast(D, Date ~ CellID, value.var = "DOAVEG_avg")
Dwide_mat <- as.matrix(Dwide[, -1])
Cells <- colnames(Dwide_mat)
# dim(Dwide_mat)


# Remove seasonality ---
EXO <- Dwide[, 1] %>%
    mutate(day_year = as.numeric(format(Date, "%j"))) %>%
    mutate(sin2 = sin(2 * pi * day_year/365.25),
           cos2 = cos(2 * pi * day_year/365.25)) %>%
    dplyr::select(sin2, cos2) %>%
    as.matrix()
# x = Dwide_mat[, 3000]; plot.ts(x)
Dwide_mat_deseas <- apply(Dwide_mat, 2, function(x) {
    m0 = lm(x ~ EXO)
    x0 = residuals(m0)
    # Replace outliers with interpolated values.
    iqr = IQR(x0)
    isoutlier = abs(x0) > 5 * iqr
    if (sum(isoutlier) > 0) {
        # print(sum(isoutlier))
        x[isoutlier] = NA
        x = zoo::na.spline(x)
        m0 = lm(x ~ EXO)
    }
    residuals(m0)
    # plot.ts(residuals(m0))
})
saveRDS(Dwide_mat_deseas,
        file = paste0("./dataderived/Dwide_mat_deseas_", year, ".rds"))
# dim(Dwide_mat_deseas)
# plot.ts(Dwide_mat_deseas[, 3300])
rm(D, Dwide, Dwide_mat, EXO, AtlanticCells, isoutlier)

# Pairwise testing ---
set.seed(123456789)

# Create cluster
cl <- parallel::makeCluster(parallel::detectCores())
print(cl)
parallel::clusterExport(cl,
                        varlist = c("Dwide_mat_deseas", "lagmax", "Cells"),
                        envir = environment())

# Run over the upper triangle # i = 1; j = 2
M3 <- parallel::parLapply(cl, 1:(ncol(Dwide_mat_deseas) - 1), fun = function(i) {
    library(vars)
    Js <- (i + 1):ncol(Dwide_mat_deseas)
    Alag_ij <- Alag_ji <- Acoef_ij <- Acoef_ji <- Gp_ij <- Gp_ji <- numeric(length(Js))
    for (j in Js) {
        ## Model, remove intercept estimate
        mod <- BigVAR::BigVAR.fit(Dwide_mat_deseas[,c(i, j)],
                                  p = lagmax,
                                  struct = "Basic",
                                  lambda = 8,
                                  intercept = TRUE)[,-1,1]
        ## Lagged coeffs of the i-th variable affecting the j-th variable
        ind <- 1:lagmax
        x <- abs(mod[2, ind])
        ## Check if not zero, then update
        if (any(x != 0)) {
            phat <- max(which(x != 0))
            # Alag[i, j] <- phat
            # For parallel version
            Alag_ij[j - i] <- phat
            # Remove last 0s
            x <- x[1:phat]
            # Weighted average, inverse proportional to the lag
            # Acoef[i, j] <- sum(x / 1:phat) / sum(1 / 1:phat)
            # For parallel version
            Acoef_ij[j - i] <- sum(x / 1:phat) / sum(1 / 1:phat)
        }
        ## Lagged coeffs of the j-th variable affecting the i-th variable
        x <- abs(mod[1, ind + lagmax])
        # Check if not zero, then update
        if (any(x != 0)) {
            phat <- max(which(x != 0))
            # Alag[j, i] <- phat
            # For parallel version
            Alag_ji[j - i] <- phat
            # Remove last 0s
            x <- x[1:phat]
            # Weighted average, inverse proportional to the lag
            # Acoef[j, i] <- sum(x / 1:phat) / sum(1 / 1:phat)
            # For parallel version
            Acoef_ji[j - i] <- sum(x / 1:phat) / sum(1 / 1:phat)
        }
        mod_var <- VAR(Dwide_mat_deseas[, c(i, j)],
                       p = max(c(1, Alag_ij[j - i], Alag_ji[j - i])),
                       type = "const")
        Gp_ij[j - i] <- causality(mod_var, cause = Cells[i])$Granger$p.value
        Gp_ji[j - i] <- causality(mod_var, cause = Cells[j])$Granger$p.value
    }
    list(i = i, Js = Js,
         Alag_ij = Alag_ij, Acoef_ij = Acoef_ij,
         Gp_ij = Gp_ij,
         Alag_ji = Alag_ji, Acoef_ji = Acoef_ji,
         Gp_ji = Gp_ji)

})
saveRDS(M3,
        file = paste0("./dataderived/CAUSpairs_", year, ".rds"))


## CBL cluster ----
# Use the code below to prepare cluster and run jobs

# In Putty, run
R --vanilla
install.packages("data.table")
install.packages("dplyr")
install.packages("BigVAR")
install.packages("zoo")
install.packages("vars")
# BigVAR deps
install.packages("abind")
# devtools::install_github("vlyubchich/funtimes")
q()

cd /local/users/lyubchich/hynet
sbatch --nodes=1 --mem=0 --time=1-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1pair_1986.R
# To force node-0
sbatch -w statcluster-n0 --mem=0 --time=1-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1pair_2015.R
squeue -u lyubchich
# scancel -u lyubchich


# 2. Visuals ----
# See "hynet_1_vis.qmd" with visualizations for 2002 and 2011.


# 3. Embedding ----

## Adjacency / Laplacian Spectral Embedding (ASE/LSE) ----

rm(list = ls())
# library(data.table)
library(dplyr)
library(igraph)
source("code/fun_create_adj.R")

rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1) %>%
    mutate(CellID = paste0("x", CellID)) %>%
    mutate(Atlantic = (LAT < 37.22 & LON > -75.96) |
               (LAT < 37.5 & LON > -75.8) |
               (LAT < 38.0 & LON > -75.6) |
               (LAT < 36.96 & LON > -75.99)
    ) %>%
    filter(!Atlantic)

YEARS = 1986L:2015L
alpha = 0.05
ndims = 5

edims <- numeric()
E_ASE <- E_LSE <- numeric()
# Full loop on weighted net is about 5 minutes.
for (year in YEARS) { # year = 2002

    # Load data obtained on cluster
    Mcoef <- readRDS(paste0("dataderived/CAUSpairs_", year, ".rds"))
    proc <- create_adj_pairs(Mcoef)

    # Thresholding / adjusting
    # proc$Acoef[proc$Acoef < 0.01] <- 0
    proc$Gp <- apply(proc$Gp, 2, p.adjust, method = "BY")
    diag(proc$Gp) <- 1

    # Unweighted directed adjacency matrix
    AM <- proc$Gp < alpha

    # Weighted directed adjacency matrix
    AMw <- proc$Acoef
    AMw[AM == 0] <- 0

    # Create graph, add attributes
    graphGranger <- igraph::graph_from_adjacency_matrix(AMw, weighted = TRUE
                                                        #AM, weighted = NULL
                                                        ,mode = "directed"
                                                        ,diag = FALSE)
    # igraph::is_weighted(graphGranger)

    # Embed
    # https://bdpedigo.github.io/networks-course/embedding.html#spectral-methods
    E1 <- igraph::embed_adjacency_matrix(graphGranger, no = ndims)
    # edims <- c(edims, igraph::dim_select(E1$D))
    # print(c(year, igraph::dim_select(E1$D)))
    if (FALSE) {
        plot.ts(E1$D)
        plot(E1$X[, 1], E1$Y[, 1])
        plot(E1$X[, 2], E1$Y[, 2])

        plot(E1$X[, 1], E1$X[, 2])
        plot(E1$Y[, 1], E1$Y[, 2])

        plot(E1$X[, 3], E1$Y[, 3])
    }
    E_ASE <- rbind(E_ASE, cbind(rca_cells$CellID, year, E1$X, E1$Y))
    E2 <- igraph::embed_laplacian_matrix(graphGranger, no = ndims, type = "OAP")
    E_LSE <- rbind(E_LSE, cbind(rca_cells$CellID, year, E2$X, E2$Y))
}
colnames(E_ASE) <- c("CellID", "Year",
                     paste0("ASE_X", 1:ndims),
                     paste0("ASE_Y", 1:ndims)
)
colnames(E_LSE) <- c("CellID", "Year",
                     paste0("LSE_X", 1:ndims),
                     paste0("LSE_Y", 1:ndims)
)
write.csv(E_ASE,
          file = "dataderived/embedding_ASE.csv",
          row.names = FALSE)
write.csv(E_LSE,
          file = "dataderived/embedding_LSE.csv",
          row.names = FALSE)

# After running
# E1 <- igraph::embed_adjacency_matrix(graphGranger, no = 50)
# i.e., up to 50 dimensions, see the summary
summary(edims)
# For weighted graphs (AMw, graph_from_adjacency_matrix):
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.0     4.0     4.0     4.2     5.0     7.0
# For unweighted graphs (AM, graph_from_adjacency_matrix):
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.000   1.250   2.000   2.067   2.750   4.000
# In most cases, 5 dimensions is enough, so change the above code
# to 5 dims and concatenate results.


## Manual ----
# Get (weighted) averages of lagged O2 from cells that influence
# the cells with benthic data, over 1-ndims neighborhoods.

rm(list = ls())
# library(data.table)
library(dplyr)
library(igraph)
source("code/fun_create_adj.R")

rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1) %>%
    mutate(CellID = paste0("x", CellID)) %>%
    mutate(Atlantic = (LAT < 37.22 & LON > -75.96) |
               (LAT < 37.5 & LON > -75.8) |
               (LAT < 38.0 & LON > -75.6) |
               (LAT < 36.96 & LON > -75.99)
    ) %>%
    filter(!Atlantic)

# Match with benthic biomass data
sta_cells <- readr::read_csv("data_benthos/stations_cells.csv") %>%
    mutate(CellID = paste0("x", CellID))
BB <- readr::read_csv("data_benthos/benthos_biomass.csv") %>%
    filter(SITE_TYPE == "RANDOM") %>%
    rename(LAT = LATITUDE,
           LON = LONGITUDE) %>%
    mutate(Atlantic = (LAT < 37.22 & LON > -75.96) |
               (LAT < 37.5 & LON > -75.8) |
               (LAT < 38.0 & LON > -75.6) |
               (LAT < 36.96 & LON > -75.99)
    ) %>%
    filter(!Atlantic) %>%
    select(STATION, SAMPLE_DATE, Year)
BB <- BB %>%
    left_join(sta_cells, by = c("STATION", "SAMPLE_DATE", "Year")) %>%
    filter(!is.na(CellID)) %>%
    group_by(CellID, SAMPLE_DATE) %>%
    summarise(Year = mean(Year)) %>%
    ungroup()
# plot(x = BB$LONGITUDE, y = BB$LATITUDE)

YEARS = 1986L:2015L
YEARS = YEARS[YEARS >= min(BB$Year) & YEARS <= max(BB$Year)]
alpha = 0.05
ndims = 2

edims <- numeric()
E_MAN <- numeric()
# Full loop on weighted net is about 5 minutes.
for (year in YEARS) { # year = 1995

    # Cells numbers 1-nrow(rca_cells) for which to get embeddings,
    # not CellID but matching the order in the adjacency matrix.
    bb <- BB %>%
        filter(Year == year) %>%
        mutate(Cell = sapply(CellID, function(i) which(i == rca_cells$CellID)))

    # Load data obtained on cluster
    O2 <- readRDS(paste0("dataderived/Dwide_mat_deseas_", year, ".rds"))
    Mcoef <- readRDS(paste0("dataderived/CAUSpairs_", year, ".rds"))
    proc <- create_adj_pairs(Mcoef)

    # Thresholding / adjusting
    # proc$Acoef[proc$Acoef < 0.01] <- 0
    proc$Gp <- apply(proc$Gp, 2, p.adjust, method = "BY")
    diag(proc$Gp) <- 1

    # Unweighted directed adjacency matrix
    AM <- proc$Gp < alpha

    # Weighted directed adjacency matrix
    AMw <- proc$Acoef
    AMw[AM == 0] <- 0

    # Lags of dependence
    Alag <- proc$Alag
    Alag[AM == 0] <- 0

    # In-degrees
    INDegrees <- apply(Alag, 2, function(x) sum(x > 0))

    # Get neighborhood information
    sapply(1:nrow(bb), function(x) { # x = w = 1
        seed <- bb$Cell[x]
        seed_date <- bb$SAMPLE_DATE[x]
        seed_day <- as.numeric(format(seed_date, "%j"))
        # Select in-waves
        waves <- as.list(rep(NA, ndims + 1))
        waves[[1]] <- data.frame(ids = seed, lags = NA, degrees = INDegrees[x])
        for (w in 1:ndims) {
            # Select influencers of waves[[w]]
            ids <- unlist(lapply(waves[[w]]$ids, function(nd) which(Alag[, nd] > 0)))
            lags <- unlist(lapply(waves[[w]]$ids, function(nd) Alag[which(Alag[, nd] > 0), nd]))
            D <- data.frame(ids = ids,
                            lags = lags,
                            degees = INDegrees[ids])
            # Get averages
            D$o2 <- sapply(1:nrow(D), function(rd) {
                mean(O2[, D$ids[rd]])
            })


            waves[[w + 1]] <- D
        }

    })


    E_LSE <- rbind(E_LSE, cbind(rca_cells$CellID, year, E2$X, E2$Y))
}








colnames(E_ASE) <- c("CellID", "Year",
                     paste0("ASE_X", 1:ndims),
                     paste0("ASE_Y", 1:ndims)
)
colnames(E_LSE) <- c("CellID", "Year",
                     paste0("LSE_X", 1:ndims),
                     paste0("LSE_Y", 1:ndims)
)
write.csv(E_ASE,
          file = "dataderived/embedding_ASE.csv",
          row.names = FALSE)
write.csv(E_LSE,
          file = "dataderived/embedding_LSE.csv",
          row.names = FALSE)

# After running
# E1 <- igraph::embed_adjacency_matrix(graphGranger, no = 50)
# i.e., up to 50 dimensions, see the summary
summary(edims)



