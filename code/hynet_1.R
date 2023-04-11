# Hypoxia network for 1 year

# CBL cluster ----
# In Putty, run
R --vanilla
install.packages("data.table")
install.packages("dplyr")
install.packages("BigVAR")
# BigVAR deps
install.packages("abind")
# devtools::install_github("vlyubchich/funtimes")
q()

cd /local/users/lyubchich/hynet
# srun -l R CMD BATCH "--vanilla" test_parallel.R # prints 32 cores
sbatch --nodes=2 --mem=0 --time=5-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1.R
sbatch --nodes=2 --mem=0 --time=5-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1_2002.R
squeue -u lyubchich
# scancel -u lyubchich


# Packages ----
rm(list = ls())
library(data.table)
library(dplyr)
library(BigVAR)


# Data ----

# rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
#     filter(FSM == 1) %>%
#     mutate(CellID = paste0("x", CellID))
D <- data.table::fread("./data_rca/rca_ts_2011_2.csv",
                       select = c("DOAVEG_avg", "CellID", "Date")) %>%
    mutate(CellID = paste0("x", CellID))
# head(D)
Dwide <- data.table::dcast(D, Date ~ CellID, value.var = "DOAVEG_avg")
Dwide_mat <- as.matrix(Dwide[, -1])


# VAR ----
set.seed(123456789)
mod1 <- constructModel(Dwide_mat[,c(1, 100, 200)*10],
                       p = 30,
                       struct = "Basic",
                       gran = c(5000, 50),
                       h = 1,
                       cv = "Rolling",
                       verbose = FALSE,
                       IC = TRUE,
                       model.controls = list(intercept = TRUE))
# mod1
results <- cv.BigVAR(mod1)
# results
plot(results)
# extract optimal lambda by running the above code with *1, *2,..., *5, *10
results@OptimalLambda
# [1] 1.583282 1.896021 1.935361 2.254627 6.748183 10.38759

# Use fixed lambda to save time
mod2 <- constructModel(Dwide_mat[,c(1, 100, 200)*10],
                       p = 30,
                       struct = "Basic",
                       ownlambdas = TRUE,
                       gran = c(1, 2, 4, 8, 16),
                       h = 1,
                       cv = "Rolling",
                       verbose = FALSE,
                       IC = TRUE,
                       model.controls = list(intercept = TRUE))
results <- cv.BigVAR(mod2)
plot(results)
SparsityPlot.BigVAR.results(results)
M <- as.matrix(coef(results)[, -1])
print(results@OptimalLambda)
saveRDS(M, file = "./dataderived/BigVAR_2011_coef_test.rds")


# Load VAR ----

## Load and extract coefficients ----
if (FALSE) {
    gc()
    library(BigVAR)
    library(lattice)
    M <- readRDS("dataderived/BigVAR_2002.rds")
    # plot(M)
    lagmax <- M@lagmax
    # VAR coefficients, without intercepts
    Mcoef <- as.matrix(coef(M)[, -1])
    # Mcoef[1:9, 1:7]
    saveRDS(Mcoef, file = paste0("./dataderived/BigVAR_2002_coef_lagmax",
                                 lagmax, ".rds"))
}

## Process coefficients ----
lagmax = 3
Mcoef <- readRDS(paste0("./dataderived/BigVAR_2011_coef_lagmax",
                        lagmax, ".rds"))
Mcoef <- as.matrix(Mcoef)
# For each spatial cell (row),
Alag <- Acoef <- matrix(0, nrow(Mcoef), nrow(Mcoef))
for (i in 1:nrow(Mcoef)) { # i = j = 1
    for (j in 1:nrow(Mcoef)) {
        # Lagged coeffs of the j's variable
        js <- (j - 1) * lagmax + 1
        js <- js:(js + lagmax - 1)
        x <- abs(Mcoef[i, js])
        # Check if not zero, then update
        if (any(x != 0)) {
            Alag[i, j] <- which.max(x)
            # Weighted average, inverse proportional to the lag
            # but do not divide because of too small values (for numerical stability)
            Acoef[i, j] <- sum(x / 1:lagmax) #/ sum(1 / 1:lagmax)
        }
    }
    if (i %% 100 == 0) {
        print(paste(i, Sys.time()))
    }
}
Alag[1:5, 1:5]
Acoef[1:5, 1:5]

## Network ----
library(dplyr)
library(igraph)
# library(network)
# If subsample
ss <- 1:500

rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1) %>%
    mutate(CellID = paste0("x", CellID))

# NET <- network::network(Acoef,
#                         matrix.type = "adjacency",
#                         directed = TRUE,
#                         loops = FALSE)
g <- igraph::graph_from_adjacency_matrix(Acoef#[ss, ss]
                                         ,weighted = TRUE
                                         ,mode = "directed"
                                         ,diag = FALSE)
coords <- layout_(g, as_star())
summary(coords)
# from -1 to 1
# library(scales)
coords <- rca_cells[ss, ] %>%
    mutate(lon = scales::rescale(LON, to = c(-1, 1)),
           lat = scales::rescale(LAT, to = c(-1, 1)),
           .keep = "none") %>%
    as.matrix()

# devtools::install_github("wjrl/RBioFabric")
# library(RBioFabric)
#
# height <- vcount(g)
# width <- ecount(g)
# aspect <- height / width
# plotWidth <- 10
# plotHeight <- 2 #plotWidth * (aspect * 1.2)
#
# pdf("images/BioFabric_2011.pdf", width = plotWidth, height = plotHeight)
# bioFabric(g, dropNodeLabels = TRUE, dropZoneLabels = TRUE)
# dev.off()

g_dd_in <- degree(g, mode = "in")
hist(g_dd_in)
g_dd_out <- degree(g, mode = "out")
hist(g_dd_out)

# map of degrees, centralities


kc <- igraph::cluster_edge_betweenness(g)
kc
