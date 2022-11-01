
rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)


# Data ----
rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
    filter(FSM == 1) %>%
    mutate(CellID = paste0("x", CellID))
D <- data.table::fread("./data_rca/rca_ts_1986_2.csv",
                       select = c("DOAVEG_avg", "CellID", "Date")) %>%
    mutate(CellID = paste0("x", CellID))
head(D)
nrow(D) / 3386

nodes <- unique(D$CellID)
# Plot some time series
tmp <- c("x332", sample(nodes, 5))
D %>%
    filter(is.element(CellID, tmp)) %>%
    mutate(CellID = as.factor(CellID)) %>%
    ggplot(aes(x = Date, y = DOAVEG_avg, col = CellID)) +
    geom_line() +
    theme_minimal()

d <- D %>%
    filter(is.element(CellID, "x332")) %>%
    mutate(day_year = as.numeric(format(Date, "%j"))) %>%
    mutate(sin2 = sin(2 * pi * day_year/365.25),
           cos2 = cos(2 * pi * day_year/365.25))

# Remove seasonality
mod0 <- lm(DOAVEG_avg ~ sin2 + cos2, data = d)
d %>%
    mutate(mod0_resid = residuals(mod0)) %>%
    ggplot(aes(x = Date, y = mod0_resid, col = CellID)) +
    geom_line() +
    theme_minimal()

# Network ----
library(vars)

Dwide <- data.table::dcast(D, Date ~ CellID, value.var = "DOAVEG_avg")
exo <- Dwide[, 1] %>%
    mutate(day_year = as.numeric(format(Date, "%j"))) %>%
    mutate(sin2 = sin(2 * pi * day_year/365.25),
           cos2 = cos(2 * pi * day_year/365.25)) %>%
    dplyr::select(sin2, cos2) %>%
    as.matrix()
# Remove the column Date so keep only the data
Dwide <- Dwide[, Date := NULL]
Cells <- colnames(Dwide)

M <- matrix(data = NA,
            nrow = length(Cells),
            ncol = length(Cells),
            dimnames = list(Cells, Cells))
M_low <- which(lower.tri(M), arr.ind = TRUE)

## sequential ----
system.time(
    for (v in 1:1000) { # v = 5
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = p_opt, type = "const", exogen = exo)
        M[i, j] <- causality(mod_var, cause = Cells[i])$Granger$p.value
        M[j, i] <- causality(mod_var, cause = Cells[j])$Granger$p.value
    }
)
12.01 * nrow(M_low) / 1000 / 3600 # 19.1186 hours

## foreach ----
library(foreach)
library(doParallel)
# Create cluster
cl <- parallel::makeCluster(parallel::detectCores())
# Register it for the foreach loop
# doParallel::registerDoParallel(cl)
#Export the dataset (could be done directly in the foreach, but this is more explicit)
parallel::clusterExport(cl,
                        varlist = c("Dwide", "M_low", "exo", "Cells"),
                        envir = environment())
# parallel::clusterExport(cl, "VARselect")

system.time(
    M2 <- foreach(v = 1:1000, .combine = rbind, .inorder = FALSE) %dopar% {
        library(data.table)
        library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = p_opt, type = "const", exogen = exo)
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    }
)
5.15 * nrow(M_low) / 1000 / 3600 # 8.198235 hours
stopCluster(cl)


## parSapply ----
library(parallel)
# Create cluster
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl,
                        varlist = c("Dwide", "M_low", "exo", "Cells"),
                        envir = environment())

system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        library(data.table)
        library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = p_opt, type = "const", exogen = exo)
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    })
)
2.61 * nrow(M_low) / 1000 / 3600 # 4.154834 hours
stopCluster(cl)
M3 <- matrix(runif(2*nrow(M_low)), nrow = 2)

## Put back into adjacency matrix ----
for (v in 1:nrow(M_low)) {
    i <- M_low[v, 1]
    j <- M_low[v, 2]
    M[i, j] <- M3[1, v]
    M[j, i] <- M3[2, v]
}

# Plot ----
library(network)
# library(randomcoloR)

# Significance level alpha
level_sig = 0.05

# Random subset of nodes
ss <- sample(1:nrow(M), 50)

N <- network(M[ss, ss] < level_sig,
             matrix.type = "adjacency",
             directed = TRUE, loops = FALSE)
# Add lon and lat information
N %v% 'lon' <- sapply(network.vertex.names(N), function(i) {
    rca_cells$LON[rca_cells$CellID == i]
})
N %v% 'lat' <- sapply(network.vertex.names(N), function(i) {
    rca_cells$LAT[rca_cells$CellID == i]
})

# col = colorRampPalette(c("blue", "orange", "red"))(3)
par(mar = c(0, 0, 0, 0))
plot.network(N, new = TRUE, coord = cbind(N%v%'lon', N%v%'lat'),
             #edge.col='grey',
             vertex.cex = 0.51,
             arrowhead.cex = 0.5,
             edge.col = rgb(0, 58, 76, alpha = 30, maxColorValue = 100), #"steelblue",
             # label = ISO[match(network.vertex.names(trade_Net_map),ii.countries)],
             label.col = "black",
             label.cex = 0.5,
             vertex.col = "grey",
             #edge.lwd=(trade_Net_map%e%'Trdvalue'),
             #arrowhead.cex = 0.5,
             vertex.border = "black", jitter = FALSE,
             # usecurve = 1, edge.curve = 0.2,
             usearrows = TRUE)
# dev.off()

# Embedding ----
library(igraph)
G <- graph_from_adjacency_matrix(M < level_sig,
                                 mode = "directed",
                                 diag = FALSE)
summary(G)

## v1 ----
E1 <- embed_adjacency_matrix(G, no = 9)
plot(E1$X[, 1], E1$Y[, 1])
plot(E1$X[, 2], E1$Y[, 2])

## v2 ----
# very slow
library(node2vec)
E2 <- node2vecR(as_edgelist(G), dim = 9, directed = TRUE)


# unused ----
output <- foreach(i = 1:(nrow(dat)/30), .combine = rbind, .inorder = FALSE) %:%
    foreach(j = 1:(ncol(dat)/30), .combine = rbind, .inorder = FALSE) %dopar%{
        row <- 30 * (i - 1) + 1
        col <- 30 * (j - 1) + 1
        c(x = max(dat[row:(row + 29), col:(col + 29)]), i = i, j = j)
    }

data.table::dcast(Date ~ CellID, value.var = "DOAVEG_avg")
matplot(d, type = "l")

p_opt <- VARselect(Dwide[, 1:2], lag.max = 10, type = "const", exogen = exo)
p_opt <- p_opt$selection["AIC(n)"]
mod_var <- VAR(Dwide[, 1:2], p = p_opt, type = "const", exogen = exo)
causality(mod_var)
