# Test different options for creating the hypoxia network

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

## Plot some time series ----
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

## Remove 1 seasonality ----
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
# 2do+ Create a matrix with selected AR(p) for p.
# 2do+ Try larger lags to justify lower limit.


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

# without lag selection
system.time(
    for (v in 1:1000) { # v = 5
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        # p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, type = "const", exogen = exo)
        # p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = 10, type = "const", exogen = exo)
        M[i, j] <- causality(mod_var, cause = Cells[i])$Granger$p.value
        M[j, i] <- causality(mod_var, cause = Cells[j])$Granger$p.value
    }
)
19.11 * nrow(M_low) / 1000 / 3600 # 30 hours for p = 10
36.65 * nrow(M_low) / 1000 / 3600 # 58 hours for p = 30


## foreach ----

library(foreach)
library(doParallel)
# Create cluster
cl <- parallel::makeCluster(parallel::detectCores())
# Register it for the foreach loop
doParallel::registerDoParallel(cl)
# Export the dataset (could be done directly in the foreach, but this is more explicit)
# parallel::clusterExport(cl,
#                         varlist = c("Dwide", "M_low", "exo", "Cells"),
#                         envir = environment())
# parallel::clusterExport(cl, "VARselect")

system.time(
    M2 <- foreach(v = 1:1000, .combine = rbind, .inorder = FALSE,
                  .packages = c("data.table", "vars")) %dopar% {
        # library(data.table)
        # library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 30, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = p_opt, type = "const", exogen = exo)
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    }
)
9.94 * nrow(M_low) / 1000 / 3600 # 15.82 hours
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
        p_opt <- VARselect(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 30, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide[, .SD, .SDcols = c(i, j)], p = p_opt, type = "const", exogen = exo)
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    })
)
7.07 * nrow(M_low) / 1000 / 3600 # 11.25 hours
stopCluster(cl)
M3 <- matrix(runif(2*nrow(M_low)), nrow = 2)

# without data.table
Dwide_mat <- as.matrix(Dwide)
parallel::clusterExport(cl,
                        varlist = c("Dwide_mat", "M_low", "exo", "Cells"),
                        envir = environment())
system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        # library(data.table)
        library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide_mat[, c(i, j)], lag.max = 30, type = "const", exogen = exo)
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide_mat[, c(i, j)], p = p_opt, type = "const", exogen = exo)
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    })
)
6.56 * nrow(M_low) / 1000 / 3600 # 10.44 hours
stopCluster(cl)

# without data.table and exo
cl <- parallel::makeCluster(parallel::detectCores())
Dwide_mat <- as.matrix(Dwide)
parallel::clusterExport(cl,
                        varlist = c("Dwide_mat", "M_low", "exo", "Cells"),
                        envir = environment())
system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        # library(data.table)
        library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        p_opt <- VARselect(Dwide_mat[, c(i, j)], lag.max = 30, type = "const")
        p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide_mat[, c(i, j)], p = p_opt, type = "const")
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    })
)
6.31 * nrow(M_low) / 1000 / 3600 # 10.04 hours
stopCluster(cl)

# without data.table and exo and lag selection
cl <- parallel::makeCluster(parallel::detectCores())
Dwide_mat <- as.matrix(Dwide)
parallel::clusterExport(cl,
                        varlist = c("Dwide_mat", "M_low", "Cells"),
                        envir = environment())
system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        # library(data.table)
        library(vars)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        # p_opt <- VARselect(Dwide_mat[, c(i, j)], lag.max = 30, type = "const")
        # p_opt <- p_opt$selection["AIC(n)"]
        mod_var <- VAR(Dwide_mat[, c(i, j)], p = 30, type = "const")
        c(causality(mod_var, cause = Cells[i])$Granger$p.value,
          causality(mod_var, cause = Cells[j])$Granger$p.value)
    })
)
8.89 * nrow(M_low) / 1000 / 3600 # 14.15 hours
stopCluster(cl)


# lmtest
system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        library(data.table)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        test1 <- lmtest::grangertest(Dwide[, .SD, .SDcols = c(i)], Dwide[, .SD, .SDcols = c(j)], order = 30)$`Pr(>F)`[2]
        test2 <- lmtest::grangertest(Dwide[, .SD, .SDcols = c(j)], Dwide[, .SD, .SDcols = c(i)], order = 30)$`Pr(>F)`[2]
        c(test1,
          test2)
    })
)
80.09 * nrow(M_low) / 1000 / 3600 # 127.4945 hours
stopCluster(cl)
M3 <- matrix(runif(2*nrow(M_low)), nrow = 2)


system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        library(data.table)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        res1 <- funtimes::causality_pred(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, p.free = TRUE, B = 0, test = 1)
        res2 <- funtimes::causality_pred(Dwide[, .SD, .SDcols = c(j, i)], lag.max = 10, p.free = TRUE, B = 0, test = 1)
        # Selected lags of the cause variables
        c(res1$FullH0$p[2],
          res2$FullH0$p[2])
    })
)
422.17 * nrow(M_low) / 1000 / 3600 # 672.0483 hours with lag.max = 30
58.63 * nrow(M_low) / 1000 / 3600 # 93.33253 hours with lag.max = 10
stopCluster(cl)

summary(as.vector(M3))
mean(M3 == 10)*100
# 2do: run with smaller lag.max or code if statements to run with bigger lag
# Use p.free, strip down the function of causality_pred.

## Put back into adjacency matrix ----
for (v in 1:nrow(M_low)) {
    i <- M_low[v, 1]
    j <- M_low[v, 2]
    M[i, j] <- M3[1, v]
    M[j, i] <- M3[2, v]
}

source("./code/lagcause.R")
parallel::clusterExport(cl,
                        varlist = c("Dwide", "M_low", "lagcause"),
                        envir = environment())
system.time(
    M3 <- parallel::parSapply(cl, 1:1000, FUN = function(v) {
        library(data.table)
        i <- M_low[v, 1]
        j <- M_low[v, 2]
        res1 <- lagcause(Dwide[, .SD, .SDcols = c(i, j)], lag.max = 10, p.free = TRUE)
        res2 <- lagcause(Dwide[, .SD, .SDcols = c(j, i)], lag.max = 10, p.free = TRUE)
        # Selected lags of the cause variables
        c(res1$FullH0$p[2],
          res2$FullH0$p[2])
    })
)

stopCluster(cl)

# 2do correct for multiple test, by row. (parApply?)
# p.adjust(, method = "BH") # or Holm

# Plot ----
library(igraph)
library(network)
# library(randomcoloR)

# Significance level alpha
level_sig = 0.05

# Random subset of nodes because the plot is too slow
ss <- sort(sample(1:nrow(M), 100))

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

#https://stackoverflow.com/questions/22453273/how-to-visualize-a-large-network-in-r
g <- graph_from_adjacency_matrix(M[ss, ss] < level_sig,
                                 mode = "directed",
                                 diag = FALSE)
dev.off()
# ?plot.igraph
coords <- layout_(g, as_star())
summary(coords)
# from -1 to 1
# library(scales)
coords <- rca_cells[ss, ] %>%
    mutate(lon = scales::rescale(LON, to = c(-1, 1)),
           lat = scales::rescale(LAT, to = c(-1, 1)),
           .keep = "none") %>%
    as.matrix()

# 2do squeeze to different scale for lon, eg -0.3 to 0.3

plot(simplify(g),
     vertex.size = 0.01 ,
     vertex.label.cex = 0.75,
     vertex.label = NA,
     vertex.label.color = "black",
     # vertex.frame.color = adjustcolor("white", alpha.f = 0),
     vertex.color = adjustcolor("white", alpha.f = 0),
     edge.arrow.size=0.001,
     edge.color = rgb(0, 58, 76, alpha = 30, maxColorValue = 100),
     # display.isolates=TRUE
     # ,
     # vertex.label=ifelse(page_rank(g)$vector > 0.1 ,
     #                     "important nodes", NA)
     layout = coords
)

# 2do try BioFabric if it is useful. Clusters?
# 2do M as clusters?

# Embedding ----
library(igraph)
G <- graph_from_adjacency_matrix(M < level_sig,
                                 mode = "directed",
                                 diag = FALSE)
summary(G)

## v1 ----
E1 <- igraph::embed_adjacency_matrix(G, no = 9)
plot(E1$X[, 1], E1$Y[, 1])
plot(E1$X[, 2], E1$Y[, 2])

## v2 ----
# very slow
library(node2vec)
E2 <- node2vecR(as_edgelist(G), dim = 9,
                walk_length = 3,
                directed = TRUE)

# 2do: differences between embeddings?
# 2do: Embedding with node information, so the DO info is included.
# (Create separate embeddings for snapshots when need to predict?)

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


# library(BigVAR)


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

# Without running cross-validation, with fixed lambda
mod3 <- BigVAR.fit(Dwide_mat[,c(1, 100, 200)*10],
                   p = 30,
                   struct = "Basic",
                   lambda = 2,
                   intercept = TRUE)


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


