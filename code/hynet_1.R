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
# sbatch --nodes=2 --mem=0 --time=5-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1.R
sbatch --nodes=2 --mem=0 --time=7-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1_2002.R
sbatch --nodes=1 --mem=0 --time=7-23 R CMD BATCH "--vanilla --no-save --no-restore" ./hynet_1pair_2002.R
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

