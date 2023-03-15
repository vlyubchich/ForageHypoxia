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
mod1 <- constructModel(Dwide_mat[,1:3],
                       p = 90,
                       struct = "Basic",
                       gran = c(50, 10),
                       h = 1,
                       cv = "Rolling",
                       verbose = FALSE,
                       IC = TRUE,
                       model.controls = list(intercept = TRUE))
# mod1
results <- cv.BigVAR(mod1)
# results
# plot(results)
# SparsityPlot.BigVAR.results(results)
# coef(results)
saveRDS(results, file = "./dataderived/BigVAR_2011.rds")
