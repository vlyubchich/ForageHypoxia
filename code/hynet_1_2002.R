# Hypoxia network for 1 year

# Packages ----
rm(list = ls())
library(data.table)
library(dplyr)
library(BigVAR)


# Data ----
# rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
#     filter(FSM == 1) %>%
#     mutate(CellID = paste0("x", CellID))
D <- data.table::fread("./data_rca/rca_ts_2002_2.csv",
                       select = c("DOAVEG_avg", "CellID", "Date")) %>%
    mutate(CellID = paste0("x", CellID))
head(D)
Dwide <- data.table::dcast(D, Date ~ CellID, value.var = "DOAVEG_avg")
Dwide_mat <- as.matrix(Dwide[, -1])


# VAR ----
set.seed(123456789)
# mod1 <- constructModel(Dwide_mat,
#                        p = 30,
#                        struct = "Basic",
#                        ownlambdas = TRUE,
#                        gran = c(1, 2, 4, 8, 16),
#                        h = 1,
#                        cv = "Rolling",
#                        verbose = FALSE,
#                        IC = TRUE,
#                        model.controls = list(intercept = TRUE))
# # mod1
# results <- cv.BigVAR(mod1)
# # plot(results)
# # SparsityPlot.BigVAR.results(results)
# M <- as.matrix(coef(results)[, -1])
# print(results@OptimalLambda)
# saveRDS(results@OptimalLambda,
#         file = "./dataderived/BigVAR_2002_lambda.rds")

M <- BigVAR.fit(Dwide_mat,
                p = 30,
                struct = "Basic",
                lambda = 2,
                intercept = TRUE)
saveRDS(M,
        file = "./dataderived/BigVAR_2002_coef.rds")
