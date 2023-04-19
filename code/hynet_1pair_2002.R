# Hypoxia network for 1 year

# Packages ----
rm(list = ls())
library(data.table)
library(dplyr)
library(BigVAR)

year = 2002


# Data ----
# rca_cells <- data.table::fread("./data_rca/rca_cells_2.csv") %>%
#     filter(FSM == 1) %>%
#     mutate(CellID = paste0("x", CellID))
D <- data.table::fread(paste0("./data_rca/rca_ts_", year, "_2.csv"),
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
lagmax = 30L
if (FALSE) {
    lmtest::grangertest(x1000 ~ x1315
                        ,order = lagmax
                        ,data = Dwide)
    tmp <- sample(3000, 10)
    Dwide_mat <- Dwide_mat[,tmp]
    summary(as.vector(Alag[lower.tri(Alag)]))
    summary(as.vector(Alag[upper.tri(Alag)]))
    Dwide_mat <- as.matrix(Dwide[, -1])
    lmtest::grangertest(Dwide_mat[,2325] ~ Dwide_mat[,151]
                        ,order = lagmax)
}
Alag <- Acoef <- Gp <- matrix(0, ncol(Dwide_mat), ncol(Dwide_mat))
# Run over the upper triangle # i = 1; j = 2
for (i in 1:(ncol(Dwide_mat) - 1)) {
    if (i %% 10 == 0) {
        print(paste(i, Sys.time()))
    }
    for (j in (i + 1):ncol(Dwide_mat)) {
        # Model, remove intercept estimate
        mod <- BigVAR.fit(Dwide_mat[,c(i, j)],
                          p = lagmax,
                          struct = "Basic",
                          lambda = 10,
                          intercept = TRUE)[,-1,1]
        ## Lagged coeffs of the i-th variable affecting the j-th variable
        ind <- 1:lagmax
        x <- abs(mod[2, ind])
        # Check if not zero, then update
        if (any(x != 0)) {
            phat <- max(which(x != 0))
            Alag[i, j] <- phat
            # Remove last 0s
            x <- x[1:phat]
            # Weighted average, inverse proportional to the lag
            Acoef[i, j] <- sum(x / 1:phat) / sum(1 / 1:phat)
        }
        ## Lagged coeffs of the j-th variable affecting the i-th variable
        x <- abs(mod[1, ind + lagmax])
        # Check if not zero, then update
        if (any(x != 0)) {
            phat <- max(which(x != 0))
            Alag[j, i] <- phat
            # Remove last 0s
            x <- x[1:phat]
            # Weighted average, inverse proportional to the lag
            Acoef[j, i] <- sum(x / 1:phat) / sum(1 / 1:phat)
        }
        Gt <- lmtest::grangertest(Dwide_mat[,j] ~ Dwide_mat[,i]
                                  ,order = lagmax)
        Gp[i, j] <- Gt$`Pr(>F)`[2]
        Gt <- lmtest::grangertest(Dwide_mat[,i] ~ Dwide_mat[,j]
                                  ,order = lagmax)
        Gp[j, i] <- Gt$`Pr(>F)`[2]
    }
}
saveRDS(Gp,
        file = paste0("./dataderived/BigVAR_", year, "_Gp.rds"))
saveRDS(Alag,
        file = paste0("./dataderived/BigVAR_", year, "_Alag.rds"))
saveRDS(Acoef,
        file = paste0("./dataderived/BigVAR_", year, "_Acoef.rds"))
