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
dim(Dwide_mat)

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
# Alag <- Acoef <- Gp <- matrix(0, ncol(Dwide_mat), ncol(Dwide_mat))


library(parallel)
# Create cluster
cl <- parallel::makeCluster(parallel::detectCores())
parallel::clusterExport(cl,
                        varlist = c("Dwide_mat", "lagmax"),
                        envir = environment())

# Run over the upper triangle # i = 1; j = 2
M3 <- parallel::parLapply(cl, 1:(ncol(Dwide_mat) - 1), fun = function(i) {
    # for (i in 1:(ncol(Dwide_mat) - 1)) {
    # if (i %% 10 == 0) {
    #     print(paste(i, Sys.time()))
    # }
    # For parallel version
    Js <- (i + 1):ncol(Dwide_mat)
    Alag_ij <- Acoef_ij <- Gp_ij <- Alag_ji <- Acoef_ji <- Gp_ji <- numeric(length(Js))
    for (j in Js) {
        # Model, remove intercept estimate
        mod <- BigVAR::BigVAR.fit(Dwide_mat[,c(i, j)],
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
        Gt <- lmtest::grangertest(Dwide_mat[,j] ~ Dwide_mat[,i]
                                  ,order = lagmax)
        # Gp[i, j] <- Gt$`Pr(>F)`[2]
        # For parallel version
        Gp_ij[j - i] <- Gt$`Pr(>F)`[2]
        Gt <- lmtest::grangertest(Dwide_mat[,i] ~ Dwide_mat[,j]
                                  ,order = lagmax)
        # Gp[j, i] <- Gt$`Pr(>F)`[2]
        # For parallel version
        Gp_ji[j - i] <- Gt$`Pr(>F)`[2]
    }
    list(i = i, Js = Js,
         Alag_ij = Alag_ij, Acoef_ij = Acoef_ij, Gp_ij = Gp_ij,
         Alag_ji = Alag_ji, Acoef_ji = Acoef_ji, Gp_ji = Gp_ji)

})
# saveRDS(Gp,
#         file = paste0("./dataderived/BigVAR_", year, "_Gp.rds"))
# saveRDS(Alag,
#         file = paste0("./dataderived/BigVAR_", year, "_Alag.rds"))
# saveRDS(Acoef,
#         file = paste0("./dataderived/BigVAR_", year, "_Acoef.rds"))
saveRDS(M3,
        file = paste0("./dataderived/VARpairs_", year, ".rds"))
