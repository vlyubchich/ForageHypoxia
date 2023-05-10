# Hypoxia network for 1 year

# Packages ----
rm(list = ls())
library(data.table)
library(dplyr)
library(BigVAR)
library(parallel)

year = 1990
lagmax = 30L


# Data ----
D <- data.table::fread(paste0("./data_rca/rca_ts_", year, "_2.csv"),
                       select = c("DOAVEG_avg", "CellID", "Date")) %>%
    mutate(CellID = paste0("x", CellID))
# head(D)
Dwide <- data.table::dcast(D, Date ~ CellID, value.var = "DOAVEG_avg")
Dwide_mat <- as.matrix(Dwide[, -1])
Cells <- colnames(Dwide_mat)
# dim(Dwide_mat)


# Remove seasonality ----
EXO <- Dwide[, 1] %>%
    mutate(day_year = as.numeric(format(Date, "%j"))) %>%
    mutate(sin2 = sin(2 * pi * day_year/365.25),
           cos2 = cos(2 * pi * day_year/365.25)) %>%
    dplyr::select(sin2, cos2) %>%
    as.matrix()
Dwide_mat_deseas <- apply(Dwide_mat, 2, function(x) {
    m0 = lm(x ~ EXO)
    residuals(m0)
})
saveRDS(Dwide_mat_deseas,
        file = paste0("./dataderived/Dwide_mat_deseas_", year, ".rds"))
# dim(Dwide_mat_deseas)
# plot.ts(Dwide_mat_deseas[, 500])
rm(D, Dwide, Dwide_mat, EXO)

# Pairwise testing ----
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
