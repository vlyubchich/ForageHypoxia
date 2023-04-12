# Create adjacency matrices with lags and averaged coefficients
# from rows to columns
create_adj <- function(X,
                       lagmax = NULL,
                       verbose = FALSE) {
    # X is a matrix of BigVAR coefficients, n * (n*lagmax)
    # lagmax is the maximal order p in VAR(p)
    # verbose prints out the progress after each 100 rows processed
    # Guess lagmax if it is not given
    if (is.null(lagmax)) {
        lagmax_tmp <- ncol(X) / nrow(X)
        if (lagmax_tmp %% 1 == 0) {
            lagmax <- lagmax_tmp
        }
    }
    Alag <- Acoef <- matrix(0, nrow(X), nrow(X))
    weigts_sum <- sum(1 / 1:lagmax)
    for (i in 1:nrow(X)) { # i = j = 1
        js <- 1:lagmax
        for (j in 1:nrow(X)) {
            # Lagged coeffs of the j's variable
            x <- abs(X[i, js])
            # Check if not zero, then update
            if (any(x != 0)) {
                Alag[j, i] <- which.max(x)
                # Weighted average, inverse proportional to the lag
                Acoef[j, i] <- sum(x / 1:lagmax) / weigts_sum
            }
            js <- js + lagmax
        }
        if (verbose) {
            if (i %% 100 == 0) {
                print(paste(i, Sys.time()))
            }
        }
    }
    return(list(Alag = Alag,
                Acoef = Acoef))
}
