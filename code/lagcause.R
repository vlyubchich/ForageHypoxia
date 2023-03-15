if (FALSE) {
    Canada <- vars::Canada
    y = Canada[,1:2]; cause = "e"; lag.max = 5; p.free = TRUE
    p = NULL; lag.restrict = 0L; k = 2
    # lagcause(Canada[,1:2], cause = "e", lag.max = 5, p.free = TRUE)
}
lagcause <- function(y,
                     cause = NULL,
                     p = NULL,
                     p.free = FALSE,
                     lag.restrict = 0L,
                     lag.max = NULL,
                     k = 2)
{
    if (!is.null(lag.max) && any(lag.max < 1L)) {
        stop("lag.max must be positive integers.")
    }
    if (!is.null(p) && any(p < 1L)) {
        stop("p must be positive integers.")
    }
    varnames <- colnames(y)
    if (length(varnames) != 2)
        stop("y must have 2 columns")
    if (is.null(cause)) {
        cause <- varnames[1]
    }
    dep <- setdiff(varnames, cause)
    n <- nrow(y)
    if (is.null(p) && is.null(lag.max)) {
        stop("Please specify p or lag.max.")
    }
    if (!is.null(lag.max)) {
        if (length(lag.max) == 2) {
            p.free = TRUE
        } else {
            lag.max <- c(lag.max, lag.max)
        }
        maxl <- max(lag.max)
        if (lag.restrict >= lag.max[2]) {
            warning("lag.restrict >= lag.max. Using lag.restrict = 0 instead.")
            lag.restrict <- 0
        }
        if (lag.restrict > 0) {
            lagX <- embed(y[, cause], maxl + 1)[, -c(1:(lag.restrict + 1)), drop = FALSE]
        } else {
            lagX <- embed(y[, cause], maxl + 1)[, -1, drop = FALSE]
        }
        lagY <- embed(y[, dep], maxl + 1)
        if (p.free) {
            best.ic <- Inf
            for (p1 in 1:lag.max[1]) {
                for (p2 in (lag.restrict + 1):lag.max[2]) {
                    fit <- stats::lm.fit(x = cbind(1,
                                                   lagY[, 2:(p1 + 1)],
                                                   lagX[, 1:(p2 - lag.restrict), drop = FALSE]),
                                         y = lagY[, 1])
                    nfit <- length(fit$residuals)
                    edf <- nfit - fit$df.residual
                    RSS <- sum(fit$residuals^2, na.rm = TRUE)
                    dev <- nfit * log(RSS/nfit)
                    fit.ic <- dev + k * edf
                    if (fit.ic < best.ic) {
                        best.ic <- fit.ic
                        p <- c(p1, p2)
                    }
                }
            }
        } else {
            IC <- sapply((lag.restrict + 1):lag.max[2], function(s) {
                fit <- stats::lm.fit(x = cbind(1,
                                               lagY[, 2:(s + 1)],
                                               lagX[, 1:(s - lag.restrict), drop = FALSE]),
                                     y = lagY[, 1])
                nfit <- length(fit$residuals)
                edf <- nfit - fit$df.residual
                RSS <- sum(fit$residuals^2, na.rm = TRUE)
                dev <- nfit * log(RSS/nfit)
                dev + k * edf
            })
            p <- which.min(IC) + lag.restrict
            p <- c(p, p)
        }
    }
    maxp <- max(p)
    if (lag.restrict >= p[2]) {
        warning("lag.restrict >= p. Using lag.restrict = 0 instead.")
        lag.restrict <- 0
    }
    if (lag.restrict > 0) {
        lagX <- embed(y[, cause], maxp + 1)[, -c(1:(lag.restrict + 1)), drop = FALSE]
    } else {
        lagX <- embed(y[, cause], maxp + 1)[, -1, drop = FALSE]
    }
    lagY <- embed(y[, dep], maxp + 1)
    m_yx <- stats::lm.fit(x = cbind(1,
                                    lagY[, 2:(p[1] + 1), drop = FALSE],
                                    lagX[, 1:(p[2] - lag.restrict), drop = FALSE]),
                          y = lagY[, 1, drop = FALSE])
    m_y <- stats::lm.fit(x = cbind(1,
                                   lagY[, 2:(p[1] + 1), drop = FALSE]),
                         y = lagY[, 1, drop = FALSE])
    # ANOVA F-test for nested models
    SSEyx <- sum(m_yx$residuals^2)
    SSEy <- sum(m_y$residuals^2)
    dfyx <- m_yx$df.residual
    dfy <- m_y$df.residual
    Fobs <- ((SSEy - SSEyx) / (dfy - dfyx)) / (SSEyx/dfyx)
    pv <- pf(Fobs, dfy - dfyx, dfyx, lower.tail = FALSE)
    # Verify on slower estimation
    # tmp <- data.frame(cbind(lagY[, 1, drop = FALSE],
    #              lagY[, 2:(p[1] + 1), drop = FALSE],
    #              lagX[, 1:(p[2] - lag.restrict), drop = FALSE]))
    # m_yx <- stats::lm(tmp[,1] ~ ., data = tmp)
    # tmp <- data.frame(cbind(lagY[, 1, drop = FALSE],
    #                         lagY[, 2:(p[1] + 1), drop = FALSE]))
    # m_y <- stats::lm(tmp[,1] ~ ., data = tmp)
    # anova(m_yx, m_y)

    list(p = p,
         cause = cause,
         pv = pv)


}
