# ==============================================================================
# functions/econometrics.R
#
# OLS with heteroskedasticity- and cluster-robust standard errors, and Wald
# tests of linear restrictions. Used by 02-06.
#
# The sandwich estimators are written out rather than taken from the sandwich
# package, so the project takes no direct dependency beyond readxl, dplyr,
# tidyr and ggplot2, and the algebra is short enough to check directly. HC0
# carries no df correction, which is what makes the numbers comparable with the
# timing exercise in 02.
#
# p-values use the normal distribution and Wald statistics the chi-squared,
# including for the clustered fits. With 135 event clusters this is the usual
# asymptotic choice; in a much smaller cluster count it would not be.
# ==============================================================================


# OLS with a robust or cluster-robust covariance matrix.
#
#   HC0   (X'X)^-1 sum_i u_i^2 x_i x_i' (X'X)^-1
#   HC1   as HC0, scaled by n / (n - k)
#   HC3   u_i^2 / (1 - h_i)^2, the leverage correction that matters at this
#         sample size
#
# `cluster` names a column in `data`; when supplied it overrides `vcov` and the
# meat is summed over clusters with the G/(G-1) * (n-1)/(n-k) factor.
fit_ols_robust <- function(formula, data, vcov = "HC0", cluster = NULL) {
  stopifnot(vcov %in% c("HC0", "HC1", "HC3"))

  fit <- stats::lm(formula, data = data)
  X   <- stats::model.matrix(fit)
  u   <- as.vector(stats::residuals(fit))
  h   <- as.vector(stats::hatvalues(fit))
  n   <- nrow(X)
  k   <- ncol(X)

  # rows lm kept, so the cluster variable lines up with the estimation sample
  rows <- seq_len(nrow(data))
  if (!is.null(fit$na.action)) rows <- rows[-unname(fit$na.action)]

  XtXi <- solve(crossprod(X))

  if (is.null(cluster)) {
    w <- switch(vcov,
      HC0 = u^2,
      HC1 = u^2 * n / (n - k),
      HC3 = u^2 / (1 - h)^2
    )
    meat <- crossprod(X * sqrt(w))
  } else {
    grp  <- split(seq_len(n), data[[cluster]][rows])
    meat <- Reduce(`+`, lapply(grp, function(ix)
      tcrossprod(colSums(X[ix, , drop = FALSE] * u[ix]))))
    G    <- length(grp)
    meat <- meat * (G / (G - 1)) * ((n - 1) / (n - k))
  }

  V  <- XtXi %*% meat %*% XtXi
  b  <- stats::coef(fit)
  se <- sqrt(diag(V))

  structure(
    list(
      b = b, se = se, t = b / se, p = 2 * stats::pnorm(-abs(b / se)),
      V = V, n = n, k = k, r2 = summary(fit)$r.squared,
      u = u, h = h, XtXi = XtXi, X = X, rows = rows,
      vcov = vcov, formula = formula, n_clusters = if (is.null(cluster)) NA_integer_ else G
    ),
    class = "ols_robust"
  )
}


# Coefficient table; `terms` selects and orders a subset.
tidy_ols <- function(fit, terms = NULL) {
  out <- data.frame(
    term = names(fit$b), estimate = unname(fit$b), se = unname(fit$se),
    t = unname(fit$t), p = unname(fit$p), n = fit$n, r2 = fit$r2,
    row.names = NULL, stringsAsFactors = FALSE
  )
  if (!is.null(terms)) out <- out[match(terms, out$term), , drop = FALSE]
  rownames(out) <- NULL
  out
}

coef_of <- function(fit, term) unname(fit$b[term])
t_of    <- function(fit, term) unname(fit$t[term])


# Restriction matrix from named weight vectors, so that a hypothesis reads as
# the algebra it is:  restrictions(names(fit$b), c(jk_mp = 1, jk_cbi = -1))
restrictions <- function(term_names, ...) {
  rows <- list(...)
  R <- matrix(0, length(rows), length(term_names),
              dimnames = list(NULL, term_names))
  for (i in seq_along(rows)) {
    stopifnot(all(names(rows[[i]]) %in% term_names))
    R[i, names(rows[[i]])] <- rows[[i]]
  }
  R
}

# Wald test of R b = r. Uses the covariance matrix carried by `fit`, so
# switching a specification to HC3 or clustered errors also switches the test.
wald_test <- function(fit, R, r = rep(0, nrow(R))) {
  R <- as.matrix(R)
  if (!is.null(colnames(R))) R <- R[, names(fit$b), drop = FALSE]
  d <- as.vector(R %*% fit$b) - r
  W <- as.numeric(t(d) %*% solve(R %*% fit$V %*% t(R)) %*% d)
  c(chi2 = W, df = nrow(R), p = stats::pchisq(W, nrow(R), lower.tail = FALSE))
}
