# functions/shock_helpers.R
#
# EA-MPD loading and the two shock decompositions.
#
# The Altavilla et al. (2019) dimensions ask which part of ECB communication
# moved markets, through Target, Timing, Forward Guidance and QE. The
# Jarocinski-Karadi (2020) decomposition asks what kind of news arrived,
# separating a monetary policy shock from information the ECB revealed about
# the economy. The two answer separate questions and each is reported on its
# own.
#
# The internally estimated rotations use the largest samples their inputs
# allow before the resulting series are aligned with the Romanian event
# window, which keeps the factor structure independent of the 135 events the
# regressions use.
#
# The published JK series is the primary measure. The internal implementation
# independently replicates the published JK construction, and 05 verifies the
# resulting series against the published file.

SHORT_TENORS <- c("OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y")
TERM_TENORS  <- c("OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y", "OIS_2Y", "OIS_5Y", "OIS_10Y")


# EA-MPD

# The date column arrives in two formats. Older events carry Excel serial
# numbers and recent ones DD/MM/YYYY text, so each group is parsed on its own
# terms; a single col_type would silently drop one of them.
parse_eampd_date <- function(x) {
  x   <- as.character(x)
  out <- as.Date(rep(NA_character_, length(x)))

  is_serial <- grepl("^\\d+(\\.0+)?$", x)
  out[is_serial] <- as.Date(as.numeric(x[is_serial]), origin = "1899-12-30")

  for (fmt in c("%d/%m/%Y", "%Y-%m-%d", "%Y-%m-%d %H:%M:%S")) {
    miss <- is.na(out)
    if (!any(miss)) break
    out[miss] <- as.Date(x[miss], format = fmt)
  }
  out
}

# Each sheet carries about 1027 columns and only the first 46 (A:AT) hold data.
# The explicit range keeps col_types aligned with the real columns; readxl
# otherwise reads the trailing padding and shifts them silently.
load_eampd_raw <- function(path = "data/raw/ECB_surprise_shocks.xlsx") {
  read_window <- function(sheet) {
    readxl::read_excel(
      path, sheet = sheet,
      range = readxl::cell_cols("A:AT"),
      col_types = c("text", rep("numeric", 45))
    ) |>
      dplyr::mutate(date = parse_eampd_date(date)) |>
      dplyr::filter(!is.na(date)) |>
      dplyr::arrange(date)
  }
  list(
    prw = read_window("Press Release Window"),
    pcw = read_window("Press Conference Window"),
    mew = read_window("Monetary Event Window")
  )
}

# OIS_* carries the risk-free surprise term structure and DE* the Bund
# reaction used for the intraday timing check in 02. STOXX50 enters the JK sign
# restriction, EURUSD covers the exchange-rate channel, and IT/ES serve as a
# periphery benchmark.
EVENT_COLS <- c(
  "date",
  "OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y", "OIS_2Y", "OIS_3Y", "OIS_5Y", "OIS_10Y",
  "DE2Y", "DE5Y", "DE10Y", "IT10Y", "ES10Y", "STOXX50", "EURUSD"
)

# ois_5y and ois_10y are missing for the early events, since the long end of
# the intraday surprise curve begins in August 2011. They stay NA and 03 drops
# the affected events from the family that needs them. A missing surprise says
# nothing about the size of the shock, so a filled zero would manufacture
# information.
event_surprises <- function(mew, start_date) {
  mew |>
    dplyr::filter(date >= start_date) |>
    dplyr::select(dplyr::all_of(EVENT_COLS)) |>
    dplyr::rename_with(tolower)
}


# Factor machinery

# Principal components of the covariance matrix. Every column already sits in
# basis points, so the covariance form preserves the relative volatility of
# each tenor, which standardising would flatten.
#
# Scores are whitened to unit variance and the scale absorbed into the
# loadings. PCA scores have covariance diag(eigenvalues), and an orthogonal
# rotation of them stays orthogonal only when that covariance is proportional
# to the identity; rotating raw scores leaves forward guidance and QE
# correlated at 0.73.
pca_cov <- function(X, k) {
  Xc   <- scale(X, center = TRUE, scale = FALSE)
  ev   <- eigen(cov(Xc), symmetric = TRUE)
  sd_j <- sqrt(ev$values[1:k])
  list(
    scores    = Xc %*% ev$vectors[, 1:k, drop = FALSE] %*% diag(1 / sd_j, k, k),
    loadings  = ev$vectors[, 1:k, drop = FALSE] %*% diag(sd_j, k, k),
    var_share = ev$values[1:k] / sum(ev$values)
  )
}

# With unit-variance orthogonal factors this returns cov(ref, f), which is the
# structural loading itself.
loading_on <- function(f, ref) as.numeric(stats::coef(stats::lm(ref ~ f))[2])

# Altavilla normalisation. Rescaling by the loading on `ref` sets that loading
# to exactly 1, so a coefficient on the factor reads as bp per 1 bp of the
# reference maturity. The direction is a multiplication. Given ref = b * f the
# wanted factor is b * f, and dividing would leave a loading of b^2, which
# rescales every reported magnitude while each t-statistic stays put.
to_bp <- function(f, ref) f * loading_on(f, ref)


# Altavilla-style dimensions

# Target comes from the Press Release window, the decision itself. Timing,
# Forward Guidance and QE come from the Press Conference window, rotated under
# two restrictions.
#
#   (i)  FG and QE have zero loading on the 1M OIS, so both leave the rate over
#        the current maintenance period alone. This pins Timing to the 1M
#        loading vector and leaves a plane for the other two.
#   (ii) within that plane, QE is the direction maximising the 10Y loading and
#        FG is its orthogonal complement, which places FG's 10Y loading at zero
#        by construction.
#
# Three restrictions for the three degrees of freedom of a 3x3 rotation, so the
# rotation is exactly identified.
#
# This is an Altavilla-style rotation, and the block carries that label
# throughout. Restriction (i) and the four normalisations below follow
# Altavilla et al. Restriction (ii) is this project's own. They identify QE as
# the direction with the smallest variance from 2 January 2002 to 7 August
# 2008, before the balance sheet became a policy instrument, which requires
# pre-crisis events. Before August 2011 the euro area long end exists only as
# German Bund yields, and this project keeps the OIS term structure pure, so
# its sample opens in August 2011 and the pre-crisis restriction lies outside
# its reach. The two rules diverge materially. Applied to the same plane, the
# published rule turns the QE direction substantially and gives forward
# guidance a sizeable 10Y loading in place of the zero imposed here.
#
# So the factors below answer which part of the term structure moved, under a
# stated rotation. The headline shocks throughout the project remain the
# published Jarocinski-Karadi series, which these factors leave untouched.
#
# On samples, the short-end OIS is complete for all 315 EA-MPD events back to
# 1999, so Target uses the full history. OIS_5Y and OIS_10Y begin on 4 August
# 2011, which confines any factor defined by the long end to the 128
# complete-term-structure events.
construct_altavilla_factors <- function(prw, pcw) {
  prw_cc <- prw[stats::complete.cases(prw[, SHORT_TENORS]), ]
  pt     <- pca_cov(as.matrix(prw_cc[, SHORT_TENORS]), 1)

  target <- as.vector(pt$scores[, 1])
  if (cor(target, prw_cc$OIS_1M) < 0) target <- -target
  target <- to_bp(target, prw_cc$OIS_1M)          # 1 unit = 1 bp on the 1M OIS

  pcw_cc <- pcw[stats::complete.cases(pcw[, TERM_TENORS]), ]
  pp     <- pca_cov(as.matrix(pcw_cc[, TERM_TENORS]), 3)
  L      <- pp$loadings
  Fm     <- pp$scores

  i1m  <- match("OIS_1M",  TERM_TENORS)
  i2y  <- match("OIS_2Y",  TERM_TENORS)
  i10y <- match("OIS_10Y", TERM_TENORS)

  u1 <- L[i1m, ] / sqrt(sum(L[i1m, ]^2))                  # Timing direction
  Qf <- qr.Q(qr(cbind(u1, diag(3))))[, 1:3]
  if (sum(Qf[, 1] * u1) < 0) Qf <- -Qf
  comp <- Qf[, 2:3, drop = FALSE]                          # plane orthogonal to u1

  w  <- as.vector(t(comp) %*% L[i10y, ])
  w  <- w / sqrt(sum(w^2))
  u3 <- as.vector(comp %*% w)                              # QE, max 10Y loading
  u2 <- as.vector(comp %*% c(-w[2], w[1]))                 # forward guidance

  U <- cbind(u1, u2, u3)
  stopifnot(max(abs(crossprod(U) - diag(3))) < 1e-10)

  Fr <- Fm %*% U
  Lr <- L  %*% U

  # a positive factor raises its own reference rate
  for (j in seq_len(3)) {
    ref <- c(i1m, i2y, i10y)[j]
    if (Lr[ref, j] < 0) { Fr[, j] <- -Fr[, j]; Lr[, j] <- -Lr[, j] }
  }

  # normalisation maturities follow Altavilla et al., with Target on 1M, Timing
  # on 6M, FG on 2Y and QE on 10Y
  refs     <- list(pcw_cc$OIS_6M, pcw_cc$OIS_2Y, pcw_cc$OIS_10Y)
  scale_bp <- vapply(seq_len(3), function(j) loading_on(Fr[, j], refs[[j]]), numeric(1))

  factors <- dplyr::full_join(
    data.frame(date = prw_cc$date, target = target),
    data.frame(
      date             = pcw_cc$date,
      timing           = Fr[, 1] * scale_bp[1],
      forward_guidance = Fr[, 2] * scale_bp[2],
      qe               = Fr[, 3] * scale_bp[3]
    ),
    by = "date"
  )

  # loadings on the normalised factors, so each reference cell is 1
  loadings <- sweep(Lr, 2, scale_bp, "/")
  dimnames(loadings) <- list(TERM_TENORS, c("timing", "forward_guidance", "qe"))

  list(
    factors    = dplyr::arrange(factors, date),
    unscaled   = data.frame(
      date = pcw_cc$date,
      timing_std = Fr[, 1], forward_guidance_std = Fr[, 2], qe_std = Fr[, 3]
    ),
    target_std = data.frame(date = prw_cc$date, target_std = as.vector(pt$scores[, 1])),
    loadings   = loadings,
    var_share  = list(target = pt$var_share, conference = pp$var_share),
    n          = c(target = nrow(prw_cc), conference = nrow(pcw_cc))
  )
}


# Jarocinski-Karadi

# Independent replication of the published Jarocinski-Karadi pipeline, from
# the raw EA-MPD workbook. Everything it returns carries the _own suffix and
# feeds the validation exercise in 05.
#
# The two helpers below mirror the authors' MATLAB algorithm, code/mypc.m and
# code/signrestr_median.m, so that there is no deliberate algorithmic
# difference between the replication and the authors' implementation. The one
# addition is an explicit PC1 sign orientation, which resolves the arbitrary
# sign an SVD returns so that the R implementation follows the authors'
# orientation.
#
# Both implementations round pc1 to eight decimals before the rotation, as
# main.m does. The _own implementation keeps the resulting median-rotation
# outputs at full precision, while the published file rounds MP and CBI to
# eight decimals at export. That final export rounding produces the residual
# deviation of roughly 5e-07 bp reported in 05.
#
# jk_mypc() follows mypc.m. Each tenor is divided by its own standard
# deviation, the principal components are taken without centering, rows missing
# every tenor return NA, and PC1 is rescaled to the standard deviation of the
# reference series. Note the standardisation, which the covariance-PCA used for
# the Altavilla-style factors does not apply.
jk_mypc <- function(tab, varlist, var2scale) {
  X     <- as.matrix(tab[, varlist])
  imiss <- rowSums(is.na(X)) == ncol(X)
  X[is.na(X)] <- 0
  X     <- sweep(X, 2, apply(X, 2, stats::sd), "/")
  sv    <- svd(X)
  v1    <- sv$v[, 1]
  if (sum(v1) < 0) v1 <- -v1              # orient PC1 to load positively
  score1 <- as.vector(X %*% v1)
  score1[imiss] <- NA
  list(
    pc1       = score1 / stats::sd(score1, na.rm = TRUE) *
                stats::sd(tab[[var2scale]], na.rm = TRUE),
    var_share = sv$d[1]^2 / sum(sv$d^2)
  )
}

# jk_signrestr_median() follows signrestr_median.m. M = [pc1, equity] is
# decomposed into two orthogonal shocks summing to pc1, the first lowering
# equities and the second raising them. The sign restrictions leave an interval
# of admissible rotation angles and the midpoint is taken, which is the
# authors' default w = 0.5. Both the data matrix and the QR factorisation are
# used undemeaned, as in the original.
jk_signrestr_median <- function(M, w = 0.5) {
  ok <- !is.na(rowSums(M))
  qd <- qr(M[ok, , drop = FALSE])
  Q  <- qr.Q(qd); R <- qr.R(qd)
  S  <- diag(sign(diag(R)), 2, 2)
  Q  <- Q %*% S; R <- S %*% R

  lo <- if (R[1, 2] > 0) atan(R[1, 2] / R[2, 2]) else 0
  hi <- if (R[1, 2] > 0) pi / 2 else atan(-R[2, 2] / R[1, 2])
  a  <- (1 - w) * lo + w * hi

  P <- matrix(c(cos(a), -sin(a), sin(a), cos(a)), 2, 2)
  D <- diag(c(R[1, 1] * cos(a), R[1, 1] * sin(a)))
  U <- matrix(NA_real_, nrow(M), 2)
  U[ok, ] <- Q %*% P %*% D
  list(U = U, arc = c(lo, hi), angle = a)
}

# main.m drops three announcements that the ECB made jointly with the Federal
# Reserve, before the component is extracted. They fall outside the Romanian
# window but inside the estimation sample, so the exclusion belongs here.
JK_JOINT_FED_DATES <- as.Date(c("2001-09-13", "2001-09-17", "2008-10-08"))

# PC1 of the short-end Monetary Event window surprises, scaled to the
# volatility of the 1Y OIS, combined with the equity surprise. Both shocks
# raise the policy rate; a monetary policy shock lowers equities, an
# information shock raises them. The median rotation lets both shocks occur in
# one event, while the poor man's rule assigns each event wholly to one and
# travels as a robustness check.
#
# main.m works in percentage points per annum, so the rotation runs on pc1/100
# beside the untouched equity surprise, and the results are returned in basis
# points to match the rest of the panel.
reconstruct_jk <- function(mew) {
  tab <- mew[!mew$date %in% JK_JOINT_FED_DATES, ]

  pcx <- jk_mypc(tab, SHORT_TENORS, "OIS_1Y")
  pc1 <- round(pcx$pc1 / 100, 8)                  # % p.a., as main.m exports
  sr  <- jk_signrestr_median(cbind(pc1, tab$STOXX50))

  keep <- !is.na(sr$U[, 1])
  list(
    shocks = data.frame(
      date          = tab$date[keep],
      jk_pc1_own    = 100 * pc1[keep],
      jk_mp_own     = 100 * sr$U[keep, 1],
      jk_cbi_own    = 100 * sr$U[keep, 2],
      jk_mp_pm_own  = 100 * ifelse(pc1[keep] * tab$STOXX50[keep] <  0, pc1[keep], 0),
      jk_cbi_pm_own = 100 * ifelse(pc1[keep] * tab$STOXX50[keep] >= 0, pc1[keep], 0)
    ),
    diagnostics = list(
      n              = sum(keep),
      pc1_var_share  = pcx$var_share,
      arc_degrees    = 180 / pi * sr$arc,
      median_degrees = 180 / pi * sr$angle,
      adds_up        = max(abs(sr$U[keep, 1] + sr$U[keep, 2] - pc1[keep]))
    )
  )
}

# The published decomposition, converted from percentage points to basis
# points. This function reads the published file alone, which keeps the primary
# measure independent of reconstruct_jk().
#
# The conversion is a unit change and the authors' own normalisation is
# retained. mypc.m rescales the first principal component to the standard
# deviation of the OIS_1Y surprise in the EA-MPD input, which is in basis
# points, and main.m then divides by 100.
#
#     pc1 = mypc(tab, irnames, "OIS_1Y")/100;
#
# The published pc1, MP and CBI therefore arrive in percentage points per
# annum, which the repository README states as pc1 "scaled to have the same
# standard deviation as the OIS1Y Monetary Event-window change (in % p.a.)".
# Multiplying by 100 recovers basis points.
#
# The operation to keep clear of this one is sd(OIS_1Y)/sd(pc1) computed on the
# events this project uses. That renormalises an already-normalised series to
# the volatility of a subsample, yields 102.11 here in place of 100, and
# inflates every reported magnitude by 2.1%. The equity column separates the
# two cases, since main.m leaves STOXX50 in EA-MPD units and it matches the
# workbook exactly, while pc1 sits a factor of 100 below.
#
# `since` keeps the diagnostics deterministic. The published file on GitHub
# runs from 1999 and this project uses the events from 2011 onwards. With the
# conversion now a constant, the filter governs the reported counts alone.
#
# Two counts answer separate questions. n_official is how many rows the file
# holds after the filter, and n_overlap how many of those the EA-MPD vintage
# also carries, which is the sample the equity check runs on.
#
# equity_max_abs_dev checks the STOXX50 input in the published file against
# the workbook read here. Agreement to machine precision confirms consistent
# parsing of that equity input, which is one of the two series the
# decomposition runs on; 05 compares the shocks themselves against
# reconstruct_jk().
load_jk_official <- function(mew, path = "data/raw/jk_shocks_official.csv",
                             since = as.Date("2011-01-01")) {
  if (!file.exists(path))
    stop("Missing ", path, ". The primary shock measure is the published ",
         "series. Download shocks_ecb_mpd_me_d.csv from ",
         "github.com/marekjarocinski/jkshocks_update_ecb and save it there.")

  off <- utils::read.csv(path)
  names(off) <- tolower(names(off))
  off$date <- as.Date(off$date)
  off <- off[off$date >= since, ]

  cc <- mew[stats::complete.cases(mew[, c(SHORT_TENORS, "STOXX50")]), ]
  m  <- merge(off, cc[, c("date", "OIS_1Y", "STOXX50")], by = "date")
  k  <- 100  # percentage points per annum -> basis points

  list(
    shocks = data.frame(
      date      = off$date,
      jk_pc1    = off$pc1        * k,
      jk_mp     = off$mp_median  * k,
      jk_cbi    = off$cbi_median * k,
      jk_mp_pm  = off$mp_pm      * k,
      jk_cbi_pm = off$cbi_pm     * k
    ),
    scale              = k,
    n_official         = nrow(off),
    n_overlap          = nrow(m),
    equity_max_abs_dev = max(abs(m$stoxx50 - m$STOXX50))
  )
}
