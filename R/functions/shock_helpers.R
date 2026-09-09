# ==============================================================================
# functions/shock_helpers.R
#
# EA-MPD loading and the two shock decompositions.
#
#   Altavilla et al. (2019)   which dimension of ECB communication moved
#                             markets: Target, Timing, Forward Guidance, QE
#   Jarocinski-Karadi (2020)  what kind of news it was: monetary policy, or
#                             information the ECB revealed about the economy
#
# Different questions, never nested or combined. Both rotations are estimated
# on the largest sample the data allow and only then subset to the Romanian
# window - never re-estimated on 135 events.
#
# The published JK series is the primary measure. The reconstruction here is an
# audit of the method, not a substitute for the original; 05 compares them.
# ==============================================================================

SHORT_TENORS <- c("OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y")
TERM_TENORS  <- c("OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y", "OIS_2Y", "OIS_5Y", "OIS_10Y")


# ------------------------------------------------------------------------------
# EA-MPD
# ------------------------------------------------------------------------------

# Mixed-format date column: older events are Excel serial numbers, recent ones
# DD/MM/YYYY text. Reading it as a single type silently loses one group.
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

# Each sheet carries about 1027 columns of which only the first 46 (A:AT) hold
# data. Without an explicit range, readxl reads the padding and col_types
# misaligns without erroring.
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

# OIS_* is the risk-free surprise term structure; DE* the Bund reaction, used
# to validate the shock series; STOXX50 is required by the JK sign restriction;
# EURUSD is the exchange-rate channel; IT/ES an optional periphery benchmark.
EVENT_COLS <- c(
  "date",
  "OIS_1M", "OIS_3M", "OIS_6M", "OIS_1Y", "OIS_2Y", "OIS_3Y", "OIS_5Y", "OIS_10Y",
  "DE2Y", "DE5Y", "DE10Y", "IT10Y", "ES10Y", "STOXX50", "EURUSD"
)

# ois_5y and ois_10y are missing for the early events - the long end of the
# intraday surprise curve does not exist before August 2011. Left as NA rather
# than filled: a missing surprise is not a zero shock, and 03 drops the
# affected events instead of imputing them.
event_surprises <- function(mew, start_date) {
  mew |>
    dplyr::filter(date >= start_date) |>
    dplyr::select(dplyr::all_of(EVENT_COLS)) |>
    dplyr::rename_with(tolower)
}


# ------------------------------------------------------------------------------
# Factor machinery
# ------------------------------------------------------------------------------

# Principal components of the covariance matrix - covariance and not
# correlation, since every column is already in basis points and standardising
# would discard the relative volatility of each tenor.
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

# With unit-variance orthogonal factors this is cov(ref, f), the structural
# loading, not merely a fitted slope.
loading_on <- function(f, ref) as.numeric(stats::coef(stats::lm(ref ~ f))[2])

# Altavilla normalisation: rescale so the loading on `ref` is exactly 1, and a
# coefficient on the factor reads as bp per 1 bp of that reference maturity.
# Note the direction - multiply. If ref = b * f the factor wanted is b * f, so
# dividing leaves a loading of b^2, which corrupts every reported magnitude
# while leaving every t-statistic untouched.
to_bp <- function(f, ref) f * loading_on(f, ref)


# ------------------------------------------------------------------------------
# Altavilla-style dimensions
# ------------------------------------------------------------------------------

# Target from the Press Release window, the decision itself. Timing, Forward
# Guidance and QE from the Press Conference window, rotated under:
#
#   (i)  FG and QE have zero loading on the 1M OIS - neither moves the rate
#        over the current maintenance period. This pins Timing to the 1M
#        loading vector and leaves a plane for the other two.
#   (ii) within that plane, QE maximises the 10Y loading; FG is its orthogonal
#        complement, which mechanically sets FG's 10Y loading to zero.
#
# Three restrictions for the three degrees of freedom of a 3x3 rotation, so the
# rotation is exactly identified.
#
# THIS IS NOT THE ALTAVILLA ET AL. IDENTIFICATION, and the block is labelled
# Altavilla-style throughout for that reason. Restriction (i) and the four
# normalisations below are theirs. Restriction (ii) is not: they split the same
# plane by requiring QE to be the direction with the smallest variance from
# 2 January 2002 to 7 August 2008, before the balance sheet became a policy
# instrument. That restriction cannot be imposed on this sample. It needs
# pre-crisis events, and before August 2011 the euro area long end exists only
# as German Bund yields, which this project does not substitute into the OIS
# term structure. The two rules are not interchangeable: applied to the same
# plane, the published rule turns the QE direction substantially and gives
# forward guidance a sizeable 10Y loading instead of the zero imposed here.
#
# What follows from that: the factors below answer "which part of the term
# structure moved" under a stated rotation, not "what Altavilla et al. call
# Timing, FG and QE". Nothing else in the project depends on them - the
# headline shocks are the published Jarocinski-Karadi series.
#
# Sample: short-end OIS is complete for all 315 EA-MPD events back to 1999, so
# Target uses the full history. OIS_5Y and OIS_10Y do not exist before 4 August
# 2011, so any factor defined by the long end is confined to the 128
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

  # normalisation maturities follow Altavilla et al.: Target 1M, Timing 6M,
  # FG 2Y, QE 10Y
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


# ------------------------------------------------------------------------------
# Jarocinski-Karadi
# ------------------------------------------------------------------------------

# PC1 of the short-end Monetary Event window surprises, scaled to the
# volatility of the 1Y OIS, combined with the equity surprise. Both shocks
# raise the policy rate; a monetary policy shock lowers equities, an
# information shock raises them.
#
# The sign restrictions leave an arc of admissible rotations rather than a
# point - Sigma has three distinct elements and B four parameters, so one
# degree of freedom survives. The median rotation is the representative taken
# here and allows both shocks in one event; the poor man's rule assigns each
# event wholly to one and is carried as a robustness check.
#
# Everything this returns carries the _own suffix: it is the audit measure,
# never the baseline.
reconstruct_jk <- function(mew) {
  cc <- mew[stats::complete.cases(mew[, c(SHORT_TENORS, "STOXX50")]), ]
  pj <- pca_cov(as.matrix(cc[, SHORT_TENORS]), 1)

  pc1 <- as.vector(pj$scores[, 1])
  if (cor(pc1, cc$OIS_1Y) < 0) pc1 <- -pc1
  pc1 <- pc1 / sd(pc1) * sd(cc$OIS_1Y)

  p1 <- pc1 - mean(pc1)
  sx <- cc$STOXX50 - mean(cc$STOXX50)

  P  <- t(chol(cov(cbind(p1, sx))))
  th <- seq(0, 2 * pi, length.out = 200001)
  ct <- cos(th); st <- sin(th)
  adm <- (P[1, 1] * ct > 0) &
         (P[2, 1] * ct + P[2, 2] * st < 0) &
         (-P[1, 1] * st > 0) &
         (-P[2, 1] * st + P[2, 2] * ct > 0)
  stopifnot(any(adm))

  th_med <- median(th[adm])
  Q <- matrix(c(cos(th_med), sin(th_med), -sin(th_med), cos(th_med)), 2, 2)
  B <- P %*% Q
  e <- solve(B) %*% t(cbind(p1, sx))

  list(
    shocks = data.frame(
      date          = cc$date,
      jk_pc1_own    = pc1,
      jk_mp_own     = B[1, 1] * e[1, ],
      jk_cbi_own    = B[1, 2] * e[2, ],
      jk_mp_pm_own  = ifelse(pc1 * cc$STOXX50 <  0, pc1, 0),
      jk_cbi_pm_own = ifelse(pc1 * cc$STOXX50 >= 0, pc1, 0)
    ),
    diagnostics = list(
      n              = nrow(cc),
      pc1_var_share  = pj$var_share[1],
      arc_degrees    = 180 / pi * range(th[adm]),
      median_degrees = 180 / pi * th_med,
      adds_up        = max(abs(B[1, 1] * e[1, ] + B[1, 2] * e[2, ] - p1))
    )
  )
}

# The published decomposition, converted from percentage points to basis
# points. Makes no reference to reconstruct_jk(), so the primary measure is
# independent of the reconstruction rather than calibrated to it.
#
# THE CONVERSION IS A UNIT CHANGE, NOT A NORMALISATION. The authors have
# already normalised: mypc.m rescales the first principal component to the
# standard deviation of the OIS_1Y surprise in the EA-MPD input, which is in
# basis points, and main.m then divides by 100:
#
#     pc1 = mypc(tab, irnames, "OIS_1Y")/100;
#
# so the published pc1, MP and CBI are in percentage points per annum - the
# README describes pc1 as "scaled to have the same standard deviation as the
# OIS1Y Monetary Event-window change (in % p.a.)". Multiplying by 100 returns
# basis points and nothing else.
#
# What this must NOT be is sd(OIS_1Y)/sd(pc1) computed on the events this
# project happens to use. That looks like the same operation and is not: it
# renormalises an already-normalised series to the volatility of a subsample,
# which on this sample gives 102.11 rather than 100 and inflates every reported
# magnitude by 2.1%. The equity column is the check that the two are different
# operations - main.m leaves STOXX50 alone, and it matches the workbook exactly,
# while pc1 does not.
#
# `since` keeps the diagnostics deterministic: the published file on GitHub runs
# from 1999 and this project uses the events from 2011 onwards. It no longer
# affects any coefficient, since the conversion is now a constant.
#
# Two counts are returned because they answer different questions: n_official
# is how many rows the file holds after the filter, n_overlap how many of them
# the EA-MPD vintage also carries, which is the sample the equity check runs on.
#
# equity_max_abs_dev compares the STOXX50 surprise the published file carries
# against the same column in the workbook read here. It is the check that the
# workbook is being read correctly: if the two files agree on the equity input
# to machine precision, a discrepancy in the shocks is about the decomposition,
# not about parsing.
load_jk_official <- function(mew, path = "data/raw/jk_shocks_official.csv",
                             since = as.Date("2011-01-01")) {
  if (!file.exists(path))
    stop("Missing ", path, ". The primary shock measure is the published ",
         "series: take shocks_ecb_mpd_me_d.csv from ",
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
