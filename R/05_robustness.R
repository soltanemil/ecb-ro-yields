# ==============================================================================
# 05_robustness.R
#
# Audit of the information-shock result. Everything here asks whether it
# survives, not whether it can be made larger.
#
# In:   data/processed/analysis_panel.rds, data/processed/shock_validation.rds
# Out:  tables/robustness_*.csv, figures/robustness_loo.png
#
# The four quoted maturities are one family of correlated outcomes, not four
# independent experiments: the stacked omnibus is the joint test of the family
# and the per-maturity tables are description. Level, Slope and Curvature are
# transformations of the same yields and are not counted again.
#
# HC0 stays the headline; HC3 appears where leverage is material, above all in
# the regime block. No observation is dropped on statistical grounds - the
# exclusions in section 3 are chosen ex ante by historical episode, and
# influence is diagnosed and reported rather than acted on.
#
# Regime and Late indicators are defined here because they are analysis
# choices; everything factual about the data was settled in 01.
# ==============================================================================

library(dplyr)
library(ggplot2)

source("R/functions/econometrics.R")
source("R/functions/data_helpers.R")

dat  <- readRDS("data/processed/analysis_panel.rds")
shv  <- readRDS("data/processed/shock_validation.rds")
JK   <- c("jk_mp", "jk_cbi")
BASE <- reformulate(JK, response = "dy_10y")


# ------------------------------------------------------------------------------
# Leave one event out
# ------------------------------------------------------------------------------

# What matters is whether the coefficient stays negative and of similar size
# across all 135 refits, not whether every refit clears a p-value threshold.
cat("=== Leave-one-event-out, information-shock coefficient ===\n")

loo <- list()
for (o in c("dy_10y", "dy_level")) {
  f    <- reformulate(JK, response = o)
  full <- fit_ols_robust(f, dat)
  est  <- t(vapply(seq_len(nrow(dat)), function(i) {
    fi <- fit_ols_robust(f, dat[-i, ])
    c(beta = coef_of(fi, "jk_cbi"), t = t_of(fi, "jk_cbi"))
  }, numeric(2)))
  loo[[o]] <- data.frame(date = dat$date, beta = est[, 1], t = est[, 2])

  lab <- if (o == "dy_10y") "10Y" else "Level"
  cat(sprintf("\n%-6s full sample: beta = %.3f (t = %.2f)\n",
              lab, coef_of(full, "jk_cbi"), t_of(full, "jk_cbi")))
  cat(sprintf("  beta : min %.3f | p25 %.3f | median %.3f | p75 %.3f | max %.3f\n",
              min(est[, 1]), quantile(est[, 1], .25), median(est[, 1]),
              quantile(est[, 1], .75), max(est[, 1])))
  cat(sprintf("  t    : min %.2f | median %.2f | max %.2f | refits with t < -1.96: %d of %d\n",
              min(est[, 2]), median(est[, 2]), max(est[, 2]),
              sum(est[, 2] < -1.96), nrow(est)))
  cat(sprintf("  most influential drop: %s (beta -> %.3f)\n",
              format(dat$date[which.max(abs(est[, 1] - coef_of(full, "jk_cbi")))]),
              est[which.max(abs(est[, 1] - coef_of(full, "jk_cbi"))), 1]))
}

f10  <- fit_ols_robust(BASE, dat)
s2   <- sum(f10$u^2) / (f10$n - f10$k)
cook <- (f10$u^2 / (f10$k * s2)) * f10$h / (1 - f10$h)^2
dfb  <- (f10$XtXi %*% t(f10$X))["jk_cbi", ] * f10$u / (1 - f10$h)

cat("\n--- Influence diagnostics, 10Y ---\n")
print(data.frame(date = dat$date[f10$rows], leverage = f10$h,
                 cooks_d = cook, dfbeta_cbi = dfb) |>
        arrange(desc(abs(dfbeta_cbi))) |> head(6),
      row.names = FALSE, digits = 3)
cat(sprintf("Cook's distance above 4/n = %.3f: %d observations\n",
            4 / f10$n, sum(cook > 4 / f10$n)))

write.csv(loo[["dy_10y"]], "tables/robustness_loo.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# Historical exclusions, chosen ex ante
# ------------------------------------------------------------------------------

# Not outlier trimming: whether the result is a general phenomenon or an
# artefact of two identifiable crisis episodes.
cat("\n=== Historical exclusions ===\n\n")

yr      <- format(dat$date, "%Y")
subsets <- list(
  "full sample"              = rep(TRUE, nrow(dat)),
  "excl. 2011-2012"          = !(yr %in% c("2011", "2012")),
  "excl. Mar-Apr 2020"       = !(dat$date >= as.Date("2020-03-01") &
                                 dat$date <= as.Date("2020-04-30")),
  "excl. all of 2020"        = yr != "2020",
  "excl. 2011-2012 and 2020" = !(yr %in% c("2011", "2012", "2020"))
)

exclusions <- do.call(rbind, lapply(names(subsets), function(nm) {
  d <- dat[subsets[[nm]], ]
  do.call(rbind, lapply(c("dy_10y", "dy_level"), function(o) {
    fit <- fit_ols_robust(reformulate(JK, response = o), d)
    data.frame(subset = nm, outcome = if (o == "dy_10y") "10Y" else "Level",
               n = fit$n,
               beta_cbi = coef_of(fit, "jk_cbi"), t_cbi = t_of(fit, "jk_cbi"),
               beta_mp  = coef_of(fit, "jk_mp"),  t_mp  = t_of(fit, "jk_mp"),
               row.names = NULL)
  }))
}))
print(exclusions, row.names = FALSE, digits = 3)
write.csv(exclusions, "tables/robustness_exclusions.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# Regimes
# ------------------------------------------------------------------------------

# One interacted model with a Wald test of equality, rather than a comparison
# of significance stars across four underpowered subsamples.
cat("\n=== Regimes ===\n\n")

REGIMES <- c(y2011_2014 = "2011-2014", y2015_2019 = "2015-2019",
             y2020_2021 = "2020-2021", y2022_2025 = "2022-2025")
brk <- as.Date(c("2011-01-01", "2015-01-01", "2020-01-01", "2022-01-01", "2026-12-31"))
dat$regime <- cut(dat$date, breaks = brk, labels = names(REGIMES), right = FALSE)

# interaction columns built explicitly, so every coefficient in the Wald
# restrictions below has a name that says what it is
for (r in names(REGIMES)) {
  d <- as.numeric(dat$regime == r)
  dat[[paste0("a_",   r)]] <- d
  dat[[paste0("mp_",  r)]] <- d * dat$jk_mp
  dat[[paste0("cbi_", r)]] <- d * dat$jk_cbi
}
reg_terms <- c(paste0("a_", names(REGIMES)),
               paste0("mp_", names(REGIMES)),
               paste0("cbi_", names(REGIMES)))
reg_f  <- as.formula(paste("dy_10y ~ 0 +", paste(reg_terms, collapse = " + ")))
reg_h0 <- fit_ols_robust(reg_f, dat, vcov = "HC0")
reg_h3 <- fit_ols_robust(reg_f, dat, vcov = "HC3")

regimes <- data.frame(
  regime    = unname(REGIMES),
  n         = as.vector(table(dat$regime)),
  beta_cbi  = unname(reg_h0$b[paste0("cbi_", names(REGIMES))]),
  t_cbi     = unname(reg_h0$t[paste0("cbi_", names(REGIMES))]),
  t_cbi_hc3 = unname(reg_h3$t[paste0("cbi_", names(REGIMES))]),
  beta_mp   = unname(reg_h0$b[paste0("mp_",  names(REGIMES))]),
  t_mp      = unname(reg_h0$t[paste0("mp_",  names(REGIMES))])
)
print(regimes, row.names = FALSE, digits = 3)

# HC3 is shown here and not elsewhere because the 2020-2021 window holds 16
# events and three parameters, so leverage is extreme exactly where the
# coefficient is most extreme.
R_eq <- restrictions(
  names(reg_h0$b),
  setNames(c(1, -1), c("cbi_y2011_2014", "cbi_y2015_2019")),
  setNames(c(1, -1), c("cbi_y2011_2014", "cbi_y2020_2021")),
  setNames(c(1, -1), c("cbi_y2011_2014", "cbi_y2022_2025"))
)
w0 <- wald_test(reg_h0, R_eq)
w3 <- wald_test(reg_h3, R_eq)
cat(sprintf("\nWald, CBI equal across regimes: chi2(%d) = %.2f, p = %.4f (HC0)\n",
            w0["df"], w0["chi2"], w0["p"]))
cat(sprintf("                                chi2(%d) = %.2f, p = %.4f (HC3)\n",
            w3["df"], w3["chi2"], w3["p"]))
cat("\nThe pooled negative result survives the leave-one-out and ex ante\n",
    "exclusion checks in sections 1 and 2, but the regime estimates are\n",
    "heterogeneous and should not be read as a stable parameter. The\n",
    "coefficient is negative in the three regimes estimated with useful\n",
    "precision. The positive 2020-21 estimate is leverage-sensitive and\n",
    "statistically uninformative under HC3 - a statement about precision, not\n",
    "about what happened.\n", sep = "")
write.csv(regimes, "tables/robustness_regimes.csv", row.names = FALSE)

# The estimand here is the difference between the 2011-14 and 2022-25 slopes.
# 2015-2021 is dropped from this regression, so delta is a contrast between the
# two endpoint regimes and NOT a full-sample post-2022 interaction, which would
# need every event and a date >= 2022 indicator. Reported and labelled as the
# endpoint contrast everywhere, in RESULTS.md and the README too.
cat("\n--- Endpoint regimes: 2011-2014 against 2022-2025 ---\n")
late <- dat |>
  filter(regime %in% c("y2011_2014", "y2022_2025")) |>
  mutate(late = as.numeric(regime == "y2022_2025"),
         mp_late = jk_mp * late, cbi_late = jk_cbi * late)

for (o in c("dy_10y", "dy_level")) {
  f  <- reformulate(c("late", JK, "mp_late", "cbi_late"), response = o)
  h0 <- fit_ols_robust(f, late, vcov = "HC0")
  h3 <- fit_ols_robust(f, late, vcov = "HC3")
  cat(sprintf("%-6s gamma_CBI (2011-14) = %6.3f | delta = %6.3f  t(HC0) = %5.2f  t(HC3) = %5.2f  p(HC3) = %.4f\n",
              if (o == "dy_10y") "10Y" else "Level",
              coef_of(h0, "jk_cbi"), coef_of(h0, "cbi_late"),
              t_of(h0, "cbi_late"), t_of(h3, "cbi_late"),
              unname(h3$p["cbi_late"])))
}
cat("N =", nrow(late), "events, 2011-2014 and 2022-2025 only. The 10Y and the\n")
cat("Level disagree, so this stays heterogeneity rather than a temporal\n")
cat("narrative.\n")


# ------------------------------------------------------------------------------
# The reconstruction against the published series
# ------------------------------------------------------------------------------

# The published series is primary. This asks whether the independent
# reconstruction agrees with it.
cat("\n=== Reconstruction against the published series ===\n\n")

pairs <- list(c("jk_pc1_own", "jk_pc1"), c("jk_mp_own", "jk_mp"),
               c("jk_cbi_own", "jk_cbi"), c("jk_mp_pm_own", "jk_mp_pm"),
               c("jk_cbi_pm_own", "jk_cbi_pm"))
validation <- do.call(rbind, lapply(pairs, function(pr) {
  x <- dat[[pr[1]]]; z <- dat[[pr[2]]]
  sl <- as.numeric(coef(lm(z ~ x + 0)))
  d  <- z - sl * x
  data.frame(reconstruction = pr[1], published = pr[2], corr = cor(x, z),
             scale = sl, max_abs_dev = max(abs(d)),
             rmse_rel_sd = sqrt(mean(d^2)) / sd(z), row.names = NULL)
}))
print(validation, row.names = FALSE, digits = 4)

dis <- (dat$jk_cbi_pm_own != 0) != (dat$jk_cbi_pm != 0)
cat(sprintf("\nPoor man's classification agreement: %.1f%% of %d events\n",
            100 * mean(!dis), nrow(dat)))
cat(sprintf("Median |pc1| at the %d disagreements: %.3f bp, against %.3f bp where\n",
            sum(dis), median(abs(dat$jk_pc1_own[dis])),
            median(abs(dat$jk_pc1_own[!dis]))))
cat("the two agree - the rule turns on the sign of a product, so it flips on\n")
cat("economically empty surprises.\n")
cat(sprintf(paste0(
  "\nThe equity surprise in the published file and in the workbook read here\n",
  "agree to %.1e over the %d shared events, so the workbook is parsed\n",
  "correctly and the remaining gap is about the decomposition, not the input.\n"),
  shv$jk$equity_max_abs_dev, shv$jk$overlap_n))
cat("That gap is consistent with the principal component being estimated on a\n",
    "different event set - this EA-MPD vintage carries ", shv$jk$diagnostics$n,
    " events, the published\nfile its own - but nothing here isolates that as ",
    "the cause. The published\nseries is primary for that reason.\n", sep = "")
write.csv(validation, "tables/robustness_jk_validation.csv", row.names = FALSE)

cat("\n--- The headline on each measure ---\n")
cat(sprintf("%-7s %26s %26s\n", "", "published (primary)", "reconstruction"))
for (o in c("dy_12m", "dy_3y", "dy_5y", "dy_10y", "dy_level")) {
  a <- fit_ols_robust(reformulate(JK, response = o), dat)
  b <- fit_ols_robust(reformulate(c("jk_mp_own", "jk_cbi_own"), response = o), dat)
  cat(sprintf("%-7s CBI %7.3f (t = %5.2f)     CBI %7.3f (t = %5.2f)\n",
              toupper(sub("dy_", "", o)),
              coef_of(a, "jk_cbi"), t_of(a, "jk_cbi"),
              coef_of(b, "jk_cbi_own"), t_of(b, "jk_cbi_own")))
}


# ------------------------------------------------------------------------------
# Median rotation against the poor man's restriction
# ------------------------------------------------------------------------------

cat("\n=== Identification scheme ===\n")
scheme <- do.call(rbind, lapply(names(MATURITIES), function(o) {
  m <- fit_ols_robust(reformulate(JK, response = o), dat)
  p <- fit_ols_robust(reformulate(c("jk_mp_pm", "jk_cbi_pm"), response = o), dat)
  data.frame(maturity = MATURITIES[[o]],
             cbi_median = coef_of(m, "jk_cbi"),   t_median = t_of(m, "jk_cbi"),
             cbi_pm     = coef_of(p, "jk_cbi_pm"), t_pm     = t_of(p, "jk_cbi_pm"),
             row.names = NULL)
}))
print(scheme, row.names = FALSE, digits = 3)
cat("These are two selection rules inside the same JK sign-restricted set, not\n")
cat("two independent identifications: the median rotation splits every event\n")
cat("between the two shocks, the poor man's rule assigns each event wholly to\n")
cat("one. Agreement means the result does not hinge on which rule is used.\n")


# ------------------------------------------------------------------------------
# Are the two shocks actually different?
# ------------------------------------------------------------------------------

# If they produced the same Romanian response the decomposition would add
# nothing - one could simply sum them.
cat("\n=== Wald test of gamma_MP = gamma_CBI ===\n\n")

mp_vs_cbi <- do.call(rbind, lapply(names(MATURITIES), function(o) {
  fit <- fit_ols_robust(reformulate(JK, response = o), dat)
  w   <- wald_test(fit, restrictions(names(fit$b), c(jk_mp = 1, jk_cbi = -1)))
  data.frame(maturity = MATURITIES[[o]],
             mp = coef_of(fit, "jk_mp"), cbi = coef_of(fit, "jk_cbi"),
             difference = coef_of(fit, "jk_mp") - coef_of(fit, "jk_cbi"),
             chi2 = w[["chi2"]], p = w[["p"]], row.names = NULL)
}))
print(mp_vs_cbi, row.names = FALSE, digits = 3)
write.csv(mp_vs_cbi, "tables/robustness_mp_vs_cbi.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# Omnibus test across the term structure
# ------------------------------------------------------------------------------

# A reshape for a joint test, not a data step: the four maturities are one
# family of correlated outcomes, estimated together with errors clustered on
# the ECB event.
cat("\n=== Omnibus test across the term structure ===\n")

long <- do.call(rbind, lapply(names(MATURITIES), function(o)
  data.frame(event = seq_len(nrow(dat)), maturity = MATURITIES[[o]],
             dy = dat[[o]], jk_mp = dat$jk_mp, jk_cbi = dat$jk_cbi)))
for (m in unname(MATURITIES)) {
  d <- as.numeric(long$maturity == m)
  long[[paste0("a_",   m)]] <- d
  long[[paste0("mp_",  m)]] <- d * long$jk_mp
  long[[paste0("cbi_", m)]] <- d * long$jk_cbi
}
mats     <- unname(MATURITIES)
stack_f  <- as.formula(paste("dy ~ 0 +", paste(c(
  paste0("a_", mats), paste0("mp_", mats), paste0("cbi_", mats)), collapse = " + ")))
stacked  <- fit_ols_robust(stack_f, long, vcov = "HC0", cluster = "event")

cat("Event-clustered standard errors,", stacked$n, "observations in",
    stacked$n_clusters, "clusters\n\n")
print(tidy_ols(stacked, c(paste0("mp_", mats), paste0("cbi_", mats)))[, 1:4],
      row.names = FALSE, digits = 3)

cbi_terms <- paste0("cbi_", mats)
w_zero <- wald_test(stacked, restrictions(names(stacked$b),
  setNames(1, cbi_terms[1]), setNames(1, cbi_terms[2]),
  setNames(1, cbi_terms[3]), setNames(1, cbi_terms[4])))
w_flat <- wald_test(stacked, restrictions(names(stacked$b),
  setNames(c(1, -1), cbi_terms[c(1, 2)]),
  setNames(c(1, -1), cbi_terms[c(1, 3)]),
  setNames(c(1, -1), cbi_terms[c(1, 4)])))

cat(sprintf("\nH0: all four CBI coefficients are zero  -> chi2(%d) = %.2f, p = %.4f\n",
            w_zero["df"], w_zero["chi2"], w_zero["p"]))
cat(sprintf("H0: the four CBI coefficients are equal -> chi2(%d) = %.2f, p = %.4f\n",
            w_flat["df"], w_flat["chi2"], w_flat["p"]))
cat("The first is the joint test of the maturity family; the second asks\n")
cat("whether the effect strengthens along the curve.\n")


# ------------------------------------------------------------------------------
# Economic size
# ------------------------------------------------------------------------------

# A coefficient in bp per bp is not interpretable without the shock's spread.
cat("\n=== How large is a typical shock? ===\n\n")
for (v in JK) {
  x <- dat[[v]]
  cat(sprintf("%-8s sd = %5.2f bp | p5 %6.2f  p25 %6.2f  p75 %6.2f  p95 %6.2f\n",
              v, sd(x), quantile(x, .05), quantile(x, .25),
              quantile(x, .75), quantile(x, .95)))
}
cat("\nOne standard deviation effects:\n")
for (o in c("dy_10y", "dy_level")) {
  fit <- fit_ols_robust(reformulate(JK, response = o), dat)
  cat(sprintf("  %-6s MP: %+6.2f bp   CBI: %+6.2f bp\n",
              if (o == "dy_10y") "10Y" else "Level",
              coef_of(fit, "jk_mp") * sd(dat$jk_mp),
              coef_of(fit, "jk_cbi") * sd(dat$jk_cbi)))
}


# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

p <- ggplot(loo[["dy_10y"]], aes(beta)) +
  geom_histogram(bins = 30, fill = "#1f4e79", colour = "white") +
  geom_vline(xintercept = coef_of(f10, "jk_cbi"), colour = "#c0762f", linewidth = 0.8) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  labs(
    title    = "Leave-one-event-out: the information-shock coefficient on the 10Y",
    subtitle = "135 refits, each dropping one Governing Council event; orange line is the full-sample estimate",
    caption  = "A result that depends on no single event has a tight distribution well away from zero.",
    x = "Coefficient on the information shock (bp per bp)", y = "Refits"
  ) +
  theme_minimal(base_size = 11)

ggsave("figures/robustness_loo.png", p, width = 9, height = 5, dpi = 150)
