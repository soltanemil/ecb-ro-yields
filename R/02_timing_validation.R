# ==============================================================================
# 02_timing_validation.R
#
# Establishes the response window used by everything downstream, by
# falsification rather than assertion. Four conditions have to hold jointly:
# the Bund reacts inside the EA-MPD intraday window; the response on the
# Romanian fixing struck before the announcement is small and statistically
# indistinguishable from zero; the first fixing able to incorporate the
# announcement does respond; and the estimate for the window after that is
# small and statistically indistinguishable from zero. Non-rejection is not a demonstrated zero, and the placebo windows
# are read that way throughout.
#
# In:   data/processed/analysis_panel.rds
# Out:  tables/timing_validation.csv, figures/timing_validation.png
#
# Runs on the raw OIS 2Y surprise, not a decomposed shock - the point is to
# test the measurement convention before any identification. HC0 throughout;
# this exercise is the reference point every later number is compared against.
#
# Reference values, 135 events. A deviation means the pipeline changed and the
# write-up needs re-checking:
#   pre       12M  0.031  3Y  0.039  5Y  0.081  10Y  0.104   (all |t| < 1)
#   response  12M  0.193  3Y  0.245  5Y  0.273  10Y  0.394 (t = 2.38)
#   post      12M  0.098  3Y  0.126  5Y  0.136  10Y  0.091   (all |t| < 1.3)
#   Bund:     DE10Y on OIS_2Y, beta = 0.625, t = 8.88, R2 = 0.56
#   Gradient: beta_10Y - beta_12M = 0.201, se = 0.112, t = 1.79, p = 0.074
# ==============================================================================

library(ggplot2)

source("R/functions/econometrics.R")
source("R/functions/data_helpers.R")

dat   <- readRDS("data/processed/analysis_panel.rds")
short <- sub("^dy_", "", names(MATURITIES))


# ------------------------------------------------------------------------------
# Response and placebo windows
# ------------------------------------------------------------------------------

WINDOWS <- c(pre = "pre_", response = "dy_", post = "post_")

results <- do.call(rbind, lapply(names(WINDOWS), function(w) {
  do.call(rbind, lapply(short, function(m) {
    fit <- fit_ols_robust(
      reformulate("ois_2y", response = paste0(WINDOWS[[w]], m)),
      data = dat, vcov = "HC0"
    )
    data.frame(
      window   = w,
      maturity = toupper(m),
      beta     = coef_of(fit, "ois_2y"),
      se       = unname(fit$se["ois_2y"]),
      t        = t_of(fit, "ois_2y")
    )
  }))
}))
results$window   <- factor(results$window, levels = names(WINDOWS))
results$maturity <- factor(results$maturity, levels = unname(MATURITIES))

cat("=== Response to the OIS 2Y surprise, bp per bp, HC0 ===\n\n")
print(results, digits = 3, row.names = FALSE)


# ------------------------------------------------------------------------------
# Does the surprise measure work where it must?
# ------------------------------------------------------------------------------

# A large, precise Bund coefficient inside the intraday window, alongside an
# imprecise near-zero on the predetermined Romanian fixing, is what separates
# "the shock is measured correctly and Romania is slow" from "the shock measure
# is broken".
bund <- fit_ols_robust(de10y ~ ois_2y, data = dat, vcov = "HC0")

cat(sprintf("\nDE10Y on OIS_2Y: beta = %.3f, t = %.2f, R2 = %.2f\n",
            coef_of(bund, "ois_2y"), t_of(bund, "ois_2y"), bund$r2))


# ------------------------------------------------------------------------------
# Maturity gradient
# ------------------------------------------------------------------------------

# Same regressor on both sides, so regressing the slope contrast returns
# beta_10Y - beta_12M together with the standard error of the difference.
# Comparing two separately significant coefficients would not.
slope_fit <- fit_ols_robust(dy_slope ~ ois_2y, data = dat, vcov = "HC0")
b10 <- coef_of(fit_ols_robust(dy_10y ~ ois_2y, data = dat), "ois_2y")
b12 <- coef_of(fit_ols_robust(dy_12m ~ ois_2y, data = dat), "ois_2y")

cat(sprintf("\nGradient: beta_10Y - beta_12M = %.3f (slope regression %.3f)\n",
            b10 - b12, coef_of(slope_fit, "ois_2y")))
cat(sprintf("          se = %.3f, t = %.2f, p = %.3f - suggestive, not established\n",
            unname(slope_fit$se["ois_2y"]), t_of(slope_fit, "ois_2y"),
            unname(slope_fit$p["ois_2y"])))


# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

write.csv(results, "tables/timing_validation.csv", row.names = FALSE)

p <- ggplot(results, aes(maturity, beta, colour = window, group = window)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
                  position = position_dodge(width = 0.35), size = 0.4) +
  scale_colour_manual(
    values = c(pre = "grey55", response = "#1f4e79", post = "grey75"),
    labels = c(pre = "Pre (predetermined)",
               response = "First fixing after the announcement",
               post = "Post (placebo)")
  ) +
  labs(
    title    = "The response is concentrated in the first Romanian fixing after the announcement",
    subtitle = "Response to a 1 bp ECB OIS 2Y surprise, 135 Governing Council events, 2011-2025",
    caption  = "BNR fixing struck 12:00 Bucharest; ECB announcement 13:45 Frankfurt to Jul 2022, 14:15 thereafter.\nBars are 95% intervals, HC0 robust.",
    x = "Romanian maturity", y = "Basis points per basis point", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figures/timing_validation.png", p, width = 9, height = 5.5, dpi = 150)
