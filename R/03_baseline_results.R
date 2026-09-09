# ==============================================================================
# 03_baseline_results.R
#
# Baseline responses of Romanian sovereign yields to ECB shocks, under two
# identifications reported in two separate tables.
#
#   Altavilla-style   Target / Timing / Forward Guidance / QE, 128 events
#   JK                MP / CBI, published series, 135 events
#
# The first family uses the short-end zero restrictions and the normalisations
# of Altavilla et al. but not their third restriction, so it is labelled
# Altavilla-style rather than Altavilla; construct_altavilla_factors() states
# the difference.
#
# In:   data/processed/analysis_panel.rds
# Out:  tables/baseline_altavilla.csv, tables/baseline_jk.csv,
#       figures/baseline_results.png
#
# HC0 is the headline, for continuity with 02; HC1 and HC3 sit alongside as
# small-sample checks rather than replacements.
#
# Each family is estimated jointly rather than one regressor at a time: the
# responding Romanian fixing comes after both ECB windows, so a single change
# can embed the press-release Target and the press-conference dimensions at
# once. The three conference factors are orthogonal by construction; Target
# comes from a different window and is not assumed orthogonal to them.
#
# The two families run on different samples on purpose. The seven January-July
# 2011 events lack a long-end surprise, so they are dropped from the first
# family rather than imputed - a missing surprise is not a zero shock - while
# the JK decomposition exists for all 135. Section 4 asks what that sample
# restriction does to the univariate Target relation, which is the only part of
# the family the seven events can speak to.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/functions/econometrics.R")
source("R/functions/data_helpers.R")

dat <- readRDS("data/processed/analysis_panel.rds")

OUTCOMES <- c(MATURITIES,
              dy_level = "Level", dy_slope = "Slope", dy_curvature = "Curvature")

estimate_family <- function(rhs, label) {
  do.call(rbind, lapply(names(OUTCOMES), function(o) {
    f   <- reformulate(rhs, response = o)
    h0  <- fit_ols_robust(f, dat, vcov = "HC0")
    h1  <- fit_ols_robust(f, dat, vcov = "HC1")
    h3  <- fit_ols_robust(f, dat, vcov = "HC3")
    keep <- rhs
    data.frame(
      family = label, outcome = OUTCOMES[[o]], term = keep,
      estimate = unname(h0$b[keep]),
      se_hc0 = unname(h0$se[keep]), t_hc0 = unname(h0$t[keep]),
      p_hc0  = unname(h0$p[keep]),
      t_hc1  = unname(h1$t[keep]), t_hc3 = unname(h3$t[keep]),
      n = h0$n, r2 = h0$r2, row.names = NULL
    )
  }))
}

show_family <- function(res, title, units) {
  cat("\n===", title, "===\n")
  cat(units, "\n\n")
  wide <- res |>
    mutate(cell = sprintf("%6.3f (%5.2f)", estimate, t_hc0)) |>
    select(outcome, term, cell) |>
    pivot_wider(names_from = term, values_from = cell)
  print(as.data.frame(wide), row.names = FALSE)
  cat("\nN =", unique(res$n), "\n")
}


# ------------------------------------------------------------------------------
# Which dimension of ECB communication moved the curve?
# ------------------------------------------------------------------------------

ALTAVILLA <- c("target", "timing", "forward_guidance", "qe")

res_altavilla <- estimate_family(ALTAVILLA, "Altavilla-style")

show_family(
  res_altavilla,
  "A. Dimensions of ECB communication, Altavilla-style rotation",
  paste("beta = bp of Romanian yield per 1 bp of the factor's own reference",
        "maturity\n(target 1M, timing 6M, forward guidance 2Y, qe 10Y)")
)

cat("\nFactor correlations on the estimation sample:\n")
print(round(cor(dat[complete.cases(dat[, ALTAVILLA]), ALTAVILLA]), 2))
cat("The three conference factors are orthogonal by construction; Target\n",
    "against the other three is the entry that matters.\n", sep = "")


# ------------------------------------------------------------------------------
# Policy news, or information about the economy?
# ------------------------------------------------------------------------------

res_jk    <- estimate_family(c("jk_mp", "jk_cbi"),       "JK median rotation")
res_jk_pm <- estimate_family(c("jk_mp_pm", "jk_cbi_pm"), "JK poor man's")

jk_units <- paste("beta = bp of Romanian yield per 1 bp of the policy-surprise",
                  "component\n(the published series, converted from % p.a. to",
                  "basis points)")

show_family(res_jk,    "B. Policy news versus central bank information", jk_units)
show_family(res_jk_pm, "B'. Same, poor man's sign restriction (robustness)", jk_units)


# ------------------------------------------------------------------------------
# What the common-sample restriction does to Target
# ------------------------------------------------------------------------------

# Target is the one factor defined on all 135 events, so it is the only one
# whose sample sensitivity can be examined at all. The comparison below is
# univariate on both sides: it asks whether the Target-yield relation shifts
# when the seven events without a long-end surprise are dropped. It does not
# and cannot speak to the Target coefficient of the joint model, which is not
# estimable on 135 events because the conference factors are missing there.
cat("\n=== Target alone, univariate, full versus common sample ===\n")
common <- complete.cases(dat[, ALTAVILLA])
for (o in c("dy_12m", "dy_10y", "dy_slope")) {
  f    <- reformulate("target", response = o)
  full <- fit_ols_robust(f, dat)
  sub  <- fit_ols_robust(f, dat[common, ])
  cat(sprintf("%-10s N=%3d beta=%6.3f (t=%5.2f)   |   N=%3d beta=%6.3f (t=%5.2f)\n",
              OUTCOMES[[o]], full$n, coef_of(full, "target"), t_of(full, "target"),
              sub$n,  coef_of(sub,  "target"), t_of(sub,  "target")))
}


# ------------------------------------------------------------------------------
# Export
# ------------------------------------------------------------------------------

write.csv(res_altavilla, "tables/baseline_altavilla.csv", row.names = FALSE)
write.csv(rbind(res_jk, res_jk_pm), "tables/baseline_jk.csv", row.names = FALSE)

LABELS <- c(
  target = "Target", timing = "Timing",
  forward_guidance = "Forward guidance", qe = "QE",
  jk_mp = "Monetary policy shock", jk_cbi = "Information shock"
)

plot_df <- rbind(res_altavilla, res_jk) |>
  filter(outcome %in% unname(MATURITIES)) |>
  mutate(
    outcome = factor(outcome, levels = unname(MATURITIES)),
    term    = factor(LABELS[term], levels = unname(LABELS))
  )

p <- ggplot(plot_df, aes(outcome, estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = estimate - 1.96 * se_hc0,
                      ymax = estimate + 1.96 * se_hc0),
                  colour = "#1f4e79", size = 0.35) +
  facet_wrap(~term, nrow = 2, scales = "free_y") +
  labs(
    title    = "Romanian sovereign yield response, by shock dimension",
    subtitle = "First fixing after the announcement; Altavilla-style factors N = 128, Jarocinski-Karadi N = 135",
    caption  = "Bars are 95% intervals, HC0 robust. Units differ by family - see the table headers.",
    x = "Romanian maturity", y = "Basis points"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

ggsave("figures/baseline_results.png", p, width = 10, height = 6.5, dpi = 150)
