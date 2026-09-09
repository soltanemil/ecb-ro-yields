# ==============================================================================
# 01_build_analysis_data.R
#
# Builds the analysis dataset. All data operations live here: loading,
# cleaning, shock construction, event matching, curve factors. Nothing
# downstream cleans, matches, joins or renames.
#
# In:   data/raw/titluri_de_stat_ro.xlsx    BNR government securities fixing
#       data/raw/ECB_surprise_shocks.xlsx   EA-MPD intraday windows
#       data/raw/jk_shocks_official.csv     published Jarocinski-Karadi shocks
# Out:  data/processed/analysis_panel.rds, data/processed/shock_validation.rds,
#       data/processed/analysis_panel.csv, tables/factor_loadings.csv,
#       figures/bnr_yields_sanity.png, figures/eampd_ois3m_hist.png,
#       figures/factor_loadings.png
#
# Romanian yields are matched to the first BNR fixing able to incorporate each
# announcement, not to calendar t+1; see first_available_fixing(). The
# published Jarocinski-Karadi shocks are the primary measure and carry no
# suffix, jk_mp and jk_cbi; the internal reconstruction carries _own and exists
# so that 05 can audit it. The Altavilla-style term-structure dimensions and
# the JK shocks are separate exercises: they sit side by side in the panel and
# are never combined. The rotation behind the first is not the one Altavilla
# et al. impose - see construct_altavilla_factors().
#
# The assertions in section 5 are the contract 02-06 rely on. If this script
# finishes without error, the panel has the shape they expect and none of them
# needs to check again.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/functions/data_helpers.R")
source("R/functions/shock_helpers.R")


# ------------------------------------------------------------------------------
# 1. Romanian yields
# ------------------------------------------------------------------------------

bnr <- load_bnr_raw() |> clean_bnr()

# first_available_fixing() indexes into this sorted date vector, so one row per
# date is a precondition of the matching, not a stylistic preference. Asserted
# on the input because the panel contract cannot reach it: that contract checks
# analysis_panel$date, which holds ECB event dates.
stopifnot(!anyDuplicated(bnr$date))

cat("=== BNR fixing ===\n")
cat("Range :", format(min(bnr$date)), "to", format(max(bnr$date)), "\n")
cat("Obs   :", nrow(bnr), "\n")
cat("Mean 10Y bid-ask spread:", round(mean(bnr$spread_10y_bps, na.rm = TRUE), 1),
    "bps\n")


# ------------------------------------------------------------------------------
# 2. ECB surprises
# ------------------------------------------------------------------------------

eampd  <- load_eampd_raw()
events <- event_surprises(eampd$mew, start_date = min(bnr$date))

cat("\n=== EA-MPD ===\n")
cat("Full history  :", nrow(eampd$mew), "events,",
    format(min(eampd$mew$date)), "to", format(max(eampd$mew$date)), "\n")
cat("Overlap sample:", nrow(events), "events from", format(min(bnr$date)), "\n")
cat("Missing in the retained columns:\n")
print(colSums(is.na(events)))
cat("OIS_5Y and OIS_10Y are absent for early events; the long end of the\n",
    "intraday surprise curve does not exist before August 2011.\n", sep = "")


# ------------------------------------------------------------------------------
# 3. Shocks
# ------------------------------------------------------------------------------

altavilla <- construct_altavilla_factors(eampd$prw, eampd$pcw)
jk_own    <- reconstruct_jk(eampd$mew)
jk        <- load_jk_official(eampd$mew)

cat("\n=== Altavilla-style dimensions ===\n")
cat("Target      :", altavilla$n[["target"]], "events, PC1 explains",
    round(100 * altavilla$var_share$target[1]), "% of short-end variance\n")
cat("Conference  :", altavilla$n[["conference"]], "events, 3 PCs explain",
    round(100 * sum(altavilla$var_share$conference)), "% of variance\n")
cat("\nRotated loadings, reference cell is 1 by construction",
    "(timing 6M, fg 2Y, qe 10Y):\n")
print(round(altavilla$loadings, 3))

cat("\n=== Jarocinski-Karadi ===\n")
cat("Reconstruction:", jk_own$diagnostics$n, "events, PC1 explains",
    round(100 * jk_own$diagnostics$pc1_var_share), "% of short-end variance\n")
cat("Admissible rotation arc:",
    paste(round(jk_own$diagnostics$arc_degrees, 1), collapse = " to "),
    "degrees; median", round(jk_own$diagnostics$median_degrees, 1), "\n")
cat("Decomposition adds up, max abs deviation:",
    signif(jk_own$diagnostics$adds_up, 3), "\n")
cat("Published series:", jk$n_official, "events from 2011,",
    jk$n_overlap, "shared with this EA-MPD vintage;\n")
cat("  converted from % p.a. to basis points by x", jk$scale,
    "- a unit change, not a renormalisation\n")
cat("  equity surprise, published file against the workbook: max abs deviation",
    signif(jk$equity_max_abs_dev, 3), "\n")


# ------------------------------------------------------------------------------
# 4. Event matching and assembly
# ------------------------------------------------------------------------------

analysis_panel <- build_event_panel(bnr, events) |>
  add_curve_factors() |>
  left_join(altavilla$factors, by = "date") |>
  left_join(jk$shocks,         by = "date") |>
  left_join(jk_own$shocks,     by = "date") |>
  arrange(date)

cat("\n=== Analysis panel ===\n")
cat("Events                       :", nrow(analysis_panel), "\n")
cat("Events on a non-fixing day   :", sum(!analysis_panel$on_fixing_day), "\n")
cat("Calendar days to the responding fixing:\n")
print(table(analysis_panel$days_to_fixing))
cat("Non-missing shock columns:\n")
print(colSums(!is.na(analysis_panel[, c(
  "target", "timing", "forward_guidance", "qe",
  "jk_mp", "jk_cbi", "jk_mp_own", "jk_cbi_own"
)])))
cat("The events short of Timing/FG/QE are the January-July 2011 meetings,\n",
    "before long-end OIS coverage begins. They are not imputed.\n", sep = "")

shock_validation <- list(
  altavilla = altavilla[c("loadings", "var_share", "n")],
  jk        = list(diagnostics = jk_own$diagnostics,
                   official_scale = jk$scale,
                   official_n = jk$n_official, overlap_n = jk$n_overlap,
                   equity_max_abs_dev = jk$equity_max_abs_dev),
  unscaled  = altavilla$target_std |>
    full_join(altavilla$unscaled, by = "date") |>
    filter(date %in% analysis_panel$date) |>
    arrange(date)
)


# ------------------------------------------------------------------------------
# 5. Data contract
#
# Nothing is written to data/processed/ until every assertion below has passed.
# A failing contract must not be able to leave a well-named but invalid panel on
# disk for a later session to pick up.
# ------------------------------------------------------------------------------

stopifnot(nrow(analysis_panel) == 135)
stopifnot(!anyDuplicated(analysis_panel$date))
stopifnot(all(analysis_panel$fixing_date > analysis_panel$date))
stopifnot(all(c(
  names(MATURITIES), rownames(CURVE_WEIGHTS),
  "date", "fixing_date", "days_to_fixing",
  "jk_pc1", "jk_mp", "jk_cbi", "jk_mp_pm", "jk_cbi_pm",
  "jk_pc1_own", "jk_mp_own", "jk_cbi_own",
  "target", "timing", "forward_guidance", "qe",
  "ois_2y", "de10y"
) %in% names(analysis_panel)))

# Both decompositions are exact by construction, so this catches a broken
# rotation or a misaligned join rather than a rounding difference.
stopifnot(max(abs(analysis_panel$jk_mp     + analysis_panel$jk_cbi     - analysis_panel$jk_pc1))     < 1e-8)
stopifnot(max(abs(analysis_panel$jk_mp_own + analysis_panel$jk_cbi_own - analysis_panel$jk_pc1_own)) < 1e-8)
stopifnot(max(abs(analysis_panel$jk_mp_pm  + analysis_panel$jk_cbi_pm  - analysis_panel$jk_pc1))     < 1e-8)

# The curve factors must be exactly the contrasts CURVE_WEIGHTS declares.
stopifnot(max(abs(
  as.matrix(analysis_panel[, colnames(CURVE_WEIGHTS)]) %*% t(CURVE_WEIGHTS) -
  as.matrix(analysis_panel[, rownames(CURVE_WEIGHTS)])
)) < 1e-9)

# lm() drops incomplete cases silently, so a single NA in a headline column
# would shrink the estimation sample without anything downstream noticing.
stopifnot(all(stats::complete.cases(
  analysis_panel[, c(names(MATURITIES), "jk_mp", "jk_cbi")])))
stopifnot(sum(stats::complete.cases(
  analysis_panel[, c("target", "timing", "forward_guidance", "qe")])) == 128)

saveRDS(analysis_panel,   "data/processed/analysis_panel.rds")
saveRDS(shock_validation, "data/processed/shock_validation.rds")
write.csv(analysis_panel, "data/processed/analysis_panel.csv", row.names = FALSE)

cat("\nData contract satisfied.\n")


# ------------------------------------------------------------------------------
# 6. Sanity figures
# ------------------------------------------------------------------------------

p_yields <- bnr |>
  select(date, mid_12m, mid_3y, mid_5y, mid_10y) |>
  pivot_longer(-date, names_to = "maturity", values_to = "yield") |>
  mutate(maturity = factor(maturity,
    levels = c("mid_12m", "mid_3y", "mid_5y", "mid_10y"),
    labels = c("12M", "3Y", "5Y", "10Y"))) |>
  ggplot(aes(date, yield)) +
  geom_line(colour = "steelblue", linewidth = 0.4) +
  facet_wrap(~maturity, ncol = 2, scales = "free_y") +
  labs(title = "Romanian government bond yields - BNR mid-fixing, daily",
       caption = "Source: National Bank of Romania", x = NULL, y = "Yield (%)") +
  theme_minimal(base_size = 11)
ggsave("figures/bnr_yields_sanity.png", p_yields, width = 10, height = 7, dpi = 150)

p_shocks <- ggplot(events, aes(ois_3m)) +
  geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  labs(title = "ECB monetary policy surprises, 3-month OIS",
       subtitle = "Monetary Event Window, overlap sample",
       caption = "Source: EA-MPD, Altavilla et al. (2019)",
       x = "Surprise (bps)", y = "Events") +
  theme_minimal(base_size = 11)
ggsave("figures/eampd_ois3m_hist.png", p_shocks, width = 8, height = 5, dpi = 150)

# The rotated loadings are exported as a table as well: they are the evidence
# that the rotation recovered an economically interpretable curve shape.
loadings_tab <- data.frame(maturity = rownames(altavilla$loadings),
                           round(altavilla$loadings, 3), row.names = NULL)
write.csv(loadings_tab, "tables/factor_loadings.csv", row.names = FALSE)

p_loadings <- loadings_tab |>
  pivot_longer(-maturity, names_to = "factor", values_to = "loading") |>
  mutate(
    maturity = factor(maturity, levels = TERM_TENORS),
    factor   = factor(factor,
      levels = c("timing", "forward_guidance", "qe"),
      labels = c("Timing", "Forward guidance", "QE"))
  ) |>
  ggplot(aes(maturity, loading, group = factor, colour = factor)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = c("#1f4e79", "#c0762f", "#4f7a52")) +
  labs(
    title    = "Rotated press-conference factors and the euro area OIS curve",
    subtitle = "Altavilla-style rotation under the Altavilla normalisation, 128 events, Aug 2011 - Oct 2025",
    caption  = "Zero loading on the 1M OIS is imposed for forward guidance and QE, and defining QE as the direction\nthat maximises the 10Y loading forces the 10Y loading of forward guidance to zero. Altavilla et al. split\nthis plane by a pre-2008 minimum-variance rule instead; the shape in between is recovered, not imposed.",
    x = "OIS maturity", y = "Loading", colour = NULL
  ) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/factor_loadings.png", p_loadings, width = 9, height = 5.5, dpi = 150)
