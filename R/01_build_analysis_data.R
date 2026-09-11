# 01_build_analysis_data.R
#
# Builds the analysis dataset. Every data operation lives here, covering
# loading, cleaning, shock construction, event matching and curve factors.
# Downstream scripts read the finished panel and analyse it.
#
# Reads the BNR government securities fixing from
# data/raw/titluri_de_stat_ro.xlsx, the EA-MPD intraday windows from
# data/raw/ECB_surprise_shocks.xlsx and the published Jarocinski-Karadi shocks
# from data/raw/jk_shocks_official.csv.
#
# Writes data/processed/analysis_panel.rds, data/processed/shock_validation.rds
# and data/processed/analysis_panel.csv, together with
# tables/factor_loadings.csv and the figures bnr_yields_sanity.png,
# eampd_ois3m_hist.png and factor_loadings.png.
#
# Romanian yields are matched to the first BNR fixing able to incorporate each
# announcement; see first_available_fixing(). The published Jarocinski-Karadi
# shocks are the primary measure and appear as jk_mp and jk_cbi, while the
# internal reconstruction carries the _own suffix and feeds the validation
# exercise in 05. The Altavilla-style term-structure dimensions and the JK shocks are
# separate exercises, sitting side by side in the panel and reported apart. The
# rotation behind the first departs from Altavilla et al.; see
# construct_altavilla_factors().
#
# The assertions in section 5 are the contract 02-06 rely on. Once this script
# runs to completion, the panel has the shape they expect and each of them can
# take it as given.

library(dplyr)
library(tidyr)
library(ggplot2)

source("R/functions/data_helpers.R")
source("R/functions/shock_helpers.R")


# 1. Romanian yields

bnr <- load_bnr_raw() |> clean_bnr()

# first_available_fixing() indexes into this sorted date vector, so one row per
# date is a precondition of the matching. The assertion sits on the input here
# because the panel contract examines analysis_panel$date, which holds ECB
# event dates.
stopifnot(!anyDuplicated(bnr$date))

cat("=== BNR fixing ===\n")
cat("The fixing runs from", format(min(bnr$date)), "to", format(max(bnr$date)), "\n")
cat("The sample holds", nrow(bnr), "daily observations\n")
cat("The 10Y bid-ask spread averages", round(mean(bnr$spread_10y_bps, na.rm = TRUE), 1),
    "bps\n")


# 2. ECB surprises

eampd  <- load_eampd_raw()
events <- event_surprises(eampd$mew, start_date = min(bnr$date))

cat("\n=== EA-MPD ===\n")
cat("The full history covers", nrow(eampd$mew), "events,",
    format(min(eampd$mew$date)), "to", format(max(eampd$mew$date)), "\n")
cat("The overlap sample holds", nrow(events), "events from", format(min(bnr$date)), "\n")
cat("Missing values in the retained columns\n")
print(colSums(is.na(events)))
cat("OIS_5Y and OIS_10Y are absent for early events, since the long end of\n",
    "the intraday surprise curve begins in August 2011.\n", sep = "")


# 3. Shocks

altavilla <- construct_altavilla_factors(eampd$prw, eampd$pcw)
jk_own    <- reconstruct_jk(eampd$mew)
jk        <- load_jk_official(eampd$mew)

cat("\n=== Altavilla-style dimensions ===\n")
cat("Target uses", altavilla$n[["target"]], "events, with PC1 explaining",
    round(100 * altavilla$var_share$target[1]), "% of short-end variance\n")
cat("The conference block uses", altavilla$n[["conference"]], "events, with 3 PCs explaining",
    round(100 * sum(altavilla$var_share$conference)), "% of variance\n")
cat("\nRotated loadings, reference cell 1 by construction",
    "(timing 6M, fg 2Y, qe 10Y)\n")
print(round(altavilla$loadings, 3))

cat("\n=== Jarocinski-Karadi ===\n")
cat("The reconstruction uses", jk_own$diagnostics$n, "events, with PC1 explaining",
    round(100 * jk_own$diagnostics$pc1_var_share), "% of short-end variance\n")
cat("The admissible rotation arc spans",
    paste(round(jk_own$diagnostics$arc_degrees, 1), collapse = " to "),
    "degrees, with median", round(jk_own$diagnostics$median_degrees, 1), "\n")
cat("The decomposition adds up to a maximum absolute deviation of",
    signif(jk_own$diagnostics$adds_up, 3), "\n")
cat("The published file holds", jk$n_official, "events from 2011, of which",
    jk$n_overlap, "are shared with this EA-MPD vintage\n")
cat("  Conversion to basis points multiplies by", jk$scale,
    "and retains the authors' own normalisation\n")
cat("  The published equity surprise matches the workbook to within",
    signif(jk$equity_max_abs_dev, 3), "\n")


# 4. Event matching and assembly

analysis_panel <- build_event_panel(bnr, events) |>
  add_curve_factors() |>
  left_join(altavilla$factors, by = "date") |>
  left_join(jk$shocks,         by = "date") |>
  left_join(jk_own$shocks,     by = "date") |>
  arrange(date)

cat("\n=== Analysis panel ===\n")
cat("The panel holds", nrow(analysis_panel), "ECB events\n")
cat("Of these,", sum(!analysis_panel$on_fixing_day), "fall on a closed Romanian day\n")
cat("Calendar days to the responding fixing\n")
print(table(analysis_panel$days_to_fixing))
cat("Non-missing observations by shock column\n")
print(colSums(!is.na(analysis_panel[, c(
  "target", "timing", "forward_guidance", "qe",
  "jk_mp", "jk_cbi", "jk_mp_own", "jk_cbi_own"
)])))
cat("The events short of Timing/FG/QE are the January-July 2011 meetings,\n",
    "which precede long-end OIS coverage, and they stay as NA.\n", sep = "")

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


# 5. Data contract
#
# Nothing reaches data/processed/ until every assertion below has passed, so a
# failing contract leaves the directory as it found it and a later session
# always picks up a panel that satisfied the contract.

stopifnot(nrow(analysis_panel) == 135)
stopifnot(!anyDuplicated(analysis_panel$date))
stopifnot(all(analysis_panel$fixing_date > analysis_panel$date))
stopifnot(all(c(
  names(MATURITIES), rownames(CURVE_WEIGHTS),
  "date", "fixing_date", "days_to_fixing",
  "jk_pc1", "jk_mp", "jk_cbi", "jk_mp_pm", "jk_cbi_pm",
  "jk_pc1_own", "jk_mp_own", "jk_cbi_own", "jk_mp_pm_own", "jk_cbi_pm_own",
  "target", "timing", "forward_guidance", "qe",
  "ois_2y", "de10y"
) %in% names(analysis_panel)))

# Both decompositions are exact by construction, so a failure here points to a
# broken rotation or a misaligned join.
stopifnot(max(abs(analysis_panel$jk_mp     + analysis_panel$jk_cbi     - analysis_panel$jk_pc1))     < 1e-8)
stopifnot(max(abs(analysis_panel$jk_mp_own + analysis_panel$jk_cbi_own - analysis_panel$jk_pc1_own)) < 1e-8)
stopifnot(max(abs(analysis_panel$jk_mp_pm  + analysis_panel$jk_cbi_pm  - analysis_panel$jk_pc1))     < 1e-8)
stopifnot(max(abs(analysis_panel$jk_mp_pm_own + analysis_panel$jk_cbi_pm_own - analysis_panel$jk_pc1_own)) < 1e-8)

# Replication regression test. reconstruct_jk() implements the authors' own
# pipeline, so each _own series has to reproduce its published counterpart.
# The tolerance follows the published file's rounding: main.m exports pc1 at
# eight decimals in percentage points, which bounds the deviation at 5e-07 bp,
# and the observed maxima are exactly zero on pc1 and on both poor man's
# series, and 4.951e-07 bp on the two median shocks. Any change to the units, the event filtering, the PCA or the
# rotation moves these by orders of magnitude and halts the pipeline here.
REPLICATION_PAIRS <- list(
  c("jk_pc1_own",    "jk_pc1"),
  c("jk_mp_own",     "jk_mp"),
  c("jk_cbi_own",    "jk_cbi"),
  c("jk_mp_pm_own",  "jk_mp_pm"),
  c("jk_cbi_pm_own", "jk_cbi_pm")
)
stopifnot(all(vapply(
  REPLICATION_PAIRS,
  function(pr) max(abs(analysis_panel[[pr[1]]] - analysis_panel[[pr[2]]])) < 1e-6,
  logical(1)
)))

# The curve factors have to reproduce the contrasts CURVE_WEIGHTS declares.
stopifnot(max(abs(
  as.matrix(analysis_panel[, colnames(CURVE_WEIGHTS)]) %*% t(CURVE_WEIGHTS) -
  as.matrix(analysis_panel[, rownames(CURVE_WEIGHTS)])
)) < 1e-9)

# lm() drops incomplete cases silently, so a single NA in a headline column
# would quietly shrink the estimation sample. These two assertions hold the
# headline sample at 135 events and the Altavilla-style sample at 128.
stopifnot(all(stats::complete.cases(
  analysis_panel[, c(names(MATURITIES), "jk_mp", "jk_cbi")])))
stopifnot(sum(stats::complete.cases(
  analysis_panel[, c("target", "timing", "forward_guidance", "qe")])) == 128)

saveRDS(analysis_panel,   "data/processed/analysis_panel.rds")
saveRDS(shock_validation, "data/processed/shock_validation.rds")
write.csv(analysis_panel, "data/processed/analysis_panel.csv", row.names = FALSE)

cat("\nData contract satisfied.\n")


# 6. Sanity figures

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

# The rotated loadings are exported as a table as well, to document the
# resulting curve shape.
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
    caption  = "Zero loading on the 1M OIS is imposed for forward guidance and QE, and defining QE as the direction\nthat maximises the 10Y loading places the 10Y loading of forward guidance at zero. Altavilla et al. split\nthis plane by a pre-2008 minimum-variance rule. The shape in between is recovered from the data.",
    x = "OIS maturity", y = "Loading", colour = NULL
  ) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/factor_loadings.png", p_loadings, width = 9, height = 5.5, dpi = 150)
