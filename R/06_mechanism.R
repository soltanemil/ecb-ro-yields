# ==============================================================================
# 06_mechanism.R
#
# Asks why a positive ECB information shock lowers Romanian sovereign yields.
# The hypothesis under test is narrow - compression of Romanian sovereign risk
# premia - rather than a generic appeal to risk appetite. The same two shocks,
# MP and CBI, are regressed on three sets of outcomes already carried in the
# panel: the euro area core term structure (DE 2Y / 5Y / 10Y), the periphery
# (IT 10Y, ES 10Y), and the differential Romanian-German yield response. HC0
# baseline, with HC3 alongside on the differential, which is the only
# specification here carrying an interpretive claim.
#
# In:   data/processed/analysis_panel.rds
# Out:  tables/mechanism.csv, figures/mechanism.png
#
# No new data: every series used here was brought into the panel by 01.
#
# EUR/RON, equity and CDS are not added as controls. If an ECB information
# shock moves the exchange rate, the exchange rate is a potential mediator, and
# conditioning on it would put a post-treatment variable on the right hand side
# and make the CBI coefficient uninterpretable. Such variables belong on the
# left hand side, as separate outcomes; EURUSD is available and is reported
# that way.
#
# STOXX50 and jk_pc1 are excluded as outcomes. MP and CBI are constructed from
# those two series, so the pair spans the same space and a regression of either
# on MP + CBI is an identity: R-squared of 1 and meaningless t-statistics.
#
# Window asymmetry, relevant to every magnitude below. The Romanian leg is
# measured over the F(i) window, roughly 24 hours, because the fixing precedes
# the announcement. Every euro area leg is measured over the EA-MPD intraday
# window, roughly 30 minutes, because those markets trade through the
# announcement. Both are the first observation of their own series able to
# incorporate the shock, but they are not a common window. Section 3 therefore
# estimates a DIFFERENTIAL RESPONSE, not the change in an observable spread:
# signs and orders of magnitude carry, a precise decomposition of a spread does
# not. Call it the differential Romanian-German yield response, or the
# difference between the Romanian first-fixing response and the German intraday
# response - not "the change in the RO-DE sovereign spread".
#
# Romania is also not euro area periphery. Italy and Spain are the periphery;
# Romania is a non-euro EU member whose response resembles theirs.
# ==============================================================================

library(dplyr)
library(ggplot2)

source("R/functions/econometrics.R")


# ------------------------------------------------------------------------------
# 1. Load analysis panel
# ------------------------------------------------------------------------------

dat <- readRDS("data/processed/analysis_panel.rds")
JK  <- c("jk_mp", "jk_cbi")

report <- function(outcomes, label) {
  do.call(rbind, lapply(names(outcomes), function(v) {
    h0 <- fit_ols_robust(reformulate(JK, response = v), dat, vcov = "HC0")
    h3 <- fit_ols_robust(reformulate(JK, response = v), dat, vcov = "HC3")
    data.frame(
      block = label, series = outcomes[[v]],
      mp = coef_of(h0, "jk_mp"),  t_mp  = t_of(h0, "jk_mp"),
      cbi = coef_of(h0, "jk_cbi"), t_cbi = t_of(h0, "jk_cbi"),
      se_cbi = unname(h0$se["jk_cbi"]),
      t_cbi_hc3 = t_of(h3, "jk_cbi"), r2 = h0$r2, n = h0$n, row.names = NULL
    )
  }))
}


# ------------------------------------------------------------------------------
# 2. The euro area, core and periphery
# ------------------------------------------------------------------------------

core <- report(c(de2y = "DE 2Y", de5y = "DE 5Y", de10y = "DE 10Y"), "core")
peri <- report(c(it10y = "IT 10Y", es10y = "ES 10Y"), "periphery")
fx   <- report(c(eurusd = "EUR/USD"), "fx")

cat("=== Euro area response to the same two shocks, intraday window ===\n\n")
print(rbind(core, peri, fx), row.names = FALSE, digits = 3)

cat("\nIn the CBI column, a positive information shock raises the German short\n")
cat("end, leaves the Bund 10Y estimate close to zero and imprecise, and\n")
cat("compresses Italian and Spanish 10Y yields. Policy-shock point estimates\n")
cat("are positive throughout. This is measured on the same shocks, in the same\n")
cat("database, at the frequency where the shock itself is defined.\n")

cat("\nEUR/USD is reported as an outcome, not a control, for the reason given\n")
cat("in the header. It responds to the policy shock and not to the\n")
cat("information shock.\n")


# ------------------------------------------------------------------------------
# 3. Differential Romanian-German yield response
# ------------------------------------------------------------------------------

# By linearity the coefficient on the difference is the difference of the two
# coefficients, and its standard error accounts for their covariance - which is
# why it is estimated directly rather than read off two separate tables.
dat$ro_minus_de <- dat$dy_10y - dat$de10y

diff_h0 <- fit_ols_robust(ro_minus_de ~ jk_mp + jk_cbi, dat, vcov = "HC0")
diff_h3 <- fit_ols_robust(ro_minus_de ~ jk_mp + jk_cbi, dat, vcov = "HC3")

cat("\n=== Differential Romanian-German 10Y yield response ===\n\n")
cat(sprintf("MP  %+7.3f  (t = %5.2f, HC3 %5.2f)\n",
            coef_of(diff_h0, "jk_mp"), t_of(diff_h0, "jk_mp"), t_of(diff_h3, "jk_mp")))
cat(sprintf("CBI %+7.3f  (t = %5.2f, HC3 %5.2f)\n",
            coef_of(diff_h0, "jk_cbi"), t_of(diff_h0, "jk_cbi"), t_of(diff_h3, "jk_cbi")))

w <- wald_test(diff_h0, restrictions(names(diff_h0$b), c(jk_mp = 1, jk_cbi = -1)))
cat(sprintf("\nEquality of the two differentials: chi2(%d) = %.2f, p = %.4f\n",
            w[["df"]], w[["chi2"]], w[["p"]]))
cat("Three separate statements, in decreasing order of what the data support:\n")
cat("  CBI against zero  the differential is negative and precisely estimated\n")
cat("  MP against zero   the differential is not distinguishable from zero\n")
cat("  CBI against MP    the Wald test above does not reject equality at 5%\n")
cat("The first two do not add up to the third. One coefficient clearing a\n")
cat("threshold while the other does not is not evidence that they differ; the\n")
cat("test of the difference is the line above, and it is suggestive at 10%\n")
cat("rather than conclusive. Given the window asymmetry in the header, read\n")
cat("all of this as a differential response, not as a change in an observable\n")
cat("spread.\n")

mechanism <- rbind(
  core, peri, fx,
  data.frame(block = "differential", series = "RO-DE differential",
             mp = coef_of(diff_h0, "jk_mp"),  t_mp = t_of(diff_h0, "jk_mp"),
             cbi = coef_of(diff_h0, "jk_cbi"), t_cbi = t_of(diff_h0, "jk_cbi"),
             se_cbi = unname(diff_h0$se["jk_cbi"]),
             t_cbi_hc3 = t_of(diff_h3, "jk_cbi"), r2 = diff_h0$r2, n = diff_h0$n,
             row.names = NULL)
)
write.csv(mechanism, "tables/mechanism.csv", row.names = FALSE)


# ------------------------------------------------------------------------------
# 4. Scope of the evidence
# ------------------------------------------------------------------------------

cat("\n=== Interpretation, at the level the evidence supports ===\n")
cat("Supported: positive ECB information shocks are associated with lower\n")
cat("  Romanian sovereign yields, alongside compression of euro area periphery\n")
cat("  yields, while the Bund 10Y estimate is close to zero and imprecise.\n")
cat("  That pattern is consistent with a sovereign risk-premium channel.\n")
cat("NOT supported: that the risk-premium channel has been identified. Doing\n")
cat("  so would need Romanian-specific pricing - EUR/RON, a local equity index,\n")
cat("  sovereign CDS - and a design that does not condition on post-treatment\n")
cat("  variables. Neither is claimed here.\n")


# ------------------------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------------------------

plot_df <- mechanism |>
  filter(series != "EUR/USD") |>
  mutate(series = factor(series, levels = c(
    "DE 2Y", "DE 5Y", "DE 10Y", "IT 10Y", "ES 10Y", "RO-DE differential")))

p <- ggplot(plot_df, aes(series, cbi)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = cbi - 1.96 * se_cbi, ymax = cbi + 1.96 * se_cbi),
                  colour = "#1f4e79", size = 0.4) +
  labs(
    title    = paste0("Positive ECB information shocks are associated with lower Romanian and euro area\n",
                      "periphery yields; the Bund 10Y estimate is near zero"),
    subtitle = "Response to a 1 bp information shock, 135 Governing Council events",
    caption  = "Euro area legs are EA-MPD intraday windows; the Romanian leg is the first fixing after the announcement.\nThe last column is the difference between those two responses, not a change in an observable spread. Bars are 95% intervals, HC0 robust.",
    x = NULL, y = "Basis points per basis point"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/mechanism.png", p, width = 9, height = 5, dpi = 150)
