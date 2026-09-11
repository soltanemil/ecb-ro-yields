# 06_mechanism.R
#
# Examines whether the cross-market pattern is consistent with compression of
# Romanian sovereign risk premia following positive ECB information shocks.
# Three sets of outcomes are regressed on the same two shocks, MP and CBI,
# covering the euro area core term structure (DE 2Y / 5Y / 10Y), the periphery
# (IT 10Y, ES 10Y) and the differential Romanian-German yield response. HC0 is
# the baseline covariance estimator; HC3 is recorded alongside for the CBI
# coefficient throughout and reported explicitly for the differential.
#
# Reads data/processed/analysis_panel.rds and writes tables/mechanism.csv and
# figures/mechanism.png.
#
# Every series used here was brought into the panel by 01, so this block reads
# the panel and estimates.
#
# EUR/RON, equity and CDS belong on the left hand side, as separate outcomes.
# An ECB information shock that moves the exchange rate makes the exchange rate
# a potential mediator, and conditioning on a potential mediator in the
# headline regression would change the estimand and can introduce
# post-treatment bias. EURUSD is available and is reported as an outcome for
# that reason.
#
# The outcomes exclude STOXX50 and jk_pc1. MP and CBI are constructed from
# those two series, so the pair spans the same space and a regression of either
# on MP + CBI returns an identity, with an R-squared of 1 and meaningless
# t-statistics.
#
# Window asymmetry, relevant to every magnitude below. The Romanian leg is
# measured over the F(i) window, because the fixing precedes the announcement;
# that window usually spans about 24 hours, with weekend and holiday
# extensions. Every euro area leg is measured over the full EA-MPD Monetary
# Event Window, spanning the press release and the press conference, because
# those markets trade through the announcement. Each leg is the first
# observation of its own series able to incorporate the shock, and the two
# windows differ in length. Section 3
# therefore estimates the difference between the Romanian first-fixing response
# and the German intraday response, a DIFFERENTIAL RESPONSE whose signs and
# orders of magnitude carry information. An observed change in the RO-DE
# sovereign spread would require a common window and is a separate object.
#
# On terminology, Italy and Spain are the euro area periphery. Romania is a
# non-euro EU member whose response resembles theirs.

library(dplyr)
library(ggplot2)

source("R/functions/econometrics.R")


# 1. Load analysis panel

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


# 2. The euro area, core and periphery

core <- report(c(de2y = "DE 2Y", de5y = "DE 5Y", de10y = "DE 10Y"), "core")
peri <- report(c(it10y = "IT 10Y", es10y = "ES 10Y"), "periphery")
fx   <- report(c(eurusd = "EUR/USD"), "fx")

cat("=== Euro area response to the same two shocks, Monetary Event Window ===\n\n")
print(rbind(core, peri, fx), row.names = FALSE, digits = 3)

cat("\nIn the CBI column, a positive information shock raises the German short\n")
cat("end, leaves the Bund 10Y estimate close to zero and imprecise, and is\n")
cat("associated with lower Italian and Spanish 10Y yields. Policy-shock point\n")
cat("estimates are positive throughout. This is measured on the same shocks, in\n")
cat("the same database, at the frequency where the shock is defined.\n")

cat("\nEUR/USD is reported as an outcome, for the reason given in the header.\n")
cat("It responds to the policy shock, while its response to the information\n")
cat("shock is indistinguishable from zero.\n")


# 3. Differential Romanian-German yield response

# By linearity the coefficient on the difference equals the difference of the
# two coefficients, and estimating it directly delivers a standard error that
# carries their covariance.
dat$ro_minus_de <- dat$dy_10y - dat$de10y

diff_h0 <- fit_ols_robust(ro_minus_de ~ jk_mp + jk_cbi, dat, vcov = "HC0")
diff_h3 <- fit_ols_robust(ro_minus_de ~ jk_mp + jk_cbi, dat, vcov = "HC3")

cat("\n=== Differential Romanian-German 10Y yield response ===\n\n")
cat(sprintf("MP  %+7.3f  (t = %5.2f, HC3 %5.2f)\n",
            coef_of(diff_h0, "jk_mp"), t_of(diff_h0, "jk_mp"), t_of(diff_h3, "jk_mp")))
cat(sprintf("CBI %+7.3f  (t = %5.2f, HC3 %5.2f)\n",
            coef_of(diff_h0, "jk_cbi"), t_of(diff_h0, "jk_cbi"), t_of(diff_h3, "jk_cbi")))

w <- wald_test(diff_h0, restrictions(names(diff_h0$b), c(jk_mp = 1, jk_cbi = -1)))
cat(sprintf("\nEquality of the two differentials, chi2(%d) = %.2f, p = %.4f\n",
            w[["df"]], w[["chi2"]], w[["p"]]))
cat("Three separate statements, in decreasing order of evidential support\n")
cat("  CBI against zero  the differential is negative and precisely estimated\n")
cat("  MP against zero   the differential is indistinguishable from zero\n")
cat("  CBI against MP    suggestive at 10%, above the 5% threshold\n")
cat("The third statement rests on the Wald test alone. One coefficient clearing\n")
cat("a threshold while the other stays below it speaks to each coefficient\n")
cat("against zero, and the test of the difference is the line above, which is\n")
cat("suggestive at 10%. Given the window asymmetry in the header, all of this\n")
cat("is a differential response, and an observed spread change would need a\n")
cat("common window.\n")

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


# 4. Scope of the evidence

cat("\n=== Interpretation, at the level the evidence supports ===\n")
cat("The evidence supports an association between positive ECB information\n")
cat("  shocks and lower Romanian sovereign yields, alongside lower euro area\n")
cat("  periphery yield estimates, with the Bund 10Y estimate close to zero\n")
cat("  and imprecise, a pattern consistent with a sovereign risk-premium channel.\n")
cat("Identifying that channel separately would require Romanian-specific\n")
cat("  pricing measures - EUR/RON, a local equity index, sovereign CDS - and a\n")
cat("  design that explicitly handles their post-treatment status. This block\n")
cat("  reports the association.\n")


# 5. Export

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
    caption  = "Euro area legs use the EA-MPD Monetary Event Window; the Romanian leg is the first fixing after the announcement.\nThe last column is the difference between those two responses, measured over different windows. Bars are 95% intervals, HC0 robust.",
    x = NULL, y = "Basis points per basis point"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/mechanism.png", p, width = 9, height = 5, dpi = 150)
