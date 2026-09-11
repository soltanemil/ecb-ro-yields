# 04_yield_curve_results.R
#
# The curve factors, the agreement between the imposed contrasts and an
# unrestricted decomposition, and the baseline JK specification run on them.
# HC0 with HC3 alongside.
#
# Reads data/processed/analysis_panel.rds and writes tables/curve_factors.csv
# and figures/curve_factors.png.
#
# Level, Slope and Curvature are fixed contrasts of the four quoted maturities
# (CURVE_WEIGHTS), so all three summarise the same maturity family. Section 2
# checks that they track the principal components of the same yields; any
# divergence belongs in the write-up, since these are imposed contrasts rather
# than estimated factors.
#
# The omnibus tests in 05 supply the reading rule. The joint null that the
# four CBI maturity coefficients are zero is rejected. The equality test
# leaves the data compatible with a common response across maturities. Level
# therefore provides a useful summary, while a maturity gradient remains
# unestablished.

library(dplyr)
library(ggplot2)

source("R/functions/econometrics.R")
source("R/functions/data_helpers.R")

dat     <- readRDS("data/processed/analysis_panel.rds")
FACTORS <- c(dy_level = "Level", dy_slope = "Slope", dy_curvature = "Curvature")


# What the fixed contrasts capture

cat("=== Contrasts ===\n\n")
print(CURVE_WEIGHTS)
cat("\nA rise in Level lifts the whole curve, a rise in Slope lifts the long\n",
    "end relative to the short, and a rise in Curvature lifts the belly.\n",
    sep = "")

X  <- as.matrix(dat[, names(MATURITIES)])
Xc <- scale(X, center = TRUE, scale = FALSE)
ev <- eigen(cov(Xc), symmetric = TRUE)

cat("\n=== Against an unrestricted decomposition ===\n")
cat("Variance shares of the principal components of Romanian yield changes\n")
cat(paste(sprintf("  PC%d  %4.1f%%", 1:4, 100 * ev$values / sum(ev$values)),
          collapse = "\n"), "\n\n")

pcs     <- Xc %*% ev$vectors[, 1:3]
cor_tab <- round(cor(pcs, dat[, names(FACTORS)]), 3)
dimnames(cor_tab) <- list(paste0("PC", 1:3), unname(FACTORS))
print(cor_tab)

cat("\nPC loadings\n")
lt <- round(ev$vectors[, 1:3], 3)
dimnames(lt) <- list(unname(MATURITIES), paste0("PC", 1:3))
print(lt)
cat("\nThe correlation table shows how closely the three fixed contrasts track\n",
    "the leading principal components. PC signs are arbitrary, so absolute\n",
    "correlations carry the comparison.\n", sep = "")


# Baseline specification on the three factors

out <- do.call(rbind, lapply(names(FACTORS), function(o) {
  f  <- reformulate(c("jk_mp", "jk_cbi"), response = o)
  h0 <- fit_ols_robust(f, dat, vcov = "HC0")
  h3 <- fit_ols_robust(f, dat, vcov = "HC3")
  data.frame(
    factor = FACTORS[[o]], shock = c("MP", "CBI"),
    beta   = unname(h0$b[c("jk_mp", "jk_cbi")]),
    se_hc0 = unname(h0$se[c("jk_mp", "jk_cbi")]),
    t_hc0  = unname(h0$t[c("jk_mp", "jk_cbi")]),
    t_hc3  = unname(h3$t[c("jk_mp", "jk_cbi")]),
    n      = h0$n, row.names = NULL
  )
}))

cat("\n=== JK specification on the curve factors ===\n")
print(out, row.names = FALSE, digits = 3)

cat("\nEffect of a one standard deviation information shock\n")
sd_cbi <- sd(dat$jk_cbi)
for (i in which(out$shock == "CBI"))
  cat(sprintf("  %-10s %+6.2f bp\n", out$factor[i], out$beta[i] * sd_cbi))


# Export

write.csv(out, "tables/curve_factors.csv", row.names = FALSE)

p <- out |>
  mutate(
    factor = factor(factor, levels = unname(FACTORS)),
    shock  = factor(shock, levels = c("MP", "CBI"),
                    labels = c("Monetary policy", "Information"))
  ) |>
  ggplot(aes(factor, beta, colour = shock)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_pointrange(aes(ymin = beta - 1.96 * se_hc0, ymax = beta + 1.96 * se_hc0),
                  position = position_dodge(width = 0.4), size = 0.4) +
  scale_colour_manual(values = c("#c0762f", "#1f4e79")) +
  labs(
    title    = "Romanian yield-curve responses to monetary-policy and information shocks",
    subtitle = "Fixed contrasts of the four quoted maturities, 135 events",
    caption  = "Basis points per 1 bp of shock. Bars are 95% intervals, HC0 robust.\nLevel, slope and curvature are transformations of the same four yields and summarise one maturity family.",
    x = NULL, y = "Basis points", colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figures/curve_factors.png", p, width = 8, height = 5, dpi = 150)
