# Canonical results

The numbers this project reports, where each comes from, and how each is worded
here. Its purpose is to keep three distinctions apart: the reconstruction
against the published shock series, HC0 against HC3, and the Romanian window
against the euro area window.

Regenerate with `R/01` through `R/06` in order. `R/02` carries the twelve timing
coefficients in its own header; if those do not reproduce, nothing below should
be quoted until the discrepancy is explained.

---

## Standing conventions

1. **The published Jarocinski-Karadi series is the reported measure.** In the
   code it is `jk_mp` and `jk_cbi`, no suffix. Anything ending `_own` is the
   internal reconstruction and exists solely for the validation in `05`. Values
   from the reconstruction differ in the third decimal and are not quoted
   outside that comparison.
2. **HC0 is the headline covariance estimator for the single-equation
   specifications.** The stacked cross-maturity omnibus tests (#25) cluster by
   event instead. HC3 is a small-sample check, reported alongside, and quoted
   where a result is being stressed. A coefficient and its t-statistic always
   come from the same estimator.
3. **The Romanian response is measured over the F(i) window**, roughly 24 hours,
   because the fixing precedes the announcement. **Every euro area leg is an
   EA-MPD intraday window**, roughly 30 minutes. Both are the first observation
   of their own series able to incorporate the shock, but they are not a common
   window.
4. **Coefficients are basis points per basis point.** The published shocks are
   converted from % p.a. by ×100, a unit change. They are never renormalised to
   the volatility of this project's subsample; doing so gave ×102.11 and
   inflated every magnitude by 2.1% until it was corrected.
5. **Romania is not euro area periphery.** Italy and Spain are. Romania is a
   non-euro EU member whose response resembles theirs.
6. **The Target/Timing/FG/QE block is Altavilla-style, not Altavilla.** The
   short-end restriction and the four normalisations follow the paper; the
   third restriction — QE as the minimum-variance direction over 2 January 2002
   to 7 August 2008 — is not imposed, because the sample begins in August 2011.
   See `construct_altavilla_factors()`.

---

## The register

Basis points per basis point throughout. HC0 t-statistics in brackets unless
stated.

### Measurement — `02_timing_validation.R`

| # | Result | Value | |
|---|---|---|---|
| 1 | Response, first fixing after the announcement, 10Y | **0.394 (2.38)** | 12M 0.193, 3Y 0.245, 5Y 0.273 |
| 2 | Predetermined fixing, 10Y | 0.104 (0.99) | all four \|t\| < 1 |
| 3 | Following window, 10Y | 0.091 (0.62) | all four \|t\| < 1.3 |
| 4 | Bund reaction inside the intraday window | **0.625 (8.88)**, R² 0.56 | the surprise measure moves the Bund |
|   | Rows 2 and 3 are small estimates that are not statistically distinguishable from zero, not demonstrated zeros. | | |
| 5 | Maturity gradient, raw surprise | 0.201, se 0.112, t 1.79, **p 0.074** | not significant |

### Construction — `01_build_analysis_data.R`

| # | Result | Value |
|---|---|---|
| 6 | Events in the panel | **135**, of which 1 on a day the Romanian market was shut |
| 7 | Days to the responding fixing | 132 at one day, 2 at four, 1 at five |
| 8 | Altavilla-style estimation samples | Target 315 events (PC1 89%); Timing/FG/QE 128 events (3 PCs 98%) |
| 9 | Rotated loadings | reference cells exactly 1.000; QE 1.054 at 5Y; Timing 0.854 at 5Y |
| 10 | Published series unit conversion | ×100, percentage points per annum to basis points. Not a renormalisation: `mypc.m` already scales pc1 to the OIS_1Y standard deviation and `main.m` divides by 100 |

### Baseline — `03_baseline_results.R`, N = 135

| # | Result | 12M | 3Y | 5Y | 10Y |
|---|---|---|---|---|---|
| 11 | Information shock | −0.413 (−2.56) | −0.566 (−2.71) | −0.625 (−2.74) | **−0.817 (−3.51)** |
| 12 | Monetary policy shock | +0.565 (1.57) | +0.705 (1.58) | +0.691 (1.55) | +0.714 (1.81) |

| # | Result | Value |
|---|---|---|
| 13 | Level | CBI **−0.605 (−3.12)**; MP +0.669 (1.64) |
| 14 | Altavilla-style family, N = 128 | no coefficient reaches \|t\| = 1.6 at any maturity; QE rises 0.148 → 0.472 along the curve |

### Curve — `04_yield_curve_results.R`

| # | Result | Value |
|---|---|---|
| 15 | Contrasts against an unrestricted decomposition | PC1 91.7% of variance, corr with Level **−1.000**; PC2 with Slope −0.928 |
| 16 | Slope | CBI −0.404, HC0 −2.28, **HC3 −1.98** — loses the threshold under HC3 |
| 17 | Curvature | CBI −0.021 (−0.10) — nil |

### Audit — `05_robustness.R`

| # | Result | Value |
|---|---|---|
| 18 | Leave-one-event-out, 10Y | **−0.876 to −0.677**, all 135 refits t < −1.96; most influential 12 March 2020 |
| 19 | Ex ante exclusions, 10Y CBI | 2011-12 −1.137; all of 2020 −0.705; both −1.085; every one t < −3.6 |
| 20 | Same exclusions, MP | falls **0.714 → 0.247** — the fragile half |
| 21 | Regimes, 10Y CBI | −0.189 / −1.308 / +5.097 / −1.734; HC3 t −2.04 / −1.96 / **0.73** / −4.30 |
| 22 | Equality across regimes | χ²(3) 41.51 HC0, **16.90 HC3, p 0.0007** — rejected |
| 23 | Endpoint contrast, 2011-14 against 2022-25 (2015-2021 excluded from this regression, so not a full-sample post-2022 interaction) | 10Y δ **−1.544** (HC3 −3.73, p 0.0002); Level δ −0.685 (HC3 −1.33, **p 0.182**) |
| 24 | γ_MP = γ_CBI | p 0.043 / 0.031 / 0.024 / **0.004** — rejected at every maturity |
| 25 | Omnibus, event-clustered | all four CBI zero: χ²(4) **13.31, p 0.0099**; all four equal: χ²(3) 5.42, **p 0.1432** |
| 26 | Reconstruction against published | corr 0.9973 to 0.9988; poor man's agreement 87.4%; equity input agrees to 4.9e-15 over 135 events |
| 27 | Economic size | σ(CBI) 3.11 bp; 10Y **−2.54 bp**, Level −1.88 bp |

### Mechanism — `06_mechanism.R`

| # | Result | MP | CBI |
|---|---|---|---|
| 28 | DE 2Y | +1.019 (6.14) | +0.508 (3.05) |
| 29 | DE 10Y | +0.511 (3.45) | **+0.042 (0.34)** |
| 30 | IT 10Y | +1.510 (8.64) | **−0.937 (−6.69)** |
| 31 | ES 10Y | +1.025 (8.25) | **−0.756 (−6.21)** |
| 32 | RO−DE differential | +0.203 (0.42) | **−0.858 (−4.07)**, HC3 −3.50 |
| 33 | EUR/USD, as an outcome | +0.071 (3.17) | +0.007 (0.46) |

34. Equality of the MP and CBI differential responses: χ²(1) = 2.73, p = 0.098.
    A single Wald test comparing the two coefficients in #32, so it belongs to
    neither column above. It does not reject equality at 5%.

---

## Interpretation notes

How each result is worded here, and the reading it does not support.

| Claim | Reported as | Not as |
|---|---|---|
| The response window | "the first Romanian fixing able to incorporate the announcement" | "the day after", "t+1" |
| #11 | "positive information shocks are associated with lower Romanian yields" | "information shocks cause Romanian yields to fall" |
| #11 vs #24 | "the two components generate statistically different responses" | "monetary tightening raises Romanian yields" as a headline — #12 never reaches significance and #20 shows it is fragile |
| #25 | "information shocks move the curve broadly" | "the effect is concentrated at the long end" — equality across maturities is not rejected |
| #2, #3 | "small and statistically indistinguishable from zero" | "does not react", "no response" — non-rejection is not a demonstrated zero |
| #29 | "the Bund 10Y estimate is close to zero and imprecise" | "the Bund 10Y is untouched" |
| #32 vs #34 | "the CBI differential is negative and precise; the MP differential is not distinguishable from zero; the difference between them is suggestive at 10%" | that one being significant and the other not shows the two differ |
| #14 | "no individually significant coefficient; the estimates are comparatively imprecise" | "splitting the sample four ways exhausts the available power" — no power calculation supports that |
| #21, #23 | "statistically distinguishable from zero" | "identified" — a parameter can be identified and imprecisely estimated |
| #5, #16 | omit, or report as suggestive | any steepening claim; both sit above 0.05 |
| #18, #19 | "the pooled negative coefficient is robust to leave-one-out and the ex ante exclusions" | "the sign is robust across all specifications" — #21 includes a positive regime estimate |
| #21, #22 | "negative in the regimes estimated with useful precision, not a stable magnitude" | quoting the pooled coefficient alone as if it were a parameter |
| #21, 2020-21 | "highly leverage-sensitive and statistically uninformative under HC3" | "the sign reverses during the pandemic" — that is a claim about precision, not about the world |
| #23 | "distinguishable from zero on the 10Y but not on the Level, so recorded as heterogeneity" | "transmission has strengthened since 2022"; "post-2022 interaction" — 2015-2021 is not in that regression |
| #26 | "correlations above 0.997" | "above 0.998" — one pair is 0.9973 |
| #32 | "the differential Romanian-German yield response" | "the change in the RO-DE sovereign spread" — different windows |
| #30, #31, #32 | "consistent with a sovereign risk-premium channel" | "identifies the risk-premium channel"; "Romania is euro area periphery" |
| #33 | "EUR/USD is reported as an outcome, not a control" | anything implying the FX channel has been controlled for |
| #8, #9, #14 | "Altavilla-style rotation" | "the Altavilla factors" — the QE restriction is not the paper's |
| Any | values from the published series | values from the `_own` reconstruction |
| The naive finding | "the aggregate, undecomposed surprise measure produces a weak response" | "ECB policy barely transmits to Romania" — the project shows that reading is an artefact of aggregation |
| #12, #20 | "positive in the pooled sample, roughly +0.70, but not robust" | "insignificant" alone — the coefficient also collapses to 0.24 once the crisis periods are removed, which is the more informative failure |
| #11, #30–32 | "yields carry, above the risk-free component, premia associated with sovereign risk, liquidity, inflation and currency exposure; the cross-market pattern is consistent with compression of the risk component" | narrowing it to "credit and liquidity risk" — a RON 10Y yield contains BNR policy expectations, inflation and a term premium too |
| #30–32 | "the Romanian response RESEMBLES that of the euro area periphery" | "Romania behaves like the periphery" — Italy and Spain are comparators under the same identification, but Romania differs in monetary, currency and institutional structure |
