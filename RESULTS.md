# Canonical results

The numbers this project reports, where each comes from, and how each is worded
here. Three distinctions stay apart throughout. The reconstruction stands
against the published shock series, HC0 against HC3, and the Romanian window
against the euro area window.

Regenerate with `R/01` through `R/06` in order. `R/02` carries the twelve timing
coefficients in its own header. Should those fail to reproduce, everything
below waits until the discrepancy has been explained.

---

## Standing conventions

1. **The published Jarocinski-Karadi series is the reported measure.** In the
   code it is `jk_mp` and `jk_cbi`, no suffix. Anything ending `_own` is the
   internal reconstruction and exists solely for the validation in `05`. The
   reconstruction reproduces the published series to numerical precision, and
   reported values are quoted from the published series.
2. **HC0 is the headline covariance estimator for the single-equation
   specifications.** The stacked cross-maturity omnibus tests (#25) cluster by
   event instead. HC3 is a small-sample check, reported alongside, and quoted
   where a result is being stressed. A coefficient and its t-statistic always
   come from the same estimator.
3. **The Romanian response is measured over the F(i) window**, usually about
   24 hours and longer across weekends and holidays, because the fixing
   precedes the announcement. **Every euro area leg is measured over the full
   EA-MPD Monetary Event Window**, spanning the press release and the press
   conference. The two responses are measured over different windows, each
   chosen to capture the first market response available in its data source.
4. **Coefficients are basis points per basis point.** The published shocks are
   converted from % p.a. by ×100, a unit change. The factor stays at 100 for
   every reported figure. Re-estimating it from the volatility of this
   project's subsample gave ×102.11 and inflated every magnitude by 2.1% until
   it was corrected.
5. **Italy and Spain are euro area periphery. Romania is a non-euro EU
   member** whose response resembles theirs.
6. **The Target/Timing/FG/QE block is Altavilla-style. It departs from
   Altavilla et al.** The short-end restriction and the four normalisations
   follow the paper. Their third restriction, QE as the minimum-variance
   direction over 2 January 2002 to 7 August 2008, is left out, because the
   sample begins in August 2011. See `construct_altavilla_factors()`.

---

## The register

Basis points per basis point throughout. HC0 t-statistics in brackets unless
stated.

### Measurement (`02_timing_validation.R`)

| # | Result | Value | |
|---|---|---|---|
| 1 | Response, first fixing after the announcement, 10Y | **0.394 (2.38)** | 12M 0.193, 3Y 0.245, 5Y 0.273 |
| 2 | Predetermined fixing, 10Y | 0.104 (0.99) | all four \|t\| < 1 |
| 3 | Following window, 10Y | 0.091 (0.62) | all four \|t\| < 1.3 |
| 4 | Bund reaction inside the intraday window | **0.625 (8.88)**, R² 0.56 | the surprise measure moves the Bund |
|   | Rows 2 and 3 report small estimates whose distance from zero the data cannot resolve. They leave a true zero undemonstrated. | | |
| 5 | Maturity gradient, raw surprise | 0.201, se 0.112, t 1.79, **p 0.074** | above the 5% threshold |

### Construction (`01_build_analysis_data.R`)

| # | Result | Value |
|---|---|---|
| 6 | Events in the panel | **135**, of which 1 on a day the Romanian market was shut |
| 7 | Days to the responding fixing | 132 at one day, 2 at four, 1 at five |
| 8 | Altavilla-style estimation samples | Target 315 events (PC1 89%); Timing/FG/QE 128 events (3 PCs 98%) |
| 9 | Rotated loadings | reference cells exactly 1.000; QE 1.054 at 5Y; Timing 0.854 at 5Y |
| 10 | Published series unit conversion | ×100, percentage points per annum to basis points. `mypc.m` already scales pc1 to the OIS_1Y standard deviation and `main.m` divides by 100, so the factor carries units alone |

### Baseline (`03_baseline_results.R`), N = 135

| # | Result | 12M | 3Y | 5Y | 10Y |
|---|---|---|---|---|---|
| 11 | Information shock | −0.413 (−2.56) | −0.566 (−2.71) | −0.625 (−2.74) | **−0.817 (−3.51)** |
| 12 | Monetary policy shock | +0.565 (1.57) | +0.705 (1.58) | +0.691 (1.55) | +0.714 (1.81) |

| # | Result | Value |
|---|---|---|
| 13 | Level | CBI **−0.605 (−3.12)**; MP +0.669 (1.64) |
| 14 | Altavilla-style family, N = 128 | no coefficient reaches \|t\| = 1.6 at any maturity; QE rises 0.148 → 0.472 along the curve |

### Curve (`04_yield_curve_results.R`)

| # | Result | Value |
|---|---|---|
| 15 | Contrasts against an unrestricted decomposition | PC1 91.7% of variance, corr with Level **−1.000**; PC2 with Slope −0.928 |
| 16 | Slope | CBI −0.404, HC0 −2.28 (p 0.022), HC3 **−1.98 (p 0.048)**. Borderline under HC3, still just inside the 5% normal threshold |
| 17 | Curvature | CBI −0.021 (−0.10), near zero |

### Audit (`05_robustness.R`)

| # | Result | Value |
|---|---|---|
| 18 | Leave-one-event-out, 10Y | **−0.876 to −0.677**, all 135 refits t < −1.96; most influential 12 March 2020 |
| 19 | Episode-based exclusions, 10Y CBI | 2011-12 −1.137; all of 2020 −0.705; both −1.085; every one t < −3.6 |
| 20 | Same exclusions, MP | falls **0.714 → 0.247**, the fragile half |
| 21 | Regimes, 10Y CBI | −0.189 / −1.308 / +5.097 / −1.734; HC3 t −2.04 / −1.96 / **0.73** / −4.30 |
| 22 | Equality across regimes | χ²(3) 41.51 HC0, **16.90 HC3, p 0.0007**, rejected |
| 23 | Endpoint contrast, 2011-14 against 2022-25 (2015-2021 excluded from this regression, so not a full-sample post-2022 interaction) | 10Y δ **−1.544** (HC3 −3.73, p 0.0002); Level δ −0.685 (HC3 −1.33, **p 0.182**) |
| 24 | γ_MP = γ_CBI | p 0.043 / 0.031 / 0.024 / **0.004**, rejected at every maturity |
| 25 | Omnibus, event-clustered | all four CBI zero, χ²(4) **13.31, p 0.0099**; all four equal, χ²(3) 5.42, **p 0.1432** |
| 26 | Independent replication of published JK construction | exact replication of the authors' PCA, event filtering and median rotation, with corr **1.000**, slope 1.000, deviation against the published values **exactly zero** on pc1 and on both poor man's series and 5.0e-07 bp on the median shocks, which is the published file's export rounding; poor man's classification agreement **100%**; equity input agrees to 4.9e-15 |
| 27 | Economic size | σ(CBI) 3.11 bp; 10Y **−2.54 bp**, Level −1.88 bp |

### Mechanism (`06_mechanism.R`)

| # | Result | MP | CBI |
|---|---|---|---|
| 28 | DE 2Y | +1.019 (6.14) | +0.508 (3.05) |
| 29 | DE 10Y | +0.511 (3.45) | **+0.042 (0.34)** |
| 30 | IT 10Y | +1.510 (8.64) | **−0.937 (−6.69)** |
| 31 | ES 10Y | +1.025 (8.25) | **−0.756 (−6.21)** |
| 32 | RO−DE differential | +0.203 (0.42) | **−0.858 (−4.07)**, HC3 −3.50 |
| 33 | EUR/USD, as an outcome | +0.071 (3.17) | +0.007 (0.46) |

34. Equality of the MP and CBI differential responses gives χ²(1) = 2.73,
    p = 0.098. A single Wald test comparing the two coefficients in #32 sits
    outside both columns above. Equality survives at the 5% level.

---
