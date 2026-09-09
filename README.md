# Monetary Policy or Information? ECB Spillovers to Romanian Sovereign Yields

## Research question

Do ECB surprises transmit to the Romanian sovereign yield curve, and does it
matter whether a surprise is monetary policy news or information the ECB
reveals about the economy?

Romania sits outside the euro area and runs its own monetary policy, but the
two markets are linked through euro funding, non-resident participation in the
government bond market, and regional risk appetite. The answer to the second
question matters: the two components move Romanian yields in opposite
directions, which is why an undecomposed surprise measure appears to do so
little.

## Identification and timing

The identifying assumption rests on a timing asymmetry, which determines the
entire specification.

| | Time (Frankfurt) |
|---|---|
| BNR government securities fixing | 11:00 (12:00 Bucharest) |
| ECB announcement, to 20 July 2022 | 13:45 |
| ECB announcement, from 21 July 2022 | 14:15 |

In both regimes the Romanian fixing on the announcement day is struck
**before** the announcement, and is therefore predetermined with respect to the
surprise. The response is measured on the first fixing able to incorporate it:

```
dy_i = y_{F(i)} - y_{F(i)-1}
```

where `F(i)` indexes the first BNR fixing struck strictly after ECB event `i`.
This is not `lead(date, 1)`: it is defined when an event falls on a Romanian
public holiday, it spans weekends and holiday bridges correctly, and it
discards no events.

`R/02_timing_validation.R` treats this as a falsification design rather than an
assertion. Four conditions must hold jointly, and do: the Bund reacts strongly
inside the EA-MPD intraday window; the estimate on the predetermined Romanian
fixing is small and statistically indistinguishable from zero; the first
subsequent fixing does respond; and the estimate for the window after that is
again small and indistinguishable from zero.

## Data

| Dataset | Source | Coverage |
|---|---|---|
| Romanian government bond yields (6M, 12M, 3Y, 5Y, 10Y) | National Bank of Romania, daily bid-ask fixing | 2011-01-10 to 2026-04-27, 3,844 observations |
| ECB intraday surprises | EA-MPD, Altavilla et al. (2019) | 315 events from 1999; 135 in the Romanian overlap |
| Decomposed policy and information shocks | Jarocinski and Karadi, published series | the same 135 events |

None of the three inputs is redistributed here. To reproduce the project, place
in `data/raw/`:

- `titluri_de_stat_ro.xlsx` — the BNR government securities fixing
- `ECB_surprise_shocks.xlsx` — the EA-MPD workbook
- `jk_shocks_official.csv` — `shocks_ecb_mpd_me_d.csv` from
  [github.com/marekjarocinski/jkshocks_update_ecb](https://github.com/marekjarocinski/jkshocks_update_ecb),
  renamed. Only the rows from 2011 onwards are used, and `load_jk_official()`
  enforces that with an explicit filter so the reported event counts do not
  depend on how the file was trimmed before saving.

## Shock measures

The primary measure is the published Jarocinski-Karadi decomposition, converted
from percentage points to basis points by multiplying by 100. It carries no
suffix in the code: `jk_mp`, `jk_cbi`.

The conversion is a unit change and nothing more. The authors have already
normalised the series — `mypc.m` rescales the first principal component to the
standard deviation of the OIS_1Y surprise in the EA-MPD input, which is in
basis points, and `main.m` then divides by 100, which is why their README calls
`pc1` *"scaled to have the same standard deviation as the OIS1Y Monetary
Event-window change (in % p.a.)"*. Re-estimating a scale factor as
`sd(OIS_1Y)/sd(pc1)` on this project's 135 events would renormalise an
already-normalised series to the volatility of a subsample; on this sample that
gives 102.11 instead of 100 and inflates every reported magnitude by 2.1%.

The decomposition is also reconstructed from the raw EA-MPD windows, under the
`_own` suffix. On the 135 overlapping events the reconstruction correlates with
the published series at 0.997 to 0.999, agrees with its poor man's
classification on 87% of events — every disagreement sitting on a near-zero
surprise — and gives the same headline. The equity surprise, which is the
other input to the decomposition, agrees between the published file and the
workbook read here to machine precision, so the remaining gap is about the
principal component rather than about parsing; it is consistent with the PCA
being estimated on a different event set, though nothing here isolates that as
the cause. **No headline result is estimated on the reconstruction**; it exists so that the method can
be audited, and it is reported in `05`.

The Target / Timing / Forward Guidance / QE dimensions are a *separate*
exercise, answering which part of the term structure moved rather than what
kind of news it was. The two are never nested or combined.

That second block is **Altavilla-style, not a replication of Altavilla et al.**
Their short-end restriction — forward guidance and QE do not load on the 1M OIS
— and their four normalisations (Target 1M, Timing 6M, FG 2Y, QE 10Y) are
imposed here. Their third restriction is not: they identify QE as the direction
with the smallest variance between 2 January 2002 and 7 August 2008, which
needs pre-crisis events, and before August 2011 the euro area long end exists
only as German Bund yields that this project does not substitute into the OIS
curve. The rotation used instead defines QE as the direction maximising the 10Y
loading, which mechanically sets forward guidance's 10Y loading to zero. The
two rules are not interchangeable — applied to the same plane the published
rule turns the QE direction substantially — so the factors below should be read
as a stated term-structure rotation, not as the paper's Timing, FG and QE.

## Repository structure

```
R/
  functions/
    data_helpers.R           BNR cleaning, event matching, curve factors
    shock_helpers.R          EA-MPD loading, both shock decompositions
    econometrics.R           robust and clustered OLS, Wald tests
  01_build_analysis_data.R   every data operation in the project
  02_timing_validation.R     is the measurement convention right?
  03_baseline_results.R      what are the baseline coefficients?
  04_yield_curve_results.R   what happens to the curve?
  05_robustness.R            do the results survive?
  06_mechanism.R             cross-market evidence on the mechanism
data/raw/                    inputs (gitignored, see above)
data/processed/              analysis_panel.rds (gitignored, regenerated by 01)
figures/  tables/            generated outputs (gitignored)
```

The architecture is deliberate, and it rests on two rules.

**`01` may transform data. `02` to `06` may not.** All loading, cleaning, shock
construction, event matching and curve-factor construction happens in `01`,
which ends with assertions that constitute a contract: if it finishes without
error, the panel has the shape the analysis expects and nothing downstream
checks again. `02` to `06` read `analysis_panel.rds` and only analyse it. None
of them cleans, matches, joins or renames anything. If a new variable is needed
in three months, there is exactly one place it goes.

**`functions/` may hold reusable operations, never a hidden specification.**
`fit_ols_robust()`, `wald_test()` and `first_available_fixing()` belong there.
`run_headline_model()` never would, because it would conceal which regression
is being estimated. A little repetition in `03` to `06` is worth more than an
abstraction that hides the model — every specification in this project should
be readable as the algebra it is, at the point where it is estimated.

## Running it

Tested with R 4.5.2. Package versions are recorded in `renv.lock`; the
lockfile records package versions and the R version it was created under, and
does not install or pin the R executable itself.

```r
renv::restore()
source("R/01_build_analysis_data.R")   # writes data/processed/analysis_panel.rds
source("R/02_timing_validation.R")
source("R/03_baseline_results.R")
source("R/04_yield_curve_results.R")
source("R/05_robustness.R")
source("R/06_mechanism.R")
```

`02` carries its reference values in its header. If they do not reproduce, the
pipeline has changed and the numbers below need re-checking.

## Results

| Maturity | Monetary policy | Information |
|---|---|---|
| 12M | +0.565 (1.57) | -0.413 (-2.56) |
| 3Y | +0.705 (1.58) | -0.566 (-2.71) |
| 5Y | +0.691 (1.55) | -0.625 (-2.74) |
| 10Y | +0.714 (1.81) | **-0.817 (-3.51)** |
| Level | +0.669 (1.64) | -0.605 (-3.12) |

Basis points per 1 bp of shock, HC0 robust t-statistics in brackets, 135
events. Equality of the two coefficients is rejected at every maturity, most
sharply at the 10Y (p = 0.004), so ECB surprises are not homogeneous from
Romania's point of view.

What is established, and what is not:

* The **pooled negative coefficient is robust** to leave-one-out and to the
  ex ante sample exclusions. All 135 leave-one-out refits lie between -0.88 and
  -0.68 with t below -1.96; excluding 2011-2012, or 2020, or both, strengthens
  rather than weakens it; the median rotation and the poor man's rule give the
  same answer, as do the published and reconstructed shock series. Regime
  heterogeneity is a separate question, treated below.
* The effect is a **level shift**, not an established gradient. The joint null
  that the four maturity coefficients are zero is rejected (p = 0.010); the
  null that they are equal is not (p = 0.143).
* The **magnitude is not a stable parameter**. Equality across policy regimes
  is rejected (p = 0.0007, HC3). The coefficient is negative in the three
  regimes estimated with useful precision, from -0.19 in 2011-2014 to -1.73 in
  2022-2025. The positive estimate in 2020-2021 is highly leverage-sensitive
  and statistically uninformative under HC3 — a statement about precision, not
  about what happened.
* Comparing the **two endpoint regimes** — 2011-2014 against 2022-2025, with
  2015-2021 dropped from that regression — the contrast is statistically
  distinguishable from zero on the 10Y (-1.544, t = -3.73 under HC3) but not on
  the Level (-0.685, t = -1.33). This is a contrast between two endpoint
  samples, not a full-sample post-2022 interaction. Because the two outcomes do
  not tell the same story, it is recorded as heterogeneity rather than built
  into a temporal narrative.
* The **policy-shock side is the weaker half**. Its coefficient is positive
  pooled but falls from 0.71 to 0.25 once the sovereign crisis and 2020 are
  removed.

Economic size: a one standard deviation information shock, 3.11 bp, is
associated with a fall of 2.5 bp in the Romanian 10Y.

The Altavilla-style decomposition produces no individually significant
coefficient at any maturity in the 128-event sample; the estimates are
comparatively imprecise. QE has an increasing profile along the curve.

## Mechanism

The same two shocks, applied to euro area yields already carried in the panel,
place the Romanian result in a recognisable pattern. Responses to a 1 bp shock,
HC0 t-statistics in brackets. Euro area outcomes use EA-MPD intraday windows;
the Romanian leg uses the F(i) fixing window:

| | Monetary policy | Information |
|---|---|---|
| DE 2Y | +1.019 (6.14) | +0.508 (3.05) |
| DE 10Y | +0.511 (3.45) | +0.042 (0.34) |
| IT 10Y | +1.510 (8.64) | **-0.937 (-6.69)** |
| ES 10Y | +1.025 (8.25) | **-0.756 (-6.21)** |
| RO 10Y net of DE 10Y | +0.203 (0.42) | **-0.858 (-4.07)** |

A positive information shock raises the German short end, leaves the Bund 10Y
estimate close to zero and imprecise (+0.042, t = 0.34), and compresses Italian
and Spanish 10Y yields. The Romanian differential is negative, as the Italian
and Spanish responses are; the magnitudes are not compared across the two,
because the windows differ. Policy-shock point estimates are positive across
the comparison, although the Romanian differential is imprecise (+0.203,
t = 0.42).

This is consistent with a sovereign risk-premium channel, and it is as far as
the evidence goes. Identifying that channel would need Romanian-specific
pricing — EUR/RON, a local equity index, sovereign CDS — and a design that does
not condition on post-treatment variables. EUR/USD is reported as an outcome
rather than a control for that reason, and responds to the policy shock but not
to the information shock.

Two caveats are recorded in `06`. The Romanian leg is measured over
the roughly 24-hour F(i) window while every euro area leg is an intraday
window, so these are differential responses rather than changes in an
observable spread. And STOXX50 and `jk_pc1` are excluded as outcomes by
construction: the two shocks are built from precisely those series, so
regressing either on them is an identity.

## Canonical results

`RESULTS.md` lists the project's reported numbers, their source script, and how
each is worded here. It is the reference for any sentence containing a
figure, and it separates the three distinctions that are easiest to blur: the
reconstruction against the published shock series, HC0 against HC3, and the
Romanian window against the euro area window.

## Licence

The original code in this repository is released under the MIT Licence; see
`LICENSE`. That covers the code only. The three raw datasets are not
redistributed here and are not covered by it: the BNR fixing, the EA-MPD
workbook and the published Jarocinski-Karadi series remain subject to their
providers' own terms, and each has to be obtained from the source named above.

## References

Altavilla, C., Brugnolini, L., Gurkaynak, R. S., Motto, R., & Ragusa, G. (2019).
Measuring euro area monetary policy. *Journal of Monetary Economics*, 108, 162-179.

Jarocinski, M., & Karadi, P. (2020). Deconstructing monetary policy surprises:
the role of information shocks. *American Economic Journal: Macroeconomics*, 12(2), 1-43.

McQuade, P., Falagiarda, M., & Tirpak, M. (2015). Spillovers from the ECB's
non-standard monetary policies on non-euro area EU countries: evidence from an
event-study analysis. *ECB Working Paper* No. 1869.
