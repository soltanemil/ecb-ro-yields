# Monetary Policy or Information? ECB Spillovers to Romanian Sovereign Yields

## Research question

Do ECB surprises transmit to the Romanian sovereign yield curve, and does it
matter whether a surprise is monetary policy news or information the ECB
reveals about the economy?

Romania sits outside the euro area and runs its own monetary policy. The two
markets are linked all the same, through euro funding, through non-resident
participation in the government bond market, and through regional risk
appetite. Whether the distinction between the two kinds of surprise matters is
an empirical question, and here it does. The estimated Romanian responses to
the two components carry opposite signs, so an aggregate surprise measure
averages over heterogeneity that the decomposition brings out.

## Identification and timing

Everything in the specification follows from a single timing asymmetry.

| | Time (Frankfurt) |
|---|---|
| BNR government securities fixing | 11:00 (12:00 Bucharest) |
| ECB announcement, to 20 July 2022 | 13:45 |
| ECB announcement, from 21 July 2022 | 14:15 |

In both regimes the Romanian fixing on the announcement day is struck
**before** the announcement, and is therefore predetermined with respect to the
surprise. The response has to be read off the first fixing able to incorporate
it.

```
dy_i = y_{F(i)} - y_{F(i)-1}
```

where `F(i)` indexes the first BNR fixing struck strictly after ECB event `i`.
Indexing the fixings that actually exist handles an event falling on a Romanian
public holiday, carries weekends and holiday bridges correctly, and discards no
events. A calendar lead such as `lead(date, 1)` would break on all three
counts.

`R/02_timing_validation.R` puts the convention through a falsification test
instead of asserting it. Four conditions have to hold jointly, and all four do.
The Bund reacts strongly inside the EA-MPD intraday window. The estimate on the
predetermined Romanian fixing is small and statistically indistinguishable from
zero. The first subsequent fixing shows the largest response, and its 10Y
estimate is statistically distinguishable from zero. The estimate for the
window after that is again small and indistinguishable from zero.

## Data

| Dataset | Source | Coverage |
|---|---|---|
| Romanian government bond yields (6M, 12M, 3Y, 5Y, 10Y) | National Bank of Romania, daily bid-ask fixing | 2011-01-10 to 2026-04-27, 3,844 observations |
| ECB intraday surprises | EA-MPD, Altavilla et al. (2019) | 315 events from 1999; 135 in the Romanian overlap |
| Decomposed policy and information shocks | Jarocinski and Karadi, published series | the same 135 events |

None of the three inputs is redistributed here. Reproducing the project means
placing three files in `data/raw/`.

- `titluri_de_stat_ro.xlsx`, the BNR government securities fixing
- `ECB_surprise_shocks.xlsx`, the EA-MPD workbook
- `jk_shocks_official.csv`, which is `shocks_ecb_mpd_me_d.csv` from
  [github.com/marekjarocinski/jkshocks_update_ecb](https://github.com/marekjarocinski/jkshocks_update_ecb),
  renamed. Only the rows from 2011 onwards are used. `load_jk_official()`
  enforces that with an explicit filter, so the reported event counts do not
  depend on how the file was trimmed before saving.

## Shock measures

The primary measure is the published Jarocinski-Karadi decomposition, converted
from percentage points to basis points by multiplying by 100. In the code it
carries no suffix and appears as `jk_mp` and `jk_cbi`.

That conversion changes units and does nothing else. Normalisation has already
been carried out by the authors. `mypc.m` rescales the first principal
component to the standard deviation of the OIS_1Y surprise in the EA-MPD input,
which is measured in basis points, and `main.m` then divides by 100. Their
README accordingly describes `pc1` as *"scaled to have the same standard
deviation as the OIS1Y Monetary Event-window change (in % p.a.)"*. Re-estimating
a scale factor as `sd(OIS_1Y)/sd(pc1)` over this project's 135 events would
renormalise an already normalised series to the volatility of a subsample. On
this sample it returns 102.11 in place of 100, inflating every reported
magnitude by 2.1%.

A second version of the decomposition is built here from the raw EA-MPD
workbook and carries the `_own` suffix. It transcribes the authors' own MATLAB,
`mypc.m` for the principal component and `signrestr_median.m` for the
sign-restricted rotation, and it applies the same event filtering, excluding
the three joint Fed and ECB announcements that `main.m` drops before
extraction. The result reproduces the published series exactly. Compared
against the published values with no rescaling, `pc1` and both poor man's
series are bit-identical, and the two median shocks deviate by at most
5.0e-07 bp, which is the eight-decimal rounding the published file carries at
export. Correlations are 1.000 with a unit slope, the poor man's classification
agrees on all 135 events, and regressions on the two series return the same
coefficients to three decimals.

What that agreement establishes is the implementation, and nothing beyond it.
The two series are one estimator implemented twice, so their coincidence says
that the decomposition described in the paper is the decomposition in the
published file, and that the workbook is read correctly here. It carries no
independent empirical content. **No headline result is estimated on the
reconstruction.** The published series stays primary, being the authors' own
output, and the replication check is reported in `05`.

The Target, Timing, Forward Guidance and QE dimensions form a *separate*
exercise. They answer which part of the term structure moved, where the
decomposition above answers what kind of news moved it. The two exercises are
never nested or combined.

**That second block is Altavilla-style. It does not replicate Altavilla et
al.** Their short-end restriction, under which forward guidance and QE do not
load on the 1M OIS, is imposed here, as are their four normalisations (Target
1M, Timing 6M, FG 2Y, QE 10Y). Their third restriction is left out. They
identify QE as the direction with the smallest variance between 2 January 2002
and 7 August 2008, which requires pre-crisis events, and before August 2011 the
euro area long end exists only as German Bund yields, which this project does
not substitute into the OIS curve. The rotation used in their place defines QE
as the direction maximising the 10Y loading, and that mechanically sets the 10Y
loading of forward guidance to zero. The two rules are not interchangeable.
Applied to the same plane, the published rule turns the QE direction
substantially, so the factors below should be read as a stated term-structure
rotation and not as the paper's Timing, FG and QE.

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

Two rules govern the layout, and both are deliberate.

**`01` may transform data. `02` to `06` may not.** All loading, cleaning, shock
construction, event matching and curve-factor construction happens in `01`.
That script ends with assertions which constitute a contract. If it finishes
without error, the panel has the shape the analysis expects, and nothing
downstream checks again. `02` to `06` read `analysis_panel.rds` and analyse it.
None of them cleans, matches, joins or renames anything. A new variable needed
in three months has exactly one place to go.

**`functions/` may hold reusable operations, never a hidden specification.**
`fit_ols_robust()`, `wald_test()` and `first_available_fixing()` belong there.
`run_headline_model()` never would, because it would conceal which regression
is being estimated. A little repetition across `03` to `06` is worth more than
an abstraction that hides the model. Every specification in this project should
be readable as the algebra it is, at the point where it is estimated.

## Running it

Tested with R 4.5.2. `renv.lock` records the package versions and the R version
under which it was created. It does not install or pin the R executable itself.

```r
renv::restore()
source("R/01_build_analysis_data.R")   # writes data/processed/analysis_panel.rds
source("R/02_timing_validation.R")
source("R/03_baseline_results.R")
source("R/04_yield_curve_results.R")
source("R/05_robustness.R")
source("R/06_mechanism.R")
```

`02` carries its reference values in its header. Should they fail to reproduce,
the pipeline has changed and the numbers below need re-checking.

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
sharply at the 10Y (p = 0.004). Seen from Romania, then, ECB surprises are not
homogeneous.

Five points set out what the estimates support and where they stop short.

* The **pooled negative coefficient is robust** to leave-one-out and to the
  episode-based sample exclusions. All 135 leave-one-out refits lie between
  -0.88 and -0.68 with t below -1.96. Under each episode-based exclusion the
  coefficient remains negative and precisely estimated, though its magnitude
  varies across the exclusion samples. The median rotation and the poor man's
  decomposition give the same qualitative result, though the poor man's
  decomposition imposes the additional restriction that only one shock is
  present in each announcement. Regime heterogeneity is a separate question,
  treated below.
* The response is **broad across maturities**. The joint null that the four
  CBI maturity coefficients are zero is rejected (p = 0.010). The joint
  equality test leaves the data compatible with a common response (p = 0.143),
  while the endpoint slope contrast is negative and borderline under HC3
  (p = 0.048). Evidence for a maturity gradient is therefore sensitive to the
  test used.
* The **magnitude is not a stable parameter**. Equality across policy regimes
  is rejected (p = 0.0007, HC3). In the three regimes estimated with useful
  precision the coefficient is negative, running from -0.19 in 2011-2014 to
  -1.73 in 2022-2025. The positive estimate in 2020-2021 is highly
  leverage-sensitive and statistically uninformative under HC3, which is a
  statement about precision and not about what happened.
* Comparing the **two endpoint regimes**, 2011-2014 against 2022-2025 with
  2015-2021 dropped from that regression, the contrast is precise on the 10Y
  (-1.544, t = -3.73 under HC3) and imprecise on the Level (-0.685, t = -1.33).
  The evidence from this endpoint comparison is therefore outcome-specific. It
  does not by itself establish that the two outcome contrasts differ from one
  another. It is a contrast between two endpoint samples, and a full-sample
  post-2022 interaction would be a separate estimand.
* The **policy-shock side is the weaker half**. Its pooled coefficient is
  positive, and it falls from 0.71 to 0.25 once the sovereign crisis and 2020
  are removed.

In economic terms, a one standard deviation information shock of 3.11 bp is
associated with a fall of 2.5 bp in the Romanian 10Y.

The Altavilla-style decomposition produces no individually significant
coefficient at any maturity in the 128-event sample, and the estimates are
comparatively imprecise. QE has an increasing profile along the curve.

## Mechanism

Applied to the euro area yields already carried in the panel, the same two
shocks place the Romanian result in a recognisable pattern. The table reports
responses to a 1 bp shock with HC0 t-statistics in brackets. Euro area outcomes
use the EA-MPD Monetary Event Window, while the Romanian leg uses the F(i)
fixing window.

| | Monetary policy | Information |
|---|---|---|
| DE 2Y | +1.019 (6.14) | +0.508 (3.05) |
| DE 10Y | +0.511 (3.45) | +0.042 (0.34) |
| IT 10Y | +1.510 (8.64) | **-0.937 (-6.69)** |
| ES 10Y | +1.025 (8.25) | **-0.756 (-6.21)** |
| RO 10Y net of DE 10Y | +0.203 (0.42) | **-0.858 (-4.07)** |

A positive information shock raises the German short end, leaves the Bund 10Y
estimate close to zero and imprecise (+0.042, t = 0.34), and is associated with
lower Italian and Spanish 10Y yields. The Romanian differential is negative, as
the Italian and Spanish responses are. Magnitudes are not compared across the
two, because the windows differ. Policy-shock point estimates are positive
across the comparison, although the Romanian differential is imprecise
(+0.203, t = 0.42).

All of this is consistent with a sovereign risk-premium channel, and that is as
far as the evidence goes. Identifying the channel would call for
Romanian-specific pricing (EUR/RON, a local equity index, sovereign CDS)
together with a design that handles their post-treatment status explicitly.
Conditioning on a potential mediator in the headline regression would change
the estimand and can introduce post-treatment bias, which is why EUR/USD enters
as an outcome. It responds to the policy shock, while its response to the
information shock is indistinguishable from zero.

Two caveats are recorded in `06`. The Romanian leg is measured over the F(i)
window, usually about 24 hours and longer across weekends and holidays, while
every euro area leg is measured over the full EA-MPD Monetary Event Window,
which spans the press release and the press conference. The responses are
therefore differential, and reading an observed spread change off them would
require a common window. STOXX50 and `jk_pc1` are excluded as outcomes by
construction, since the two shocks are built from precisely those series and
regressing either on them is an identity.

## Canonical results

`RESULTS.md` lists the project's reported numbers, their source script, and how
each is worded here. Any sentence containing a figure can be checked against
it. It also keeps apart the three distinctions that are easiest to blur. These
are the reconstruction against the published shock series, HC0 against HC3, and
the Romanian window against the euro area window.

## Licence

The repository is released under the MIT Licence, reproduced in `LICENSE`. It
covers the code and the documentation here. The three raw datasets form no part
of this repository and are not redistributed. The BNR fixing, the EA-MPD
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
