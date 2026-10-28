# Pre-registration

Written 2026-10-28, after ingest and QC (row counts, censoring pattern, series lengths) and **before any concentration statistic, trend test, or compliance average was computed**. Committed before `R/stats.R` exists.

## Data

Michigan EGLE public water supply PFAS results (`R/ingest.R`): 23,951 sample rows, 2020-08-18 to 2025-12-24, 7 regulated PFAS. After QC (`R/clean.R`, `data/qc_log.csv`): 156,223 analyte results at 5,435 sampling points in 2,795 systems. Every nondetect is reported as `<2` ng/L, so there is a single reporting limit (RL = 2 ppt).

What is known before analysis: ~95% of results are nondetects; 674 points have at least one detect; 208 points have ≥ 8 sampling dates and at least one detect.

## Known biases, stated in advance

- **Informative sampling.** Michigan's rule (R 325.10717d) moves points with detects to quarterly monitoring, so contaminated points have more samples. Statewide summaries therefore use **one value per sampling point** (its maximum), never pooled samples.
- **Drinking water, not monitoring wells.** These are treated or untreated entry points of public supplies. The methods are the ones used for groundwater detection monitoring (EPA 2009 Unified Guidance), but the results describe delivered drinking water.
- **Results are screening estimates, not regulatory determinations.** EGLE's own compliance calculation uses data and rules this repo doesn't have.

## Analyses

**A. Occurrence.** Per analyte: share of points with ≥ 1 detect; 90th/95th/99th percentiles of point maxima (nonparametric). Data sets more than 80% censored get **no mean**: with one RL and that much censoring, any mean is mostly an assumption about the unobserved tail (Helsel 2012). This applies statewide.

**B. Screening running annual average (RAA).** For each point and analyte, the RAA at each sample date is the mean of samples in the trailing 365 days. Nondetects enter as 0, the federal 2024 rule's convention; Michigan's rule text on this was not verified. A sensitivity analysis sets nondetects to the RL (2 ppt), giving an upper bound. Each point's maximum RAA is compared with the Michigan MCLs, and for PFOA/PFOS with the federal 4 ppt MCL. Reported: points over each limit under each convention, and how many points change status between conventions.

**C. Trends.** For each point × analyte with ≥ 8 sampling dates spanning ≥ 2 years and ≥ 3 detects, run the Akritas–Theil–Sen test for censored data (NADA2 `ATSmini`; slope in ppt/yr). Benjamini–Hochberg FDR at q = 0.10 across all tests. Then a two-sided sign test on increasing vs decreasing counts among FDR-significant trends.
- *Falsifier for "levels are falling":* if FDR-significant increases are not outnumbered by decreases (sign test p ≥ 0.05), I won't claim a statewide decline.

**D. Intrawell detection monitoring (Unified Guidance Ch. 6, 19).** For each point × analyte with ≥ 8 dates, the first 8 dates are the baseline.
- Baseline all nondetect: double quantification rule, so an SSI (statistically significant increase) is two consecutive detects after baseline.
- Otherwise: nonparametric prediction limit = baseline maximum, with a 1-of-2 verification resample; an SSI is a later result above the limit that is confirmed by the next result.

Reported: SSI counts by analyte, plus an example series.

**E. Method check (simulation).** Lognormal data (σ_log = 1) with censoring at 50%, 70%, and 85% at a single RL, for n = 8, 20, 50, with 2,000 replicates each. Compare five mean estimators: substitute 0, RL/2, or RL; Kaplan–Meier; lognormal ROS. Report bias and RMSE.
- *Expectation (Helsel 2012):* ROS has the lowest RMSE at high censoring, and RL/2 substitution is biased.
- *Falsifier:* if RL/2 substitution matches ROS RMSE (within 10%) across settings, I report that substitution is adequate here and drop the claim that it's wrong.

**F. Map.** Hexbins (EGLE's location unit) shaded by share of points with any detect, and by any screening RAA over a Michigan MCL.

## Amendments (before any statistic was computed)

1. **QC rule.** 43 detects are ≤ 0.01 ppt and 6 more are under 0.1 ppt, which is below any plausible EPA 537.1 detection limit. They look like µg/L values entered as ng/L. All results under 0.1 ppt are excluded, not converted (`R/clean.R` step 3b).
2. **B, full-year RAA.** A one-sample "annual average" isn't an annual average. B is reported two ways: as pre-registered (every trailing-365-day window), and restricted to windows holding ≥ 4 samples. Separately, points with any single result over an MCL are counted.
3. **D, eligibility.** A point needs > 8 dates (8 baseline plus at least one later result). An exceedance on a point's last sample can't be verified, so it's reported as pending.

## Not doing

No causal claims about sources, no health-risk claims, no individual system named as "in violation".
