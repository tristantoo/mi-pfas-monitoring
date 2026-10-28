# PFAS in Michigan public water supplies, 2020–2025

This project runs censored-data environmental statistics on every PFAS result Michigan EGLE has published for public water supplies: 156,174 results at 5,434 sampling points in 2,795 systems. The methods are the ones used for groundwater detection monitoring (EPA 2009 Unified Guidance; Helsel 2012).

**Read the report:** [tristkenn.com/pfas-michigan](https://tristkenn.com/pfas-michigan/). It is also in this repo as `report.html`, rendered from [report.qmd](report.qmd). The pre-registered plan is in [PREREG.md](PREREG.md).

## Findings

- **95.5% of results are nondetects** (< 2 ng/L), so no statewide mean is reported.
- **The federal 4 ppt PFOA/PFOS limit is the binding constraint.**
  - Over a Michigan MCL in some full year: 22 points (15 systems).
  - Over federal 4 ppt: 76 points (52 systems).
  - In each point's latest full year: 15 points (10 systems) for Michigan vs 63 (45) for federal.
- **No measurable statewide decline.** There are 43 increasing vs 32 decreasing FDR-significant trends (Akritas–Theil–Sen test), sign test p = 0.25. The pre-registered falsifier fired.
- **63 verified statistically significant increases (SSIs) at 39 points.** These come from intrawell prediction limits with verification resampling, plus the double quantification rule.
- **Simulation:** with a single reporting limit, RL/2 substitution matched or beat ROS at the pre-registered spread, and Kaplan–Meier was never best. This is the second pre-registered falsifier that fired.

## Run

```
Rscript R/ingest.R   # download current EGLE data -> data/raw/ (+ MANIFEST.txt)
Rscript R/run.R      # QC, analyses A-F -> outputs/, figures/
Rscript -e 'testthat::test_dir("tests/testthat")'
quarto render report.qmd
```

Requires R ≥ 4.3 with EnvStats, NADA2, sf, ggplot2, jsonlite, testthat, and Quarto.

| Path | Contents |
|---|---|
| `R/ingest.R` | Paged ArcGIS download, row-count check, manifest |
| `R/clean.R` | Censored-result parsing, duplicate handling, unit-error rule, QC log |
| `R/stats.R` | Occurrence, running annual average, ATS trends, SSI rules, estimator simulation |
| `R/run.R` | Pipeline, figures |
| `outputs/` | Every table behind the report, as CSV |

Data: State of Michigan / EGLE open data. These are screening estimates, not compliance determinations.

Built with AI coding assistance (Claude Code); I reviewed the methods, code, and results.
