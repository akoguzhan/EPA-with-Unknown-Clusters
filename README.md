# Testing Clustered Equal Predictive Ability with Unknown Clusters: code

RStudio project: open `EPA with Unknown Clusters.Rproj`. All scripts locate files with
`here::here()` relative to this folder, so no path needs editing.

## Structure

```
R/PanelKmeansInference   package: Panel Kmeans, information criterion, selective inference
R/clusteredEPA           package: EPA tests (O-EPA, predetermined, naive, Split Sample,
                         selective, cluster-specific, unit-level benchmarks), simulation design
Simulation/              run_simulations.R (all Monte Carlo experiments), results in Results/
Application/Exchange Rates/  run_application.R (the whole empirical application),
                         robustness_samples.R (other treatments of extreme returns)
```

## 1. Setup (once, in a fresh R session)

```r
source("Simulation/00_setup.R")
```

It installs the CRAN dependencies and the two packages, and compiles the C++ routines
(truncation sets of the selective tests in `PanelKmeansInference/inst/cpp/walk.cpp`, simulated
panels in `clusteredEPA/inst/cpp/gen.cpp`; Rtools is needed on Windows). Without a compiler the
packages use R versions with the same results, many times slower. Re-run it, after restarting R,
whenever the package code changes. All functions of both packages have roxygen documentation
(`?selective_test`, `?epa_clustered_selective`, ...).

## 2. Simulations

```r
source("Simulation/run_simulations.R")
```

N = 80 and 1000 replications in every design. Choose the experiments with the vector
`experiments` at the top of the script; each finished design is saved immediately, so the script
can be interrupted and restarted. Output in `Simulation/Results`: `results_<experiment>.rds`
(p-values of every replication) and `results_<experiment>.xlsx` (rejection rates at 5%; suffix
`_cond` for the conditional tests). The break experiment also writes the two figures of the paper
(`*_CEPA_size_sbreaks.pdf`; copy them to the paper's `Figures/` folder).

| Experiment | Content | Referee comment |
|---|---|---|
| `main` | Size (T = 20, 50, 100, 200) and power (Cases 1 and 2, T = 50, 200, psi = 0.125..0.5): Predetermined, Naive, Split Sample (gamma = 0.1, 0.2, 0.3), selective F_SI,r, joint selective C-EPA, homogeneity (merged, joint), weighted Bonferroni, O-EPA; r in {-2, -5, -10, -20, -Inf}; unit-level tests of joint EPA and power-enhanced selective tests | R1 8(b), 9(f), 10(b); Co-Editor 6 |
| `variants` | Conditioning on partition vs all paths; Xi vs Euclidean vs shrunk metric; chi vs F reference; c_B in {0.5, 1, 2}; Forgy vs random-partition initialization; S = 1, 10, 100; K = 3 fixed; kappa = 3 | R1 2(b), 7, 8, 9(d), 10(a); R2 1; R3 3 |
| `split` | A true cluster split by Panel Kmeans (K0 = 2, K = 3) | R2 2 |
| `breaks` | Breaks in the cluster means, null holding on average; global vs local demeaning | R1 1, 9(e) |
| `ic` | K by the IC (kappa = 1.5, 3) and by cross-validation | R1 10(a) |
| `within` | Heterogeneity within the true clusters | Co-Editor 5 |
| `jepa_het` | All units heterogeneous: clustering as dimension reduction for joint EPA vs unit-level tests | Co-Editor 6 |

Running time per replication (CPU, both test types): about 3-5 s for `main` and 1-3 s for the
smaller experiments; about 40 s for `variants`, dominated by S = 100 initializations. The full
study needs roughly 200 CPU hours, mostly `variants`; use a server and set `n_cores`.

## 3. Application

```r
source("Application/Exchange Rates/run_application.R")
```

Main sample: the standardized realized returns winsorized at +-5 and one currency of each group
of duplicated series (monetary unions, hard pegs) kept, N = 86 (`sample_rule <- "main"`);
`"all"` (all 122 currencies), `"exclude"` and `"none"` give the robustness runs, with a suffix in
the output file names. For the post-crisis and post-COVID dummies the Split Sample test is
computed with the training proportion 0.1 only. Steps (see the header of the script): optional regeneration of the
forecasts (slow, off by default); data preparation, summary statistics and scatter plots; main results with K chosen by
the IC and for each fixed K = 2..5, Split Sample with gamma in {0.1, 0.2, 0.3}, Holm adjustment
across methods; robustness of the selective tests (r, c_B, F reference, path conditioning,
Euclidean metric, kappa = 3, 10 initializations, local demeaning); unit-level tests of joint EPA;
cluster-specific selective tests that each selected cluster mean is zero, with the Holm
adjustment across clusters (Section 5.3 of the paper), which identify the clusters for which EPA
fails; K by the IC and by cross-validation; cluster-specific regressions. Output in
`Application/Exchange Rates/Tables` (`application_main.xlsx` with sheets `IC`, `fixed_K`,
`unit_level` and `cluster_specific`, `application_robustness.xlsx`, `Khat_CV_IC.xlsx`,
`cluster_regressions.xlsx`) and `application_results.RData`.

`robustness_samples.R` runs the main tests for the other treatments of extreme returns and
subsamples reported in the supplement (output `Tables/robustness_samples.xlsx`); run it after
`run_application.R`, which prepares the data.

See `NOTES_revision.md` for the changes to the code and the correspondence with the referees'
comments.
