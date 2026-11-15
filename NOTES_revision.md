# Notes on the code revision (October 2026)

The code implements the selective tests of the revised paper. There is no bootstrap test (the
J-EPA bootstrap belongs to a separate paper); joint EPA is tested with the selective tests,
mainly the joint selective C-EPA test.

## 1. PanelKmeansInference (0.2.0)

* **Initialization (Forgy).** Algorithm 1 starts from K units drawn at random, independently
  of the data; their time averages are the first centers (`init_type = "forgy"`, default).
  The random assignment of all units (`init_type = "partition"`) is kept for comparison only.
  With it, the first centers are close to the overall mean, and with separated clusters (P = 1)
  some clusters are empty after the first assignment for almost every start.
* **Best initialization bug.** The old code chose the retained start with
  `sum(apply(objectives, 2, min))`, a minimum over units within each cluster. It now uses the
  within-cluster sum of squares of the final partition. Ties go to the first start, and starts
  with an empty cluster are discarded. If all starts for some K are discarded, further starts are
  drawn.
* **Selection outcome.** The selection outcome is (K_IC, selected partition), computed with the
  initializations held fixed. `pk_cluster()` still records every path, because they are needed
  to compute the truncation set.
* **Selective tests (`selective_test()`).**
  - `conditioning = "partition"` (default). The truncation set is the set of phi >= 0 at which
    the whole procedure, with the same initializations, selects K_IC and the observed partition
    (up to labels). Validity holds conditionally on the initializations, hence uniformly over
    them; the losing starts and the iterations are not conditioned on.

    The set is computed exactly (`R/partition_walk.R`, C++ kernel in `inst/cpp/walk.cpp`):
    1. Each start is followed along phi through the segments on which its Lloyd path is constant
       (quadratic constraints).
    2. The retained start for each K is found by WSS comparisons (quadratics).
    3. K_IC is found by the IC comparisons (polynomials of degree 2P).

    The walk proceeds in blocks downwards from D and stops once the p-value is determined.
    Brute-force checks (re-running the whole procedure on 600 to 1500 values of phi) give the
    same set in every case tried: 12 random designs re-checked on the final code (P = 1, 2; K
    fixed and IC; S = 1, 3; both initializations; partition and paths conditioning).
  - `conditioning = "paths"` gives the earlier, much smaller set. It conditions on all paths, the
    retained starts and K_IC.
  - `metric = "xi"` (default) uses v = Xi_hat L Sigma^(-1/2) dir. `"euclidean"` (Chen and
    Witten) and `"shrink"` are for comparison only.

    `"shrink"` shrinks Xi_hat L towards L(L'L)^(-1) Sigma_hat. Its estimated intensity was close
    to 1 even under cross-sectional dependence. With the wide partition-conditioning sets,
    "shrink" and "euclidean" mis-calibrated the null p-values in a pilot (Kolmogorov-Smirnov
    p = 0.01-0.03, 80 replications), while "xi" was well calibrated (KS p = 0.26-0.52).
  - `ref = "chi"` (default) or `"F"` (Hotelling law). F made the joint tests very conservative in
    a pilot, so it is an option only.
  - Speed (0.2.1): about 0.1-0.5 s per test at N = 80, Kmax = 5, S = 10. Steps 1 and 2 of the
    walk are in C++ (`winners_cpp`), and the walks for K other than K_IC are run only on the
    pieces where the observed partition is retained for K_IC, since the IC matters only there
    (the walk is exact on any subinterval). The p-values are the same as with 0.2.0 (largest
    difference 2e-12 over 231 tests on 48 designs, P = 1, 2, K fixed and IC). The R fallback
    gives the same p-values.
  - The simulated panels use C++ for the AR recursion and the rolling-window forecast errors
    (clusteredEPA 0.2.1, same random draws as the R version). A full replication of `01` at
    N = 80, T = 200 went from about 29 s to 3.6 s.
* **Restrictions.** `L_pairwise()`, `L_cluster()`, `L_joint("cepa" | "homogeneity")` (Remark on
  joint tests).
* **B.** `choose_B()`: B = ceil(c_B T^(2/3)). By default c_B comes from Sun's (2013, eq. 8) CPE
  rule, c_B = 0.42293 |b|^(-1/3), with b = pi^2 rho / (3 (1 - rho)^2) from an AR(1) fit to the
  cross-sectional average. B is bounded below by Kmax*P + 4 and is the same for all
  full-sample tests on a sample.
* Local demeaning of the long-run variance (`demean = "local"`) for the break experiment.

## 2. clusteredEPA (0.2.0)

* `epa_clustered_selective()` returns:
  - F_SI,r;
  - the joint selective C-EPA and homogeneity tests;
  - the merged homogeneity test;
  - the weighted Bonferroni combination;
  - O-EPA, with the same B.

  The old version ignored `Ninit` and `iter.max`.
* `epa_default_B()`, `generate_epa_simulation_data()` (extended design, see 3b).

## 3. Points that the revised paper now reflects

**(a) Conditioning and Proposition 8.** Proposition 8 is restored. It is true because the
selected partition has K_IC clusters. The truncation set nevertheless includes the IC
comparisons, because the event "the procedure selects C" is smaller than "Algorithm 1 with K
clusters returns C". The remark after it explains:
- why the initializations need not be conditioned on (validity uniformly over them);
- why conditioning on all paths is valid but costly.

**(b) Monte Carlo design.** With estimated forecasts, one noise variance calibrated to
E(DL) = psi leaves E(Y_{t-1} DL) far from 0 under the null (0.2 to 0.5). The generator therefore
calibrates both moments, E(H DL) = psi (1, alpha_i)'. It uses forecast noise whose variance
depends on Y_{t-1}, added to Forecaster 1 or 2 depending on the sign. Both moments are matched
within Monte Carlo error.

**(c) Application, signs.** The paper's convention is d = Loss_model - Loss_AR1 (negative =
improvement). Table S.1 and the text used it for Elastic Net, but the other four rows came from
d = Loss_AR1 - Loss_model and had the wrong sign. With the correct sign, the mean loss of every
method is at least that of AR(1):

| AR(1) | AR(p) | ElNet | XGBoost | SVM | RF |
|---|---|---|---|---|---|
| 0.979 | 1.319 | 1.013 | 1.521 | 0.981 | 1.053 |

The sentence saying that AR(p), XGBoost and RF beat AR(1) on average is corrected. In the code:
- `d_` variables are Loss_AR1 - Loss_model throughout (`prepare_data.R` fixed d_ElNet);
- the summary table is written in the paper's convention;
- the p-values do not depend on the sign.

**(d) Application settings.** 100 initializations (10000 before), Kmax = 5, K by the IC and each
fixed K = 2..5, Split Sample with gamma in {0.1, 0.2, 0.3}, Holm adjustment across methods.
The forecasts use a rolling 60-month window ("recursive" corrected in the text).

## 4. Pilot evidence so far (small numbers of replications; to be confirmed by the full runs)

* Under the null, the selective tests with partition conditioning were close to 5% in most
  cells (F_SI 0.03-0.05, merged homogeneity 0.02-0.08). The joint tests had a few cells at
  0.10-0.15 with 60 replications at T = 200, which the full runs should check. Uniformity checks
  of the joint C-EPA p-values gave KS p-values between 0.03 and 0.50.
* Partition conditioning gives truncation sets 10 to 100 times wider than paths conditioning,
  and much smaller p-values under alternatives. For example, at d = 1, K = 3, S = 10 and
  Euclidean metric, the p-values went from 0.02 and 0.65 to 0.0000 and 0.0000.

## 5. Version 0.2.2 (11 October 2026): one script per part, documentation, robustness

* **Scripts.** `Simulation/run_simulations.R` runs all Monte Carlo experiments (N = 80, 1000
  replications); `Application/Exchange Rates/run_application.R` runs the whole application. The
  previous numbered scripts were removed in 0.2.3.
* **Documentation.** Every function of both packages, exported or internal, has a roxygen
  header (title, description, parameters, value).
* **New functions (clusteredEPA).** `epa_unit_tests()`: unit-level tests of joint EPA (Bonferroni
  max test, Delta test of Pesaran and Yamagata with and without a correction for cross-sectional
  correlation, inverse normal combination, defactored tests). `power_enhance()`: weighted
  Bonferroni combination of a selective p-value with the unit-level max test (Fan, Liao and Yao,
  2015), size at most alpha in finite samples. Both are benchmarks only; the paper's tests are the
  selective C-EPA tests. Power enhancement is valid for joint EPA, not for C-EPA.
* **Option `ref`.** `epa_clustered_selective()` and `panel_selective_tests()` pass `ref = "chi"`
  (default) or `"F"` to `selective_test()`. In the pilot with 1000 replications the conditional
  joint test rejected 11% under the null at T = 50 (7% at T = 200) with the chi reference; the
  `variants` experiment reports the F reference.
* **Bug fix.** `Genmean_pcombine()` (positive-order generalized mean, not used by default) reset
  its order to 1/(r - 1) by mistake; removed.

### Referee comments addressed by the code

| Comment | Where |
|---|---|
| Co-Editor 5 (heterogeneity within true clusters) | `run_simulations.R`, experiment `within` |
| Co-Editor 6 (unit-level null, unconditional competitors) | `main` and `jepa_het` (unit-level tests, power enhancement); application step 4 |
| R1 1, 9(e) (time-averaged null, breaks) | `breaks`, local demeaning |
| R1 2(b), 7(b) (orthogonality, metric) | `variants`: Xi vs Euclidean vs shrunk |
| R1 7(a), 8 (multiple starts) | `variants`: S = 1, 10, 100, random-partition initialization, path conditioning |
| R1 8(b) (choice of r) | `main` (r grid); application step 3 |
| R1 8(c), 9(d) (choice of B) | `variants` (c_B, F reference); application step 3 |
| R1 9(a)-(c) (DGP) | `generate_epa_simulation_data()`: fixed effects, factor with heterogeneous loadings, estimated rolling-window forecasts |
| R1 9(f) (power when O-EPA holds) | `main`: joint test and weighted Bonferroni vs F_SI,r |
| R1 10(a) (IC vs CV, fixed K) | `ic`; application steps 2 and 5 |
| R1 10(b) (gamma) | `main` and application step 2: gamma in {0.1, 0.2, 0.3} |
| R2 2 (conditional limits) | `split` (selected clusters within one true cluster) |
| R3 3 (number of clusters) | conditioning on K_IC in all selective tests; `variants`: K = 3, kappa = 3 |

## 6. Version 0.2.3 (3 November 2026): cluster-specific tests, clean-up

* **Cluster-specific tests (Section 5.3 of the paper).** `panel_bycluster_tests()`
  (PanelKmeansInference) computes, for each selected cluster, the selective test of
  theta_k = 0 (restriction L_k = a_k x I_P, the joint test with R = e_k' x I_P), and returns the
  Holm-adjusted p-values, the estimated cluster means and the cluster sizes; it now passes `ref`.
  `epa_clustered_selective(..., cluster_specific = TRUE)` (clusteredEPA) adds them to its output
  (`cluster_pvals`, `cluster_holm`, `cluster_holm_pval`, `cluster_centers`, `cluster_sizes`).
  The default is FALSE, so the simulations are unchanged and do not report these tests.
* **Application.** `run_application.R` computes the cluster-specific tests for K chosen by the IC
  and writes them to the sheet `cluster_specific` of `application_main.xlsx` (and `tab_cs` in
  `application_results.RData`): cluster, size, mean loss differential in the paper's convention
  (Loss_model - Loss_AR1; second column for the interaction with the testing function),
  selective and Holm-adjusted p-values, the currencies of the cluster, and whether C-EPA is
  rejected at 5% by F_SI,r or the joint test.
* **Pilot evidence (not in the paper).** In a pilot with P = 1, K fixed, N = 80, T = 200 and 500
  replications, the cluster-specific tests had size 0.04-0.07 and identified the violating
  clusters (power 0.75-0.85 with separated clusters), but Holm's global test of C-EPA was less
  powerful than O-EPA combined with the pairwise tests in every design and had no power against a
  common shift of all units. Hence their role after a rejection only.
* **Removed files (not used by `run_simulations.R` or `run_application.R`).** The numbered
  simulation and application scripts and `sim_tools.R`, both `legacy/` folders,
  `Application/Exchange Rates/dataprep.m`, `PanelKmeansInference/R/truncated_distribution.r`
  (truncated normal routines of Chen and Witten, replaced by the truncated chi functions in
  `selective.R`), `clusteredEPA/R/forecasting.R` (ARX forecasting functions not used by the
  simulations or the application) and their help pages, and the RStudio session files.
* **Main sample of the application.** `sample_rule <- "winsorize"` (all 122 currencies, realized
  returns winsorized at +-5) is now the default and writes the outputs without suffix;
  `"exclude"` (cut-off 8, 96 currencies) and `"none"` are robustness runs (suffixes `_exclude`,
  `_none`). `prepare_data.R` saves the main panel as `forecast_panel_ready.rds` and the others as
  `forecast_panel_ready_exclude.rds` and `forecast_panel_ready_none.rds`, and computes the summary
  statistics and the scatter plots on the main panel. `regenerate_forecasts` is FALSE by default.
* **Open point.** The selective p-values depend on the draw of the Forgy initializations (16 Oct
  2026: global boosting, homogeneity 0.006 in one draw and 0.85 in another). Averaging the
  p-values over several seeds (twice the mean, valid under any dependence) is proposed but not
  implemented; `run_application.R` uses `set.seed(1)` and one draw.

## 7. Checks before the final runs (6 November 2026)

* `Simulation/00_setup.R` now installs every package the scripts use: `patchwork` and `ggrepel`
  (figures of `R/prepare_data.R`, which `run_application.R` sources), `plm`, `lmtest` and
  `sandwich` (cluster regressions, step 6, silently skipped without them), and `e1071` and
  `randomForest` (forecasts, only when they are regenerated).
* `robustness_samples.R` now includes SVM, so that it covers the 13 methods of the paper and the
  Holm adjustment is across 13 methods, as in the main tables.
* `run_simulations.R`: the `eucl_fail` experiment now runs before the figures of the break
  experiment and the final message; results and seeds are unchanged.
* Help page of `selective_test()`: the shrinkage estimator is described as the alternative
  procedure of the paper, not as a comparison only.
* Smoke tests: all simulation experiments not yet in `Results/` (split, breaks, ic, within,
  jepa_het, eucl_fail) ran with 2 replications per design without errors; the application ran
  end to end with two methods and 5 initializations.
* **Superseded by the next two items.** Known limitation found in the smoke test. The main sample contains 28 currencies whose
  returns duplicate another one exactly (monetary unions and hard pegs), so their time averages
  tie. The default partition conditioning handles the ties, but path conditioning (robustness
  column `paths`) can return an empty truncation set and NA, as for SVM, unconditional.
* **Duplicated currencies dropped from the main sample.** `prepare_data.R` keeps one currency of
  each group of series whose standardized returns have an absolute correlation above 0.98 (the
  anchor when it is in the sample: euro area, United Kingdom, Australia, South Africa, India,
  Singapore). This drops 36 of the 122 currencies and leaves N = 86. Panels: `forecast_panel_ready.rds`
  (main), `_all` (all currencies, winsorized), `_exclude` and `_none` (duplicates dropped),
  `_raw` (all currencies, no treatment; input of `robustness_samples.R`). `sample_rule` is now
  "main", "all", "exclude" or "none". `robustness_samples.R` applies the same rule to every sample
  except `all_winsor5`, and uses a new checkpoint file (`robustness_samples_checkpoint_v2.rds`).
  With the ties gone, path conditioning gives p-values in all cases of the smoke test.
* **Split Sample with dummies.** In the application the Split Sample test for the post-crisis and
  post-COVID dummies uses the training proportion 0.1 only (with 0.2 or 0.3 the post-crisis
  dummy is almost constant in the testing sample). `epa_clustered_split()` now ignores the
  columns of Z that are identically zero in the training sample when estimating the clusters (the
  product with a dummy that is zero there), which otherwise make the IC determinant zero. This
  does not affect the simulations, where no column is zero.
