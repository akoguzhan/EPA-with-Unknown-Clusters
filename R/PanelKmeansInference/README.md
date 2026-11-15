# PanelKmeansInference

Panel Kmeans clustering and selective inference for cluster means in panel data.

## Installation

```r
devtools::install_github("akoguzhan/PanelKmeansInference")
```

or, from the RStudio project of the paper, `source("Simulation/00_setup.R")`.

## Main functions (version 0.2.3)

* `panel_kmeans_estimation()`: Panel Kmeans with `Ninit` random initializations (by default K
  randomly chosen units as initial centers), K fixed or chosen by the IC; returns the retained run
  and the full selection outcome.
* `panel_selective_fit()`: clustering plus the cosine transforms of the unit-level series.
* `selective_test(fit, L)`: selective p-value of `L' mu = 0`, conditioning on the selected partition
  and K_IC given the (data-independent) initializations (`conditioning = "partition"`, default,
  computed exactly by walking along the perturbation path), or on all paths (`"paths"`);
  `metric = "xi"` (default), `"shrink"` or `"euclidean"`.
  Restrictions: `L_pairwise()`, `L_cluster()`, `L_joint(type = "cepa" | "homogeneity")`.
* `panel_selective_tests()`, `panel_homogeneity_test()`: pairwise, merged and joint tests.
* `panel_bycluster_tests()`: selective tests that each selected cluster mean is zero, with
  Holm-adjusted p-values, to identify the clusters for which the hypothesis fails.
* `choose_B()`: number of cosines, `B = ceil(c_B T^(2/3))`.
