# =============================================================================
# Empirical application of "Testing Clustered Equal Predictive Ability with Unknown Clusters":
# forecasting monthly exchange-rate returns (Section 7 of the paper).
#
# One script for the whole application. Run from the RStudio project (paths use here::here()):
#
#   source("Application/Exchange Rates/run_application.R")
#
# Steps:
#   0. (optional, slow) regenerate the forecasts with a rolling window of 60 months
#      (R/generate_forecasts.R -> Data/forecast_panel_final.rds);
#   1. clean the panel, summary statistics and scatter plots (R/prepare_data.R ->
#      Data/forecast_panel_ready.rds, Tables/summary_stats.xlsx, Figures/);
#   2. main results: for each forecasting method (loss differential relative to AR(1)) and each
#      testing function (unconditional; lagged target; post-crisis; post-COVID):
#        - K chosen by the IC (kappa = 1.5) from 2 to 5: selective F_SI,r, joint selective
#          C-EPA, selective homogeneity (merged, joint), weighted Bonferroni, O-EPA, naive test,
#          Split Sample with training proportion 0.1, 0.2, 0.3 (R1 10(b));
#        - each fixed K = 2, ..., 5: the same tests (R1 10(a));
#        - Holm-adjusted p-values across the methods;
#        - cluster-specific selective tests that each selected cluster mean is zero, with the Holm
#          adjustment across the clusters (Section 5.3), to identify the clusters for which EPA
#          fails (reported for every method; the paper discusses those with a C-EPA rejection);
#   3. robustness of the selective tests with K chosen by the IC (Referees 1 and 2): merging
#      parameter r; constant c_B of B = ceil(c_B T^(2/3)); F reference law; conditioning on all
#      paths; Euclidean metric; shrinkage estimator of Xi L; kappa = 3; 10 instead of 100 initializations; local demeaning;
#   4. unit-level tests of EPA for every unit, as a benchmark (Co-Editor 6);
#   5. number of clusters by the IC and by cross-validation (R1 10(a));
#   6. cluster-specific regressions for the significant conditional tests.
#
# Loss differentials: d = Loss_AR1 - Loss_model in the data; the paper's convention is
# Loss_model - Loss_AR1 (negative = improvement). The p-values do not depend on the sign.
#
# Output (Tables/): application_main.xlsx (sheets IC, fixed_K, unit_level, cluster_specific),
# application_robustness.xlsx, Khat_CV_IC.xlsx, cluster_regressions.xlsx; application_results.RData.
# =============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(zoo)
  library(PanelKmeansInference)
  library(clusteredEPA)
})

# Stop early if the installed packages predate this script (e.g. no unit-level tests, no F reference)
if (packageVersion("clusteredEPA") < "0.2.3" || packageVersion("PanelKmeansInference") < "0.2.3")
  stop("Installed packages are older than 0.2.3: restart R and run source(\"Simulation/00_setup.R\").")

app_dir <- here("Application", "Exchange Rates")
tab_dir <- file.path(app_dir, "Tables")
dir.create(tab_dir, showWarnings = FALSE)

# ---- settings -------------------------------------------------------------------------------
regenerate_forecasts <- FALSE   # step 0 takes hours; Data/forecast_panel_final.rds already has
                                # all 13 methods and the window standardization
sample_rule <- "main"           # "main": realized returns winsorized, duplicated currencies dropped;
                                # robustness: "all" (all currencies), "exclude", "none"; see
                                # R/prepare_data.R. Non-main rules add a suffix to the outputs
                                # (e.g. application_main_all.xlsx)
run_prepare          <- TRUE    # step 1 (fast)
Kmax     <- 5
Ninit    <- 100                 # Forgy initializations of Panel Kmeans (per K)
iter.max <- 100
kappa    <- 1.5
props    <- c(0.1, 0.2, 0.3)    # training proportions of the Split Sample test
props_dummy <- 0.1              # for the post-crisis and post-COVID dummies: with 0.2 or 0.3 the
                                # post-crisis dummy is (almost) constant in the testing sample
r_grid   <- c(-2, -5, -10, -20, -Inf)
cB_grid  <- c(0.5, 1, 2)
alpha    <- 0.05
shrink_main <- TRUE             # also report the selective tests with the shrinkage estimator of Xi L
set.seed(1)

if (regenerate_forecasts) source(file.path(app_dir, "R", "generate_forecasts.R"), local = new.env())
if (run_prepare) source(file.path(app_dir, "R", "prepare_data.R"), local = new.env())

invisible(PanelKmeansInference:::.get_cpp_walk())       # compile/load the C++ walk once

# ---- data -----------------------------------------------------------------------------------
sample_rule <- match.arg(sample_rule, c("main", "all", "exclude", "none"))
sfx <- if (sample_rule == "main") "" else paste0("_", sample_rule)
out <- function(f) file.path(tab_dir, sub("\\.xlsx$", paste0(sfx, ".xlsx"), f))
df <- readRDS(file.path(app_dir, "Data", paste0("forecast_panel_ready", sfx, ".rds")))
cat(sprintf("Sample rule '%s': N = %d currencies, T = %d months\n", sample_rule,
            length(unique(df$Asset)), length(unique(df$Date))))
df$time <- as.integer(factor(df$Date))
df$id <- as.integer(factor(df$Asset))
df$post_crisis <- as.numeric(df$Date >= as.yearmon("2008-09"))
df$post_covid <- as.numeric(df$Date >= as.yearmon("2020-01"))
df <- df[order(df$id, df$time), ]

methods <- c(AR_BIC = "d_AR_BIC", ElNet = "d_ElNet", XGBoost = "d_XGBoost_X", SVM = "d_SVM_X",
             RF = "d_RF_X", Zero = "d_ZERO", AR1_JS = "d_AR1_JS", Pooled = "d_POOL",
             EB = "d_EB", Grouped = "d_GRP", FAR = "d_FAR", GlobalGBM = "d_GBM",
             Combination = "d_COMB")
if (any(!methods %in% names(df))) {
  warning("Skipped (not in the data; set regenerate_forecasts <- TRUE to add them): ",
          paste(names(methods)[!methods %in% names(df)], collapse = ", "), immediate. = TRUE)
  methods <- methods[methods %in% names(df)]
}
conds <- list(Unconditional = NULL, Lagged = "Actual_Lag1", Crisis = "post_crisis",
              Covid = "post_covid")

# Moments Z = d H_t for one method and testing function, on the balanced panel
make_panel <- function(v, h) {
  dd <- if (identical(h, "Actual_Lag1")) df[!is.na(df$Actual_Lag1), ] else df
  dd <- dd[!is.na(dd[[v]]), ]
  Z <- if (is.null(h)) matrix(dd[[v]], ncol = 1) else cbind(dd[[v]], dd[[v]] * dd[[h]])
  list(Z = Z, id = dd$id, time = dd$time, asset = unique(dd$Asset[order(dd$id)]))
}

safe <- function(expr) tryCatch(expr, error = function(e) {
  message("  skipped (error): ", conditionMessage(e)); NA_real_ })
pv <- function(x, el = "pval") if (is.list(x)) x[[el]] else NA_real_

# Selective tests, naive test and Split Sample tests for one fitted selection
test_row <- function(p, fit, B, props, Kfixed = NULL, use_props = props) {
  s <- epa_clustered_selective(p$Z, p$id, p$time, fit = fit)
  sp <- sapply(props, function(g) if (!g %in% use_props) NA_real_ else pv(safe(
    if (is.null(Kfixed)) epa_clustered_split(p$Z, p$id, p$time, Kmax = Kmax, prop = g, Ninit = Ninit,
                                             iter.max = iter.max)
    else epa_clustered_split(p$Z, p$id, p$time, K = Kfixed, prop = g, Ninit = Ninit,
                             iter.max = iter.max))))
  names(sp) <- paste0("p_split_", props)
  c(K = fit$K, B = B,
    p_sel = s$pval, p_sel_joint = s$joint_cepa_pval, p_hom = s$hom_test_pval,
    p_hom_joint = s$joint_hom_pval, p_wbonf = s$wbonf_pval, p_oepa = s$overall_pval,
    p_naive = pv(safe(epa_clustered_known(p$Z, p$id, p$time, as.integer(fit$clustering), lrv_par = B))),
    sp, if (is.null(Kfixed)) shrink_cols(p, fit))
}

# Selective tests with the shrinkage estimator of Xi L (alternative procedure; same fit, so the
# selection is identical and only the decomposition changes). Computed only when shrink_main = TRUE.
shrink_cols <- function(p, fit) {
  if (!shrink_main) return(NULL)
  s <- safe(epa_clustered_selective(p$Z, p$id, p$time, fit = fit, metric = "shrink"))
  if (!is.list(s)) return(c(p_sel_sh = NA, p_sel_joint_sh = NA, p_hom_sh = NA, p_hom_joint_sh = NA))
  c(p_sel_sh = s$pval, p_sel_joint_sh = s$joint_cepa_pval, p_hom_sh = s$hom_test_pval,
    p_hom_joint_sh = s$joint_hom_pval)
}

# =============================================================================
# 2. Main results: K by the IC and each fixed K
# =============================================================================
rows_ic <- list(); rows_k <- list(); rows_rob <- list(); rows_unit <- list(); rows_cs <- list()
clusters <- list(); pairwise <- list()
for (m in names(methods)) for (cn in names(conds)) {
  cat(sprintf("[%s] %s, %s\n", format(Sys.time(), "%H:%M:%S"), m, cn))
  p <- make_panel(methods[[m]], conds[[cn]])
  B <- epa_default_B(p$Z, p$id, p$time, Kmax)

  fit <- panel_selective_fit(p$Z, p$id, p$time, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                             kappa = kappa, lrv_par = B)
  rows_ic[[length(rows_ic) + 1]] <- data.frame(method = m, test_function = cn,
                                               as.list(test_row(p, fit, B, props,
                                                                use_props = if (cn %in% c("Crisis", "Covid")) props_dummy else props)))
  clusters[[paste(m, cn)]] <- data.frame(Asset = p$asset, cluster = fit$clustering)

  # Cluster-specific selective tests (Section 5.3): H0^k: mean of cluster k = 0, Holm across k.
  # Cluster means are reported in the paper's convention (Loss_model - Loss_AR1), hence the sign.
  cs <- safe(panel_bycluster_tests(fit = fit))
  if (is.list(cs)) {
    ctr <- -cs$centers
    rows_cs[[length(rows_cs) + 1]] <- data.frame(
      method = m, test_function = cn, K = fit$K, cluster = seq_len(fit$K), size = cs$sizes,
      mean_dL = ctr[, 1], mean_dL_x_h = if (ncol(ctr) > 1) ctr[, 2] else NA_real_,
      p_cluster = cs$bycluster_pvals, p_cluster_holm = cs$bycluster_holm,
      assets = vapply(seq_len(fit$K), function(k)
        paste(p$asset[fit$clustering == k], collapse = ", "), ""), row.names = NULL)
  }

  for (K in 2:Kmax) {
    fK <- panel_selective_fit(p$Z, p$id, p$time, K = K, Ninit = Ninit, iter.max = iter.max, lrv_par = B)
    rows_k[[length(rows_k) + 1]] <- data.frame(method = m, test_function = cn,
                                               as.list(test_row(p, fK, B, 0.1, Kfixed = K)))
  }

  # ---- 3. robustness of the selective tests (K by the IC) ----------------------------------
  sel <- function(f, ...) {
    s <- safe(epa_clustered_selective(p$Z, p$id, p$time, fit = f, ...))
    if (!is.list(s)) return(c(K = NA, p_sel = NA, p_sel_joint = NA, p_hom = NA, p_hom_joint = NA))
    c(K = f$K, p_sel = s$pval, p_sel_joint = s$joint_cepa_pval, p_hom = s$hom_test_pval,
      p_hom_joint = s$joint_hom_pval)
  }
  base <- epa_clustered_selective(p$Z, p$id, p$time, fit = fit)
  rob <- list(baseline = sel(fit),
              paths = sel(fit, conditioning = "paths"),
              euclidean = sel(fit, metric = "euclidean"),
              shrink = sel(fit, metric = "shrink"),
              F_reference = sel(fit, ref = "F"),
              kappa3 = sel(panel_selective_fit(p$Z, p$id, p$time, Kmax = Kmax, Ninit = Ninit,
                                               iter.max = iter.max, kappa = 3, lrv_par = B)),
              S10 = sel(panel_selective_fit(p$Z, p$id, p$time, Kmax = Kmax, Ninit = 10,
                                            iter.max = iter.max, kappa = kappa, lrv_par = B)),
              local_demeaning = sel(panel_selective_fit(p$Z, p$id, p$time, Kmax = Kmax, Ninit = Ninit,
                                                        iter.max = iter.max, kappa = kappa,
                                                        lrv_par = B, demean = "local"),
                                    demean = "local"))
  for (cb in cB_grid) {
    Bc <- epa_default_B(p$Z, p$id, p$time, Kmax, cB = cb)
    rob[[paste0("cB_", cb)]] <- c(sel(panel_selective_fit(p$Z, p$id, p$time, Kmax = Kmax, Ninit = Ninit,
                                                          iter.max = iter.max, kappa = kappa,
                                                          lrv_par = Bc)), B = Bc)
  }
  for (r in r_grid) {
    rob[[paste0("r_", r)]] <- c(K = fit$K,
      p_sel = merge_pvalues(c(base$pairwise_pvalues, base$overall_pval), "Genmean_neq", r),
      p_hom = merge_pvalues(base$pairwise_pvalues, "Genmean_neq", r))
  }
  rows_rob[[length(rows_rob) + 1]] <- bind_rows(lapply(names(rob), function(nm)
    data.frame(method = m, test_function = cn, variant = nm, as.list(rob[[nm]]))))
  pairwise[[paste(m, cn)]] <- base$pairwise_pvalues

  # ---- 4. unit-level tests of EPA for every unit (benchmark) -------------------------------
  ut <- safe(epa_unit_tests(p$Z, p$id, p$time, B = B, p_oepa = base$overall_pval))
  if (is.list(ut)) {
    rows_unit[[length(rows_unit) + 1]] <- data.frame(
      method = m, test_function = cn, p_maxt = ut$p_maxt, p_delta = ut$p_delta,
      p_delta_cc = ut$p_delta_cc, p_choi = ut$p_choi, p_cce_delta = ut$p_cce_delta,
      p_cce_max = ut$p_cce_max, p_pe_sel_joint = power_enhance(base$joint_cepa_pval, ut$p_maxt))
  }
}

holm <- function(tab, by) {
  tab %>% group_by(across(all_of(by))) %>%
    mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "holm"), .names = "{.col}_holm")) %>%
    ungroup() %>% as.data.frame()
}
tab_ic <- holm(bind_rows(rows_ic), "test_function")
tab_k <- holm(bind_rows(rows_k), c("test_function", "K"))
tab_rob <- bind_rows(rows_rob)
tab_unit <- bind_rows(rows_unit)
# cluster-specific tests, with an indicator of a C-EPA rejection (F_SI,r or joint test) at alpha
tab_cs <- bind_rows(rows_cs) %>%
  left_join(tab_ic %>% transmute(method, test_function,
                                 cepa_rejected = (p_sel <= alpha) | (p_sel_joint <= alpha)),
            by = c("method", "test_function")) %>% as.data.frame()

print(tab_ic %>% select(method, test_function, K, p_sel, p_sel_joint, p_hom_joint, p_oepa,
                        p_split_0.1) %>% mutate(across(starts_with("p_"), ~ round(.x, 3))))
print(tab_cs %>% filter(cepa_rejected) %>%
        select(method, test_function, K, cluster, size, mean_dL, p_cluster, p_cluster_holm) %>%
        mutate(across(c(mean_dL, starts_with("p_")), ~ round(.x, 3))))

write.xlsx(list(IC = tab_ic, fixed_K = tab_k, unit_level = tab_unit, cluster_specific = tab_cs),
           out("application_main.xlsx"))
write.xlsx(list(robustness = tab_rob,
                summary = tab_rob %>% select(method, test_function, variant, p_sel, p_sel_joint) %>%
                  pivot_wider(names_from = variant, values_from = c(p_sel, p_sel_joint)) %>%
                  as.data.frame()),
           out("application_robustness.xlsx"))

# =============================================================================
# 5. Number of clusters: IC and cross-validation
# =============================================================================
khat <- bind_rows(lapply(names(methods), function(m) bind_rows(lapply(names(conds), function(cn) {
  p <- make_panel(methods[[m]], conds[[cn]])
  data.frame(method = m, test_function = cn,
             K_IC = tab_ic$K[tab_ic$method == m & tab_ic$test_function == cn],   # selected in step 2
             K_CV = length(unique(panel_kmeans_cv_time(p$Z, p$id, p$time, Kmax = Kmax, nfold = 10,
                                                       Ninit = Ninit, iter.max = iter.max,
                                                       n_cores = 1)$final_cluster)))
}))))
print(khat)
write.xlsx(khat, out("Khat_CV_IC.xlsx"))

# =============================================================================
# 6. Cluster-specific regressions for the significant conditional tests
# =============================================================================
if (all(sapply(c("plm", "lmtest", "sandwich"), requireNamespace, quietly = TRUE))) {
  reg <- list()
  for (m in names(methods)) for (cn in setdiff(names(conds), "Unconditional")) {
    row <- tab_ic[tab_ic$method == m & tab_ic$test_function == cn, ]
    if (!nrow(row) || is.na(row$p_sel) || row$p_sel > alpha) next
    h <- conds[[cn]]
    dd <- if (h == "Actual_Lag1") df[!is.na(df$Actual_Lag1), ] else df
    cl <- clusters[[paste(m, cn)]]
    dd$cluster <- factor(cl$cluster[match(dd$Asset, cl$Asset)])
    X <- model.matrix(~ cluster - 1, dd)
    XH <- X * dd[[h]]
    colnames(XH) <- paste0(colnames(X), "_x_", h)
    rd <- data.frame(y = dd[[methods[[m]]]], X, XH, id = dd$id, time = dd$time)
    fml <- as.formula(paste("y ~ -1 +", paste(c(colnames(X), colnames(XH)), collapse = " + ")))
    mod <- plm::plm(fml, data = plm::pdata.frame(rd, index = c("id", "time")), model = "pooling")
    ct <- lmtest::coeftest(mod, vcov. = plm::vcovSCC(mod, type = "HC0", maxlag = 5))
    reg[[paste(m, cn)]] <- data.frame(method = m, test_function = cn, term = rownames(ct),
                                      estimate = ct[, 1], std.error = ct[, 2], p.value = ct[, 4])
  }
  reg_df <- bind_rows(reg)
  if (nrow(reg_df)) write.xlsx(reg_df, out("cluster_regressions.xlsx"))
}

save(tab_ic, tab_k, tab_rob, tab_unit, tab_cs, khat, clusters, pairwise,
     file = file.path(app_dir, paste0("application_results", sfx, ".RData")))
cat("Done. Tables in", tab_dir, "\n")
