# =============================================================================
# Monte Carlo study of "Testing Clustered Equal Predictive Ability with Unknown Clusters"
#
# One script for all experiments. Run from the RStudio project (paths use here::here()), after
# Simulation/00_setup.R has installed the two packages in a fresh R session.
#
#   source("Simulation/run_simulations.R")
#
# Every experiment uses N = 80 units and n_sim = 1000 replications. Each design is saved as soon
# as it is finished, so the script can be interrupted and restarted: finished designs are
# skipped. Results are written to Simulation/Results:
#   results_<experiment>.rds   p-values (and diagnostics) of every replication
#   results_<experiment>.xlsx  rejection rates at 5% of every column starting with "p_",
#                              averages of the other columns
# Columns with suffix "_cond" are the conditional tests (testing functions H = (1, Y_{t-1})'),
# the others the unconditional tests (H = 1).
#
# Experiments (set `experiments` below to run a subset):
#
#   main       Size and power of the C-EPA tests (Section 6 of the paper).
#              Tests: predetermined clusters; naive; Split Sample with training proportion
#              0.1, 0.2, 0.3 (R1 10(b)); selective F_SI,r; joint selective C-EPA test;
#              selective homogeneity tests (merged, joint); weighted Bonferroni; O-EPA.
#              Robustness to the merging parameter r in {-2, -5, -10, -20, -Inf} (R1 8(b)).
#              Benchmarks for the hypothesis that EPA holds for every unit (Co-Editor 6):
#              unit-level max, Delta, correlation-corrected Delta, Choi and defactored tests,
#              and the power-enhanced selective tests.
#   variants   Robustness of the selective tests (Referees 1 and 2):
#              conditioning on the selected partition vs on all paths (R1 7(a), 8);
#              Xi metric vs Euclidean metric of Chen and Witten vs shrunk Xi (R1 2(b), 7(b));
#              chi vs F reference law; the constant c_B of B = ceil(c_B T^(2/3)) (R1 8(c), 9(d));
#              Forgy vs random-partition initialization; S = 1, 10, 100 initializations
#              (R1 7(a), 8); K fixed at 3, and the IC with kappa = 1.5 or 3 (R1 10(a), R3 3).
#   split      Size of the selective test when a true cluster is split (K0 = 2, K = 3).
#   breaks     Breaks in the cluster means, C-EPA holding on average over time (R1 1, 9(e)):
#              global vs local demeaning of the long-run variance. Figures for the paper.
#   ic         Estimation of K: IC with kappa = 1.5 and 3, and cross-validation (R1 10(a)).
#   within     Heterogeneity within the true clusters (Co-Editor 5): C-EPA holds for the true
#              clusters while units within a cluster differ (variant "null"), and separated
#              true clusters with heterogeneity within them (variant "alternative").
#   eucl_fail  Euclidean decomposition not valid: metric of Xi, Euclidean, shrinkage
#   jepa_het   All units heterogeneous, no true clusters (Co-Editor 6): clustering as dimension
#              reduction for joint EPA, against unit-level tests.
#
# Running time (CPU, N = 80, one replication, both test types): about 3-5 s for "main",
# 15-40 s for "variants" (S = 100 and path conditioning dominate), 1-3 s for the others.
# Use a server for the full study; set n_cores below.
# =============================================================================
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  library(parallel)
  library(dplyr)
  library(openxlsx)
  library(PanelKmeansInference)
  library(clusteredEPA)
})

# ---- settings -------------------------------------------------------------------------------
n_sim   <- 1000                  # replications per design
N       <- 80                    # units
Kmax    <- 5                     # largest number of clusters for the IC
Ninit   <- 10                    # Forgy initializations of Panel Kmeans (per K)
n_cores <- NULL                  # NULL: all physical cores but one
seed    <- 20250708
experiments <- c("main", "variants", "split", "breaks", "ic", "within", "jepa_het", "eucl_fail")

res_dir <- here("Simulation", "Results")
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# =============================================================================
# Tools
# =============================================================================

# True clusters of Section 6: N/4, N/4 and N/2 units
true_clusters <- function(N) rep(1:3, c(floor(N / 4), floor(N / 2) - floor(N / 4), N - floor(N / 2)))

# Expected loss differentials of the three clusters
psi_values <- function(d, case) {
  switch(case,
         size = c(0, 0, 0),
         overall_fails = d / 2 + d * c(-1.2, -0.8, 1),    # Case 1: O-EPA fails
         overall_holds = d * c(-1.2, -0.8, 1),            # Case 2: O-EPA holds
         stop("unknown case"))
}

# NT x P matrix of moments (unit-major): unconditional H = 1, conditional H = (1, Y_{t-1})'
build_Z <- function(sim, conditional) {
  DL <- as.vector(sim$DL)
  if (conditional) cbind(DL, DL * as.vector(sim$Ylag)) else matrix(DL, ncol = 1)
}

safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
pv <- function(x, el = "pval") if (is.list(x)) x[[el]] else NA_real_

# Adjusted Rand index, recovery of the true partition, number of clusters
cl_stats <- function(cl, gamma) {
  if (is.null(cl) || anyNA(cl)) return(c(ari = NA, recovery = NA, Khat = NA))
  K0 <- length(unique(gamma))
  c(ari = mclust::adjustedRandIndex(cl, gamma),
    recovery = as.numeric(length(unique(cl)) == K0 && same_cl(cl, gamma, K0)),
    Khat = length(unique(cl)))
}

# Merged p-values for several r (robustness to the merging parameter)
r_grid <- c(-2, -5, -10, -20, -Inf)
merge_r <- function(p, prefix) {
  out <- sapply(r_grid, function(r) safe(merge_pvalues(p, "Genmean_neq", r)))
  names(out) <- paste0(prefix, "_r", ifelse(is.infinite(r_grid), "Inf", abs(r_grid)))
  out
}

# The selective tests of the paper on a fitted selection
sel_block <- function(Z, id, time, fit, tag = "", ...) {
  s <- safe(epa_clustered_selective(Z, id, time, fit = fit, ...))
  nm <- paste0(c("p_sel", "p_sel_joint", "p_hom", "p_hom_joint", "p_wbonf"), tag)
  if (!is.list(s)) return(setNames(rep(NA, 5), nm))
  setNames(c(s$pval, s$joint_cepa_pval, s$hom_test_pval, s$joint_hom_pval, s$wbonf_pval), nm)
}

# Run n_sim replications of rep_fun(sim, ...) in parallel (PSOCK cluster: works on Windows)
run_parallel <- function(rep_fun, n_sim, args, n_cores = NULL, seed = 1) {
  if (is.null(n_cores)) n_cores <- max(1, detectCores(logical = FALSE) - 1)
  invisible(PanelKmeansInference:::.get_cpp_walk())                  # compile the C++ code once
  invisible(clusteredEPA:::.get_cpp_gen("rolling_errors_cpp"))
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  clusterSetRNGStream(cl, iseed = seed)
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(PanelKmeansInference); library(clusteredEPA) })
    NULL
  })
  fns <- Filter(function(x) is.function(get(x, envir = globalenv())), ls(globalenv()))
  clusterExport(cl, c(fns, "r_grid"), envir = globalenv())
  res <- parLapply(cl, seq_len(n_sim), function(s, repf, repargs) {
    out <- tryCatch(do.call(repf, c(list(sim = s), repargs)),
                    error = function(e) structure(conditionMessage(e), class = "rep_error"))
    if (inherits(out, "rep_error")) return(data.frame(error = as.character(out)))
    as.data.frame(as.list(out))
  }, repf = rep_fun, repargs = args)
  res <- bind_rows(res)
  if ("error" %in% names(res)) {
    n_err <- sum(!is.na(res$error))
    if (n_err) warning(n_err, " replications failed; first error: ", res$error[!is.na(res$error)][1])
    res <- res[is.na(res$error), setdiff(names(res), "error"), drop = FALSE]
  }
  res
}

# Rejection rates at level alpha of all columns starting with "p_", averages of the others
rej_rates <- function(df, alpha = 0.05) {
  pc <- grep("^p_", names(df), value = TRUE)
  oc <- setdiff(names(df), pc)
  c(sapply(df[pc], function(x) mean(x <= alpha, na.rm = TRUE)),
    sapply(df[oc], function(x) mean(as.numeric(x), na.rm = TRUE)),
    n_ok = nrow(df))
}

# Append a design to results_<stub>.rds and rewrite results_<stub>.xlsx
save_design <- function(raw, design, stub) {
  raw <- cbind(design[rep(1, nrow(raw)), , drop = FALSE], raw)
  f_rds <- file.path(res_dir, paste0("results_", stub, ".rds"))
  all_raw <- if (file.exists(f_rds)) bind_rows(readRDS(f_rds), raw) else raw
  saveRDS(all_raw, f_rds)
  key <- names(design)
  summ <- all_raw %>% group_by(across(all_of(key))) %>%
    group_modify(~ as.data.frame(as.list(rej_rates(as.data.frame(.x))))) %>% ungroup()
  write.xlsx(as.data.frame(summ), file.path(res_dir, paste0("results_", stub, ".xlsx")))
  invisible(summ)
}

design_done <- function(design, stub) {
  f_rds <- file.path(res_dir, paste0("results_", stub, ".rds"))
  if (!file.exists(f_rds)) return(FALSE)
  old <- readRDS(f_rds)
  key <- names(design)
  keys_old <- do.call(paste, c(lapply(old[key], as.character), sep = "|"))
  any(keys_old == paste(sapply(design, as.character), collapse = "|"))
}

# Run all designs of one experiment
run_experiment <- function(stub, rep_fun, designs, extra = list()) {
  for (j in seq_len(nrow(designs))) {
    des <- designs[j, , drop = FALSE]
    if (design_done(des, stub)) next
    cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%H:%M:%S"), stub,
                paste(names(des), sapply(des, as.character), sep = " = ", collapse = ", ")))
    des_seed <- seed + sum(utf8ToInt(paste(stub, paste(sapply(des, as.character), collapse = "|"))))
    raw <- run_parallel(rep_fun, n_sim, c(as.list(des), extra), n_cores = n_cores, seed = des_seed)
    save_design(raw, des, stub)
  }
}

# =============================================================================
# Replications
# =============================================================================

# ---- main: size and power -------------------------------------------------------------------
rep_main <- function(sim, N, Tobs, d, case, Kmax, Ninit) {
  gamma <- true_clusters(N)
  dat <- generate_epa_simulation_data(N, Tobs, gamma, psi_values(d, case))
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    B <- epa_default_B(Z, id, time, Kmax)
    fit <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B)
    sel <- epa_clustered_selective(Z, id, time, fit = fit)
    spl <- lapply(c(0.1, 0.2, 0.3), function(g)
      safe(epa_clustered_split(Z, id, time, Kmax = Kmax, prop = g, Ninit = Ninit)))
    ut <- safe(epa_unit_tests(Z, id, time, B = B, p_oepa = sel$overall_pval))
    u <- function(el) if (is.list(ut)) ut[[el]] else NA
    res <- c(
      p_known = pv(safe(epa_clustered_known(Z, id, time, gamma, lrv_par = B))),
      p_naive = pv(safe(epa_clustered_known(Z, id, time, as.integer(fit$clustering), lrv_par = B))),
      p_split_0.1 = pv(spl[[1]]), p_split_0.2 = pv(spl[[2]]), p_split_0.3 = pv(spl[[3]]),
      p_sel = sel$pval, p_sel_joint = sel$joint_cepa_pval,
      p_hom = sel$hom_test_pval, p_hom_joint = sel$joint_hom_pval,
      p_wbonf = sel$wbonf_pval, p_oepa = sel$overall_pval,
      merge_r(c(sel$pairwise_pvalues, sel$overall_pval), "p_sel"),
      merge_r(sel$pairwise_pvalues, "p_hom"),
      # benchmarks for EPA for every unit
      p_unit_maxt = u("p_maxt"), p_unit_delta = u("p_delta"), p_unit_delta_cc = u("p_delta_cc"),
      p_unit_choi = u("p_choi"), p_unit_cce_delta = u("p_cce_delta"), p_unit_cce_max = u("p_cce_max"),
      p_pe_sel_joint = safe(power_enhance(sel$joint_cepa_pval, u("p_maxt"))),
      p_pe_sel = safe(power_enhance(sel$pval, u("p_maxt"))),
      cl_stats(fit$clustering, gamma),
      setNames(if (is.list(spl[[2]])) cl_stats(spl[[2]]$clustering, gamma) else rep(NA, 3),
               c("ari_split", "recovery_split", "Khat_split")),
      B = B, in_set = mean(sel$in_set))
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- variants: robustness of the selective tests --------------------------------------------
rep_variants <- function(sim, N, Tobs, d, case, Kmax, Ninit) {
  gamma <- true_clusters(N)
  dat <- generate_epa_simulation_data(N, Tobs, gamma, psi_values(d, case))
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    B <- epa_default_B(Z, id, time, Kmax)
    fit <- function(...) panel_selective_fit(Z, id, time, lrv_par = B, ...)
    f10 <- fit(Kmax = Kmax, Ninit = Ninit)
    res <- c(
      sel_block(Z, id, time, f10, "_base"),
      sel_block(Z, id, time, f10, "_paths", conditioning = "paths"),
      sel_block(Z, id, time, f10, "_eucl", metric = "euclidean"),
      sel_block(Z, id, time, f10, "_shrink", metric = "shrink"),
      sel_block(Z, id, time, f10, "_F", ref = "F"),
      sel_block(Z, id, time, fit(Kmax = Kmax, Ninit = Ninit, init_type = "partition"), "_rpart"),
      sel_block(Z, id, time, fit(Kmax = Kmax, Ninit = 1), "_S1"),
      sel_block(Z, id, time, fit(Kmax = Kmax, Ninit = 100), "_S100"),
      sel_block(Z, id, time, fit(K = 3, Ninit = Ninit), "_K3"),
      sel_block(Z, id, time, fit(Kmax = Kmax, Ninit = Ninit, kappa = 3), "_kappa3"),
      unlist(lapply(c(0.5, 1, 2), function(cb) {
        Bc <- epa_default_B(Z, id, time, Kmax, cB = cb)
        fc <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = Bc)
        c(sel_block(Z, id, time, fc, paste0("_cB", cb)), setNames(Bc, paste0("B_cB", cb)))
      })),
      Khat = f10$K, B = B)
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- split: a true cluster split by Panel Kmeans (K0 = 2, K = 3) ----------------------------
rep_split <- function(sim, N, Tobs, a) {
  gamma <- rep(1:2, each = N / 2)
  dat <- generate_epa_simulation_data(N, Tobs, gamma, psi = c(-a, a), rho_vec = c(0.1, 0.3))
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    P <- ncol(Z)
    B <- epa_default_B(Z, id, time, Kmax = 3)
    fit <- panel_selective_fit(Z, id, time, K = 3, Ninit = 10, lrv_par = B)
    cl <- fit$selection$final_cluster
    pairs <- t(combn(3, 2))
    pure <- apply(pairs, 1, function(pr) length(unique(gamma[cl %in% pr])) == 1)
    res <- c(p_sel = NA, p_naive = NA, has_pure = any(pure))
    if (any(pure)) {                     # a pair of selected clusters within one true cluster
      pr <- pairs[which(pure)[1], ]
      L <- L_pairwise(fit, pr[1], pr[2])
      res["p_sel"] <- safe(selective_test(fit, L)$pval)
      G <- fit$Lambda %*% L
      y <- crossprod(L, fit$x)
      W <- as.numeric(crossprod(y, solve(crossprod(G) / B, y))) * (B - P + 1) / (P * B)
      res["p_naive"] <- pf(W, P, B - P + 1, lower.tail = FALSE)
    }
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- breaks: C-EPA holds on average over time ------------------------------------------------
rep_breaks <- function(sim, N, Tobs, d, Kmax, Ninit) {
  gamma <- true_clusters(N)
  psi1 <- psi_values(d, "overall_fails")
  dat <- generate_epa_simulation_data(N, Tobs, gamma, psi = psi1, psi_break = -psi1)
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    B <- epa_default_B(Z, id, time, Kmax)
    fit <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B)
    fitl <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B, demean = "local")
    s <- epa_clustered_selective(Z, id, time, fit = fit)
    sl <- epa_clustered_selective(Z, id, time, fit = fitl, demean = "local")
    res <- c(p_known = pv(safe(epa_clustered_known(Z, id, time, gamma, lrv_par = B))),
             p_known_local = pv(safe(epa_clustered_known(Z, id, time, gamma, lrv_par = B, demean = "local"))),
             p_split = pv(safe(epa_clustered_split(Z, id, time, Kmax = Kmax, prop = 0.2, Ninit = Ninit))),
             p_split_0.5 = pv(safe(epa_clustered_split(Z, id, time, Kmax = Kmax, prop = 0.5, Ninit = Ninit))),
             p_sel = s$pval, p_sel_joint = s$joint_cepa_pval,
             p_sel_local = sl$pval, p_sel_joint_local = sl$joint_cepa_pval,
             p_oepa = s$overall_pval, p_oepa_local = sl$overall_pval)
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- ic: estimation of the number of clusters -----------------------------------------------
rep_ic <- function(sim, N, Tobs, d, case, Kmax) {
  gamma <- true_clusters(N)
  dat <- generate_epa_simulation_data(N, Tobs, gamma, psi_values(d, case))
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    r15 <- panel_kmeans_estimation(Z, id, time, Kmax = Kmax, kappa = 1.5)$final_cluster
    r3 <- panel_kmeans_estimation(Z, id, time, Kmax = Kmax, kappa = 3)$final_cluster
    rcv <- panel_kmeans_cv_time(Z, id, time, Kmax = Kmax, nfold = 10)$final_cluster
    res <- c(setNames(cl_stats(r15, gamma), c("ari_ic15", "rec_ic15", "K_ic15")),
             setNames(cl_stats(r3, gamma), c("ari_ic3", "rec_ic3", "K_ic3")),
             setNames(cl_stats(rcv, gamma), c("ari_cv", "rec_cv", "K_cv")),
             correct_ic15 = length(unique(r15)) == 3, correct_ic3 = length(unique(r3)) == 3,
             correct_cv = length(unique(rcv)) == 3)
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- within: heterogeneity within the true clusters -----------------------------------------
# psi_i = theta_{k_i} + s u_i, with sum of u_i = 0 within each true cluster (fixed across
# replications), so that the true cluster means are theta_k. Variant "null": theta = 0, C-EPA
# holds for the true clusters. Variant "alternative": theta = psi_values(0.5, "overall_fails").
within_params <- function(N, s, variant) {
  gamma <- true_clusters(N)
  set.seed(4242)
  u <- rnorm(N)
  u <- u - ave(u, gamma)
  u <- u / sd(u)
  theta <- if (variant == "null") c(0, 0, 0) else psi_values(0.5, "overall_fails")
  list(psi = theta[gamma] + s * u, rho = c(0.1, 0.2, 0.3)[gamma])
}
rep_within <- function(sim, N, Tobs, s, variant, psi, rho, Kmax, Ninit) {
  gamma <- true_clusters(N)
  dat <- generate_epa_simulation_data(N, Tobs, gamma = 1:N, psi = psi, rho_vec = rho)
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    B <- epa_default_B(Z, id, time, Kmax)
    fit <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B)
    sel <- epa_clustered_selective(Z, id, time, fit = fit)
    res <- c(p_known = pv(safe(epa_clustered_known(Z, id, time, gamma, lrv_par = B))),
             p_naive = pv(safe(epa_clustered_known(Z, id, time, as.integer(fit$clustering), lrv_par = B))),
             p_split = pv(safe(epa_clustered_split(Z, id, time, Kmax = Kmax, prop = 0.2, Ninit = Ninit))),
             p_sel = sel$pval, p_sel_joint = sel$joint_cepa_pval,
             p_hom = sel$hom_test_pval, p_hom_joint = sel$joint_hom_pval,
             p_wbonf = sel$wbonf_pval, p_oepa = sel$overall_pval,
             cl_stats(fit$clustering, gamma), B = B)
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- jepa_het: all units heterogeneous ------------------------------------------------------
# psi_i = d z_i with z_i fixed across replications: size (d = 0); offsetting (O-EPA holds);
# shifted (O-EPA fails); sparse (four units with |z_i| = 4, O-EPA holds).
het_params <- function(N, variant) {
  set.seed(1000 + N)
  z <- rnorm(N)
  z <- switch(variant, size = 0 * z, offsetting = z - mean(z), shifted = z + 0.5,
              sparse = c(4, 4, -4, -4, rep(0, N - 4)))
  list(z = z, rho = runif(N, 0.1, 0.3))
}
rep_jepa_het <- function(sim, N, Tobs, d, variant, z, rho, Kmax, Ninit) {
  dat <- generate_epa_simulation_data(N, Tobs, gamma = 1:N, psi = d * z, rho_vec = rho)
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  delta <- log(log(Tobs)) * sqrt(log(N))            # screening threshold of Fan, Liao and Yao
  out <- c()
  for (cond in c(FALSE, TRUE)) {
    Z <- build_Z(dat, cond)
    B <- epa_default_B(Z, id, time, Kmax)
    fit <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B)
    sel <- epa_clustered_selective(Z, id, time, fit = fit)
    ut <- epa_unit_tests(Z, id, time, B = B, p_oepa = sel$overall_pval)
    res <- c(p_sel = sel$pval, p_sel_joint = sel$joint_cepa_pval,
             p_hom = sel$hom_test_pval, p_hom_joint = sel$joint_hom_pval,
             p_wbonf = sel$wbonf_pval, p_oepa = sel$overall_pval,
             p_split = pv(safe(epa_clustered_split(Z, id, time, Kmax = Kmax, prop = 0.2, Ninit = Ninit))),
             p_unit_maxt = ut$p_maxt, p_unit_delta = ut$p_delta, p_unit_delta_cc = ut$p_delta_cc,
             p_unit_choi = ut$p_choi, p_unit_cce_delta = ut$p_cce_delta, p_unit_cce_max = ut$p_cce_max,
             p_pe_sel_joint = power_enhance(sel$joint_cepa_pval, ut$p_maxt),
             p_pe_sel = power_enhance(sel$pval, ut$p_maxt),
             p_fly_sel_joint = if (any(sqrt(ut$unit_W) > delta)) 0 else sel$joint_cepa_pval,
             Khat = fit$K, B = B, n_screened = sum(sqrt(ut$unit_W) > delta))
    names(res) <- paste0(names(res), if (cond) "_cond" else "")
    out <- c(out, res)
  }
  out
}

# ---- eucl_fail: a design in which the Euclidean decomposition is not valid ------------------
# Loss differentials generated directly, Z_it = psi_{k_i} + lambda_i F_t + sigma_i e_it, with a
# common AR(1) factor F_t, strongly heterogeneous loadings lambda_i ~ U(0, 2), long-run variances
# that differ across units (sigma_i in {0.5, 2}) and unit-specific AR(1) errors with coefficients
# in [0, 0.6]. The pairwise differences of cluster means are then far from eigenvectors of Xi, so
# that the Euclidean residual is correlated with the tested difference. Compares the metric of Xi,
# the Euclidean decomposition and the shrinkage estimator of Xi L (size: d = 0).
eucl_params <- function(N) {
  set.seed(777)
  list(lambda = runif(N, 0, 2), sigma = sample(rep(c(0.5, 2), length.out = N)), phi = runif(N, 0, 0.6))
}
gen_eucl_fail <- function(N, Tobs, d, lambda, sigma, phi, burn = 100) {
  gamma <- true_clusters(N)
  psi <- psi_values(d, if (d == 0) "size" else "overall_holds")[gamma]
  TT <- Tobs + burn
  Fc <- as.numeric(arima.sim(list(ar = 0.5), TT)) * sqrt(1 - 0.25)
  E <- sapply(seq_len(N), function(i) as.numeric(arima.sim(list(ar = phi[i]), TT)) * sqrt(1 - phi[i]^2))
  Zm <- outer(Fc, lambda) + sweep(E, 2, sigma, "*")
  Zm <- Zm[(burn + 1):TT, , drop = FALSE]
  Zm <- sweep(Zm, 2, psi, "+")
  matrix(as.vector(Zm), ncol = 1)                       # unit-major (id = rep(1:N, each = Tobs))
}
rep_eucl_fail <- function(sim, N, Tobs, d, lambda, sigma, phi, Kmax, Ninit) {
  Z <- gen_eucl_fail(N, Tobs, d, lambda, sigma, phi)
  id <- rep(1:N, each = Tobs); time <- rep(1:Tobs, N)
  B <- epa_default_B(Z, id, time, Kmax)
  fit <- panel_selective_fit(Z, id, time, Kmax = Kmax, Ninit = Ninit, lrv_par = B)
  c(sel_block(Z, id, time, fit, "_xi"),
    sel_block(Z, id, time, fit, "_eucl", metric = "euclidean"),
    sel_block(Z, id, time, fit, "_shrink", metric = "shrink"),
    p_oepa = pv(safe(overall_EPA_test(Z = Z, id = id, time = time, lrv = "EWC", lrv_par = B))),
    Khat = fit$K, B = B)
}

# =============================================================================
# Designs and runs
# =============================================================================
g <- function(...) expand.grid(..., stringsAsFactors = FALSE)

if ("main" %in% experiments) {
  run_experiment("main", rep_main,
                 bind_rows(g(N = N, Tobs = c(20, 50, 100, 200), d = 0, case = "size"),
                           g(N = N, Tobs = c(50, 200), d = c(0.125, 0.25, 0.375, 0.5),
                             case = c("overall_fails", "overall_holds"))),
                 list(Kmax = Kmax, Ninit = Ninit))
}
if ("variants" %in% experiments) {
  run_experiment("variants", rep_variants,
                 bind_rows(g(N = N, Tobs = c(50, 200), d = 0, case = "size"),
                           g(N = N, Tobs = c(50, 200), d = c(0.25, 0.5),
                             case = c("overall_fails", "overall_holds"))),
                 list(Kmax = Kmax, Ninit = Ninit))
}
if ("split" %in% experiments) {
  run_experiment("split", rep_split, g(N = N, Tobs = c(50, 200), a = c(0.25, 0.5)))
}
if ("breaks" %in% experiments) {
  run_experiment("breaks", rep_breaks, g(N = N, Tobs = c(20, 50, 100, 200), d = 0.25),
                 list(Kmax = Kmax, Ninit = Ninit))
}
if ("ic" %in% experiments) {
  run_experiment("ic", rep_ic,
                 bind_rows(g(N = N, Tobs = c(20, 50, 100, 200), d = 0, case = "size"),
                           g(N = N, Tobs = c(50, 200), d = c(0.125, 0.25, 0.375, 0.5),
                             case = c("overall_fails", "overall_holds"))),
                 list(Kmax = Kmax))
}
if ("within" %in% experiments) {
  for (variant in c("null", "alternative")) for (s in c(0, 0.125, 0.25, 0.5)) {
    wp <- within_params(N, s, variant)
    run_experiment("within", rep_within, g(N = N, Tobs = c(50, 200), s = s, variant = variant),
                   list(psi = wp$psi, rho = wp$rho, Kmax = Kmax, Ninit = Ninit))
  }
}
if ("jepa_het" %in% experiments) {
  for (variant in c("size", "offsetting", "shifted", "sparse")) {
    hp <- het_params(N, variant)
    dd <- if (variant == "size") 0 else c(0.125, 0.25, 0.5)
    run_experiment("jepa_het", rep_jepa_het, g(N = N, Tobs = c(50, 200), d = dd, variant = variant),
                   list(z = hp$z, rho = hp$rho, Kmax = Kmax, Ninit = Ninit))
  }
}

if ("eucl_fail" %in% experiments) {
  ep <- eucl_params(N)
  run_experiment("eucl_fail", rep_eucl_fail,
                 bind_rows(g(N = N, Tobs = c(20, 50, 100, 200), d = 0),
                           g(N = N, Tobs = c(50, 200), d = c(0.25, 0.5))),
                 list(lambda = ep$lambda, sigma = ep$sigma, phi = ep$phi, Kmax = Kmax, Ninit = Ninit))
}

# ---- figures of the break experiment (Figure 1 of the paper) --------------------------------
# (written to Simulation/Results; copy them to the Figures folder of the paper)
f_breaks <- file.path(res_dir, "results_breaks.xlsx")
if (file.exists(f_breaks) && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  library(tidyr)
  res <- read.xlsx(f_breaks)
  make_plot <- function(df, sfx) {
    cols <- paste0(c("p_known", "p_split", "p_split_0.5", "p_sel", "p_sel_local"), sfx)
    long <- df %>% select(Tobs, all_of(cols)) %>%
      pivot_longer(-Tobs, names_to = "Test", values_to = "RejectionRate") %>%
      mutate(Test = factor(if (nzchar(sfx)) sub(sfx, "", Test, fixed = TRUE) else Test,
                           levels = c("p_known", "p_split", "p_split_0.5", "p_sel", "p_sel_local"),
                           labels = c("Predetermined", "Split Sample (0.2)", "Split Sample (0.5)",
                                      "Selective Inference", "Selective Inference (local demeaning)")))
    ggplot(long, aes(x = Tobs, y = RejectionRate, color = Test, shape = Test)) +
      geom_line(linewidth = 1) + geom_point(size = 3) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
      labs(x = "T", y = "Rejection rate") +
      ylim(0, max(0.4, max(long$RejectionRate, na.rm = TRUE))) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      guides(color = guide_legend(nrow = 2))
  }
  ggsave(file.path(res_dir, "unconditional_CEPA_size_sbreaks.pdf"), make_plot(res, ""), width = 5.5, height = 4.5)
  ggsave(file.path(res_dir, "conditional_CEPA_size_sbreaks.pdf"), make_plot(res, "_cond"), width = 5.5, height = 4.5)
}
cat("Done. Results in", res_dir, "\n")
