# =============================================================================
# Robustness of the application to the sample and the loss function.
#
# Runs the selective tests (K by the IC) for a grid of sample definitions, forecasting methods
# and testing functions, and reports for each case the O-EPA, naive, selective C-EPA (merged and
# joint) and selective homogeneity (merged and joint) p-values, together with the cluster
# sizes, the cluster means of the loss differential and the members of the small clusters.
#
# Sample definitions (all start from Data/forecast_panel_ready_raw.rds written by
# R/prepare_data.R, i.e. all 122 currencies without treatment; forecasts are truncated at +-5,
# and, except in "all_winsor5", only one currency of each group of duplicated series is kept,
# with the rule of R/prepare_data.R):
#   winsor5         main results: realized returns winsorized at +-5
#   exclude8/12/20  currencies with |standardized realized return| > 8, 12 or 20 at some date
#                   dropped instead of winsorizing
#   winsor3         realized returns winsorized at +-3
#   none            no treatment of the realized returns
#   all_winsor5     winsor5 with all currencies, including the duplicated series
#   post2010        winsor5 restricted to the forecasts from January 2010
#   abs_winsor5     winsor5 with the absolute loss instead of the quadratic loss
#
# Output: Tables/robustness_samples.xlsx (sheet "tests": one row per case; sheet "clusters":
# cluster sizes, means and members). Each finished case is saved, so the script can be
# interrupted and restarted.
# =============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(openxlsx); library(zoo)
  library(PanelKmeansInference); library(clusteredEPA)
})

app_dir <- here("Application", "Exchange Rates")
out_file <- file.path(app_dir, "Tables", "robustness_samples.xlsx")
ckpt_file <- file.path(app_dir, "Tables", "robustness_samples_checkpoint_v2.rds")

# ---- settings -------------------------------------------------------------------------------
samples <- c("winsor5", "exclude8", "exclude12", "exclude20", "winsor3", "none",
             "all_winsor5", "post2010", "abs_winsor5")
methods <- c(Zero = "ZERO", AR1_JS = "AR1_JS", Pooled = "POOL", EB = "EB", Grouped = "GRP",
             AR_BIC = "AR_BIC", ElNet = "ElNet", XGBoost = "XGBoost_X", SVM = "SVM_X", RF = "RF_X",
             FAR = "FAR", GlobalGBM = "GBM", Combination = "COMB")   # forecast columns
conds <- list(Unconditional = NULL, Lagged = "Actual_Lag1")
Kmax <- 5; Ninit <- 100; iter.max <- 100; kappa <- 1.5
cut_forecast <- 5
set.seed(1)
invisible(PanelKmeansInference:::.get_cpp_walk())

# ---- data -----------------------------------------------------------------------------------
base <- readRDS(file.path(app_dir, "Data", "forecast_panel_ready_raw.rds"))
fc_cols <- intersect(c("AR1", methods), names(base))
missing <- setdiff(methods, fc_cols)
if (length(missing)) {
  warning("Forecasts not in the data (regenerate them): ", paste(missing, collapse = ", "), immediate. = TRUE)
  methods <- methods[methods %in% fc_cols]
}
base <- base %>% arrange(Asset, Date)
for (m in fc_cols) base[[m]] <- pmin(pmax(base[[m]], -cut_forecast), cut_forecast)

# Loss differentials d = Loss_AR1 - Loss_model (sign irrelevant for the tests) and lagged target
make_d <- function(dat, loss = c("quad", "abs")) {
  loss <- match.arg(loss)
  L <- if (loss == "quad") function(e) e^2 else abs
  for (m in setdiff(fc_cols, "AR1")) dat[[paste0("d_", m)]] <- L(dat$AR1 - dat$Actual) - L(dat[[m]] - dat$Actual)
  dat %>% group_by(Asset) %>% mutate(Actual_Lag1 = lag(Actual)) %>% ungroup()
}
max_abs <- tapply(abs(base$Actual), base$Asset, max, na.rm = TRUE)
exclude <- function(dat, cut) dat[!dat$Asset %in% names(max_abs)[max_abs > cut], ]
winsor <- function(dat, cut) mutate(dat, Actual = pmin(pmax(Actual, -cut), cut))
# One currency of each group of duplicated series, as in R/prepare_data.R (anchor first)
anchors <- c("G163", "GBR", "AUS", "ZAF", "IND", "SGP")
dedup <- function(dat, thr = 0.98) {
  W <- base %>% select(Date, Asset, Actual) %>% pivot_wider(names_from = Asset, values_from = Actual)
  C <- abs(cor(as.matrix(W[, -1]), use = "pairwise.complete.obs"))
  code <- sub("\\..*", "", colnames(C))
  ord <- order(!(code %in% anchors), match(code, anchors), colnames(C))
  keep <- character(0)
  for (a in colnames(C)[ord]) if (!length(keep) || all(C[a, keep] <= thr)) keep <- c(keep, a)
  dat[dat$Asset %in% keep, ]
}
build_sample <- function(s) switch(s,
  winsor5       = make_d(dedup(winsor(base, 5))),
  exclude8      = make_d(dedup(exclude(base, 8))),
  exclude12     = make_d(dedup(exclude(base, 12))),
  exclude20     = make_d(dedup(exclude(base, 20))),
  winsor3       = make_d(dedup(winsor(base, 3))),
  none          = make_d(dedup(base)),
  all_winsor5   = make_d(winsor(base, 5)),
  post2010      = make_d(dedup(winsor(base, 5))) %>% filter(Date >= as.yearmon("2010-01")),
  abs_winsor5   = make_d(dedup(winsor(base, 5)), loss = "abs"))

# ---- tests for one case ---------------------------------------------------------------------
safe <- function(expr) tryCatch(expr, error = function(e) { message("  error: ", conditionMessage(e)); NULL })
run_case <- function(dat, s, m, cn) {
  v <- paste0("d_", methods[[m]]); h <- conds[[cn]]
  dd <- if (identical(h, "Actual_Lag1")) dat[!is.na(dat$Actual_Lag1), ] else dat
  dd <- dd[!is.na(dd[[v]]), ]
  dd$time <- as.integer(factor(dd$Date)); dd$id <- as.integer(factor(dd$Asset))
  dd <- dd[order(dd$id, dd$time), ]
  Z <- if (is.null(h)) matrix(dd[[v]], ncol = 1) else cbind(dd[[v]], dd[[v]] * dd[[h]])
  B <- epa_default_B(Z, dd$id, dd$time, Kmax)
  fit <- panel_selective_fit(Z, dd$id, dd$time, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                             kappa = kappa, lrv_par = B)
  sel <- safe(epa_clustered_selective(Z, dd$id, dd$time, fit = fit))
  nai <- safe(epa_clustered_known(Z, dd$id, dd$time, as.integer(fit$clustering), lrv_par = B))
  assets <- unique(dd$Asset[order(dd$id)])
  cl <- data.frame(Asset = assets, cluster = fit$clustering)
  # cluster means of the loss differential in the paper convention (Loss_model - Loss_AR1)
  cm <- dd %>% mutate(cluster = cl$cluster[match(Asset, cl$Asset)]) %>% group_by(cluster) %>%
    summarise(n = n_distinct(Asset), mean = -mean(.data[[v]]), .groups = "drop") %>%
    mutate(members = sapply(cluster, function(k) {
      a <- sub("\\..*", "", cl$Asset[cl$cluster == k]); if (length(a) <= 12) paste(a, collapse = " ") else "" }))
  row <- data.frame(sample = s, method = m, test_function = cn, N = length(assets),
                    T = length(unique(dd$time)), K = fit$K, B = B,
                    sizes = paste(sort(table(fit$clustering), decreasing = TRUE), collapse = "/"),
                    p_oepa = if (is.null(sel)) NA else sel$overall_pval,
                    p_naive = if (is.null(nai)) NA else nai$pval,
                    p_sel = if (is.null(sel)) NA else sel$pval,
                    p_sel_joint = if (is.null(sel)) NA else sel$joint_cepa_pval,
                    p_hom = if (is.null(sel)) NA else sel$hom_test_pval,
                    p_hom_joint = if (is.null(sel)) NA else sel$joint_hom_pval,
                    p_pair_min = if (is.null(sel)) NA else min(sel$pairwise_pvalues))
  list(row = row, clusters = data.frame(sample = s, method = m, test_function = cn, cm))
}

# ---- loop with checkpoint -------------------------------------------------------------------
res <- if (file.exists(ckpt_file)) readRDS(ckpt_file) else list()
for (s in samples) {
  dat <- NULL
  for (m in names(methods)) for (cn in names(conds)) {
    key <- paste(s, m, cn)
    if (!is.null(res[[key]])) next
    if (is.null(dat)) { dat <- build_sample(s); message(sprintf("sample %s: N = %d", s, length(unique(dat$Asset)))) }
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), key))
    res[[key]] <- safe(run_case(dat, s, m, cn))
    saveRDS(res, ckpt_file)
  }
}

tests <- bind_rows(lapply(res, `[[`, "row"))
clus <- bind_rows(lapply(res, `[[`, "clusters"))
# Holm across methods within sample and testing function
tests <- tests %>% group_by(sample, test_function) %>%
  mutate(across(c(p_oepa, p_sel, p_sel_joint, p_hom, p_hom_joint), ~ p.adjust(.x, "holm"), .names = "{.col}_holm")) %>%
  ungroup() %>% as.data.frame()
write.xlsx(list(tests = tests, clusters = clus), out_file)

# Overview: number of cases with a selective homogeneity rejection at 5%, by sample
print(tests %>% group_by(sample) %>%
        summarise(N = first(N), hom = sum(p_hom < 0.05, na.rm = TRUE),
                  hom_joint = sum(p_hom_joint < 0.05, na.rm = TRUE),
                  hom_holm = sum(p_hom_holm < 0.05, na.rm = TRUE),
                  sel = sum(p_sel < 0.05, na.rm = TRUE), naive = sum(p_naive < 0.05, na.rm = TRUE),
                  .groups = "drop") %>% as.data.frame())
cat("Done:", out_file, "\n")
