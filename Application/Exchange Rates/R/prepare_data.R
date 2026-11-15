rm(list=ls())
# ==============================================================================
# Cleaning and summarizing data
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(openxlsx)
library(tidyr)
library(ggrepel)

# Read data
forecast_panel_final <- readRDS(here("Application", "Exchange Rates", "Data", "forecast_panel_final.rds"))

# Rename some variables
forecast_panel_final <- forecast_panel_final %>%
  rename(ElNet = LASSO,
         XGBoost_X = XGB,
         SVM_X = SVMX,
         RF_X = RF,
         Loss_ElNet = Loss_LASSO,
         Loss_XGBoost_X = Loss_XGB,
         Loss_SVM_X = Loss_SVMX,
         Loss_RF_X = Loss_RF,
         d_ElNet = d_LASSO,
         d_XGBoost_X = d_XGB,
         d_SVM_X = d_SVMX,
         d_RF_X = d_RF)

# ==============================================================================
# Outliers, NAs
# ==============================================================================
# Check extremes of AR(1) loss
boxplot(forecast_panel_final$Loss_AR1, main = "Boxplot of AR(1) loss")
forecast_panel_final %>% filter(Loss_AR1 > 600)

# Identify problematic assets
bad_assets <- forecast_panel_final %>%
  filter(Loss_AR1 > 600) %>%
  pull(Asset) %>%
  unique()

# Drop all rows corresponding to these assets
forecast_panel_cleaned <- forecast_panel_final %>%
  filter(!Asset %in% bad_assets)

# Check summaries
forecast_panel_cleaned %>%
  summarise(across(starts_with("d_"), ~ mean(.x, na.rm = TRUE)))
forecast_panel_cleaned %>%
  summarise(across(starts_with("d_"), ~ sum(is.na(.x))))

# Exclude some forecasts
# Elastic net with macro predictors (too many extremes) and Markov-switching AR (too many NAs);
# no longer computed by generate_forecasts.R, dropped here if present in an older file
forecast_panel_cleaned <- forecast_panel_cleaned %>%
  select(-any_of(c("LASSOX", "Loss_LASSOX", "d_LASSOX", "MSAR", "Loss_MSAR", "d_MSAR")))

# Replace 3 NAs of Lasso with RW
forecast_panel_cleaned <- forecast_panel_cleaned %>%
  mutate(ElNet = ifelse(is.na(ElNet), Actual_Lag1, ElNet))
forecast_panel_cleaned$Loss_ElNet <- (forecast_panel_cleaned$ElNet - forecast_panel_cleaned$Actual)^2
# Same convention as generate_forecasts.R: d_model = Loss_AR1 - Loss_model (positive: the model
# improves on the AR(1) benchmark). The previous line used Loss_ElNet - Loss_AR1, which flipped
# the sign of d_ElNet only (Wald tests are invariant to the sign; summary statistics are not).
forecast_panel_cleaned$d_ElNet <- forecast_panel_cleaned$Loss_AR1 - forecast_panel_cleaned$Loss_ElNet

# Summarize once more
forecast_panel_cleaned %>%
  summarise(across(starts_with("d_"), ~ mean(.x, na.rm = TRUE)))

# ==============================================================================
# Extreme forecasts and realized values
# ==============================================================================
# Two sources of extreme loss differentials, fixed in advance:
# (i) explosive forecasts: some estimated models (mostly AR(p) with p chosen by BIC) forecast
#     more than 5 window standard deviations, e.g. -13.6 for Bangladesh in 2012, giving single
#     loss differentials of about -200. As in the "insanity filter" of Swanson and White (1997),
#     every forecast is truncated at +-cut_forecast (in the standardized units of the window);
# (ii) devaluations and regime changes: a few currencies have realized returns of many standard
#     deviations (e.g. the floating of the Egyptian pound in 2016), which produce singleton
#     clusters. Treatments:
#   "winsorize": truncate the forecasts (i), winsorize the realized returns at +-cut_forecast;
#               losses, loss differentials and the lagged target are recomputed;
#   "exclude":   truncate the forecasts (i) and drop every currency whose standardized realized
#               return exceeds cut_exclude in absolute value at some date;
#   "none":      no treatment.
# (iii) duplicated series: currencies in a monetary union or with a hard peg to the same anchor
#     have (almost) identical returns against the dollar. Identical series make the long-run
#     covariance matrix Xi singular and tie in Panel Kmeans, so only one currency of each group
#     (absolute correlation of the standardized returns above cut_dedup) is kept, the anchor when
#     it is in the sample (euro area, United Kingdom, Australia, South Africa, India, Singapore),
#     otherwise the first in alphabetical order.
# Panels (run_application.R chooses one with `sample_rule`):
#   forecast_panel_ready.rds          "main":    winsorize, duplicates dropped
#   forecast_panel_ready_all.rds      "all":     winsorize, all currencies
#   forecast_panel_ready_exclude.rds  "exclude": exclude, duplicates dropped
#   forecast_panel_ready_none.rds     "none":    no treatment, duplicates dropped
#   forecast_panel_ready_raw.rds      no treatment, all currencies (input of robustness_samples.R)
# The summary statistics and the figures are computed on the main panel.
cut_forecast <- 5
cut_exclude  <- 8
cut_dedup    <- 0.98
anchors      <- c("G163", "GBR", "AUS", "ZAF", "IND", "SGP")   # kept first within their groups

# Currencies kept after dropping duplicated series (see (iii) above)
dedup_assets <- function(dat, thr = cut_dedup, first = anchors) {
  W <- dat %>% select(Date, Asset, Actual) %>% pivot_wider(names_from = Asset, values_from = Actual)
  C <- abs(cor(as.matrix(W[, -1]), use = "pairwise.complete.obs"))
  code <- sub("\\..*", "", colnames(C))
  ord <- order(!(code %in% first), match(code, first), colnames(C))
  keep <- character(0)
  for (a in colnames(C)[ord]) if (!length(keep) || all(C[a, keep] <= thr)) keep <- c(keep, a)
  sort(keep)
}

# Recompute Loss_m = (m - Actual)^2 and d_m = Loss_AR1 - Loss_m for every model with a forecast column
recompute_losses <- function(dat) {
  models <- sub("^Loss_", "", grep("^Loss_", names(dat), value = TRUE))
  models <- models[models %in% names(dat)]
  for (m in models) dat[[paste0("Loss_", m)]] <- (dat[[m]] - dat$Actual)^2
  for (m in setdiff(models, "AR1")) dat[[paste0("d_", m)]] <- dat$Loss_AR1 - dat[[paste0("Loss_", m)]]
  dat
}
truncate_forecasts <- function(dat, cut) {
  models <- sub("^Loss_", "", grep("^Loss_", names(dat), value = TRUE))
  models <- models[models %in% names(dat)]
  n_cut <- sum(sapply(models, function(m) sum(abs(dat[[m]]) > cut, na.rm = TRUE)))
  message(sprintf("Forecasts truncated at +-%g: %d of %d", cut, n_cut, length(models) * nrow(dat)))
  for (m in models) dat[[m]] <- pmin(pmax(dat[[m]], -cut), cut)
  recompute_losses(dat)
}

panel_none <- forecast_panel_cleaned
forecast_panel_cleaned <- truncate_forecasts(forecast_panel_cleaned, cut_forecast)

max_abs <- tapply(abs(forecast_panel_cleaned$Actual), forecast_panel_cleaned$Asset, max, na.rm = TRUE)
extreme_assets <- names(max_abs)[max_abs > cut_exclude]
message(sprintf("Currencies with |standardized return| > %g at some date (%d of %d): %s",
                cut_exclude, length(extreme_assets), length(max_abs),
                paste(sub("\\..*", "", extreme_assets), collapse = ", ")))

panel_winsorize <- forecast_panel_cleaned %>%
  mutate(Actual = pmin(pmax(Actual, -cut_forecast), cut_forecast)) %>%
  recompute_losses() %>%
  arrange(Asset, Date) %>% group_by(Asset) %>% mutate(Actual_Lag1 = lag(Actual)) %>% ungroup()
keep_assets <- dedup_assets(panel_none)
message(sprintf("Duplicated series dropped (%d of %d): %s", length(unique(panel_none$Asset)) - length(keep_assets),
                length(unique(panel_none$Asset)),
                paste(sub("\\..*", "", setdiff(unique(panel_none$Asset), keep_assets)), collapse = ", ")))
panel_raw <- panel_none
panel_all <- panel_winsorize
panel_main <- panel_winsorize %>% filter(Asset %in% keep_assets)
panel_exclude <- forecast_panel_cleaned %>% filter(!Asset %in% extreme_assets, Asset %in% keep_assets)
panel_none <- panel_none %>% filter(Asset %in% keep_assets)
forecast_panel_cleaned <- panel_main             # main sample
message(sprintf("Main sample: %d currencies; all: %d; exclude rule: %d currencies",
                length(unique(panel_main$Asset)), length(unique(panel_all$Asset)),
                length(unique(panel_exclude$Asset))))

# ==============================================================================
# Save summary stats
# ==============================================================================
# Select variables for summary statistics
vars <- c("d_AR_BIC", "d_ElNet", "d_XGBoost_X", "d_SVM_X", "d_RF_X",
          "d_ZERO", "d_AR1_JS", "d_POOL", "d_EB", "d_GRP", "d_FAR", "d_GBM", "d_COMB")
# The last eight exist only if the forecasts were regenerated with the current generate_forecasts.R
missing <- setdiff(vars, names(forecast_panel_cleaned))
if (length(missing)) message("Not in the forecast file (regenerate the forecasts to add them): ",
                             paste(missing, collapse = ", "))
vars <- intersect(vars, names(forecast_panel_cleaned))

# Compute summary statistics in tidy form
# Reported in the convention of the paper (Table S.1): Loss_model - Loss_AR1, so that negative
# values indicate an improvement over the AR(1) benchmark. The d_ variables use the opposite sign.
summary_table <- forecast_panel_cleaned |>
  select(all_of(vars)) |>
  mutate(across(everything(), ~ -.x)) |>
  pivot_longer(cols = everything(), names_to = "Variable") |>
  group_by(Variable) |>
  summarise(
    Mean = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    N = sum(!is.na(value)),
    .groups = "drop"
  )

# Export to Excel
openxlsx::write.xlsx(summary_table, here::here("Application", "Exchange Rates", "Tables", "summary_stats.xlsx"))

# ==============================================================================
# Plots
# ==============================================================================
# Create folder if missing
fig_dir <- here("Application", "Exchange Rates", "Figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Variables
d_vars <- c("Loss_AR_BIC", "Loss_ElNet", "Loss_XGBoost_X", "Loss_SVM_X", "Loss_RF_X",
            "Loss_ZERO", "Loss_AR1_JS", "Loss_POOL", "Loss_EB", "Loss_GRP",
            "Loss_FAR", "Loss_GBM", "Loss_COMB")
d_vars <- intersect(d_vars, names(forecast_panel_cleaned))
method_labels <- c(Loss_AR_BIC = "AR(p)", Loss_ElNet = "Elastic Net", Loss_XGBoost_X = "XGBoost",
                   Loss_SVM_X = "SVM", Loss_RF_X = "Random forest", Loss_ZERO = "No change",
                   Loss_AR1_JS = "Shrunk AR(1)", Loss_POOL = "Pooled AR(1)",
                   Loss_EB = "Empirical Bayes AR(1)", Loss_GRP = "Grouped AR(1)",
                   Loss_FAR = "Factor ridge", Loss_GBM = "Global boosting", Loss_COMB = "Combination")
loss_breaks <- c(0, 1, 5, 25, 100)   # losses are at most (5 + 5)^2 = 100 after truncation

# Loop through each model
for (var_name in d_vars) {
  
  df_plot <- forecast_panel_cleaned
  lab <- method_labels[[var_name]]

  p <- ggplot(df_plot, aes(x = Loss_AR1, y = .data[[var_name]])) +
    geom_point(alpha = 0.35, size = 1.6, color = "#1f77b4", stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    scale_x_continuous(trans = "log1p", breaks = loss_breaks, limits = c(0, 100)) +
    scale_y_continuous(trans = "log1p", breaks = loss_breaks, limits = c(0, 100)) +
    coord_equal() +
    labs(x = "Loss of AR(1)", y = paste0("Loss of ", lab)) +
    theme_classic(base_size = 13) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      axis.line = element_line(color = "black")
    )
  
  # Save PDF
  ggsave(filename = file.path(fig_dir, paste0("scatter_log_Loss_AR1_vs_", var_name, ".pdf")),
         plot = p, width = 4.8, height = 4.8, device = cairo_pdf)
  
  # Save EPS
  ggsave(filename = file.path(fig_dir, paste0("scatter_log_Loss_AR1_vs_", var_name, ".eps")),
         plot = p, width = 4.8, height = 4.8, device = cairo_ps)
}

# ==============================================================================
# Save the cleaned file
# ==============================================================================
saveRDS(panel_main, file = here("Application", "Exchange Rates", "Data", "forecast_panel_ready.rds"))
saveRDS(panel_all, file = here("Application", "Exchange Rates", "Data", "forecast_panel_ready_all.rds"))
saveRDS(panel_raw, file = here("Application", "Exchange Rates", "Data", "forecast_panel_ready_raw.rds"))
saveRDS(panel_exclude, file = here("Application", "Exchange Rates", "Data", "forecast_panel_ready_exclude.rds"))
saveRDS(panel_none, file = here("Application", "Exchange Rates", "Data", "forecast_panel_ready_none.rds"))