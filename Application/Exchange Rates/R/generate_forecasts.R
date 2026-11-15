rm(list = ls())
# ============================
# Forecasting Exchange Rate Returns
# ============================
library(readxl)
library(dplyr)
library(tidyr)
library(zoo)
library(glmnet)
library(xgboost)
library(here)
library(fbi)
library(foreach)

source(here("Application", "Exchange Rates", "R", "tools.r"))

# === Model Selection Inputs ===
model_flags <- list(
  AR1 = TRUE, AR2 = FALSE, AR_BIC = TRUE, RW = FALSE, KF = FALSE,
  LASSO = TRUE, XGB = TRUE, TH = FALSE, TH_SOFT = FALSE,
  LASSOX = FALSE, SVMX = TRUE, RF = TRUE, MSAR = FALSE,  # LASSOX and MSAR are dropped in
                                                         # prepare_data.R (extremes, NAs): not computed
  ZERO = TRUE, AR1_JS = TRUE,                 # no-change and shrunk AR(1) forecasts
  POOL = TRUE, EB = TRUE, GRP = TRUE,         # panel AR(1) forecasts (computed after the unit loop)
  FAR = TRUE, GBM = TRUE,                     # factor-augmented ridge, global boosting (idem)
  COMB = TRUE                                 # equal-weight combination of the ML forecasts
)
panel_models <- c("POOL", "EB", "GRP")
factor_models <- c("FAR", "GBM")
comb_members <- c("AR_BIC", "LASSO", "XGB", "SVMX", "RF")

base_model <- "AR1"
checkpoint_path <- here("Application", "Exchange Rates", "R", "forecast_checkpoint_v3.rds")

# === Load Fred-MD data ===
fred_md <- fredmd(here("Application", "Exchange Rates", "Data", "fred_md.csv"), transform = TRUE)
X = fred_md[1:nrow(fred_md),2:ncol(fred_md)]
X = tw_apc(X,kmax=8)$data
fred_md[1:nrow(fred_md),2:ncol(fred_md)] = X
fred_md <- data.frame(fred_md) %>%
  mutate(Date = as.yearmon(date)) %>%
  select(-date) %>%
  filter(Date >= as.yearmon("1999-02") & Date <= as.yearmon("2023-12")) %>%
  select(-c("EXSZUSx", "EXJPUSx", "EXUSUKx", "EXCAUSx"))
fred_md <- fred_md[, c(ncol(fred_md), 1:(ncol(fred_md) - 1))]

# === Load FX data ===
fx_df <- read_excel(here("Application", "Exchange Rates", "Data", "imf_fx.xlsx")) %>%
  mutate(Date = as.yearmon(Date, format = "%Y-M%m")) %>%
  filter(Date >= as.yearmon("1999-01") & Date <= as.yearmon("2023-12"))

active_df <- fx_df
active_df[,-1] <- lapply(active_df[,-1], function(x) c(NA, diff(log(x))))
active_df <- active_df[-1, ]
active_df <- active_df[, colSums(is.na(active_df)) == 0]
active_df[,-1] <- 100 * active_df[,-1]

# === Setup ===
# Rolling window of r = 60 months: each forecast uses the last r observations only. Within each
# window the returns are demeaned and scaled by the window's own mean and standard deviation, and
# forecasts, actual values and losses are expressed in these units, so that losses are comparable
# across currencies without using any information after the forecast date.
tau <- 1
r <- 60
Tmax <- nrow(active_df)
p <- Tmax - r - 1
dates <- active_df$Date
ret_mat <- as.matrix(active_df[,-1])
n_assets <- ncol(ret_mat)

# Filter invalid series
invalid_assets <- rep(FALSE, n_assets)
for (nn in 1:n_assets) {
  for (jj in tau:p) {
    if (invalid_assets[nn]) break
    y <- ret_mat[(jj+1):(r+jj), nn]
    if (!is_valid_series(y)) {
      invalid_assets[nn] <- TRUE
      break
    }
  }
}
ret_mat <- ret_mat[, !invalid_assets]
# (no full-sample standardization: it would use information after the forecast dates)
n_assets <- ncol(ret_mat)
forecast_dates <- dates[(r + tau + 1):Tmax]

# === Load checkpoint if available ===
if (file.exists(checkpoint_path)) {
  cat("Resuming from checkpoint...\n")
  checkpoint <- readRDS(checkpoint_path)
  forecast_matrices <- checkpoint$forecast_matrices
  actual_std <- checkpoint$actual_std
  start_asset <- checkpoint$last_asset + 1
} else {
  forecast_matrices <- list()
  for (model in names(model_flags)) {
    if (model_flags[[model]]) {
      forecast_matrices[[model]] <- matrix(NA, nrow = p, ncol = n_assets)
      colnames(forecast_matrices[[model]]) <- colnames(ret_mat)
    }
  }
  actual_std <- matrix(NA, nrow = p, ncol = n_assets, dimnames = list(NULL, colnames(ret_mat)))
  start_asset <- 1
}

# === Forecast Loop with checkpointing ===
# The checkpoint stores the unit-by-unit forecasts; delete R/forecast_checkpoint_v3.rds after
# changing a forecasting method, otherwise the stored forecasts are reused.
for (nn in seq_len(n_assets)[seq_len(n_assets) >= start_asset]) {   # empty if the checkpoint is complete
  cat(sprintf("Running asset nn = %d/%d\n", nn, n_assets))
  for (jj in tau:p) {
    yr <- ret_mat[(jj+1):(r+jj), nn]
    m_w <- mean(yr); s_w <- sd(yr)
    y <- (yr - m_w) / s_w                                    # standardized window
    actual_std[jj, nn] <- (ret_mat[r+jj+1, nn] - m_w) / s_w  # target in the same units
    X_macro <- fred_md[(jj+1):(r+jj), -1]
    
    if (model_flags$AR1)      forecast_matrices$AR1[jj, nn]      <- forecast_ar1(y)
    if (model_flags$AR2)      forecast_matrices$AR2[jj, nn]      <- forecast_ar2(y)
    if (model_flags$RW)       forecast_matrices$RW[jj, nn]       <- forecast_rw(y)
    if (model_flags$KF)       forecast_matrices$KF[jj, nn]       <- forecast_kf(y)
    if (model_flags$AR_BIC)   forecast_matrices$AR_BIC[jj, nn]   <- forecast_ar_bic(y)
    if (model_flags$LASSO)    forecast_matrices$LASSO[jj, nn]    <- forecast_lasso(y)
    if (model_flags$XGB)      forecast_matrices$XGB[jj, nn]      <- forecast_xgboost_macro(y, as.matrix(X_macro))
    if (model_flags$TH)       forecast_matrices$TH[jj, nn]       <- forecast_threshold_hard(y, as.matrix(X_macro))
    if (model_flags$TH_SOFT)  forecast_matrices$TH_SOFT[jj, nn]  <- forecast_threshold_soft(y, as.matrix(X_macro))
    if (model_flags$LASSOX)   forecast_matrices$LASSOX[jj, nn]   <- forecast_lasso_macro(y, as.matrix(X_macro))
    if (model_flags$SVMX)     forecast_matrices$SVMX[jj, nn]     <- forecast_svr_macro(y, as.matrix(X_macro))
    if (model_flags$RF)       forecast_matrices$RF[jj, nn]       <- forecast_rf(y, as.matrix(X_macro))
    if (model_flags$MSAR)     forecast_matrices$MSAR[jj, nn]     <- forecast_msar(y)
    if (model_flags$ZERO)     forecast_matrices$ZERO[jj, nn]     <- forecast_zero(m_w, s_w)
    if (model_flags$AR1_JS)   forecast_matrices$AR1_JS[jj, nn]   <- forecast_ar1_js(y)
  }

  # Save checkpoint
  saveRDS(list(forecast_matrices = forecast_matrices, actual_std = actual_std, last_asset = nn),
          file = checkpoint_path)
}

# === Panel AR(1) forecasts: pooled, empirical Bayes, grouped (need all units at each date) ===
if (any(unlist(model_flags[panel_models]))) {
  for (jj in tau:p) {
    U <- apply(ret_mat[(jj+1):(r+jj), , drop = FALSE], 2, function(yr) (yr - mean(yr)) / sd(yr))
    fp <- forecast_panel_ar1(U)
    for (model in panel_models) if (model_flags[[model]]) forecast_matrices[[model]][jj, ] <- fp[[model]]
  }
}

# === Factor-augmented ridge and global boosting (window factors shared by all units) ===
if (any(unlist(model_flags[factor_models]))) {
  for (jj in tau:p) {
    if (jj %% 20 == 0) cat(sprintf("Factor models: window %d/%d\n", jj, p))
    R_win <- ret_mat[(jj+1):(r+jj), , drop = FALSE]
    s_w <- apply(R_win, 2, sd)
    U <- sweep(sweep(R_win, 2, colMeans(R_win)), 2, s_w, "/")
    Fm <- window_factors(R_win, fred_md[(jj+1):(r+jj), -1], n_pc = 3)
    if (model_flags$FAR) forecast_matrices$FAR[jj, ] <- apply(U, 2, forecast_far, Fm = Fm)
    if (model_flags$GBM) forecast_matrices$GBM[jj, ] <- forecast_global_gbm(U, Fm, s_w, seed = jj)
  }
}

# === Equal-weight combination of the ML forecasts ===
if (model_flags$COMB) {
  members <- intersect(comb_members, names(forecast_matrices))
  forecast_matrices$COMB[] <- Reduce(`+`, lapply(members, function(m) {
    x <- forecast_matrices[[m]]; x[is.na(x)] <- 0; x })) /
    Reduce(`+`, lapply(members, function(m) !is.na(forecast_matrices[[m]])))
}

# === Reshape and combine ===
reshape_forecasts <- function(mat, model_name) {
  df_long <- as.data.frame(mat)
  df_long$Date <- forecast_dates
  pivot_longer(df_long, -Date, names_to = "Asset", values_to = model_name)
}
forecast_panel <- Reduce(function(df1, df2) left_join(df1, df2, by = c("Date", "Asset")),
                         lapply(names(forecast_matrices), function(m) reshape_forecasts(forecast_matrices[[m]], m)))

# === Actual returns and merge ===
actual_returns <- as.data.frame(actual_std)                 # standardized by the window
actual_returns$Date <- forecast_dates
actual_long <- pivot_longer(actual_returns, -Date, names_to = "Asset", values_to = "Actual")
forecast_panel <- forecast_panel %>% left_join(actual_long, by = c("Date", "Asset"))

if ("KF" %in% names(forecast_matrices)) {
  forecast_panel <- forecast_panel %>% mutate(KF = ifelse(is.na(KF), AR2, KF))
}

# === Compute losses and differentials ===
for (model in names(forecast_matrices)) {
  forecast_panel[[paste0("Loss_", model)]] <- (forecast_panel[[model]] - forecast_panel$Actual)^2
}
loss_cols <- grep("^Loss_", colnames(forecast_panel), value = TRUE)
for (col in loss_cols) {
  model_name <- sub("Loss_", "", col)
  if (model_name != base_model) {
    diff_col <- paste0("d_", model_name)
    forecast_panel[[diff_col]] <- forecast_panel[[paste0("Loss_", base_model)]] - forecast_panel[[col]]
  }
}

# Add lagged actual
forecast_panel <- forecast_panel %>%
  arrange(Asset, Date) %>%
  group_by(Asset) %>%
  mutate(Actual_Lag1 = lag(Actual)) %>%
  ungroup()

# Save final panel
saveRDS(forecast_panel, file = here("Application", "Exchange Rates", "Data", "forecast_panel_final.rds"))
