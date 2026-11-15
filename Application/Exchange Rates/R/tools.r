#' Forecast with AR(1) model using matrix operations
#'
#' @param y A numeric vector of returns
#' @return A numeric value: the 1-step-ahead forecast from AR(1)
#' @export
forecast_ar1 <- function(y) {
  T_obs <- length(y)

  # Construct lagged data manually
  y_response <- y[2:T_obs]
  y_lag1     <- y[1:(T_obs - 1)]
  
  # Design matrix with intercept and lag
  X <- y_lag1
  
  # OLS estimation via matrix algebra
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y_response
  
  # Forecast: plug in the last observed value
  y_pred <- y[T_obs] %*% beta_hat
  return(as.numeric(y_pred))
}

#' Forecast with AR(2) model using matrix operations
#'
#' @param y A numeric vector of returns
#' @return A numeric value: the 1-step-ahead forecast from AR(2)
#' @export
forecast_ar2 <- function(y) {
  T_obs <- length(y)

  # Construct lags manually
  y_response <- y[3:T_obs]
  y_lag1     <- y[2:(T_obs - 1)]
  y_lag2     <- y[1:(T_obs - 2)]
  
  # Design matrix with intercept and lags
  X <- cbind(y_lag1, y_lag2)
  
  # OLS via matrix algebra
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y_response
  
  # Forecast using last two values
  new_X <- c(y[T_obs], y[T_obs - 1])
  y_pred <- new_X %*% beta_hat
  return(as.numeric(y_pred))
}

#' Forecast with AR(p) model selected by BIC
#'
#' Selects the optimal lag length p for an autoregressive model using BIC,
#' estimates the model, and returns a one-step-ahead forecast.
#'
#' @param y A numeric vector of returns.
#' @param max_lag The maximum number of lags to consider (default: 10).
#' @return A numeric value: the 1-step-ahead forecast from AR(p).
#' @export
forecast_ar_bic <- function(y, max_lag = 6) {
  T_obs <- length(y)

  max_lag <- min(max_lag, T_obs - 2)  # ensure we have enough obs
  BIC_vals <- rep(Inf, max_lag)
  models <- vector("list", max_lag)
  
  for (p in 1:max_lag) {
    y_mat <- embed(y, p + 1)
    y_response <- y_mat[, 1]
    X <- y_mat[, -1, drop = FALSE]
    
    # add intercept
    X_design <- cbind(1, X)
    n <- nrow(X_design)
    
    # Estimate coefficients via OLS
    beta_hat <- tryCatch({
      solve(t(X_design) %*% X_design, t(X_design) %*% y_response)
    }, error = function(e) return(NULL))
    
    if (!is.null(beta_hat)) {
      residuals <- y_response - X_design %*% beta_hat
      sigma2 <- mean(residuals^2)
      k <- ncol(X_design)
      BIC_vals[p] <- log(sigma2) + k * log(n) / n
      models[[p]] <- beta_hat
    }
  }
  
  best_p <- which.min(BIC_vals)
  beta_hat_best <- models[[best_p]]
  if (is.null(beta_hat_best)) return(forecast_rw(y))
  
  # Forecast using the most recent p values
  X_new <- rev(y[(T_obs - best_p + 1):T_obs])
  y_pred <- c(1, X_new) %*% beta_hat_best
  return(as.numeric(y_pred))
}

#' Forecast with Random Walk model
#'
#' @param y A numeric vector of returns
#' @return A numeric value: the last observed return as forecast
forecast_rw <- function(y) {
  return(tail(y, 1))
}

#' Check if a time series is usable for forecasting
#'
#' Removes flat or nearly constant series (rank-deficient) before modeling.
#'
#' @param y A numeric vector (e.g., monthly returns).
#' @param threshold Standard deviation threshold (default = 1e-5).
#'
#' @return TRUE if series is acceptable; FALSE otherwise.
#' @export
is_valid_series <- function(y, threshold = 1e-5) {
  if (!is.numeric(y)) return(FALSE)
  if (length(y) < 5) return(FALSE)
  sd_y <- sd(y, na.rm = TRUE)
  sd_lag <- sd(stats::lag(y, 1), na.rm = TRUE)
  range_y <- diff(range(y, na.rm = TRUE))
  
  if (any(is.na(c(sd_y, sd_lag, range_y)))) return(FALSE)
  if (sd_y < threshold || sd_lag < threshold || range_y < threshold) return(FALSE)
  
  return(TRUE)
}

#' Forecast with Gradient Boosted Trees (XGBoost)
#'
#' Uses XGBoost to forecast one-step-ahead log returns using lagged features.
#' Handles nonlinearity and interaction effects. Simple boosting configuration.
#'
#' @param y Numeric vector of log returns.
#' @param max_lag Maximum number of lags to include. Default is 12.
#'
#' @return A numeric forecast value for the next period, or NA if model fails.
#' @export
forecast_xgboost <- function(y, max_lag = 6) {
  # Build lagged matrix
  T_len <- length(y)
  Y <- y[(max_lag + 1):T_len]
  X <- sapply(1:max_lag, function(lag) y[(max_lag + 1 - lag):(T_len - lag)])
  colnames(X) <- paste0("Lag", 1:max_lag)
  
  dtrain <- xgboost::xgb.DMatrix(data = X, label = Y)
  
  # Train with default hyperparameters
  fit <- tryCatch(
      xgboost::xgboost(
      data = dtrain,
      nrounds = 50,
      objective = "reg:squarederror",
      verbose = 0
    ),
    error = function(e) return(NULL)
  )
  if (is.null(fit)) return(forecast_rw(y))
  
  # Predict next value using most recent lags
  newx <- matrix(rev(tail(y, max_lag)), nrow = 1)
  pred <- predict(fit, newdata = newx)
  return(as.numeric(pred))
}

#' Forecast with Elastic Net (Sparse AR) Model with CV-selected Alpha
#'
#' Computes a one-step-ahead forecast using Elastic Net regression, selecting alpha and lambda via cross-validation.
#'
#' @param y Numeric vector of log returns.
#' @param max_lag Maximum number of lags to include as regressors. Default is 6.
#' @param alpha_grid A numeric vector of candidate alpha values in [0, 1]. Default: seq(0, 1, by = 0.1)
#'
#' @return A numeric value: the one-step-ahead forecast, or NA if estimation fails.
#' @export
forecast_lasso <- function(y, max_lag = 6, alpha_grid = seq(0, 1, by = 0.1)) {
  if (!is.numeric(y) || length(y) < max_lag + 1) return(NA)
  if (sd(y) == 0) return(NA)
  
  T_len <- length(y)
  Y <- y[(max_lag + 1):T_len]
  
  X <- sapply(1:max_lag, function(lag) y[(max_lag + 1 - lag):(T_len - lag)])
  colnames(X) <- paste0("Lag", 1:max_lag)
  
  best_cv <- NULL
  best_alpha <- NA
  best_error <- Inf
  
  for (a in alpha_grid) {
    fit <- tryCatch(
      cv.glmnet(X, Y, alpha = a, standardize = FALSE, nfolds = 5),
      error = function(e) NULL
    )
    if (!is.null(fit) && min(fit$cvm) < best_error) {
      best_cv <- fit
      best_alpha <- a
      best_error <- min(fit$cvm)
    }
  }
  
  if (is.null(best_cv)) return(forecast_rw(y))
  
  newx <- matrix(rev(tail(y, max_lag)), nrow = 1)
  pred <- predict(best_cv, newx = newx, s = best_cv$lambda.min)
  
  return(as.numeric(pred))
}

#' Forecast with Elastic Net using lagged y and macro variables (CV for alpha & lambda)
#'
#' Performs one-step-ahead forecast using Elastic Net regression. Selects alpha and lambda via cross-validation.
#'
#' @param y A numeric vector of returns (dependent variable)
#' @param X A matrix or data.frame of macro predictors (must match length of y)
#' @param alpha_grid A numeric vector of alpha values to try (e.g., seq(0,1,0.1))
#' @return A numeric value: the one-step-ahead forecast or NA if fitting fails
#' @export
forecast_lasso_macro <- function(y, X, alpha_grid = seq(0, 1, by = 0.1)) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Please install the glmnet package.")
  
  T_obs <- length(y)
  if (T_obs < 2 || nrow(X) != T_obs) return(NA_real_)
  
  y_response <- y[2:T_obs]
  y_lag <- y[1:(T_obs - 1)]
  X_lag <- scale(X[1:(T_obs - 1), , drop = FALSE])  # Standardize macro predictors
  
  X_design <- as.matrix(cbind(y_lag, X_lag))
  y_response <- as.numeric(y_response)
  
  best_cv <- NULL
  best_alpha <- NA
  best_error <- Inf
  
  for (a in alpha_grid) {
    fit <- tryCatch(
      glmnet::cv.glmnet(X_design, y_response, alpha = a, nfolds = 5, standardize = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit) && min(fit$cvm) < best_error) {
      best_cv <- fit
      best_alpha <- a
      best_error <- min(fit$cvm)
    }
  }
  
  if (is.null(best_cv)) return(forecast_rw(y))
  
  # Construct new data
  x_new <- as.numeric(c(y[T_obs], scale(X[T_obs, , drop = FALSE], center = attr(X_lag, "scaled:center"), scale = attr(X_lag, "scaled:scale"))))
  
  beta_hat <- as.numeric(coef(best_cv, s = "lambda.min"))[-1]
  intercept <- as.numeric(coef(best_cv, s = "lambda.min"))[1]
  
  y_pred <- intercept + sum(x_new * beta_hat)
  return(y_pred)
}

#' Design Matrix with the Own Lag and Macro Predictors
#'
#' @description Common design of the macro-based learners (SVM, random forest, boosting):
#' the response \eqn{y_t} is regressed on \eqn{y_{t-1}} and the macro predictors at
#' \eqn{t-1}. The macro predictors are standardized with the mean and standard deviation of
#' the training rows, columns that are constant in the window are dropped, and the prediction
#' row (\eqn{y_T}, \eqn{X_T}) is standardized with the same moments.
#' @param y Standardized returns of the estimation window (length R).
#' @param X R x p matrix of macro predictors aligned with \code{y}.
#' @return A list with the (R-1) x (p+1) training matrix \code{X}, the response \code{y} and
#' the 1 x (p+1) prediction row \code{x_new} (data frames with the same column names).
#' @export
macro_design <- function(y, X) {
  X <- as.matrix(X)
  T_obs <- length(y)
  Xl <- X[1:(T_obs - 1), , drop = FALSE]
  mu <- colMeans(Xl); sdv <- apply(Xl, 2, stats::sd)
  keep <- is.finite(sdv) & sdv > 1e-8
  Xs <- sweep(sweep(Xl[, keep, drop = FALSE], 2, mu[keep]), 2, sdv[keep], "/")
  xn <- (X[T_obs, keep] - mu[keep]) / sdv[keep]
  cn <- c("y_lag", paste0("x", seq_len(sum(keep))))
  D <- as.data.frame(cbind(y[1:(T_obs - 1)], Xs)); names(D) <- cn
  N <- as.data.frame(matrix(c(y[T_obs], xn), nrow = 1)); names(N) <- cn
  list(X = D, y = y[2:T_obs], x_new = N)
}

#' Gradient Boosting Forecast with the Own Lag and Macro Predictors
#'
#' @description Boosted regression trees (\pkg{xgboost}, squared error loss) on the design of
#' \code{\link{macro_design}}: 100 rounds, learning rate 0.05, depth 2, at least 5
#' observations per leaf. Small trees and a low learning rate keep the fit stable with about 60
#' observations. Uses \code{xgb.train}, whose interface is the same in \pkg{xgboost} 1.x and 3.x.
#' @param y Standardized returns of the estimation window.
#' @param X Macro predictors aligned with \code{y}.
#' @return The one-step-ahead forecast; the last value of \code{y} if fitting fails.
#' @export
forecast_xgboost_macro <- function(y, X) {
  d <- macro_design(y, X)
  fit <- tryCatch(
    xgboost::xgb.train(params = list(objective = "reg:squarederror", eta = 0.05, max_depth = 2,
                                     min_child_weight = 5, nthread = 1),
                       data = xgboost::xgb.DMatrix(as.matrix(d$X), label = d$y),
                       nrounds = 100, verbose = 0),
    error = function(e) NULL)
  if (is.null(fit)) return(forecast_rw(y))
  as.numeric(predict(fit, xgboost::xgb.DMatrix(as.matrix(d$x_new))))
}

#' Forecast using Support Vector Regression (SVR) with lagged y and X
#'
#' @param y A numeric vector of returns
#' @param X A matrix or data.frame of macro predictors
#' @param ... Additional arguments passed to e1071::svm
#'
#' @return 1-step ahead forecast (numeric)
#' @export
forecast_svr_macro <- function(y, X, ...) {
  d <- macro_design(y, X)
  model <- tryCatch(
    e1071::svm(x = d$X, y = d$y, type = "eps-regression", ...),
    error = function(e) NULL
  )
  if (is.null(model)) return(forecast_rw(y))
  as.numeric(predict(model, d$x_new))
}

#' Forecast using Random Forest with macro predictors
#'
#' @param y Dependent variable (returns)
#' @param X Matrix or data.frame of macro predictors
#' @return Numeric 1-step-ahead forecast
#' @export
forecast_rf <- function(y, X) {
  if (length(y) < 2) return(NA_real_)
  d <- macro_design(y, X)
  rf_model <- randomForest::randomForest(x = d$X, y = d$y)
  as.numeric(predict(rf_model, newdata = d$x_new))
}

#' Forecast using Markov Switching AR(1) model (manual prediction)
#'
#' @param y A numeric vector (time series)
#' @return Numeric 1-step-ahead forecast
#' @export
forecast_msar <- function(y) {
  if (length(y) < 10) return(NA_real_)
  
  T_obs <- length(y)
  y_lag <- y[1:(T_obs - 1)]
  y_resp <- y[2:T_obs]
  df <- data.frame(y_resp = y_resp, y_lag = y_lag)
  
  # Base AR(1) model
  lm_fit <- lm(y_resp ~ y_lag, data = df)
  
  # Markov Switching with 2 regimes (wrapped in try)
  ms_model <- try(
    MSwM::msmFit(lm_fit, k = 2, sw = rep(TRUE, 3), control = list(parallel = FALSE)),
    silent = TRUE
  )
  
  if (inherits(ms_model, "try-error")) {
    return(forecast_rw(y))
  }
  
  # Extract regime-specific coefficients
  coef_mat <- as.matrix(ms_model@Coef)
  intercepts <- coef_mat[, 1]
  slopes     <- coef_mat[, 2]
  
  # Extract smoothed probabilities
  prob_t <- ms_model@Fit@smoProb[T_obs - 1, ]  # last smoothed probs at time T
  
  # Last observed y
  y_last <- y[T_obs]
  pred_r <- intercepts + slopes * y_last
  
  # Weighted average prediction
  y_pred <- sum(prob_t * pred_r)
  return(as.numeric(y_pred))
}


#' Hard Thresholding Forecast with Lag Selection and BIC-based PC Selection
#'
#' Implements the hard-thresholding method for forecasting a univariate series using a large set of predictors.
#' 1. Selects predictors by regressing \eqn{y_t} on lagged \eqn{y_{t-1}} and \eqn{x_{t-1,j}} one-by-one and thresholds based on absolute t-statistics.
#' 2. Performs PCA on the selected predictors (not their lags).
#' 3. Chooses number of PCs and their lags via BIC in an ARX regression.
#'
#' @param y Numeric vector: univariate target time series (e.g., returns).
#' @param X Matrix or data.frame: predictors (e.g., FRED-MD data) with same length as \code{y}.
#' @param threshold Numeric: minimum absolute t-statistic to include a variable.
#' @param max_pcs Integer: maximum number of principal components to consider.
#' @param max_lag Integer: maximum lag order to try for PCs.
#'
#' @return Numeric value: one-step-ahead forecast using the selected model.
#' @export
forecast_threshold_hard <- function(y, X, threshold = 1.65, max_pcs = 6, max_lag = 2) {
  T_obs <- length(y)
  if (nrow(X) != T_obs) stop("Length of y must match number of rows in X.")
  if (T_obs < 20) return(NA)
  
  # Step 1: Absolute t-stats for each predictor's lag-1 reg
  t_stats <- sapply(1:ncol(X), function(j) {
    xj_lag1 <- X[1:(T_obs - 1), j]
    yt <- y[2:T_obs]
    yt_lag1 <- y[1:(T_obs - 1)]
    reg_data <- data.frame(y = yt, y_lag1 = yt_lag1, x_lag1 = xj_lag1)
    fit <- try(lm(y ~ y_lag1 + x_lag1, data = reg_data), silent = TRUE)
    if (inherits(fit, "try-error")) return(0)
    tval <- summary(fit)$coefficients["x_lag1", "t value"]
    if (is.na(tval)) 0 else abs(tval)
  })
  
  # Step 2: Select Xs with |t| > threshold
  selected_idx <- which(t_stats > threshold)
  if (length(selected_idx) == 0) return(NA)
  X_selected <- X[, selected_idx, drop = FALSE]
  
  # Step 3: PCA on selected Xs (no lagging here)
  pca_out <- prcomp(na.omit(X_selected), center = TRUE, scale. = TRUE)
  PCs <- pca_out$x  # T_obs x m
  
  # Step 4: Build BIC grid over number of PCs and lags
  y_lag1 <- y[2:T_obs]
  yt <- y[3:T_obs]
  best_bic <- Inf
  best_pred <- NA
  
  for (k in 1:min(max_pcs, ncol(PCs))) {
    PC_mat <- as.matrix(PCs[, 1:k])

    for (lag in 0:max_lag) {
      start_idx <- lag + 2
      y_dep <- y[start_idx:T_obs]
      y_lag <- y[(start_idx - 1):(T_obs - 1)]
      
      # Build regressors: PC_lag0, PC_lag1, ..., PC_lag_l
      Z_list <- list(y_lag)
      for (l in 0:lag) {
        PC_lag <- PC_mat[(start_idx - l):(T_obs - l), , drop = FALSE]
        Z_list[[length(Z_list) + 1]] <- PC_lag
      }
      Z <- do.call(cbind, Z_list)
      fit <- try(lm(y_dep ~ Z), silent = TRUE)
      if (inherits(fit, "try-error")) next
      
      n <- length(y_dep)
      rss <- sum(resid(fit)^2)
      k_params <- ncol(Z) + 1
      bic <- n * log(rss / n) + k_params * log(n)
      if (!is.na(bic) && bic < best_bic) {
        best_bic <- bic
        coef_hat <- coef(fit)
        Z_new <- c(tail(y, 1))
        for (l in 0:lag) {
          Z_new <- c(Z_new, tail(PC_mat, lag + 1 - l)[1, 1:k])
        }
        Z_new <- c(1, Z_new)  # intercept
        best_pred <- sum(coef_hat * Z_new)
      }
    }
  }
  
  return(best_pred)
}

#' Soft Thresholding Forecasting with Lasso and PCA
#'
#' Implements soft thresholding for forecasting: selects exactly 30 variables via Lasso,
#' performs PCA on them, and uses BIC to choose number and lags of PCs.
#'
#' @param y Numeric vector of target variable (e.g., returns).
#' @param X Matrix or data.frame of predictors (e.g., FRED-MD data).
#' @param max_pc Maximum number of PCs to consider (default = 5).
#' @param max_lag Maximum number of lags of PCs to try (default = 2).
#'
#' @return One-step-ahead forecast value (numeric), or NA if model fails.
#' @export
forecast_threshold_soft <- function(y, X, max_pcs = 6, max_lag = 2) {
  T_obs <- length(y)
  if (nrow(X) != T_obs) stop("Length of y must match number of rows in X.")
  if (T_obs < 20) return(NA)
  
  y_lag <- y[-length(y)]
  X_lag <- X[-nrow(X), , drop = FALSE]
  y_dep <- y[-1]
  
  # Standardize predictors
  X_scaled <- scale(X_lag, center = F)
  X_design <- cbind(y_lag, X_scaled)
  penalty <- c(0, rep(1, ncol(X_scaled)))  # no penalty on lagged y
  
  # Fit full Lasso path
  fit_full <- glmnet::glmnet(x = X_design, y = y_dep, alpha = 1, penalty.factor = penalty)
  lambda_seq <- fit_full$lambda
  lambda_chosen <- NULL
  
  for (lambda in lambda_seq) {
    coef_lasso <- coef(fit_full, s = lambda)
    nonzero <- which(coef_lasso[-c(1, 2)] != 0)  # exclude intercept and y_lag
    if (length(nonzero) == 30) {
      lambda_chosen <- lambda
      break
    }
  }
  
  if (is.null(lambda_chosen)) return(NA)
  
  # Refit at chosen lambda
  coef_final <- coef(fit_full, s = lambda_chosen)
  selected_idx <- which(coef_final[-c(1, 2)] != 0)
  X_selected <- X[, selected_idx, drop = FALSE]
  
  # Step 3: PCA on selected Xs (no lagging here)
  pca_out <- prcomp(na.omit(X_selected), center = TRUE, scale. = TRUE)
  PCs <- pca_out$x  # T_obs x m
  
  # Step 4: Build BIC grid over number of PCs and lags
  y_lag1 <- y[2:T_obs]
  yt <- y[3:T_obs]
  best_bic <- Inf
  best_pred <- NA
  
  for (k in 1:min(max_pcs, ncol(PCs))) {
    PC_mat <- as.matrix(PCs[, 1:k])
    
    for (lag in 0:max_lag) {
      start_idx <- lag + 2
      y_dep <- y[start_idx:T_obs]
      y_lag <- y[(start_idx - 1):(T_obs - 1)]
      
      # Build regressors: PC_lag0, PC_lag1, ..., PC_lag_l
      Z_list <- list(y_lag)
      for (l in 0:lag) {
        PC_lag <- PC_mat[(start_idx - l):(T_obs - l), , drop = FALSE]
        Z_list[[length(Z_list) + 1]] <- PC_lag
      }
      Z <- do.call(cbind, Z_list)
      fit <- try(lm(y_dep ~ Z), silent = TRUE)
      if (inherits(fit, "try-error")) next
      
      n <- length(y_dep)
      rss <- sum(resid(fit)^2)
      k_params <- ncol(Z) + 1
      bic <- n * log(rss / n) + k_params * log(n)
      if (!is.na(bic) && bic < best_bic) {
        best_bic <- bic
        coef_hat <- coef(fit)
        Z_new <- c(tail(y, 1))
        for (l in 0:lag) {
          Z_new <- c(Z_new, tail(PC_mat, lag + 1 - l)[1, 1:k])
        }
        Z_new <- c(1, Z_new)  # intercept
        best_pred <- sum(coef_hat * Z_new)
      }
    }
  }
  
  return(best_pred)
}
# ==============================================================================
# Forecasts added in October 2026: no-change, shrinkage and panel AR(1) forecasts.
# All AR(1) forecasts below work on a window that has been demeaned and scaled by its own
# mean and standard deviation (see generate_forecasts.R), so they have no intercept.
# ==============================================================================

#' No-Change Forecast of the Exchange Rate (Zero Return)
#'
#' @param m,s Mean and standard deviation of the raw returns in the estimation window.
#' @return The forecast of the standardized return, \eqn{(0 - m)/s}.
#' @export
forecast_zero <- function(m, s) -m / s

#' AR(1) Forecast with a James-Stein Type Shrinkage of the Coefficient
#'
#' @description The least squares coefficient \eqn{\hat\phi} of the standardized window is
#' multiplied by \eqn{\max(0, 1 - 1/t^2)}, where t is its t-statistic, so that insignificant
#' coefficients are shrunk towards zero.
#' @param y Standardized returns of the estimation window.
#' @return The one-step-ahead forecast of the standardized return.
#' @export
forecast_ar1_js <- function(y) {
  x <- y[-length(y)]; z <- y[-1]
  phi <- sum(x * z) / sum(x^2)
  e <- z - phi * x
  se <- sqrt(sum(e^2) / (length(z) - 1) / sum(x^2))
  tstat <- phi / se
  phi * max(0, 1 - 1 / tstat^2) * y[length(y)]
}

#' Panel AR(1) Forecasts: Pooled, Empirical Bayes and Grouped
#'
#' @description For an R x N matrix of standardized estimation windows, computes one-step-ahead
#' forecasts with (i) the pooled AR(1) coefficient; (ii) empirical Bayes coefficients, which
#' shrink each unit's least squares coefficient towards the pooled one with weight
#' \eqn{\tau^2/(\tau^2 + v_i)}, where \eqn{v_i} is its sampling variance and \eqn{\tau^2} the
#' estimated cross-sectional variance of the coefficients; (iii) grouped coefficients, pooled
#' within the K groups of a one-dimensional Kmeans of the unit coefficients (initial centers at
#' the 10%, 50% and 90% quantiles for K = 3, so that the result is deterministic).
#' @param U R x N matrix of standardized returns (one column per unit).
#' @param K Number of groups of the grouped estimator.
#' @return A list with the N-vectors \code{POOL}, \code{EB} and \code{GRP} of forecasts.
#' @export
forecast_panel_ar1 <- function(U, K = 3) {
  X <- U[-nrow(U), , drop = FALSE]; Z <- U[-1, , drop = FALSE]
  last <- U[nrow(U), ]
  sxx <- colSums(X^2); sxz <- colSums(X * Z)
  phi_i <- sxz / sxx
  phi_pool <- sum(sxz) / sum(sxx)
  s2_i <- colSums((Z - sweep(X, 2, phi_i, "*"))^2) / (nrow(Z) - 1)
  v_i <- s2_i / sxx
  tau2 <- max(0, stats::var(phi_i) - mean(v_i))
  phi_eb <- phi_pool + tau2 / (tau2 + v_i) * (phi_i - phi_pool)
  cl <- tryCatch(stats::kmeans(phi_i, centers = stats::quantile(phi_i, seq(0.1, 0.9, length.out = K)))$cluster,
                 error = function(e) rep(1L, length(phi_i)))   # e.g. tied quantiles: one group
  phi_grp <- tapply(sxz, cl, sum) / tapply(sxx, cl, sum)
  list(POOL = phi_pool * last, EB = phi_eb * last, GRP = as.numeric(phi_grp[as.character(cl)]) * last)
}


#' Window Factors: Dollar Factor and Macro Principal Components
#'
#' @description Cross-sectional predictors shared by all currencies in one estimation window:
#' the dollar factor (cross-sectional mean of the raw returns, the first factor of exchange rates
#' against the US dollar in Lustig, Roussanov and Verdelhan, 2011) and the first \code{n_pc}
#' principal components of the macro predictors. Both are standardized within the window, so
#' that no information after the forecast date is used.
#' @param R_win R x N matrix of raw returns of the window.
#' @param X_win R x p matrix of macro predictors of the window.
#' @param n_pc Number of principal components.
#' @return R x (1 + n_pc) matrix with columns \code{DOL}, \code{PC1}, ...
#' @export
window_factors <- function(R_win, X_win, n_pc = 3) {
  X_win <- as.matrix(X_win)
  sdv <- apply(X_win, 2, stats::sd)
  Xs <- scale(X_win[, is.finite(sdv) & sdv > 1e-8, drop = FALSE])
  pcs <- svd(Xs, nu = n_pc, nv = 0)$u
  Fm <- scale(cbind(rowMeans(R_win), pcs))
  colnames(Fm) <- c("DOL", paste0("PC", seq_len(n_pc)))
  Fm
}

#' Factor-Augmented Ridge Forecast
#'
#' @description Ridge regression of \eqn{y_t} on \eqn{y_{t-1}} and the window factors at
#' \eqn{t-1} (dollar factor and macro principal components, \code{\link{window_factors}}),
#' a diffusion-index forecast in the sense of Stock and Watson (2002) with ridge shrinkage. The
#' variables are centred (the intercept is not penalized) and the penalty is chosen by
#' generalized cross-validation on a log grid, which is fast and deterministic.
#' @param y Standardized returns of the estimation window (length R).
#' @param Fm R x q matrix of window factors.
#' @param lambdas Grid of ridge penalties.
#' @return The one-step-ahead forecast.
#' @export
forecast_far <- function(y, Fm, lambdas = 10^seq(-2, 3, length.out = 30)) {
  T_obs <- length(y)
  X <- cbind(y[1:(T_obs - 1)], Fm[1:(T_obs - 1), , drop = FALSE])
  z <- y[2:T_obs]
  mx <- colMeans(X); mz <- mean(z)
  Xc <- sweep(X, 2, mx); zc <- z - mz
  sv <- svd(Xc); dz <- crossprod(sv$u, zc)
  n <- length(z)
  gcv <- sapply(lambdas, function(l) {
    sh <- sv$d^2 / (sv$d^2 + l)
    res <- zc - sv$u %*% (sh * dz)
    (sum(res^2) / n) / (1 - sum(sh) / n)^2
  })
  l <- lambdas[which.min(gcv)]
  beta <- sv$v %*% (sv$d / (sv$d^2 + l) * dz)
  as.numeric(mz + (c(y[T_obs], Fm[T_obs, ]) - mx) %*% beta)
}

#' Global Gradient Boosting Forecasts
#'
#' @description One boosted tree model (\pkg{xgboost}) trained on all currencies of the
#' estimation window together (a global model in the sense of Montero-Manso and Hyndman, 2021),
#' so that about N x (R - L) observations are available instead of R - L per currency. Features
#' of currency i at date t: its own standardized lags \eqn{y_{i,t-1}, ..., y_{i,t-L}}, the
#' window factors at t-1 (dollar factor, macro principal components) and two characteristics of
#' the currency in the window (log standard deviation of its raw returns and its first-order
#' autocorrelation), which let the trees treat groups of currencies differently.
#' @param U R x N matrix of standardized returns of the window (one column per currency).
#' @param Fm R x q matrix of window factors.
#' @param s_w N-vector of window standard deviations of the raw returns.
#' @param L Number of own lags.
#' @param seed Seed of the row and column subsampling.
#' @return N-vector of one-step-ahead forecasts (zeros, the window mean, if fitting fails).
#' @export
forecast_global_gbm <- function(U, Fm, s_w, L = 6, seed = 1) {
  R <- nrow(U); N <- ncol(U)
  rho <- colSums(U[-1, , drop = FALSE] * U[-R, , drop = FALSE]) / colSums(U^2)
  feat <- function(i, t) c(U[t - (1:L) + 1, i], Fm[t, ], log(s_w[i]), rho[i])  # info up to t
  tt <- L:(R - 1)
  Xtr <- do.call(rbind, lapply(seq_len(N), function(i) t(sapply(tt, function(t) feat(i, t)))))
  ytr <- as.vector(U[tt + 1, ])
  Xnew <- t(sapply(seq_len(N), function(i) feat(i, R)))
  set.seed(seed)
  fit <- tryCatch(
    xgboost::xgb.train(params = list(objective = "reg:squarederror", eta = 0.03, max_depth = 3,
                                     min_child_weight = 50, subsample = 0.8,
                                     colsample_bytree = 0.8, nthread = 1),
                       data = xgboost::xgb.DMatrix(Xtr, label = ytr), nrounds = 200, verbose = 0),
    error = function(e) NULL)
  if (is.null(fit)) return(rep(0, N))
  as.numeric(predict(fit, xgboost::xgb.DMatrix(Xnew)))
}
