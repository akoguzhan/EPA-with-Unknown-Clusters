#' Hybrid Bonferroni Combination of Two P-Value Mergers
#'
#' @description Merges a set of p-values with two user-supplied merging functions and returns
#' twice the minimum of the two merged p-values, capped at 1 (Bonferroni over the two mergers).
#' This follows the hybrid strategy proposed in Vovk, Wang, and Wang (2022). The p-values are first
#' clamped to \code{[.Machine$double.eps, 1]}.
#'
#' @param pvals Numeric vector of p-values in (0, 1].
#' @param combine_fun1,combine_fun2 Functions that take \code{pvals} and return a merged p-value.
#'
#' @return A numeric scalar: the combined p-value using the hybrid Bonferroni rule.
#'
#' @examples
#' bonferroni_compound_pcombine(runif(5), Geomean_pcombine, bonferroni_pcombine)
#'
#' @export
bonferroni_compound_pcombine <- function(pvals, combine_fun1, combine_fun2) {
  # Ensure all p-values are in (0,1]
  pvals <- pmin(pmax(pvals, .Machine$double.eps), 1)

  # Compute combined p-values
  p_comb1 <- combine_fun1(pvals)
  p_comb2 <- combine_fun2(pvals)

  # Apply hybrid rule: Bonferroni over two combinations
  pval <- min(2 * min(p_comb1, p_comb2), 1)
  pval
}

#' Geometric Mean Combination Multiplied by e
#'
#' @description Merges p-values by \eqn{\min\{e\,\exp(n^{-1}\sum_i\log p_i), 1\}}, the geometric
#' mean multiplied by Euler's number. This merger is valid under arbitrary dependence (Vovk and
#' Wang, 2020). The p-values are first clamped to \code{[.Machine$double.eps, 1]}.
#'
#' @param pvals A numeric vector of p-values in (0, 1].
#'
#' @return A numeric scalar: the combined p-value.
#' @references Vovk, V., & Wang, R. (2020). Combining p-values via averaging. \emph{Biometrika}, 107(4), 791-808.
#' @export
Geomean_pcombine <- function(pvals) {
  pvals <- pmin(pmax(pvals, .Machine$double.eps), 1)  # clamp to avoid log(0)
  gmean <- exp(mean(log(pvals)))
  pval <- min(exp(1) * gmean, 1)
  pval
}

#' Generalized Mean P-Value Combination (Positive Order)
#'
#' @description Merges p-values by the generalized mean of positive order r (Vovk and Wang, 2020),
#' valid under arbitrary dependence: \eqn{\min\{\min(r+1, n)^{1/r}(n^{-1}\sum_i p_i^r)^{1/r}, 1\}},
#' with n the number of p-values. For negative orders use \code{Genmean_rneg_pcombine}.
#'
#' @param pvals Numeric vector of p-values in [0, 1].
#' @param r Positive order of the generalized mean (default 20).
#'
#' @return A numeric scalar: the combined p-value.
#'
#' @references Vovk, V., Wang, R., & Wang, R. (2022). Admissible ways of merging
#' p-values under arbitrary dependence.
#'
#' @export
Genmean_pcombine <- function(pvals, r = 20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n_p <- length(pvals)
  M_value = (mean(pvals^r))^(1/r)
  pval <- min(min(r+1,n_p)^(1/r) * M_value, 1)
  pval
}

#' Generalized Mean P-Value Combination with Negative Order
#'
#' @description Merges p-values by the M-family function with order \eqn{r < -1},
#' \eqn{\min\{\frac{r}{r+1}n^{1+1/r}M_r(p), 1\}}, where \eqn{M_r(p) = (n^{-1}\sum_i p_i^r)^{1/r}} is
#' computed on the log scale to avoid overflow. This merger is valid under arbitrary dependence.
#' With a single p-value it is returned unchanged; \code{r = -Inf} gives the Bonferroni p-value.
#'
#' @param pvals Numeric vector of p-values in [0, 1].
#' @param r Order of the generalized mean (default -20).
#'
#' @return A numeric scalar: the combined p-value.
#'
#' @references Vovk, V., Wang, R., & Wang, R. (2022). Admissible ways of merging
#' p-values under arbitrary dependence.
#'
#' @export
Genmean_rneg_pcombine <- function(pvals, r = -20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  
  n_p <- length(pvals)
  if (n_p == 1) return(pvals)
  if (is.infinite(r)) return(min(n_p * min(pvals), 1))  # r = -Inf: Bonferroni
  pvals <- pmax(pvals, .Machine$double.xmin)
  # computed on the log scale to avoid overflow of pvals^r
  lp <- r * log(pvals)
  m <- max(lp)
  M_value <- exp((m + log(mean(exp(lp - m)))) / r)
  pval <- min((r/(r+1))* n_p^(1+1/r) * M_value, 1)
  pval
}

#' IU P-value Combination Test (Spreng and Urga, 2022)
#'
#' @description Merges p-values by \eqn{\min\{\frac{r}{r-1}M_{-r}(p), 1\}}, where
#' \eqn{M_{-r}(p) = (n^{-1}\sum_i p_i^{-r})^{-1/r}} is the generalized mean of order \eqn{-r}.
#'
#' @param pvals A numeric vector of p-values (all in (0, 1]).
#' @param r Power parameter (default = 20).
#'
#' @return A numeric scalar: the combined p-value.
#' @references Spreng, L., & Urga, G. (2022). Combining p-values for Multivariate Predictive Ability Testing.
#' @export
iu_pcombine <- function(pvals, r = 20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  Prn <- (mean(pvals^(-r)))^(1 / r)
  pval <- min(r / ((r - 1) * Prn), 1)
  pval
}

#' Bonferroni Combination Method
#'
#' @description Combines multiple p-values using the Bonferroni correction:
#' \eqn{p = \min(n \min_i p_i, 1)}.
#'
#' @param pvals A numeric vector of p-values (each in [0, 1]).
#' @return A single Bonferroni-adjusted p-value.
#' @examples
#' pvals <- c(0.01, 0.20, 0.05)
#' bonferroni_pcombine(pvals)
#' @export
bonferroni_pcombine <- function(pvals) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  pval <- min(n * min(pvals), 1)
  pval
}

#' Cauchy Combination Test (CCT)
#'
#' @description Combines dependent p-values using the Cauchy combination test proposed by Liu and
#' Xie (2020): \eqn{T = \sum_i w_i\tan[(0.5 - p_i)\pi]} and the p-value is
#' \eqn{0.5 - \arctan(T)/\pi}, clipped to [0, 1].
#'
#' @param p Numeric vector of p-values (all in (0,1)).
#' @param weights Optional non-negative weights (rescaled to sum to 1); default: equal weights.
#'
#' @return Combined p-value.
#' @export
cauchy_pcombine <- function(p, weights = NULL) {
  if (any(p < 0 | p > 1)) stop("All p-values must lie strictly in [0,1].")
  
  n <- length(p)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  } else {
    if (length(weights) != n || any(weights < 0)) stop("Weights must be non-negative and match p-values in length.")
    weights <- weights / sum(weights)  # normalize
  }

  T_stat <- sum(weights * tan((0.5 - p) * pi))
  pval <- min(max(0.5 - atan(T_stat) / pi, 0), 1)
  pval
}
#' Merge p-values with a Named Merging Function
#'
#' @description Dispatches to one of the merging functions of the package by name. Returns
#' \code{NA} if any p-value is missing.
#'
#' @param pvals Numeric vector of p-values.
#' @param pcombine_fun One of "Genmean_neq" (M-family with r < -1, Vovk and Wang, 2020),
#'   "Genmean", "Geomean", "iu", "bonferroni", "cauchy".
#' @param r Order of the generalized mean (default -20; -Inf gives Bonferroni).
#' @return The merged p-value.
#' @export
merge_pvalues <- function(pvals, pcombine_fun = "Genmean_neq", r = -20) {
  if (anyNA(pvals)) return(NA_real_)
  pvals <- as.numeric(pvals)
  out <- switch(
    pcombine_fun,
    Geomean     = Geomean_pcombine(pvals),
    Genmean     = Genmean_pcombine(pvals, r = r),
    Genmean_neq = Genmean_rneg_pcombine(pvals, r = r),
    iu          = iu_pcombine(pvals, r = r),
    bonferroni  = bonferroni_pcombine(pvals),
    cauchy      = cauchy_pcombine(pvals),
    stop("Invalid pcombine_fun specified.")
  )
  unname(as.numeric(out))
}

#' Weighted Bonferroni Combination of Two p-values
#'
#' @description \eqn{\min\{p_1/w, p_2/(1-w)\}\wedge 1}, valid under arbitrary dependence
#' (Remark 9 of the paper, with \eqn{p_1} the homogeneity and \eqn{p_2} the O-EPA p-value).
#' @param p1,p2 p-values.
#' @param w Weight in (0,1) of \code{p1}.
#' @return The combined p-value.
#' @export
weighted_bonferroni <- function(p1, p2, w = 0.5) {
  min(p1 / w, p2 / (1 - w), 1)
}

#' Holm-Adjusted p-values
#'
#' @description Wrapper of \code{stats::p.adjust} with \code{method = "holm"}.
#'
#' @param pvals Numeric vector.
#' @return Holm-adjusted p-values (same order).
#' @export
holm_adjust <- function(pvals) stats::p.adjust(pvals, method = "holm")
