#' Number of Cosines of the Orthonormal Series (OS) Estimator
#'
#' @description Chooses the number of cosines \eqn{B = \lceil c_B T^{2/3} \rceil}. If \code{cB}
#' is \code{NULL}, \eqn{c_B} is given by the testing-optimal (CPE) rule of Sun (2013, eq. 8,
#' p = 1 and 5\% level), \eqn{c_B = 0.42293\,|\bar b|^{-1/3}}, where
#' \eqn{\bar b = \pi^2\rho/[3(1-\rho)^2]} is the relative bias of the cosine OS estimator under an
#' AR(1) approximation with coefficient \eqn{\rho}, averaged over the columns of \code{X}.
#' The same \eqn{B} is used for all tests computed on the same sample.
#'
#' @param X A T x P matrix used to fit the AR(1) approximation (typically the cross-sectional
#'   average of the panel, \eqn{\bar Z_t}). Ignored if \code{cB} is given.
#' @param cB Optional numeric constant \eqn{c_B}.
#' @param Bmin Lower bound on B (e.g. \eqn{K_{max}P+4}, so that all Wald matrices are invertible).
#' @param Tobs Number of time periods (defaults to \code{nrow(X)}).
#'
#' @return Integer B.
#' @references Sun, Y. (2013). A heteroskedasticity and autocorrelation robust F test using an
#'   orthonormal series variance estimator. \emph{Econometrics Journal}, 16, 1--26.
#' @export
choose_B <- function(X = NULL, cB = NULL, Bmin = 1, Tobs = NULL) {
  if (is.null(Tobs)) Tobs <- NROW(X)
  if (is.null(cB)) cB <- sun_cB(X)
  B <- if (is.finite(cB)) ceiling(cB * Tobs^(2/3)) else Tobs
  B <- max(B, Bmin)
  as.integer(min(B, Tobs))
}

#' Constant of the Testing-Optimal Rule for B under an AR(1) Approximation
#'
#' @description Computes \eqn{c_B = 0.42293\,\bar b^{-1/3}} with
#' \eqn{\bar b = \pi^2\rho/[3(1-\rho)^2]} averaged over the columns of \code{X}, where \eqn{\rho} is
#' the absolute first-order autocorrelation of each demeaned column, capped at 0.97 (see
#' \code{choose_B}).
#'
#' @param X A T x P matrix.
#' @return The constant \eqn{c_B} (Inf if the estimated autocorrelation is zero).
#' @keywords internal
sun_cB <- function(X) {
  X <- as.matrix(X)
  rho <- apply(X, 2, function(x) {
    x <- x - mean(x)
    n <- length(x)
    den <- sum(x[-n]^2)
    if (den <= 0) return(0)
    sum(x[-1] * x[-n]) / den
  })
  rho <- pmin(abs(rho), 0.97)
  bbar <- mean(pi^2 * rho / (3 * (1 - rho)^2))
  if (bbar <= 0) return(Inf)
  0.42293 * bbar^(-1/3)
}

#' Local (Nonparametric) Demeaning
#'
#' @description Subtracts a Nadaraya-Watson (Epanechnikov) estimate of a time-varying mean from
#' each column, as in Harvey, Leybourne and Zu (2024). Used only for long-run variance estimation.
#'
#' @param X A T x M matrix.
#' @param bw Half-width of the kernel window in periods; default \code{ceiling(T^(4/5)/2)}.
#' @return The T x M matrix of locally demeaned columns.
#' @keywords internal
local_demean <- function(X, bw = NULL) {
  X <- as.matrix(X)
  Tobs <- nrow(X)
  if (is.null(bw)) bw <- ceiling(Tobs^(4/5) / 2)
  tt <- 1:Tobs
  Wk <- outer(tt, tt, function(t, s) pmax(1 - ((t - s) / bw)^2, 0))
  Wk <- Wk / rowSums(Wk)
  X - Wk %*% X
}

#' Cosine Transforms of the OS Estimator
#'
#' @description Computes \eqn{\hat\Lambda_j = \sqrt{2/T}\sum_t (X_t - \bar X)\cos[\pi j(t-1/2)/T]},
#' j = 1, ..., B, for every column of \code{X}.
#'
#' @param X A T x M matrix (columns are series).
#' @param B Number of cosines.
#' @param demean "global" (full-sample mean) or "local" (see \code{local_demean}).
#' @param bw Bandwidth for local demeaning.
#' @return A B x M matrix whose j-th row is \eqn{\hat\Lambda_j'}.
#' @export
cosine_transforms <- function(X, B, demean = c("global", "local"), bw = NULL) {
  demean <- match.arg(demean)
  X <- as.matrix(X)
  Tobs <- nrow(X)
  if (B < 1 || B > Tobs) stop("B must be between 1 and T.")
  Xc <- if (demean == "global") sweep(X, 2, colMeans(X)) else local_demean(X, bw)
  Cmat <- cos(pi * outer(1:B, (1:Tobs) - 0.5) / Tobs)  # B x T
  sqrt(2 / Tobs) * (Cmat %*% Xc)
}

#' Long-run Variance Estimation via the Orthonormal Series (EWC) Estimator
#'
#' @description \eqn{\hat\Omega = B^{-1}\sum_{j=1}^B \hat\Lambda_j\hat\Lambda_j'} with cosine basis.
#'
#' @param X A T x P numeric matrix.
#' @param lrv_par Number of cosines B. If NULL, chosen by \code{choose_B(X, cB, Bmin = P + 4)}.
#' @param cB Optional constant of the rule \eqn{B = \lceil c_B T^{2/3}\rceil}.
#' @param demean "global" or "local".
#' @param bw Bandwidth for local demeaning.
#'
#' @return A list with the P x P estimate \code{S} and the \code{lrv_par} used.
#' @export
EWC <- function(X, lrv_par = NULL, cB = NULL, demean = "global", bw = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)
  if (is.null(lrv_par)) lrv_par <- choose_B(X, cB = cB, Bmin = P + 4)
  if (lrv_par < 1 || lrv_par > Tobs) stop("lrv_par must be between 1 and Tobs.")
  Lj <- cosine_transforms(X, lrv_par, demean = demean, bw = bw)  # B x P
  OmegaHat <- crossprod(Lj) / lrv_par
  return(list(S = OmegaHat, lrv_par = lrv_par))
}

#' Long-run Variance Estimation via Newey-West
#'
#' @description Estimates the long-run variance of a time series using the Newey-West estimator
#' (Bartlett weights \eqn{1 - h/(L+1)}, with L = \code{lrv_par}) on the demeaned series.
#'
#' @param X A TxK numeric matrix.
#' @param lrv_par An integer specifying the maximum lag order. If NULL, uses 0 (sample covariance).
#'
#' @return A KxK estimated long-run variance matrix.
#' @export
NeweyWest <- function(X, lrv_par = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)

  if (Tobs < 2) stop("Not enough observations.")
  if (is.null(lrv_par)) lrv_par <- 0
  if (lrv_par < 0 || lrv_par >= Tobs) stop("lrv_par must be between 0 and Tobs - 1.")

  # Center the data
  X <- scale(X, scale = FALSE)

  # Sample covariance
  OmegaHat <- t(X) %*% X / Tobs

  if (lrv_par > 0) {
    for (h in 1:lrv_par) {
      X1 <- X[(h + 1):Tobs, , drop = FALSE]
      X2 <- X[1:(Tobs - h), , drop = FALSE]
      gamma_h <- (t(X1) %*% X2 + t(X2) %*% X1) / Tobs
      weight <- 1 - h / (lrv_par + 1)
      OmegaHat <- OmegaHat + weight * gamma_h
    }
  }
  return(OmegaHat)
}
