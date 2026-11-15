# =============================================================================
# Selective inference for Panel Kmeans
#
# Notation follows Section 4 of the paper and Section S.4 of the supplement:
#   X = sqrt(T) * Zbar (NP vector, unit-major), Lambda_j the cosine transforms of the unit-demeaned
#   series (B x NP), L an NP x q restriction matrix, Sigma = L' Xi L, W = Xi L (NP x q).
#   D = || Sigma^{-1/2} L'X ||, dir = Sigma^{-1/2} L'X / D,
#   U = X - W Sigma^{-1} L'X,  V = W Sigma^{-1/2} dir,  X(phi) = U + phi V.
# The truncation set is {phi >= 0 : s_hat[X(phi)] = s_hat(X)}, where the selection outcome s_hat
# contains the paths of all initializations for all K, the retained initializations and K_IC.
# =============================================================================

#' Fit Panel Kmeans and Prepare the Quantities for Selective Inference
#'
#' @description Runs Panel Kmeans (with the selection of K by the IC if \code{Kmax} is given) and
#' computes the cosine transforms \eqn{\hat\Lambda_j} of the unit-level series, which give
#' \eqn{\hat\Xi = B^{-1}\sum_j\hat\Lambda_j\hat\Lambda_j'} without forming it.
#'
#' @param Z NT x P matrix.
#' @param id,time NT vectors of unit and time identifiers.
#' @param K Fixed number of clusters (alternative to \code{Kmax}).
#' @param Kmax Largest K; K is then chosen by the IC from \code{Kmin:Kmax}.
#' @param Kmin Smallest K (default 2).
#' @param Ninit Number of random initializations per K (S in the paper).
#' @param iter.max Maximum number of Lloyd iterations.
#' @param kappa IC tuning constant (\eqn{\varsigma}).
#' @param lrv_par Number of cosines B. If NULL, \code{choose_B} is applied to the cross-sectional
#'   average with lower bound \eqn{K_{max}P+4}.
#' @param cB Optional constant of the rule \eqn{B=\lceil c_BT^{2/3}\rceil}.
#' @param demean "global" (default) or "local" demeaning for the long-run variance.
#' @param bw Bandwidth for local demeaning.
#' @param n_cores Cores for the initializations.
#' @param init_type "forgy" (default) or "partition"; see \code{pk_cluster}.
#'
#' @return An object of class \code{pk_fit}.
#' @export
panel_selective_fit <- function(Z, id, time, K = NULL, Kmax = NULL, Kmin = 2,
                                Ninit = 10, iter.max = 100, kappa = 1.5,
                                lrv_par = NULL, cB = NULL, demean = "global", bw = NULL,
                                n_cores = 1, init_type = "forgy") {
  if (!is.null(K) && !is.null(Kmax)) stop("Specify only one of K or Kmax, not both.")
  if (is.null(K) && is.null(Kmax)) stop("Either K or Kmax must be provided.")
  Z <- as.matrix(Z)
  pd <- panel_prepare(Z, id, time)
  Kset <- if (!is.null(Kmax)) Kmin:Kmax else K
  sel <- pk_cluster(pd$xbar, Kset, Ninit = Ninit, iter.max = iter.max, Gamma_w = pd$Gamma_w,
                    Tobs = pd$Tobs, kappa = kappa, n_cores = n_cores, init_type = init_type)
  if (is.na(sel$Khat)) stop("All initializations failed for every K.")
  N <- pd$N; P <- pd$P; Tobs <- pd$Tobs
  Zo <- cross_section_average(pd)
  if (is.null(lrv_par)) lrv_par <- choose_B(Zo, cB = cB, Bmin = max(Kset) * P + 4)
  Lambda <- cosine_transforms(pd$Zt, lrv_par, demean = demean, bw = bw)
  structure(list(panel = pd, selection = sel, Lambda = Lambda, B = lrv_par,
                 x = sqrt(Tobs) * as.vector(t(pd$xbar)),
                 clustering = setNames(sel$final_cluster, as.character(pd$ids)),
                 K = sel$Khat, Kset = Kset, demean = demean),
            class = "pk_fit")
}

#' Cross-Sectional Average of a Prepared Panel
#'
#' @description Averages the P series over units at each time period, \eqn{\bar Z_t = N^{-1}\sum_i Z_{it}}.
#'
#' @param pd Output of \code{panel_prepare}.
#' @return T x P matrix.
#' @keywords internal
cross_section_average <- function(pd) {
  sapply(1:pd$P, function(p) rowMeans(pd$Zt[, (0:(pd$N - 1)) * pd$P + p, drop = FALSE]))
}

#' Weight Matrix of Cluster Means
#'
#' @description \eqn{A(\mathcal C) = [a_1\otimes I_P,\dots,a_K\otimes I_P]} with
#' \eqn{(a_k)_i = 1\{k_i=k\}/|\mathcal C_k|}, so that \eqn{A'X} stacks the K cluster means.
#'
#' @param gamma Cluster assignment (length N, labels 1..K).
#' @param K Number of clusters.
#' @param P Dimension of Z.
#' @return NP x KP matrix.
#' @export
A_matrix <- function(gamma, K, P) {
  a <- sapply(1:K, function(k) (gamma == k) / sum(gamma == k))
  if (!is.matrix(a)) a <- matrix(a, ncol = K)
  kronecker(a, diag(P))
}

#' Restriction Matrices for Selective Tests
#'
#' @description \code{L_pairwise}: \eqn{\delta_{k,g}\otimes I_P} (pairwise equality of clusters k and g);
#' \code{L_cluster}: \eqn{a_k\otimes I_P} (cluster k has zero mean);
#' \code{L_joint}: \eqn{A(\hat{\mathcal C})\mathsf R'} with \eqn{\mathsf R = I_{KP}} (\code{"cepa"},
#' C-EPA for the selected clusters) or the K-1 differences \eqn{(e_k-e_K)'\otimes I_P}
#' (\code{"homogeneity"}).
#'
#' @param fit A \code{pk_fit} object.
#' @param k,g Cluster labels.
#' @param type "cepa" or "homogeneity".
#' @return An NP x q matrix.
#' @export
L_pairwise <- function(fit, k, g) {
  P <- fit$panel$P
  A <- A_matrix(fit$selection$final_cluster, fit$K, P)
  A[, (k - 1) * P + 1:P, drop = FALSE] - A[, (g - 1) * P + 1:P, drop = FALSE]
}

#' @rdname L_pairwise
#' @export
L_cluster <- function(fit, k) {
  P <- fit$panel$P
  A <- A_matrix(fit$selection$final_cluster, fit$K, P)
  A[, (k - 1) * P + 1:P, drop = FALSE]
}

#' @rdname L_pairwise
#' @export
L_joint <- function(fit, type = c("cepa", "homogeneity")) {
  type <- match.arg(type)
  P <- fit$panel$P
  K <- fit$K
  A <- A_matrix(fit$selection$final_cluster, K, P)
  if (type == "cepa") return(A)
  do.call(cbind, lapply(1:(K - 1), function(k)
    A[, (k - 1) * P + 1:P, drop = FALSE] - A[, (K - 1) * P + 1:P, drop = FALSE]))
}

# Symmetric inverse square root
#' Symmetric Inverse Square Root
#'
#' @description Computes \eqn{S^{-1/2}} from the eigendecomposition of the symmetric matrix
#' \code{S}; stops with an error if \code{S} is not positive definite.
#'
#' @param S Symmetric matrix.
#'
#' @return A matrix of the same dimension as \code{S}.
#' @keywords internal
.inv_sqrt <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  if (any(e$values <= 0)) stop("The long-run variance of the tested quantity is not positive definite; increase B.")
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

#' Selective p-value for a Linear Hypothesis on Cluster Means
#'
#' @description Computes the selective p-value of \eqn{L'\mu^0 = 0} (Lemma 3, Remark 7 and
#' Algorithm S.1 of the paper): the truncated \eqn{\chi_q} survival function at the observed
#' statistic D, with truncation set \eqn{\mathcal T} along the perturbation path
#' \eqn{X(\phi) = u + \phi v}.
#'
#' With \code{conditioning = "partition"} (default), \eqn{\mathcal T} is the set of \eqn{\phi\ge0} at
#' which the clustering procedure, run with the same (data-independent) initializations, selects the
#' observed partition with \eqn{\hat K_{IC}} clusters. Validity holds conditionally on the
#' initializations, hence uniformly in them; the paths of the initializations are not conditioned on.
#' The set is computed exactly by walking along \eqn{\phi} through the segments on which all paths are
#' constant (each obtained from the quadratic and IC constraints of \code{truncation_set}); the walk
#' stops early, conservatively, once the p-value is determined to the requested accuracy.
#' With \code{conditioning = "paths"}, \eqn{\mathcal T} is the (smaller) set on which all paths, the
#' retained initializations and \eqn{\hat K_{IC}} are unchanged.
#'
#' @param fit A \code{pk_fit} object.
#' @param L NP x q restriction matrix (e.g. from \code{L_pairwise}).
#' @param metric How \eqn{v} is built: "xi" (default) uses \eqn{W = \hat\Xi L}, the decomposition
#'   in the metric of the long-run covariance; "euclidean" uses \eqn{L(L'L)^{-1}\hat\Sigma} (Chen and
#'   Witten, 2023; valid only if the column space of L is invariant under \eqn{\Xi}); "shrink" uses
#'   \eqn{(1-w)\hat\Xi L + w L(L'L)^{-1}\hat\Sigma} with a Ledoit-Wolf type intensity w. All satisfy
#'   \eqn{L'W = \hat\Sigma}. "xi" is the default of the paper; "shrink" is the alternative
#'   procedure reported in the Monte Carlo study and the application (it raises power, with size
#'   close to the nominal level in the designs of the paper); "euclidean" is for comparisons.
#' @param conditioning "partition" (default) or "paths".
#' @param condition_K,condition_starts Only for \code{conditioning = "paths"}: condition on
#'   \eqn{\hat K_{IC}} and the paths for all K, and on all initializations (defaults TRUE; FALSE is
#'   not valid when K or the initialization is selected, and is kept for comparisons).
#' @param ref Reference distribution of \eqn{D^2}: "chi" (default, \eqn{\chi^2_q}) or "F", the
#'   Hotelling law \eqn{qB/(B-q+1)F_{q,B-q+1}} (its finite-B analogue). In pilot simulations the F
#'   reference made the joint selective tests very conservative, so it is not the default.
#' @param eps Relative numerical tolerance.
#' @param tol Accuracy of the p-value in the partition walk: the walk stops when the p-value is below
#'   \code{tol} or its upper and lower bounds differ by less than \code{tol} in relative terms;
#'   the upper bound is reported.
#' @param max_seg Maximum number of segments in the walk.
#'
#' @return A list with \code{stat} (D), \code{df}, \code{pval}, the truncation set \code{intervals},
#'   \code{in_set}, \code{w} (shrinkage intensity) and \code{nseg} (segments visited).
#' @export
selective_test <- function(fit, L, metric = c("xi", "shrink", "euclidean"),
                           conditioning = c("partition", "paths"), ref = c("chi", "F"),
                           condition_K = TRUE, condition_starts = TRUE, eps = 1e-9, tol = 1e-3,
                           max_seg = 50000) {
  metric <- match.arg(metric)
  conditioning <- match.arg(conditioning)
  ref <- match.arg(ref)
  L <- as.matrix(L)
  q <- ncol(L)
  X <- fit$x
  Lam <- fit$Lambda
  B <- nrow(Lam)
  G <- Lam %*% L                       # B x q, cosine transforms of the tested quantity
  Sigma <- crossprod(G) / B            # q x q
  y <- as.vector(crossprod(L, X))      # L'X
  W_xi <- crossprod(Lam, G) / B
  W_eu <- L %*% solve(crossprod(L), Sigma)
  w <- switch(metric, xi = 0, euclidean = 1, shrink = {
    xi_diag <- colSums(Lam^2) / B
    num <- sum(outer(xi_diag, diag(Sigma)) + W_xi^2) / B
    den <- sum((W_xi - W_eu)^2)
    if (den > 0) min(1, num / den) else 1
  })
  W <- (1 - w) * W_xi + w * W_eu
  Bref <- if (ref == "F") B else NULL
  if (!is.null(Bref) && B - q + 1 < 1) stop("B must be at least q for the F reference.")
  Sih <- .inv_sqrt(Sigma)
  ws <- as.vector(Sih %*% y)
  D <- sqrt(sum(ws^2))
  dirv <- ws / D
  u <- as.vector(X - W %*% solve(Sigma, y))
  v <- as.vector(W %*% (Sih %*% dirv))

  nseg <- NA_integer_
  if (conditioning == "paths") {
    Tset <- truncation_set(fit, u, v, D, condition_K = condition_K,
                           condition_starts = condition_starts, eps = eps)
    pval <- truncated_chi_pvalue(D, q, Tset, Bref)
  } else {
    pw <- .partition_walk(fit, u, v, D, q, eps = eps, tol = tol, max_seg = max_seg, B = Bref)
    Tset <- pw$intervals; pval <- pw$pval; nseg <- pw$nseg
  }
  in_set <- nrow(Tset) > 0 && any(Tset[, 1] <= D * (1 + 1e-7) + 1e-12 & Tset[, 2] >= D * (1 - 1e-7) - 1e-12)
  if (!in_set) warning("The observed statistic is not in the computed truncation set (numerical tolerance).")
  list(stat = D, df = q, pval = pval, intervals = Tset, in_set = in_set, w = w, nseg = nseg)
}

# log P(a^2 <= R <= b^2), vectorized, where R ~ chi^2_q (B = NULL) or R ~ q B/(B-q+1) F(q, B-q+1)
# (Hotelling reference, B = number of cosines of the OS estimator)
#' Log Probability of an Interval on the Chi Scale
#'
#' @description Computes \eqn{\log\Pr(a^2 \le R \le b^2)}, vectorized, where \eqn{R\sim\chi^2_q}
#' (\code{B = NULL}) or R follows the Hotelling law \eqn{qB/(B-q+1)F_{q,B-q+1}}. Upper or lower tail
#' probabilities are used depending on the position of \eqn{a^2} relative to the median, for
#' accuracy in the tails. Returns \code{-Inf} where \eqn{b \le a}.
#'
#' @param a,b Lower and upper ends (chi scale).
#' @param q Degrees of freedom.
#' @param B Optional number of cosines of the OS estimator (Hotelling reference).
#'
#' @return A numeric vector of log probabilities.
#' @keywords internal
.logmass <- function(a, b, q, B = NULL) {
  if (is.null(B)) {
    cdf <- function(x, lower, log) pchisq(x, q, lower.tail = lower, log.p = log)
    med <- qchisq(0.5, q)
  } else {
    c0 <- (B - q + 1) / (q * B)
    cdf <- function(x, lower, log) pf(c0 * x, q, B - q + 1, lower.tail = lower, log.p = log)
    med <- qf(0.5, q, B - q + 1) / c0
  }
  a2 <- a^2; b2 <- b^2
  useS <- a2 >= med
  la <- cdf(a2, FALSE, TRUE); lb <- cdf(b2, FALSE, TRUE)
  Fa <- cdf(a2, TRUE, TRUE); Fb <- cdf(b2, TRUE, TRUE)
  out <- ifelse(useS, la + log1p(-exp(lb - la)), Fb + log1p(-exp(Fa - Fb)))
  out[b <= a] <- -Inf
  out
}
#' Log-Sum-Exp
#'
#' @description Computes \eqn{\log\sum_i \exp(x_i)} stably; returns \code{-Inf} if the maximum of
#' \code{x} is not finite.
#'
#' @param x Numeric vector.
#'
#' @return A numeric scalar.
#' @keywords internal
.lse <- function(x) { m <- max(x); if (!is.finite(m)) return(-Inf); m + log(sum(exp(x - m))) }

# Canonical key of a selection outcome: K_hat and the selected partition (up to relabeling)
#' Canonical Key of a Selection Outcome
#'
#' @description Pastes \eqn{\hat K} and the key of the selected partition (see
#' \code{.partition_key}), so that outcomes differing only by a relabeling share the same key.
#'
#' @param sel A \code{pk_selection} object.
#'
#' @return A character scalar.
#' @keywords internal
.sel_key <- function(sel) paste(sel$Khat, .partition_key(sel$final_cluster))

# ---- constraint coefficients -------------------------------------------------

# Assignment constraints: for each unit i assigned to k = anew[i] given the groups gprev,
#   ||x_i(phi) - xbar_k(phi)||^2 - ||x_i(phi) - xbar_g(phi)||^2 <= 0, g != k.
# Returns a matrix with columns (A, B, C) of A phi^2 + B phi + C.
#' Coefficients of the Assignment Constraints
#'
#' @description For each unit i assigned to cluster \code{anew[i]} given the groups \code{gprev},
#' and each other cluster g, returns the coefficients of the quadratic in \eqn{\phi}
#' \eqn{\|x_i(\phi) - \bar x_{k}(\phi)\|^2 - \|x_i(\phi) - \bar x_g(\phi)\|^2}, which must be
#' non-positive for the assignment to be reproduced at \eqn{X(\phi) = U + \phi V}.
#'
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param gprev Assignment used to compute the centers (labels in 0..K).
#' @param anew New assignment (labels in 1..K).
#' @param K Number of clusters.
#'
#' @return An N(K-1) x 3 matrix with columns (A, B, C) of \eqn{A\phi^2 + B\phi + C}.
#' @keywords internal
.assign_coefs <- function(U, V, gprev, anew, K) {
  N <- nrow(U)
  Mg <- .onehot(gprev, K)
  cnt <- tabulate(gprev, K)
  Ub <- crossprod(Mg, U) / cnt
  Vb <- crossprod(Mg, V) / cnt
  fA <- matrix(rowSums(Vb^2), N, K, byrow = TRUE) - 2 * V %*% t(Vb)
  fC <- matrix(rowSums(Ub^2), N, K, byrow = TRUE) - 2 * U %*% t(Ub)
  fB <- matrix(rowSums(Ub * Vb), N, K, byrow = TRUE) - U %*% t(Vb) - V %*% t(Ub)
  idx <- cbind(1:N, anew)
  mask <- col(fA) != anew
  cbind((fA[idx] - fA)[mask], 2 * (fB[idx] - fB)[mask], (fC[idx] - fC)[mask])
}

# WSS(phi; g) = c1 phi^2 + c2 phi + c3
#' Within-Cluster Sum of Squares as a Quadratic in phi
#'
#' @description Coefficients of \eqn{WSS(\phi; g) = c_1\phi^2 + c_2\phi + c_3} for the partition
#' \code{g} evaluated at \eqn{X(\phi) = U + \phi V}.
#'
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param g Assignment (labels in 1..K).
#' @param K Number of clusters.
#'
#' @return A numeric vector \code{c(c1, c2, c3)}.
#' @keywords internal
.wss_coefs <- function(U, V, g, K) {
  Mg <- .onehot(g, K)
  cnt <- tabulate(g, K)
  Ru <- U - (crossprod(Mg, U) / cnt)[g, , drop = FALSE]
  Rv <- V - (crossprod(Mg, V) / cnt)[g, , drop = FALSE]
  c(sum(Rv^2), 2 * sum(Ru * Rv), sum(Ru^2))
}

# Forbidden sets {phi >= 0: A phi^2 + B phi + C > tol} of a batch of quadratics.
# Returns a two-column matrix of intervals (lo, hi).
#' Forbidden Sets of a Batch of Quadratic Constraints
#'
#' @description For each quadratic \eqn{A\phi^2 + B\phi + C}, computes the set of \eqn{\phi \ge 0}
#' where it exceeds a tolerance relative to the scale of its coefficients at \eqn{\phi = D};
#' constant, linear and quadratic cases are handled separately.
#'
#' @param A,B,C Numeric vectors of coefficients (same length).
#' @param D Observed statistic, used for the scale of the tolerance.
#' @param eps Relative tolerance.
#'
#' @return A two-column matrix of (possibly overlapping) intervals (lo, hi).
#' @keywords internal
.quad_forbidden <- function(A, B, C, D, eps = 1e-9) {
  if (!length(A)) return(matrix(numeric(0), 0, 2))
  D2 <- max(D^2, 1); D1 <- max(D, 1)
  sc <- abs(A) * D2 + abs(B) * D1 + abs(C)
  sc[sc == 0] <- 1
  C <- C - eps * sc
  zA <- abs(A) * D2 <= 1e-12 * sc
  zB <- abs(B) * D1 <= 1e-12 * sc
  lo <- list(); hi <- list()

  # constant
  i <- zA & zB & C > 0
  if (any(i)) { lo$c <- rep(0, sum(i)); hi$c <- rep(Inf, sum(i)) }
  # linear
  i <- zA & !zB
  if (any(i)) {
    r <- -C[i] / B[i]; bp <- B[i] > 0
    lo$l1 <- pmax(r[bp], 0); hi$l1 <- rep(Inf, sum(bp))
    rn <- r[!bp]; rn <- rn[rn > 0]
    lo$l2 <- rep(0, length(rn)); hi$l2 <- rn
  }
  # quadratic
  i <- !zA
  if (any(i)) {
    a <- A[i]; b <- B[i]; c <- C[i]
    disc <- b^2 - 4 * a * c
    pos <- a > 0
    # A > 0, disc < 0: all forbidden
    j <- pos & disc < 0
    if (any(j)) { lo$q0 <- rep(0, sum(j)); hi$q0 <- rep(Inf, sum(j)) }
    j <- disc >= 0
    if (any(j)) {
      sq <- sqrt(disc[j])
      qq <- -0.5 * (b[j] + ifelse(b[j] >= 0, 1, -1) * sq)
      r_a <- qq / a[j]
      r_b <- ifelse(qq == 0, 0, c[j] / qq)
      r1 <- pmin(r_a, r_b); r2 <- pmax(r_a, r_b)
      pj <- pos[j]
      # A > 0: forbidden [0, r1) and (r2, Inf)
      k1 <- pj & r1 > 0
      lo$q1 <- rep(0, sum(k1)); hi$q1 <- r1[k1]
      lo$q2 <- pmax(r2[pj], 0); hi$q2 <- rep(Inf, sum(pj))
      # A < 0: forbidden (r1, r2)
      k2 <- !pj & r2 > 0
      lo$q3 <- pmax(r1[k2], 0); hi$q3 <- r2[k2]
    }
  }
  cbind(unlist(lo, use.names = FALSE), unlist(hi, use.names = FALSE))
}

# Union of intervals (two-column matrix) as sorted disjoint intervals
#' Union of Intervals
#'
#' @description Merges a set of intervals into sorted disjoint intervals; empty intervals
#' (upper end not above lower end) are dropped.
#'
#' @param M Two-column matrix of intervals, or \code{NULL}.
#'
#' @return A two-column matrix of sorted disjoint intervals (possibly with zero rows).
#' @keywords internal
.union_intervals <- function(M) {
  if (is.null(M) || !nrow(M)) return(matrix(numeric(0), 0, 2))
  keep <- M[, 2] > M[, 1]
  M <- M[keep, , drop = FALSE]
  if (!nrow(M)) return(matrix(numeric(0), 0, 2))
  M <- M[order(M[, 1]), , drop = FALSE]
  ch <- cummax(M[, 2])
  start <- c(TRUE, M[-1, 1] > ch[-nrow(M)])
  grp <- cumsum(start)
  cbind(M[start, 1], as.numeric(tapply(M[, 2], grp, max)))
}

# Complement in [0, Inf) of a union of sorted disjoint intervals
#' Complement of a Union of Intervals
#'
#' @description Computes the complement in \eqn{[0,\infty)} of a union of sorted disjoint intervals.
#'
#' @param F Two-column matrix of sorted disjoint intervals in \eqn{[0,\infty)}.
#'
#' @return A two-column matrix of intervals.
#' @keywords internal
.complement <- function(F) {
  starts <- c(0, F[, 2])
  ends <- c(F[, 1], Inf)
  keep <- ends > starts
  cbind(starts[keep], ends[keep])
}

# Forbidden set {phi >= 0: f(phi) > tol} for a polynomial f of degree <= deg, evaluated by fun
#' Forbidden Set of a Polynomial on the Half-Line
#'
#' @description Computes \eqn{\{\phi \ge 0: f(\phi) > tol\}} for a polynomial f of degree at most
#' \code{deg}: the coefficients are obtained by interpolation on equispaced nodes in
#' \eqn{[0, \max(2D, 1)]}, the positive real roots by \code{polyroot}, and the sign of f is checked
#' between consecutive roots. The tolerance is relative to the scale of f at the nodes and at D.
#'
#' @param fun Function evaluating f at a scalar.
#' @param deg Degree bound of f.
#' @param D Observed statistic, used for the interpolation range and the scale.
#' @param eps Relative tolerance.
#'
#' @return Two-column matrix of forbidden intervals.
#' @keywords internal
.poly_forbidden <- function(fun, deg, D, eps = 1e-9) {
  s <- max(2 * D, 1)
  nodes <- s * (0:deg) / deg
  fv <- vapply(nodes, fun, 0)
  scale <- max(abs(fv), abs(fun(D)), 1e-300)
  Vm <- outer(nodes / s, 0:deg, "^")
  cf <- tryCatch(solve(Vm, fv), error = function(e) NULL)
  roots <- numeric(0)
  if (!is.null(cf)) {
    cf[abs(cf) <= 1e-12 * max(abs(cf))] <- 0
    nz <- which(cf != 0)
    if (length(nz) && max(nz) > 1) {
      cf <- cf[1:max(nz)]
      rr <- polyroot(cf)
      rr <- Re(rr[abs(Im(rr)) <= 1e-7 * (1 + abs(Re(rr)))]) * s
      roots <- sort(unique(rr[rr > 0]))
    }
  }
  br <- c(0, roots, Inf)
  mids <- c((br[-length(br)] + br[-1]) / 2)
  mids[length(mids)] <- if (length(roots)) 2 * max(roots) + 1 else max(2 * D, 1)
  bad <- vapply(mids, fun, 0) > eps * scale
  cbind(br[-length(br)][bad], br[-1][bad])
}

# M_K(phi) for the retained partition g with K clusters
#' Matrix of the Information Criterion as a Quadratic in phi
#'
#' @description For the partition \code{g}, returns the coefficients of
#' \eqn{M_K(\phi) = \hat\Gamma_w + (NT)^{-1}R(\phi)'R(\phi) = M_0 + \phi M_1 + \phi^2 M_2}, where
#' \eqn{R(\phi)} are the within-cluster residuals of \eqn{X(\phi) = U + \phi V}; the IC is
#' \eqn{\log\det M_K(\phi)} plus the penalty.
#'
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param g Assignment (labels in 1..K).
#' @param K Number of clusters.
#' @param Gamma_w P x P within-unit second moment.
#' @param N Number of units.
#' @param Tobs Number of time periods.
#'
#' @return A list with P x P matrices \code{M0}, \code{M1}, \code{M2}.
#' @keywords internal
.M_poly <- function(U, V, g, K, Gamma_w, N, Tobs) {
  Mg <- .onehot(g, K)
  cnt <- tabulate(g, K)
  Ru <- U - (crossprod(Mg, U) / cnt)[g, , drop = FALSE]
  Rv <- V - (crossprod(Mg, V) / cnt)[g, , drop = FALSE]
  NT <- N * Tobs
  list(M0 = Gamma_w + crossprod(Ru) / NT,
       M1 = (crossprod(Ru, Rv) + crossprod(Rv, Ru)) / NT,
       M2 = crossprod(Rv) / NT)
}

#' Truncation Set of the Selective Test
#'
#' @description Intersects the sets of \eqn{\phi\ge 0} on which (i) every assignment of every
#' initialization of Algorithm 1 is reproduced (including the first assignment from the random
#' initial groups and the steps of discarded initializations), (ii) the retained initialization is
#' unchanged, and (iii) \eqn{\hat K_{IC}} is unchanged, along \eqn{X(\phi) = u + \phi v}.
#'
#' @param fit A \code{pk_fit} object.
#' @param u,v NP vectors.
#' @param D Observed statistic.
#' @param condition_K,condition_starts See \code{selective_test}.
#' @param eps Numerical tolerance.
#' @return Two-column matrix of disjoint intervals.
#' @export
truncation_set <- function(fit, u, v, D, condition_K = TRUE, condition_starts = TRUE, eps = 1e-9) {
  sel <- fit$selection
  N <- sel$N; P <- sel$P; Tobs <- sel$Tobs
  U <- matrix(u, N, P, byrow = TRUE)
  V <- matrix(v, N, P, byrow = TRUE)
  Ks <- if (condition_K) sel$Kset else sel$Khat
  forb <- matrix(numeric(0), 0, 2)

  for (K in Ks) {
    bk <- sel$byK[[as.character(K)]]
    if (condition_starts) {
      ridx <- seq_along(bk$runs)
    } else {
      if (is.na(bk$best)) next
      ridx <- bk$best
    }
    co <- vector("list", 0)
    for (s in ridx) {
      run <- bk$runs[[s]]
      prev <- run$init
      for (m in seq_along(run$path)) {
        co[[length(co) + 1]] <- .assign_coefs(U, V, prev, run$path[[m]], K)
        prev <- run$path[[m]]
      }
    }
    if (condition_starts && !is.na(bk$best) && length(bk$first_idx) > 1) {
      cb <- .wss_coefs(U, V, bk$runs[[bk$best]]$final, K)
      for (s in setdiff(bk$first_idx, bk$best)) {
        co[[length(co) + 1]] <- matrix(cb - .wss_coefs(U, V, bk$runs[[s]]$final, K), 1, 3)
      }
    }
    if (length(co)) {
      Cm <- do.call(rbind, co)
      forb <- .union_intervals(rbind(forb, .quad_forbidden(Cm[, 1], Cm[, 2], Cm[, 3], D, eps)))
    }
  }

  if (condition_K && length(sel$Kset) > 1) {
    pen <- function(K) (K * P + N) * sel$kappa * log(N * Tobs) / (N * Tobs)
    Kh <- sel$Khat
    bh <- sel$byK[[as.character(Kh)]]
    Mh <- .M_poly(U, V, bh$runs[[bh$best]]$final, Kh, sel$Gamma_w, N, Tobs)
    detM <- function(M, phi) det(M$M0 + phi * M$M1 + phi^2 * M$M2)
    f0 <- exp(pen(Kh)) * detM(Mh, D)
    for (K in setdiff(sel$Kset, Kh)) {
      bk <- sel$byK[[as.character(K)]]
      if (is.na(bk$best)) next
      MK <- .M_poly(U, V, bk$runs[[bk$best]]$final, K, sel$Gamma_w, N, Tobs)
      fun <- function(phi) (exp(pen(Kh)) * detM(Mh, phi) - exp(pen(K)) * detM(MK, phi)) / f0
      hi_r <- 4 * D + 10
      fb <- .poly_forbidden_piece(fun, 2 * P, 0, hi_r, eps)
      if (fun(2 * hi_r) > 0) fb <- rbind(fb, c(hi_r, Inf))
      forb <- .union_intervals(rbind(forb, fb))
    }
  }
  .complement(forb)
}

#' Survival Function of a Truncated Chi Distribution
#'
#' @description \eqn{\Pr(\chi_q \ge d \mid \chi_q \in \mathcal T)}, computed on the log scale.
#' With \code{B} given, the reference for \eqn{d^2} is the Hotelling law
#' \eqn{qB/(B-q+1)\,F_{q,B-q+1}}, which accounts for the estimation of the long-run variance by the OS
#' estimator with B cosines and converges to \eqn{\chi^2_q} as \eqn{B\to\infty}.
#' @param d Observed statistic (chi scale).
#' @param q Degrees of freedom.
#' @param intervals Two-column matrix of disjoint intervals (chi scale).
#' @param B Optional number of cosines (Hotelling reference).
#' @return The p-value (NA if the truncation set has zero probability).
#' @export
truncated_chi_pvalue <- function(d, q, intervals, B = NULL) {
  if (!nrow(intervals)) return(NA_real_)
  lden <- .lse(.logmass(intervals[, 1], intervals[, 2], q, B))
  up <- intervals[, 2] > d
  if (!any(up)) return(0)
  lnum <- .lse(.logmass(pmax(intervals[up, 1], d), intervals[up, 2], q, B))
  if (!is.finite(lden)) return(NA_real_)
  min(max(exp(lnum - lden), 0), 1)
}

#' Selective Tests of Pairwise Equality, Joint Homogeneity and Joint C-EPA
#'
#' @description Computes, for a fitted selection, the selective p-values of all pairwise equality
#' hypotheses, their M-family merging (\eqn{F_{homo,r}}), and the joint selective tests of
#' Remark 7 (homogeneity, \eqn{\mathsf R} = differences; C-EPA, \eqn{\mathsf R = I_{KP}}).
#'
#' @param fit A \code{pk_fit} object.
#' @param pairs Optional two-column matrix of cluster pairs (default all pairs).
#' @param pcombine_fun Merging function: "Genmean_neq" (default, r < -1), "Genmean", "Geomean",
#'   "iu", "bonferroni", "cauchy".
#' @param r Order of the generalized mean (default -20).
#' @param joint Compute the joint tests (default TRUE).
#' @param metric,conditioning,condition_K,condition_starts,ref See \code{selective_test}.
#'
#' @return A list with \code{pairwise_pvalues}, \code{pairs}, \code{pairwise_stats}, \code{hom_pval}
#'   (merged), \code{joint_hom_pval}, \code{joint_cepa_pval}, \code{clustering}, \code{K}, \code{B},
#'   and \code{in_set} (diagnostic).
#' @export
panel_selective_tests <- function(fit, pairs = NULL, pcombine_fun = "Genmean_neq", r = -20,
                                  joint = TRUE, metric = "xi", conditioning = "partition",
                                  condition_K = TRUE, condition_starts = TRUE, ref = "chi") {
  K <- fit$K
  pair_idx <- if (is.null(pairs)) t(combn(K, 2)) else as.matrix(pairs)
  np <- nrow(pair_idx)
  pv <- st <- numeric(np)
  ins <- logical(np)
  for (j in seq_len(np)) {
    res <- selective_test(fit, L_pairwise(fit, pair_idx[j, 1], pair_idx[j, 2]), metric = metric,
                          conditioning = conditioning, condition_K = condition_K,
                          condition_starts = condition_starts, ref = ref)
    pv[j] <- res$pval; st[j] <- res$stat; ins[j] <- res$in_set
  }
  names(pv) <- paste0("pair_", pair_idx[, 1], "_", pair_idx[, 2])
  out <- list(pairwise_pvalues = pv, pairs = split(pair_idx, seq_len(np)), pairwise_stats = st,
              hom_pval = merge_pvalues(pv, pcombine_fun, r),
              clustering = fit$clustering, K = K, B = fit$B, in_set = ins)
  if (joint) {
    jh <- selective_test(fit, L_joint(fit, "homogeneity"), metric = metric,
                         conditioning = conditioning, condition_K = condition_K,
                         condition_starts = condition_starts, ref = ref)
    jc <- selective_test(fit, L_joint(fit, "cepa"), metric = metric,
                         conditioning = conditioning, condition_K = condition_K,
                         condition_starts = condition_starts, ref = ref)
    out$joint_hom_pval <- jh$pval
    out$joint_cepa_pval <- jc$pval
    out$in_set <- c(out$in_set, jh$in_set, jc$in_set)
  }
  out
}
