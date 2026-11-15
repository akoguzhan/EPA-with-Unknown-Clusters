#' Prepare a Balanced Panel for Panel Kmeans and Selective Inference
#'
#' @description Sorts the panel by unit and time and computes the quantities on which Panel Kmeans
#' and the selective tests depend: the unit time averages \eqn{\bar Z_i}, the T x NP matrix of the
#' stacked series \eqn{Z_t = (Z_{1t}',\dots,Z_{Nt}')'} and the within-unit second moment
#' \eqn{\hat\Gamma_w = (NT)^{-1}\sum_{i,t}(Z_{it}-\bar Z_i)(Z_{it}-\bar Z_i)'}.
#'
#' @param Z NT x P matrix.
#' @param id NT vector of unit identifiers.
#' @param time NT vector of time identifiers.
#'
#' @return A list with \code{Zt} (T x NP, column \eqn{(i-1)P+p}), \code{xbar} (N x P matrix of time
#' averages), \code{Gamma_w}, \code{ids} (sorted unit identifiers), \code{times}, \code{N}, \code{Tobs}, \code{P}.
#' @export
panel_prepare <- function(Z, id, time) {
  Z <- as.matrix(Z)
  if (length(id) != nrow(Z) || length(time) != nrow(Z)) stop("id and time must have nrow(Z) elements.")
  if (anyNA(Z)) stop("Z contains missing values.")
  P <- ncol(Z)
  ord <- order(id, time)
  Z <- Z[ord, , drop = FALSE]
  id <- id[ord]
  time <- time[ord]
  ids <- unique(id)
  times <- sort(unique(time))
  N <- length(ids)
  Tobs <- length(times)
  if (N * Tobs != nrow(Z)) stop("The panel must be balanced.")
  if (!all(time == rep(times, times = N))) stop("The panel must be balanced (same time periods for all units).")

  Zt <- matrix(NA_real_, Tobs, N * P)
  for (p in 1:P) Zt[, (0:(N - 1)) * P + p] <- matrix(Z[, p], nrow = Tobs, ncol = N)
  mu <- colMeans(Zt)
  xbar <- matrix(mu, nrow = N, ncol = P, byrow = TRUE)
  Zc <- sweep(Zt, 2, mu)
  arr <- array(Zc, dim = c(Tobs, P, N))
  M <- matrix(aperm(arr, c(1, 3, 2)), Tobs * N, P)
  Gamma_w <- crossprod(M) / (N * Tobs)

  list(Zt = Zt, xbar = xbar, Gamma_w = Gamma_w, ids = ids, times = times,
       N = N, Tobs = Tobs, P = P)
}

# Canonical key of a partition (invariant to relabeling)
#' Canonical Key of a Partition
#'
#' @description Relabels the clusters in order of first appearance and collapses the labels into a
#' comma-separated string, so that two assignments that differ only by a permutation of the labels
#' share the same key.
#'
#' @param g Integer vector of cluster labels (length N).
#'
#' @return A character scalar.
#' @keywords internal
.partition_key <- function(g) paste(match(g, unique(g)), collapse = ",")

# Cluster centers (K x P) of the rows of x for labels g in 1..K (all present)
#' Cluster Centers
#'
#' @description Computes the mean of the rows of \code{x} within each cluster. Rows with label 0
#' (units not used as initial centers under Forgy initialization) are ignored.
#'
#' @param x N x P numeric matrix.
#' @param g Integer vector of labels in 0..K (length N).
#' @param K Number of clusters.
#'
#' @return A K x P matrix of cluster means.
#' @keywords internal
.centers <- function(x, g, K) {
  # label 0: unit not used to compute the initial centers (Forgy initialization)
  crossprod(.onehot(g, K), x) / tabulate(g, K)
}

# N x K indicator matrix of labels g in 1..K (rows with label 0 are zero)
#' Indicator Matrix of Cluster Labels
#'
#' @description Builds the N x K matrix whose (i, k) entry is 1 if unit i has label k and 0
#' otherwise; rows of units with label 0 are zero.
#'
#' @param g Integer vector of labels in 0..K (length N).
#' @param K Number of clusters.
#'
#' @return An N x K numeric matrix.
#' @keywords internal
.onehot <- function(g, K) {
  M <- matrix(0, length(g), K)
  i <- which(g > 0)
  M[cbind(i, g[i])] <- 1
  M
}

# Within-cluster sum of squares
#' Within-Cluster Sum of Squares
#'
#' @description Sum over units of the squared Euclidean distance between each row of \code{x} and
#' the center of its cluster.
#'
#' @param x N x P numeric matrix.
#' @param g Integer vector of labels in 1..K (length N).
#' @param K Number of clusters.
#'
#' @return A numeric scalar.
#' @keywords internal
.wss <- function(x, g, K) {
  sum((x - .centers(x, g, K)[g, , drop = FALSE])^2)
}

#' One Run of Panel Kmeans (Algorithm 1) on Unit Time Averages
#'
#' @description Lloyd iterations from a random initial assignment. Ties are broken in favour of the
#' smallest cluster label. The run stops when the assignment does not change, after \code{iter.max}
#' assignments, or when a cluster becomes empty (the run is then discarded, \code{valid = FALSE}).
#'
#' @param x N x P matrix of time averages.
#' @param K Number of clusters.
#' @param g0 Initial assignment (integer vector in 1..K with all clusters present).
#' @param iter.max Maximum number of assignment steps.
#'
#' @return A list with \code{init}, \code{path} (list of assignments), \code{valid}, \code{final}.
#' @keywords internal
.lloyd_path <- function(x, K, g0, iter.max) {
  x2 <- rowSums(x^2)
  path <- vector("list", iter.max)
  g <- g0
  m <- 0
  valid <- TRUE
  while (m < iter.max) {
    ctr <- .centers(x, g, K)
    d <- outer(x2, rowSums(ctr^2), "+") - 2 * x %*% t(ctr)
    gnew <- max.col(-d, ties.method = "first")
    m <- m + 1
    path[[m]] <- gnew
    if (any(tabulate(gnew, K) == 0)) { valid <- FALSE; break }
    if (all(gnew == g)) break
    g <- gnew
  }
  list(init = g0, path = path[seq_len(m)], valid = valid, final = path[[m]])
}

# Random initial assignments, all K clusters present, drawn independently of the data
# "partition": every unit assigned at random to a cluster (all clusters present);
# "forgy": K units drawn at random serve as initial centers (label 0 for the other units).
#' Random Initial Assignments for Panel Kmeans
#'
#' @description Draws \code{Ninit} initial assignments independently of the data. With
#' \code{"forgy"}, K randomly chosen units receive labels 1..K and all other units label 0 (the K
#' units serve as initial centers); with \code{"partition"}, every unit is assigned at random so
#' that all K clusters are present.
#'
#' @param N Number of units.
#' @param K Number of clusters.
#' @param Ninit Number of initial assignments.
#' @param init_type \code{"forgy"} (default) or \code{"partition"}.
#'
#' @return A list of \code{Ninit} integer vectors of length N.
#' @keywords internal
.random_inits <- function(N, K, Ninit, init_type = "forgy") {
  if (init_type == "forgy") {
    return(lapply(seq_len(Ninit), function(s) { g <- integer(N); g[sample.int(N, K)] <- 1:K; g }))
  }
  lapply(seq_len(Ninit), function(s) sample(rep_len(1:K, N)))
}

#' Panel Kmeans with Multiple Initializations and Choice of K (Selection Outcome)
#'
#' @description Runs Algorithm 1 from \code{Ninit} random initial assignments for every
#' \eqn{K \in} \code{Kset}, retains for each K the first initialization attaining the smallest
#' within-cluster sum of squares (identical partitions are compared once), and chooses
#' \eqn{\hat K_{IC}} by the information criterion
#' \deqn{IC(K) = \log\det[\hat\Gamma_w + N^{-1}\sum_i(\bar Z_i-\hat\theta_{\hat k_i})(\bar Z_i-\hat\theta_{\hat k_i})'] + (KP+N)\varsigma\log(NT)/(NT).}
#' The output is the selection outcome \eqn{\hat s} on which the selective tests condition.
#'
#' @param xbar N x P matrix of unit time averages.
#' @param Kset Integer vector of candidate numbers of clusters (a single value fixes K).
#' @param Ninit Number of random initializations per K.
#' @param iter.max Maximum number of assignment steps.
#' @param Gamma_w P x P within-unit second moment.
#' @param Tobs Number of time periods.
#' @param kappa Tuning constant \eqn{\varsigma} of the IC penalty.
#' @param inits Optional list (named by K) of lists of initial assignments.
#' @param n_cores Number of cores (forked processes; ignored on Windows).
#' @param init_type "forgy" (default: K randomly chosen units as initial centers) or "partition"
#'   (random assignment of every unit). Both are drawn independently of the data. With
#'   "partition" and well-separated clusters, the initial centers nearly coincide and initial
#'   assignments often leave a cluster empty.
#' @param max_extra If all \code{Ninit} initializations for some K produce an empty cluster,
#'   further initializations are drawn until one is valid, at most \code{max_extra}.
#'
#' @return An object of class \code{pk_selection}.
#' @export
pk_cluster <- function(xbar, Kset, Ninit = 10, iter.max = 100, Gamma_w, Tobs, kappa = 1.5,
                       inits = NULL, n_cores = 1, max_extra = 100, init_type = "forgy") {
  N <- nrow(xbar)
  P <- ncol(xbar)
  byK <- list()
  IC <- rep(NA_real_, length(Kset))
  names(IC) <- as.character(Kset)
  pen <- function(K) (K * P + N) * kappa * log(N * Tobs) / (N * Tobs)

  for (K in Kset) {
    if (K < 1 || K > N) stop("Invalid K.")
    ini_all <- if (!is.null(inits)) inits[[as.character(K)]] else .random_inits(N, K, Ninit, init_type)
    ini <- ini_all[seq_len(min(Ninit, length(ini_all)))]
    pool <- ini_all[-seq_len(min(Ninit, length(ini_all)))]   # recorded extra initializations
    runfun <- function(g0) .lloyd_path(xbar, K, g0, iter.max)
    if (n_cores > 1 && .Platform$OS.type != "windows") {
      runs <- parallel::mclapply(ini, runfun, mc.cores = n_cores)
    } else {
      runs <- lapply(ini, runfun)
    }
    valid <- vapply(runs, function(r) r$valid, TRUE)
    # If every initialization produced an empty cluster, draw further initializations until one
    # is valid (at most max_extra). The discarded runs remain part of the selection outcome.
    extra <- 0
    while (!any(valid) && extra < max_extra && (is.null(inits) || extra < length(pool))) {
      g0 <- if (!is.null(inits)) pool[[extra + 1]] else .random_inits(N, K, 1, init_type)[[1]]
      runs[[length(runs) + 1]] <- runfun(g0)
      valid <- c(valid, runs[[length(runs)]]$valid)
      extra <- extra + 1
    }
    out <- list(K = K, runs = runs, valid = valid, best = NA_integer_)
    if (any(valid)) {
      keys <- rep(NA_character_, length(runs))
      keys[valid] <- vapply(runs[valid], function(r) .partition_key(r$final), "")
      ukeys <- unique(keys[valid])
      first_idx <- match(ukeys, keys)
      uwss <- vapply(first_idx, function(s) .wss(xbar, runs[[s]]$final, K), 0)
      b <- which.min(uwss)
      out$keys <- keys
      out$first_idx <- first_idx
      out$uwss <- uwss
      out$best <- first_idx[b]
      gb <- runs[[out$best]]$final
      R <- xbar - .centers(xbar, gb, K)[gb, , drop = FALSE]
      M <- Gamma_w + crossprod(R) / N
      IC[as.character(K)] <- log(det(M)) + pen(K)
    }
    byK[[as.character(K)]] <- out
  }

  if (all(is.na(IC))) {
    Khat <- NA_integer_
    final <- NULL
  } else {
    Khat <- Kset[which.min(IC)]
    bk <- byK[[as.character(Khat)]]
    final <- bk$runs[[bk$best]]$final
  }
  structure(list(Kset = Kset, byK = byK, IC = IC, Khat = Khat, final_cluster = final,
                 Gamma_w = Gamma_w, N = N, P = P, Tobs = Tobs, kappa = kappa,
                 Ninit = Ninit, iter.max = iter.max, init_type = init_type),
            class = "pk_selection")
}

#' Panel K-means Estimation with Optional Selection of K by the Information Criterion
#'
#' @description Panel Kmeans (Algorithm 1) on an NT x P panel. If \code{Kmax} is provided, K is
#' chosen by the IC over \code{Kmin:Kmax}; otherwise \code{K} is used. The function returns the
#' retained run for the selected K and, in \code{selection}, the full selection outcome
#' (all initializations, for all K), which the selective tests condition on.
#'
#' @param Z NT x P matrix of variables for clustering.
#' @param id NT vector of unit identifiers.
#' @param time NT vector of time identifiers.
#' @param K Integer; number of clusters. Required unless \code{Kmax} is provided.
#' @param Kmax Optional; if provided, select K from \code{Kmin:Kmax} by the IC.
#' @param Ninit Number of random initializations.
#' @param iter.max Maximum number of Lloyd iterations.
#' @param kappa Tuning constant of the IC penalty (\eqn{\varsigma} in the paper).
#' @param n_cores Number of cores for the initializations.
#' @param Kmin Smallest K considered when \code{Kmax} is given (default 2).
#' @param init_type "forgy" (default) or "partition"; see \code{pk_cluster}.
#'
#' @return A list with \code{clusters} (assignments along the retained path), \code{centers}
#' (cluster means of \eqn{\bar Z_i} after each assignment), \code{objectives} (T times the
#' within-cluster sum of squares after each assignment), \code{iter}, \code{final_cluster} (named by
#' unit id, units in increasing id order), \code{BIC_selected_K} (if \code{Kmax} is used),
#' \code{IC}, \code{selection} and \code{panel}.
#' @export
panel_kmeans_estimation <- function(
    Z, id, time,
    K = NULL, Kmax = NULL,
    Ninit = 10, iter.max = 100,
    kappa = 1.5,
    n_cores = 1, Kmin = 2, init_type = "forgy") {

  if (!is.null(K) && !is.null(Kmax)) stop("Specify only one of K or Kmax, not both.")
  if (is.null(K) && is.null(Kmax)) stop("Either K or Kmax must be provided.")
  if (!is.matrix(Z)) stop("Z must be a matrix")

  pd <- panel_prepare(Z, id, time)
  Kset <- if (!is.null(Kmax)) Kmin:Kmax else K
  sel <- pk_cluster(pd$xbar, Kset, Ninit = Ninit, iter.max = iter.max, Gamma_w = pd$Gamma_w,
                    Tobs = pd$Tobs, kappa = kappa, n_cores = n_cores, init_type = init_type)
  pk_result_from_selection(sel, pd, Kmax_used = !is.null(Kmax))
}

#' Format a Selection Outcome as the Output of panel_kmeans_estimation
#'
#' @description Extracts the retained run for the selected K from a \code{pk_selection} object and
#' computes the cluster centers and objectives along its path. Returns \code{NULL} with a warning
#' if no K has a valid run.
#'
#' @param sel A \code{pk_selection} object (output of \code{pk_cluster}).
#' @param pd Output of \code{panel_prepare} for the same panel.
#' @param Kmax_used Logical; if \code{TRUE}, the selected K is stored as \code{BIC_selected_K}.
#'
#' @return A list as described in \code{panel_kmeans_estimation}, or \code{NULL}.
#' @keywords internal
pk_result_from_selection <- function(sel, pd, Kmax_used = TRUE) {
  if (is.na(sel$Khat)) {
    warning("All initializations failed for every K.")
    return(NULL)
  }
  K <- sel$Khat
  bk <- sel$byK[[as.character(K)]]
  run <- bk$runs[[bk$best]]
  centers <- lapply(run$path, function(g) .centers(pd$xbar, g, K))
  objectives <- lapply(run$path, function(g) pd$Tobs * .wss(pd$xbar, g, K))
  gamma_final <- run$final
  names(gamma_final) <- as.character(pd$ids)
  res <- list(clusters = run$path, centers = centers, objectives = objectives,
              iter = length(run$path), final_cluster = gamma_final, IC = sel$IC,
              selection = sel, panel = pd)
  if (Kmax_used) res$BIC_selected_K <- K
  res
}

#' Panel K-means with K Chosen by Cross-Validation over Time Folds
#'
#' @description For each K, the time periods are split into \code{nfold} contiguous folds; the
#' clusters are estimated on the remaining periods and the squared prediction error of the cluster
#' means is evaluated on the fold. Cross-validation is not covered by the selective conditioning
#' event and is used only as a robustness check.
#'
#' @inheritParams panel_kmeans_estimation
#' @param nfold Number of folds.
#'
#' @return The output of \code{panel_kmeans_estimation} on the full sample with the selected K,
#'   augmented with \code{CV_selected_K} and \code{cv_errors} (mean fold error for each K).
#' @export
panel_kmeans_cv_time <- function(
    Z, id, time,
    Kmax,
    nfold = 5,
    Ninit = 10, iter.max = 100,
    n_cores = 1, Kmin = 2) {

  if (!is.matrix(Z)) stop("Z must be a matrix")
  T_all <- sort(unique(time))
  time_folds <- split(T_all, cut(seq_along(T_all), breaks = nfold, labels = FALSE))
  ids_sorted <- sort(unique(id))

  cv_errors <- numeric(Kmax - Kmin + 1)
  names(cv_errors) <- paste0("K=", Kmin:Kmax)

  for (K in Kmin:Kmax) {
    fold_errors <- numeric(nfold)
    for (f in 1:nfold) {
      test_idx <- time %in% time_folds[[f]]
      train_idx <- !test_idx
      res <- panel_kmeans_estimation(Z[train_idx, , drop = FALSE], id[train_idx], time[train_idx],
                                     K = K, Ninit = Ninit, iter.max = iter.max, n_cores = n_cores)
      if (is.null(res)) { fold_errors[f] <- Inf; next }
      theta_hat <- res$centers[[res$iter]]
      gamma_tot <- res$final_cluster[match(id[test_idx], ids_sorted)]
      U_test <- Z[test_idx, , drop = FALSE] - theta_hat[gamma_tot, , drop = FALSE]
      fold_errors[f] <- mean(rowSums(U_test^2))
    }
    cv_errors[K - Kmin + 1] <- mean(fold_errors)
  }

  best_K <- which.min(cv_errors) + Kmin - 1
  final_result <- panel_kmeans_estimation(Z, id, time, K = best_K, Ninit = Ninit,
                                          iter.max = iter.max, n_cores = n_cores)
  final_result$CV_selected_K <- best_K
  final_result$cv_errors <- cv_errors
  final_result
}
