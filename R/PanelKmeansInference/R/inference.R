#' Selective Test of Equality of Two Estimated Cluster Means
#'
#' @description Pairwise selective test of \eqn{\theta^0_{k}(\hat{\mathcal C}) = \theta^0_{g}(\hat{\mathcal C})}.
#'
#' @param fit A \code{pk_fit} object (see \code{panel_selective_fit}).
#' @param k1,k2 Cluster labels.
#' @param ... Passed to \code{selective_test} (\code{metric}, \code{condition_K}, \code{condition_starts}).
#'
#' @return A list with \code{final_interval} (truncation set), \code{test_stat} and \code{pval}.
#' @export
panel_kmeans_inference <- function(fit, k1, k2, ...) {
  res <- selective_test(fit, L_pairwise(fit, k1, k2), ...)
  list(final_interval = res$intervals, test_stat = res$stat, pval = res$pval)
}

#' Selective Homogeneity Test for Panel Kmeans
#'
#' @description Estimates the clusters by Panel Kmeans (K fixed, or chosen by the IC from 2 to
#' \code{Kmax}) and computes the selective p-values of all pairwise equality hypotheses, their
#' merging by the M-family function \eqn{F_{homo,r}} and the joint selective homogeneity test.
#'
#' @param Z NT x P matrix.
#' @param id,time NT vectors of unit and time identifiers.
#' @param K Number of clusters (alternative to \code{Kmax}).
#' @param Kmax Largest number of clusters for the IC.
#' @param Ninit Number of random initializations.
#' @param iter.max Maximum number of Lloyd iterations.
#' @param lrv Only "EWC" (the OS estimator) is supported for selective tests.
#' @param lrv_par Number of cosines B (default: \code{choose_B}).
#' @param pairs Optional two-column matrix of pairs.
#' @param pcombine_fun Merging function (default "Genmean_neq").
#' @param r Order of the generalized mean (default -20).
#' @param n_cores Cores for the initializations.
#' @param kappa IC tuning constant.
#' @param cB Constant of the rule for B.
#' @param metric "xi" (default), "shrink" or "euclidean"; see \code{selective_test}.
#' @param conditioning "partition" (default) or "paths"; see \code{selective_test}.
#' @param condition_K,condition_starts See \code{selective_test}.
#' @param init_type "forgy" (default) or "partition"; see \code{pk_cluster}.
#' @param fit Optional precomputed \code{pk_fit} object (then Z, id, time and the clustering
#'   arguments are ignored).
#'
#' @return A list with \code{pairwise_pvalues}, \code{pairs}, \code{pval} (merged),
#'   \code{joint_pval} (joint selective homogeneity test), \code{joint_cepa_pval} (joint selective
#'   C-EPA test), \code{clustering}, \code{K}, \code{B} and \code{fit}.
#' @export
panel_homogeneity_test <- function(
    Z = NULL, id = NULL, time = NULL,
    K = NULL,
    Kmax = NULL,
    Ninit = 10,
    iter.max = 100,
    lrv = "EWC",
    lrv_par = NULL,
    pairs = NULL,
    pcombine_fun = "Genmean_neq",
    r = -20,
    n_cores = 1,
    kappa = 1.5,
    cB = NULL,
    metric = "xi",
    conditioning = "partition",
    condition_K = TRUE,
    condition_starts = TRUE,
    fit = NULL,
    init_type = "forgy"
) {
  if (lrv != "EWC") stop("Selective tests are implemented with the OS (EWC) estimator only.")
  if (is.null(fit)) {
    fit <- panel_selective_fit(Z, id, time, K = K, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                               kappa = kappa, lrv_par = lrv_par, cB = cB, n_cores = n_cores,
                               init_type = init_type)
  }
  res <- panel_selective_tests(fit, pairs = pairs, pcombine_fun = pcombine_fun, r = r,
                               joint = TRUE, metric = metric, conditioning = conditioning, condition_K = condition_K,
                               condition_starts = condition_starts)
  list(
    pairwise_pvalues = res$pairwise_pvalues,
    pairs = res$pairs,
    pval = res$hom_pval,
    joint_pval = res$joint_hom_pval,
    joint_cepa_pval = res$joint_cepa_pval,
    clustering = fit$clustering,
    K = fit$K,
    B = fit$B,
    fit = fit
  )
}

#' Selective Tests that Each Estimated Cluster Mean is Zero
#'
#' @description For each estimated cluster k, selective test of \eqn{\theta^0_k(\hat{\mathcal C})=0}
#' (Section 5.3 of the paper), that is the joint test with \eqn{\mathsf R = e_k'\otimes I_P}. The
#' p-values are combined by the step-down procedure of Holm (1979), which controls the family-wise
#' error rate conditionally on the selection outcome under arbitrary dependence between them; the
#' clusters with a Holm-adjusted p-value below the level are those for which EPA fails.
#'
#' @inheritParams panel_homogeneity_test
#' @param ref Reference law: "chi" (default) or "F"; see \code{selective_test}.
#' @return A list with \code{bycluster_pvals} (selective p-values, one per cluster),
#'   \code{bycluster_holm} (Holm-adjusted p-values), \code{holm_pval} (the smallest adjusted
#'   p-value, a conditional test that all cluster means are zero), \code{bycluster_stats} (the
#'   statistics \eqn{D_k}), \code{centers} (K x P matrix of estimated cluster means of Z),
#'   \code{sizes} (cluster sizes), \code{in_set} (diagnostic), \code{clustering} and \code{fit}.
#' @export
panel_bycluster_tests <- function(
  Z = NULL, id = NULL, time = NULL,
  K = NULL,
  Kmax = NULL,
  Ninit = 10,
  iter.max = 100,
  lrv = "EWC",
  lrv_par = NULL,
  n_cores = 1,
  kappa = 1.5,
  cB = NULL,
  metric = "xi",
  conditioning = "partition",
  condition_K = TRUE,
  condition_starts = TRUE,
  ref = "chi",
  fit = NULL,
  init_type = "forgy"
) {
  if (lrv != "EWC") stop("Selective tests are implemented with the OS (EWC) estimator only.")
  if (is.null(fit)) {
    fit <- panel_selective_fit(Z, id, time, K = K, Kmax = Kmax, Ninit = Ninit, iter.max = iter.max,
                               kappa = kappa, lrv_par = lrv_par, cB = cB, n_cores = n_cores,
                               init_type = init_type)
  }
  res <- lapply(1:fit$K, function(k)
    selective_test(fit, L_cluster(fit, k), metric = metric, conditioning = conditioning,
                   condition_K = condition_K, condition_starts = condition_starts, ref = ref))
  pvals <- vapply(res, function(r) r$pval, 0)
  names(pvals) <- paste0("cluster_", 1:fit$K)
  g <- fit$selection$final_cluster
  centers <- .centers(fit$panel$xbar, g, fit$K)
  list(bycluster_pvals = pvals,
       bycluster_holm = holm_adjust(pvals),
       holm_pval = min(holm_adjust(pvals)),
       bycluster_stats = vapply(res, function(r) r$stat, 0),
       centers = centers,
       sizes = tabulate(g, fit$K),
       in_set = vapply(res, function(r) r$in_set, TRUE),
       clustering = fit$clustering, fit = fit)
}
