# =============================================================================
# Truncation set for conditioning on the selected partition and K_IC only
# (selective_test(conditioning = "partition")).
#
# Along X(phi) = u + phi v, with the recorded data-independent initializations:
#   1. for every initialization (K, s), walk along phi through the segments on which its Lloyd
#      path is constant (exact, from its quadratic assignment constraints) and record its final
#      partition, merging adjacent segments with the same final partition;
#   2. for every K, on each piece where all final partitions are fixed, find where each partition
#      is retained (smallest WSS, ties to the first initialization): quadratic comparisons;
#   3. on each piece where the retained partitions of all K are fixed and the retained partition for
#      the observed K_IC is the observed one, keep the set where K_IC is selected by the IC
#      (polynomials of degree 2P).
# The result is exact up to numerical tolerance on [0, phi_max], where the chi mass beyond phi_max
# is negligible (below 1e-10 times the mass near the observed statistic).
# =============================================================================

# Interval containing phi of the complement of a forbidden set
#' Allowed Interval Containing a Point
#'
#' @description Returns the interval of the complement in \eqn{[0,\infty)} of the union of the
#' forbidden intervals that contains \code{phi} (up to a relative tolerance), with its lower end
#' truncated at 0. Returns \code{c(phi, phi)} if no such interval exists.
#'
#' @param forb Two-column matrix of forbidden intervals.
#' @param phi Point on the perturbation path.
#'
#' @return A numeric vector of length 2 (lower and upper end).
#' @keywords internal
.containing <- function(forb, phi) {
  Ts <- .complement(.union_intervals(forb))
  tl <- 1e-9 * max(1, phi)
  j <- which(Ts[, 1] <= phi + tl & Ts[, 2] >= phi - tl)
  if (!length(j)) return(c(phi, phi))
  c(max(Ts[j[1], 1], 0), Ts[j[1], 2])
}

# All real roots of the quadratics A x^2 + B x + C - tol (points where a constraint can stop
# holding); the tolerance is relative to the scale of the coefficients at phi_ref
#' Positive Roots of a Batch of Quadratics
#'
#' @description Returns the positive real roots of \eqn{A\phi^2 + B\phi + C - \epsilon s}, where
#' \eqn{s = |A|\bar\phi^2 + |B|\bar\phi + |C|} with \eqn{\bar\phi = \max(\phi_{ref}, 1)}. Quadratics
#' with negligible leading coefficient are treated as linear; roots are computed with the
#' numerically stable formula.
#'
#' @param A,B,C Numeric vectors of coefficients (same length).
#' @param phi_ref Reference value of \eqn{\phi} for the scale of the tolerance.
#' @param eps Relative tolerance.
#'
#' @return A numeric vector of the finite positive roots (unsorted, possibly empty).
#' @keywords internal
.quad_roots <- function(A, B, C, phi_ref, eps = 1e-9) {
  P1 <- max(phi_ref, 1)
  sc <- abs(A) * P1^2 + abs(B) * P1 + abs(C)
  sc[sc == 0] <- 1
  C <- C - eps * sc
  zA <- abs(A) * P1^2 <= 1e-12 * sc
  r <- numeric(0)
  if (any(zA)) {
    b <- B[zA]; c <- C[zA]; nz <- abs(b) * P1 > 1e-12 * sc[zA]
    r <- c(r, -c[nz] / b[nz])
  }
  if (any(!zA)) {
    a <- A[!zA]; b <- B[!zA]; c <- C[!zA]
    disc <- b^2 - 4 * a * c
    ok <- disc >= 0
    if (any(ok)) {
      sq <- sqrt(disc[ok])
      qq <- -0.5 * (b[ok] + ifelse(b[ok] >= 0, 1, -1) * sq)
      r <- c(r, qq / a[ok], ifelse(qq == 0, 0, c[ok] / qq))
    }
  }
  r[is.finite(r) & r > 0]
}

# Lloyd iterations continued from assignment g (steps counted from 1), as in .lloyd_path
#' Continue Lloyd Iterations from a Given Assignment
#'
#' @description Runs the assignment steps of \code{.lloyd_path} starting from \code{g}, stopping
#' when the assignment does not change, after \code{nmax} steps, or when a cluster becomes empty.
#'
#' @param x N x P matrix of time averages.
#' @param K Number of clusters.
#' @param g Starting assignment (labels in 0..K).
#' @param nmax Maximum number of assignment steps.
#'
#' @return A list with \code{path} (list of assignments) and \code{valid} (FALSE if a cluster
#'   became empty).
#' @keywords internal
.lloyd_continue <- function(x, K, g, nmax) {
  x2 <- rowSums(x^2)
  path <- vector("list", nmax)
  m <- 0
  valid <- TRUE
  while (m < nmax) {
    ctr <- .centers(x, g, K)
    d <- outer(x2, rowSums(ctr^2), "+") - 2 * x %*% t(ctr)
    gnew <- max.col(-d, ties.method = "first")
    m <- m + 1
    path[[m]] <- gnew
    if (any(tabulate(gnew, K) == 0)) { valid <- FALSE; break }
    if (all(gnew == g)) break
    g <- gnew
  }
  list(path = path[seq_len(m)], valid = valid)
}

# Segments of one initialization over [0, phimax]: data.frame(lo, hi, key) (contiguous), with the
# final partitions in attr "parts". The path is recomputed only from the step whose constraint
# binds, since the constraint coefficients of the earlier steps do not depend on phi.
#' Segments of One Initialization Along the Perturbation Path
#'
#' @description Walks along \eqn{\phi} from \code{phimin} to \code{phimax} for the initialization
#' \code{g0}: on each segment the Lloyd path at \eqn{X(\phi) = U + \phi V} is constant, and the
#' segment ends at the first root of its quadratic assignment constraints. Adjacent segments with
#' the same final partition are merged.
#'
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param K Number of clusters.
#' @param g0 Initial assignment.
#' @param iter.max Maximum number of assignment steps.
#' @param phimin,phimax Range of \eqn{\phi}.
#' @param eps Relative numerical tolerance.
#' @param phiref Reference value of \eqn{\phi} for the scale of the tolerance.
#'
#' @return A data.frame with columns \code{lo}, \code{hi}, \code{key} (partition key, or
#'   \code{"invalid"}) and attribute \code{"parts"}, the list of final partitions named by key.
#' @keywords internal
.run_segments <- function(U, V, K, g0, iter.max, phimin, phimax, eps, phiref = phimax) {
  lo <- hi <- numeric(0); key <- character(0); parts <- list()
  path <- list(); roots <- list(); valid <- TRUE
  phi <- phimin; mstar <- 1
  nstep <- 0
  while (phi < phimax && nstep < 1e5) {
    nstep <- nstep + 1
    x <- U + phi * V
    gstart <- if (mstar == 1) g0 else path[[mstar - 1]]
    cont <- .lloyd_continue(x, K, gstart, iter.max - (mstar - 1))
    path <- c(path[seq_len(mstar - 1)], cont$path)
    valid <- cont$valid
    roots <- roots[seq_len(mstar - 1)]
    prev <- gstart
    for (m in mstar:length(path)) {
      cf <- .assign_coefs(U, V, prev, path[[m]], K)
      roots[[m]] <- .quad_roots(cf[, 1], cf[, 2], cf[, 3], max(phiref, 1), eps)
      prev <- path[[m]]
    }
    tl <- phi + 1e-10 * max(1, phi)
    nb <- vapply(roots, function(r) { r <- r[r > tl]; if (length(r)) min(r) else Inf }, 0)
    mstar <- which.min(nb)
    hi_s <- nb[mstar]
    k <- if (valid) .partition_key(path[[length(path)]]) else "invalid"
    if (k != "invalid" && is.null(parts[[k]])) parts[[k]] <- path[[length(path)]]
    n <- length(lo)
    if (n > 0 && key[n] == k) hi[n] <- hi_s else { lo <- c(lo, phi); hi <- c(hi, hi_s); key <- c(key, k) }
    if (!is.finite(hi_s)) break
    phi <- hi_s + 1e-8 * max(1, hi_s)
  }
  structure(data.frame(lo = lo, hi = hi, key = key, stringsAsFactors = FALSE), parts = parts)
}

# Optional C++ version of .run_segments (inst/cpp/walk.cpp), compiled once with Rcpp and cached;
# falls back to the R version if Rcpp or a compiler is not available.
.cpp_env <- new.env()
#' Compiled Partition Walk
#'
#' @description Compiles \code{inst/cpp/walk.cpp} with Rcpp once per session (cached on disk) and
#' returns its \code{run_segments_cpp} routine. Returns \code{NULL} if the option
#' \code{PanelKmeansInference.use_cpp} is FALSE or if Rcpp or a compiler is not available.
#'
#' @return A function, or \code{NULL}.
#' @keywords internal
.get_cpp_walk <- function() {
  if (!isTRUE(getOption("PanelKmeansInference.use_cpp", TRUE))) return(NULL)
  if (!is.null(.cpp_env$f)) return(if (isFALSE(.cpp_env$f)) NULL else .cpp_env$f$run_segments_cpp)
  f <- tryCatch({
    src <- system.file("cpp", "walk.cpp", package = "PanelKmeansInference")
    if (!nzchar(src) || !requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp not available")
    cache <- file.path(tools::R_user_dir("PanelKmeansInference", "cache"), "cpp")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    e <- new.env()
    ok <- tryCatch({
      suppressWarnings(suppressMessages(Rcpp::sourceCpp(src, env = e, cacheDir = cache, rebuild = FALSE, showOutput = FALSE)))
      TRUE
    }, error = function(err) FALSE)
    if (!ok) {   # e.g. concurrent compilation by several processes: compile privately
      e <- new.env()
      suppressWarnings(suppressMessages(Rcpp::sourceCpp(src, env = e, cacheDir = tempfile("pki_cpp"), showOutput = FALSE)))
    }
    e
  }, error = function(err) {
    message("PanelKmeansInference: C++ walk not available (", conditionMessage(err),
            "); using the slower R version.")
    FALSE
  })
  .cpp_env$f <- f
  if (isFALSE(f)) NULL else f$run_segments_cpp
}

# Other compiled routines of walk.cpp (NULL if the C++ walk is not available)
#' Compiled Routine of walk.cpp by Name
#'
#' @description Returns a routine compiled from \code{inst/cpp/walk.cpp}, compiling it first if
#' needed (see \code{.get_cpp_walk}).
#'
#' @param name Name of the compiled routine.
#'
#' @return A function, or \code{NULL} if the C++ code is not available.
#' @keywords internal
.get_cpp_fun <- function(name) {
  if (is.null(.get_cpp_walk())) return(NULL)
  .cpp_env$f[[name]]
}

# Segments of one initialization, C++ if available
#' Segments of One Initialization, Compiled if Available
#'
#' @description Same result as \code{.run_segments}, computed by the compiled routine when it is
#' available (with adjacent segments of the same partition merged), and by \code{.run_segments}
#' otherwise.
#'
#' @inheritParams .run_segments
#'
#' @return A data.frame with columns \code{lo}, \code{hi}, \code{key} and attribute
#'   \code{"parts"}, as in \code{.run_segments}.
#' @keywords internal
.run_segments_any <- function(U, V, K, g0, iter.max, phimin, phimax, eps, phiref = phimax) {
  f <- .get_cpp_walk()
  if (is.null(f)) return(.run_segments(U, V, K, g0, iter.max, phimin, phimax, eps, phiref))
  r <- f(U, V, as.integer(K), as.integer(g0), as.integer(iter.max), phimin, phimax, eps)
  ns <- length(r$lo)
  key <- character(ns); parts <- list()
  for (j in seq_len(ns)) {
    if (r$valid[j] == 1) {
      g <- r$parts[, j]
      k <- .partition_key(g)
      key[j] <- k
      if (is.null(parts[[k]])) parts[[k]] <- g
    } else key[j] <- "invalid"
  }
  # merge adjacent segments with the same key (different labelings of the same partition)
  if (ns > 1) {
    keep <- c(TRUE, key[-1] != key[-ns])
    grp <- cumsum(keep)
    lo <- r$lo[keep]; hi <- as.numeric(tapply(r$hi, grp, max)); key <- key[keep]
  } else { lo <- r$lo; hi <- r$hi }
  structure(data.frame(lo = lo, hi = hi, key = key, stringsAsFactors = FALSE), parts = parts)
}

# Forbidden set {phi in [a, b]: f(phi) > tol} of a polynomial f of degree <= deg, robust version:
# roots from a local interpolation on [a, b] plus a scan for sign changes refined by bisection
#' Forbidden Set of a Polynomial on an Interval
#'
#' @description Computes \eqn{\{\phi\in[a,b]: f(\phi) > tol\}} for a polynomial f of degree at most
#' \code{deg}, with \code{tol} relative to the scale of f on \eqn{[a,b]}. Candidate roots come from an
#' interpolation at Chebyshev nodes and from sign changes on a grid refined by \code{uniroot}; the
#' sign of f is then checked at the midpoint of each piece. An infinite \code{b} is replaced by
#' \code{a + 10}.
#'
#' @param fun Function evaluating f.
#' @param deg Degree bound of f.
#' @param a,b Interval ends.
#' @param eps Relative tolerance.
#' @param ngrid Number of grid points for the sign-change scan.
#' @param vectorized Logical; whether \code{fun} accepts a vector argument.
#'
#' @return Two-column matrix of forbidden intervals.
#' @keywords internal
.poly_forbidden_piece <- function(fun, deg, a, b, eps = 1e-9, ngrid = 64, vectorized = FALSE) {
  ev <- if (vectorized) fun else function(x) vapply(x, fun, 0)
  if (!is.finite(b)) b <- a + 10
  if (b <= a) return(matrix(numeric(0), 0, 2))
  tt <- (a + b) / 2 + (b - a) / 2 * cos(pi * (0:deg + 0.5) / (deg + 1))   # Chebyshev nodes
  fv <- ev(tt)
  grid <- seq(a, b, length.out = ngrid)
  fg <- ev(grid)
  scale <- max(abs(c(fv, fg)), 1e-300)
  tol <- eps * scale
  roots <- numeric(0)
  z <- (tt - a) / (b - a)
  cf <- tryCatch(solve(outer(z, 0:deg, "^"), fv), error = function(e) NULL)
  if (!is.null(cf)) {
    cf[abs(cf) <= 1e-13 * max(abs(cf))] <- 0
    nz <- which(cf != 0)
    if (length(nz) && max(nz) > 1) {
      rr <- polyroot(cf[1:max(nz)])
      rr <- Re(rr[abs(Im(rr)) <= 1e-4 * (1 + abs(Re(rr)))])
      roots <- a + (b - a) * rr
    }
  }
  g <- fg - tol
  sc <- which(g[-1] * g[-ngrid] < 0)
  for (j in sc) roots <- c(roots, stats::uniroot(function(x) fun(x) - tol, c(grid[j], grid[j + 1]),
                                                  tol = 1e-12 * max(1, b))$root)
  roots <- sort(unique(roots[roots > a & roots < b]))
  br <- c(a, roots, b)
  mids <- (br[-1] + br[-length(br)]) / 2
  bad <- ev(mids) > tol
  cbind(br[-length(br)][bad], br[-1][bad])
}

# Value of a contiguous segment table at points m
#' Key of a Segment Table at Given Points
#'
#' @description Looks up, for each point in \code{m}, the key of the segment of \code{tab} that
#' contains it; \code{NA} if the point lies outside all segments.
#'
#' @param tab A data.frame with columns \code{lo}, \code{hi}, \code{key}, sorted by \code{lo}.
#' @param m Numeric vector of points.
#'
#' @return A character vector of the same length as \code{m}.
#' @keywords internal
.key_at <- function(tab, m) {
  j <- findInterval(m, tab$lo)
  out <- rep(NA_character_, length(m))
  ok <- j > 0
  out[ok] <- tab$key[j[ok]]
  out[ok][m[ok] > tab$hi[j[ok]] + 1e-7 * pmax(1, m[ok])] <- NA_character_
  out
}

# Retained partition for K, piecewise in phi on [lo_b, hi_b]: data.frame(lo, hi, key) with the
# partitions in attr "parts" (key "none" where no initialization gives a valid partition).
# C++ (winners_cpp in walk.cpp) if available, otherwise the R version below; both give the same pieces.
#' Retained Partition for a Given K Along the Perturbation Path
#'
#' @description For every recorded initialization of K, computes its final partition piecewise in
#' \eqn{\phi} on \eqn{[lo_b, hi_b]}, and on each piece determines the retained partition (smallest
#' within-cluster sum of squares, ties to the first initialization; recorded extra initializations
#' are used only if no main initialization is valid). Uses the compiled \code{winners_cpp} if
#' available.
#'
#' @param sel A \code{pk_selection} object.
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param K Number of clusters.
#' @param lo_b,hi_b Range of \eqn{\phi}.
#' @param eps Relative numerical tolerance.
#' @param phiref Reference value of \eqn{\phi} for the scale of the tolerance.
#'
#' @return A list with \code{w}, a data.frame (\code{lo}, \code{hi}, \code{key}; key \code{"none"}
#'   where no initialization is valid) with attribute \code{"parts"}, and \code{nseg}, the number of
#'   segments visited.
#' @keywords internal
.winners_K <- function(sel, U, V, K, lo_b, hi_b, eps, phiref = hi_b) {
  bk <- sel$byK[[as.character(K)]]
  nmain <- min(sel$Ninit, length(bk$runs))
  f <- .get_cpp_fun("winners_cpp")
  if (!is.null(f)) {
    r <- f(U, V, as.integer(K), lapply(bk$runs, function(x) as.integer(x$init)), as.integer(nmain),
           as.integer(sel$iter.max), lo_b, hi_b, eps, phiref)
    parts <- lapply(seq_len(ncol(r$parts)), function(j) r$parts[, j])
    names(parts) <- vapply(parts, .partition_key, "")
    key <- ifelse(r$key < 0, "none", names(parts)[pmax(r$key, 0) + 1])
    w <- data.frame(lo = r$lo, hi = r$hi, key = key, stringsAsFactors = FALSE)
    attr(w, "parts") <- parts
    return(list(w = w, nseg = r$nseg))
  }
  tabs <- lapply(bk$runs, function(r) .run_segments(U, V, K, r$init, sel$iter.max, lo_b, hi_b, eps, phiref))
  nseg <- sum(vapply(tabs, nrow, 0L))
  parts <- do.call(c, lapply(tabs, function(t) attr(t, "parts")))
  parts <- parts[!duplicated(names(parts))]
  br <- sort(unique(c(lo_b, unlist(lapply(tabs, function(t) c(t$lo, t$hi))), hi_b)))
  br <- br[br >= lo_b & br <= hi_b]
  mids <- (br[-1] + br[-length(br)]) / 2
  KS <- vapply(tabs, .key_at, character(length(mids)), m = mids)
  if (!is.matrix(KS)) KS <- matrix(KS, nrow = length(mids))
  # candidate partitions (in order of first initialization) on each elementary piece
  cand_of <- function(ks) {
    main <- ks[seq_len(nmain)]
    cand <- main[!is.na(main) & main != "invalid"]
    if (!length(cand) && length(ks) > nmain) {           # recorded extra initializations
      ex <- ks[-seq_len(nmain)]
      ex <- ex[!is.na(ex) & ex != "invalid"]
      cand <- if (length(ex)) ex[1] else character(0)
    }
    unique(cand)
  }
  sig <- apply(KS, 1, function(ks) paste(cand_of(ks), collapse = "|"))
  # collapse consecutive pieces with the same candidate set
  st <- c(TRUE, sig[-1] != sig[-length(sig)])
  gi <- cumsum(st)
  A0 <- br[-length(br)][st]
  B0 <- as.numeric(tapply(br[-1], gi, max))
  S0 <- sig[st]
  wcf <- list()
  olo <- list(); ohi <- list(); okey <- list()
  for (j in seq_along(A0)) {
    a <- A0[j]; b <- B0[j]
    if (b <= a) next
    if (!nzchar(S0[j])) { olo[[j]] <- a; ohi[[j]] <- b; okey[[j]] <- "none"; next }
    uk <- strsplit(S0[j], "|", fixed = TRUE)[[1]]
    if (length(uk) == 1) { olo[[j]] <- a; ohi[[j]] <- b; okey[[j]] <- uk; next }
    for (k in uk) if (is.null(wcf[[k]])) wcf[[k]] <- .wss_coefs(U, V, parts[[k]], K)
    cf <- do.call(rbind, wcf[uk])
    # winner changes only at roots of pairwise differences
    pr <- utils::combn(length(uk), 2)
    dA <- cf[pr[1, ], 1] - cf[pr[2, ], 1]; dB <- cf[pr[1, ], 2] - cf[pr[2, ], 2]
    dC <- cf[pr[1, ], 3] - cf[pr[2, ], 3]
    rr <- .quad_roots(dA, dB, dC + 0, max(b, 1), 0)
    cuts <- sort(unique(c(a, rr[rr > a & rr < b], b)))
    mm <- (cuts[-1] + cuts[-length(cuts)]) / 2
    wss <- outer(mm^2, cf[, 1]) + outer(mm, cf[, 2]) + matrix(cf[, 3], length(mm), length(uk), byrow = TRUE)
    win <- max.col(-wss, ties.method = "first")
    olo[[j]] <- cuts[-length(cuts)]; ohi[[j]] <- cuts[-1]; okey[[j]] <- uk[win]
  }
  w <- data.frame(lo = unlist(olo), hi = unlist(ohi), key = unlist(okey), stringsAsFactors = FALSE)
  w <- w[order(w$lo), ]
  attr(w, "parts") <- parts
  list(w = w, nseg = nseg)
}

# det(M0 + phi M1 + phi^2 M2) for a vector phi (closed form for P <= 2)
#' Determinant of a Matrix Quadratic in phi
#'
#' @description Evaluates \eqn{\det(M_0 + \phi M_1 + \phi^2 M_2)} for a vector of \eqn{\phi}, in
#' closed form for P of 1 or 2 and with \code{det} otherwise.
#'
#' @param M A list with P x P matrices \code{M0}, \code{M1}, \code{M2} (output of \code{.M_poly}).
#' @param phi Numeric vector.
#'
#' @return A numeric vector of the same length as \code{phi}.
#' @keywords internal
.detM_vec <- function(M, phi) {
  P <- nrow(M$M0)
  e <- function(a, b) M$M0[a, b] + phi * M$M1[a, b] + phi^2 * M$M2[a, b]
  if (P == 1) return(e(1, 1))
  if (P == 2) return(e(1, 1) * e(2, 2) - e(1, 2) * e(2, 1))
  vapply(phi, function(x) det(M$M0 + x * M$M1 + x^2 * M$M2), 0)
}

#' Truncation Set on One Block of the Perturbation Path
#'
#' @description Computes the subset of \eqn{[lo_b, hi_b]} on which the observed partition is retained
#' for \eqn{\hat K} and \eqn{\hat K_{IC}} is selected by the information criterion: retained
#' partitions are obtained by \code{.winners_K} (for the other K only where the observed partition is
#' retained), and the IC comparisons are polynomial inequalities of degree 2P solved by
#' \code{.poly_forbidden_piece}.
#'
#' @param fit A \code{pk_fit} object.
#' @param U,V N x P matrices defining \eqn{X(\phi) = U + \phi V}.
#' @param q Degrees of freedom of the test (not used in the computation).
#' @param lo_b,hi_b Ends of the block.
#' @param eps Relative numerical tolerance.
#'
#' @return A list with \code{intervals} (two-column matrix of disjoint intervals) and \code{nseg}
#'   (segments visited).
#' @keywords internal
.good_block <- function(fit, U, V, q, lo_b, hi_b, eps) {
  sel <- fit$selection
  N <- sel$N; P <- sel$P; Tobs <- sel$Tobs
  Kset <- sel$Kset
  Kh <- sel$Khat
  target <- .partition_key(sel$final_cluster)
  phimax <- hi_b
  pen <- function(K) (K * P + N) * sel$kappa * log(N * Tobs) / (N * Tobs)
  nseg <- 0

  # 1-2. retained partition of each K, piecewise in phi. For K_hat on the whole block; for the
  # other K only on the pieces where the observed partition is retained for K_hat, since the
  # information criterion matters only there (the walk is exact on any subinterval).
  winners <- list()
  wk <- .winners_K(sel, U, V, Kh, lo_b, hi_b, eps, phiref = hi_b)
  nseg <- nseg + wk$nseg
  winners[[as.character(Kh)]] <- wk$w
  on_t <- wk$w$key == target
  if (!any(on_t)) return(list(intervals = matrix(numeric(0), 0, 2), nseg = nseg))
  tiv <- .union_intervals(cbind(wk$w$lo[on_t], wk$w$hi[on_t]))
  for (K in setdiff(Kset, Kh)) {
    ws <- lapply(seq_len(nrow(tiv)), function(r) .winners_K(sel, U, V, K, tiv[r, 1], tiv[r, 2], eps, phiref = hi_b))
    nseg <- nseg + sum(vapply(ws, function(x) x$nseg, 0))
    w <- do.call(rbind, lapply(ws, function(x) x$w))
    parts <- do.call(c, lapply(ws, function(x) attr(x$w, "parts")))
    attr(w, "parts") <- parts[!duplicated(names(parts))]
    winners[[as.character(K)]] <- w
  }

  # 3. selection of K_IC on pieces where the retained partition for K_hat is the observed one
  wh <- winners[[as.character(Kh)]]
  ph <- attr(wh, "parts")
  cand <- wh[wh$key == target, , drop = FALSE]
  good <- list()
  others <- setdiff(Kset, Kh)
  detM <- .detM_vec
  Mh <- if (nrow(cand)) .M_poly(U, V, ph[[target]], Kh, sel$Gamma_w, N, Tobs) else NULL
  Mcache <- list()
  getM <- function(K, k) {
    nm <- paste(K, k)
    if (is.null(Mcache[[nm]]))
      Mcache[[nm]] <<- .M_poly(U, V, attr(winners[[as.character(K)]], "parts")[[k]], K, sel$Gamma_w, N, Tobs)
    Mcache[[nm]]
  }
  for (j in seq_len(nrow(cand))) {
    a <- cand$lo[j]; b <- cand$hi[j]
    if (!length(others)) { good[[length(good) + 1]] <- c(a, b); next }
    br <- sort(unique(c(a, b, unlist(lapply(others, function(K) {
      w <- winners[[as.character(K)]]; c(w$lo, w$hi) })))))
    br <- br[br >= a & br <= b]
    for (i in seq_len(length(br) - 1)) {
      a2 <- br[i]; b2 <- br[i + 1]
      if (b2 <= a2) next
      m <- (a2 + b2) / 2
      f0 <- exp(pen(Kh)) * detM(Mh, m)
      forb <- rbind(c(0, a2), c(b2, Inf))
      for (K in others) {
        w <- winners[[as.character(K)]]
        k <- .key_at(w, m)
        if (is.na(k) || k == "none") next
        MK <- getM(K, k)
        fun <- function(phi) (exp(pen(Kh)) * detM(Mh, phi) - exp(pen(K)) * detM(MK, phi)) / f0
        forb <- rbind(forb, .poly_forbidden_piece(fun, 2 * P, a2, b2, eps, vectorized = TRUE))
      }
      Ts <- .complement(.union_intervals(forb))
      if (nrow(Ts)) good[[length(good) + 1]] <- Ts
    }
  }
  Tset <- if (length(good)) .union_intervals(do.call(rbind, good)) else matrix(numeric(0), 0, 2)
  list(intervals = Tset, nseg = nseg)
}

#' Selective p-value by the Partition Walk
#'
#' @description Computes the truncation set for \code{conditioning = "partition"} block by block:
#' first on \eqn{[\max(0, D - 0.5), \phi_{max}]}, then on blocks of doubling length towards 0, until
#' the p-value is below \code{tol}, its upper and lower bounds agree to relative accuracy \code{tol},
#' 0 is reached, or \code{max_seg} segments have been visited. \eqn{\phi_{max}} is extended until the
#' reference mass beyond it is negligible.
#'
#' @param fit A \code{pk_fit} object.
#' @param u,v NP vectors defining \eqn{X(\phi) = u + \phi v}.
#' @param D Observed statistic.
#' @param q Degrees of freedom.
#' @param eps Relative numerical tolerance.
#' @param tol Accuracy of the p-value.
#' @param max_seg Maximum number of segments.
#' @param B Optional number of cosines (Hotelling reference); \code{NULL} for the chi reference.
#'
#' @return A list with \code{pval} (upper bound if the walk stopped early), \code{intervals} and
#'   \code{nseg}.
#' @keywords internal
.partition_walk <- function(fit, u, v, D, q, eps = 1e-9, tol = 1e-3, max_seg = 50000, B = NULL) {
  sel <- fit$selection
  N <- sel$N; P <- sel$P
  U <- matrix(u, N, P, byrow = TRUE)
  V <- matrix(v, N, P, byrow = TRUE)
  phimax <- D + max(3, 10 / max(D, 1e-3))
  if (!is.null(B)) {   # heavier tails of the F reference: extend until the tail is negligible
    while (.logmass(phimax, Inf, q, B) > .logmass(max(D - 0.5, 0), D + 0.5, q, B) + log(1e-10)) phimax <- 2 * phimax
  }
  # first block around and above the observed statistic, then blocks of doubling length downwards
  a <- max(0, D - 0.5)
  gb <- .good_block(fit, U, V, q, a, phimax, eps)
  good <- gb$intervals; nseg <- gb$nseg
  delta <- 1
  repeat {
    if (nrow(good)) {
      up <- good[, 2] > D
      lnum <- if (any(up)) .lse(.logmass(pmax(good[up, 1], D), good[up, 2], q, B)) else -Inf
      lden <- .lse(.logmass(good[, 1], good[, 2], q, B))
      p_up <- exp(lnum - lden)
      p_lo <- exp(lnum - .lse(c(lden, .logmass(0, a, q, B))))
      if (a <= 0 || p_up < tol || (p_up - p_lo) <= tol * p_up) break
    } else if (a <= 0) break
    if (nseg >= max_seg) { warning("max_seg reached in the partition walk; the p-value is an upper bound."); break }
    a_new <- max(0, a - delta)
    gb <- .good_block(fit, U, V, q, a_new, a, eps)
    good <- .union_intervals(rbind(good, gb$intervals)); nseg <- nseg + gb$nseg
    a <- a_new
    delta <- 2 * delta
  }
  if (nrow(good) && good[nrow(good), 2] >= phimax * (1 - 1e-12)) good[nrow(good), 2] <- Inf
  pval <- if (nrow(good)) truncated_chi_pvalue(D, q, good, B) else NA_real_
  list(pval = pval, intervals = good, nseg = nseg)
}
