.onAttach <- function(libname, pkgname) {
  pkg_date <- utils::packageDescription("nmfkc", fields = "Date")
  formatted_date <- format(as.Date(pkg_date), "%d %b %Y")
  packageStartupMessage(
    paste("Package: nmfkc (Version", as.character(utils::packageVersion("nmfkc")),
          ", released on", formatted_date, ")")
  )
  packageStartupMessage("https://ksatohds.github.io/nmfkc/")
}







#' @title Create a kernel matrix from covariates
#' @description
#' \code{nmfkc.kernel} constructs a kernel matrix from covariate matrices.
#' It supports Gaussian, Exponential, Periodic, Linear, Normalized Linear, and Polynomial kernels.
#'
#' @param U Covariate matrix \eqn{U(K,N) = (u_1, \dots, u_N)}. Each row may be normalized in advance.
#' @param V Covariate matrix \eqn{V(K,M) = (v_1, \dots, v_M)}, typically used for prediction. If \code{NULL}, the default is \code{U}.
#' @param kernel Kernel function to use. Default is \code{"Gaussian"}. Options are \code{"Gaussian"}, \code{"Exponential"}, \code{"Periodic"}, \code{"Linear"}, \code{"NormalizedLinear"}, and \code{"Polynomial"}.
#' @param ... Additional arguments passed to the specific kernel function (e.g., \code{beta}, \code{degree}).
#'
#' @return Kernel matrix \eqn{A(N,M)}.
#' @seealso \code{\link{nmfkc.kernel}}, \code{\link{nmfkc.cv}}
#' @export
#' @source Satoh, K. (2024). Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model.
#'   arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' dim(A)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel <- function(U, V = NULL,
                         kernel = c("Gaussian","Exponential","Periodic",
                                    "Linear","NormalizedLinear","Polynomial"),...){
  k_params <- list(...)
  U <- as.matrix(U); storage.mode(U) <- "double"
  if (is.null(V)) V <- U else V <- as.matrix(V)
  storage.mode(V) <- "double"
  if (nrow(U) == 0) stop("'U' must have at least one row (feature).")
  if (nrow(U) != nrow(V)) stop("'U' and 'V' must have the same number of rows (features).")
  kernel <- match.arg(kernel)
  G <- crossprod(U, V)

  # Determine the specific parameter used for the kernel (e.g., beta or degree)
  # This section extracts the effective parameter value to store in attributes.
  effective_param <- NULL
  if (kernel %in% c("Gaussian", "Exponential")) {
    effective_param <- if (!is.null(k_params$beta)) k_params$beta else 0.5
  } else if (kernel == "Periodic") {
    effective_param <- if (!is.null(k_params$beta)) k_params$beta else c(1,1)
  } else if (kernel == "Polynomial") {
    # For Polynomial, the effective param is (beta, degree) pair
    effective_param <- list(beta = if (!is.null(k_params$beta)) k_params$beta else 0,
                            degree = if (!is.null(k_params$degree)) k_params$degree else 2)
  }

  if (kernel %in% c("Gaussian","Exponential","Periodic")) {
    u2 <- colSums(U * U)
    v2 <- colSums(V * V)
    D2 <- outer(u2, v2, "+") - 2 * G
    D2 <- pmax(D2, 0)
  }

  K <- switch(kernel,
              Gaussian ={
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0.5
                exp(-beta * D2)
              },
              Exponential = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0.5
                d <- sqrt(D2)
                exp(-beta * d)
              },
              Periodic = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else c(1,1)
                if (length(beta) < 2L)
                  stop("For 'Periodic', set beta as c(beta1, beta2).")
                d <- sqrt(D2)
                exp(-beta[1] * (sin(beta[2] * d)^2))
              },
              Linear = G,
              NormalizedLinear = {
                nu <- sqrt(colSums(U * U))
                nv <- sqrt(colSums(V * V))
                nu[nu == 0] <- .Machine$double.eps
                nv[nv == 0] <- .Machine$double.eps
                G / outer(nu, nv)
              },
              Polynomial = {
                beta <- if (!is.null(k_params$beta)) k_params$beta else 0
                degree <- if (!is.null(k_params$degree)) k_params$degree else 2
                (G + beta)^degree
              }
  )
  K[K < 0] <- 0
  dimnames(K) <- list(colnames(U), colnames(V))

  # --- Store Kernel Metadata as Attributes (NEW BLOCK) ---
  # These attributes allow nmfkc to identify the kernel type and parameters used.
  attr(K, "params") <- effective_param
  attr(K, "kernel") <- kernel
  attr(K, "function.name") <- "nmfkc.kernel"
  # ------------------------------------------------------

  return(K)
}



#' @title Estimate Gaussian/RBF kernel parameter beta from covariates (supports landmarks)
#'
#' @description
#' Computes a data-driven reference scale for the Gaussian/RBF kernel from covariates
#' using a robust "median nearest-neighbor (or nearest-landmark) distance" heuristic,
#' and returns the corresponding kernel parameter \eqn{\beta}.
#'
#' The Gaussian/RBF kernel is assumed to be written in the form
#' \deqn{k(u,v) = \exp\{-\beta \|u-v\|^2\} = \exp\{-\|u-v\|^2/(2\sigma^2)\},}
#' hence \eqn{\beta = 1/(2\sigma^2)}. This function first estimates a typical distance
#' scale \eqn{\sigma_0} by the median of distances, then sets \eqn{\beta_0 = 1/(2\sigma_0^2)}.
#'
#' If \code{Uk} is \code{NULL}, \eqn{\sigma_0} is estimated as the median of
#' nearest-neighbor distances within \code{U} (excluding self-distance).
#' If \code{Uk} is provided, \eqn{\sigma_0} is estimated as the median of
#' nearest-landmark distances from each sample in \code{U} to its closest landmark in \code{Uk}.
#'
#' To control memory usage for large \code{N} (and \code{M}), distances are computed in blocks.
#' Optionally, columns of \code{U} can be randomly subsampled via \code{sample_size} to reduce cost.
#'
#' @details
#' \strong{Candidate grid:}
#' Along with \code{beta}, the function returns \code{beta_candidates}, a small logarithmic grid
#' suitable for cross-validation.
#'
#' In the landmark case (\code{Uk} provided), the grid is designed to be symmetric on the
#' bandwidth scale \eqn{\sigma} around \eqn{\sigma_0} over one decade:
#' \deqn{\sigma = \sigma_0 \times 10^{t}, \quad t \in \{-1,-2/3,-1/3,0,1/3,2/3,1\}.}
#' Using \eqn{\beta = 1/(2\sigma^2)}, this corresponds to
#' \deqn{\beta = \beta_0 \times 10^{-2t}.}
#'
#' When \code{Uk} is \code{NULL}, a shorter coarse grid may be returned (see \code{Value}).
#'
#' \strong{Notes:}
#' \itemize{
#'   \item When \code{Uk} is identical to \code{U}, the function detects this case and excludes
#'         self-distances (distance 0) to avoid \eqn{\sigma_0=0}.
#'   \item \code{sample_size} performs random subsampling without setting a seed. For reproducible
#'         results, set \code{set.seed()} before calling this function.
#' }
#'
#' @param U A numeric matrix of covariates (\code{K x N}); columns are samples.
#' @param Uk An optional numeric matrix of landmarks (\code{K x M}); columns are landmark points.
#'   If provided, distances are computed from samples in \code{U} to landmarks in \code{Uk}.
#' @param block_size Integer. Number of columns of \code{U} processed per block when computing
#'   distances (controls memory usage). If \code{N <= 1000}, it is automatically set to \code{N}.
#' @param block_size_Uk Integer. Number of columns of \code{Uk} processed per block when \code{Uk}
#'   is not \code{NULL} (controls memory usage). If \code{M <= 2000}, it is automatically set to \code{M}.
#' @param sample_size Integer or \code{NULL}. If not \code{NULL}, randomly subsamples this many columns
#'   of \code{U} (without replacement) before computing distances, to reduce computational cost.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{beta}: Estimated kernel parameter \eqn{\beta_0 = 1/(2\sigma_0^2)}.
#'   \item \code{beta_candidates}: Numeric vector of candidate \eqn{\beta} values (logarithmic grid)
#'         intended for cross-validation.
#'   \item \code{dist_median}: The estimated distance scale \eqn{\sigma_0} (median of nearest-neighbor
#'         or nearest-landmark distances).
#'   \item \code{block_size_used}: The effective block size(s) used. Either a scalar (no \code{Uk}) or
#'         a named vector \code{c(U=..., Uk=...)} when \code{Uk} is provided.
#'   \item \code{sample_size_used}: The number of columns of \code{U} actually used (after subsampling).
#'   \item \code{uk_is_u}: Logical flag indicating whether \code{Uk} was detected as identical to \code{U}
#'         (only returned when \code{Uk} is provided).
#' }
#'
#' @examples
#' # Basic (nearest-neighbor within U)
#' # beta_info <- nmfkc.kernel.beta.nearest.med(U)
#' # beta0 <- beta_info$beta
#' # betas <- beta_info$beta_candidates
#'
#' # With landmarks (nearest-landmark distances)
#' # beta_info <- nmfkc.kernel.beta.nearest.med(U, Uk)
#' # beta0 <- beta_info$beta
#' # betas <- beta_info$beta_candidates
#'
#' @export
nmfkc.kernel.beta.nearest.med <- function(
    U,
    Uk = NULL,
    block_size = 1000,
    block_size_Uk = 2000,
    sample_size = NULL
){
  U <- as.matrix(U)
  storage.mode(U) <- "double"
  N <- ncol(U)
  if (N < 2) stop("U must have at least 2 columns.")

  # ---- optional subsampling over U columns (for speed) ----
  if (!is.null(sample_size)) {
    sample_size <- as.integer(sample_size)
    if (sample_size <= 1) stop("sample_size must be >= 2")
    if (sample_size < N) {
      idxU <- sample.int(N, sample_size)
      U <- U[, idxU, drop = FALSE]
      N <- ncol(U)
    }
  }
  sample_size_used <- N

  # If Uk is NULL, behave like the original function: NN within U (exclude self).
  if (is.null(Uk)) {
    X  <- t(U)                    # N x K
    if (N <= 1000) block_size <- N
    XX <- rowSums(X * X)          # length N
    min_d2 <- rep(Inf, N)

    for (i in seq(1, N, by = block_size)) {
      i2 <- min(i + block_size - 1, N)
      Xi <- X[i:i2, , drop = FALSE]
      Xi_norm <- rowSums(Xi * Xi)

      dist2 <- outer(Xi_norm, XX, "+") - 2 * Xi %*% t(X)
      idx <- i:i2
      dist2[cbind(seq_along(idx), idx)] <- Inf  # exclude self
      dist2[dist2 < 0] <- 0

      nn_local <- apply(dist2, 1, min)
      min_d2[idx] <- pmin(min_d2[idx], nn_local)

      rm(Xi, Xi_norm, dist2); gc(FALSE)
    }

    d_med <- stats::median(sqrt(min_d2))
    if (d_med <= 0) {
      warning("Median nearest-neighbor distance is 0 (identical points exist). Using machine epsilon.")
      d_med <- sqrt(.Machine$double.eps)
    }
    beta  <- 1 / (2 * d_med^2)
    return(list(
      beta = beta,
      beta_candidates = beta * 10^c(-2:1),
      dist_median = d_med,
      block_size_used = block_size,
      sample_size_used = sample_size_used
    ))
  }

  # ---- Uk provided: nearest landmark distance (U vs Uk) ----
  Uk <- as.matrix(Uk)
  storage.mode(Uk) <- "double"
  if (nrow(Uk) != nrow(U)) stop("nrow(Uk) must equal nrow(U).")

  M <- ncol(Uk)
  if (M < 1) stop("Uk must have at least 1 column.")

  # We'll compute, for each column of U, min_j ||u_i - uk_j||^2.
  # Block over U (columns) and optionally over Uk to control memory.
  if (N <= 1000) block_size <- N
  if (M <= 2000) block_size_Uk <- M

  U2  <- colSums(U * U)           # length N
  Uk2 <- colSums(Uk * Uk)         # length M
  min_d2 <- rep(Inf, N)

  # Detect "Uk is U" case (same object / identical values) to exclude self-distance.
  uk_is_u <- (M == N) && isTRUE(all.equal(Uk, U, check.attributes = FALSE))

  for (i in seq(1, N, by = block_size)) {
    i2 <- min(i + block_size - 1, N)
    Ui <- U[, i:i2, drop = FALSE]
    Ui2 <- U2[i:i2]

    # we may need to accumulate min over Uk blocks
    cur_min <- rep(Inf, ncol(Ui))

    for (j in seq(1, M, by = block_size_Uk)) {
      j2 <- min(j + block_size_Uk - 1, M)
      Ukj <- Uk[, j:j2, drop = FALSE]
      Ukj2 <- Uk2[j:j2]

      # dist2: (block_U) x (block_Uk)
      # dist2 = ||u||^2 + ||uk||^2 - 2 uk^T u
      G <- crossprod(Ukj, Ui)  # (block_Uk) x (block_U)
      dist2 <- outer(Ui2, Ukj2, "+") - 2 * t(G)
      dist2[dist2 < 0] <- 0

      # exclude self-distances only when Uk == U
      if (uk_is_u) {
        # global indices in U: col(Ui) = (i..i2), col(Ukj) = (j..j2)
        # where they overlap, set that cell to Inf
        gi <- i:i2
        gj <- j:j2
        common <- intersect(gi, gj)
        if (length(common) > 0) {
          # map global index -> local positions
          li <- match(common, gi)
          lj <- match(common, gj)
          dist2[cbind(li, lj)] <- Inf
        }
      }

      cur_min <- pmin(cur_min, apply(dist2, 1, min))

      rm(Ukj, Ukj2, G, dist2); gc(FALSE)
    }

    min_d2[i:i2] <- pmin(min_d2[i:i2], cur_min)
    rm(Ui, Ui2, cur_min); gc(FALSE)
  }

  d_med <- stats::median(sqrt(min_d2))
  if (d_med <= 0) {
    warning("Median nearest-neighbor distance is 0 (identical points exist). Using machine epsilon.")
    d_med <- sqrt(.Machine$double.eps)
  }
  beta  <- 1 / (2 * d_med^2)

  t <- c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1)
  betas <- beta * 10^(-2*t)

  list(
    beta = beta,
    beta_candidates = betas,
    dist_median = d_med,
    block_size_used = c(U = block_size, Uk = block_size_Uk),
    sample_size_used = sample_size_used,
    uk_is_u = uk_is_u
  )
}





#' @title Optimize beta of the Gaussian kernel function by cross-validation
#' @description
#' \code{nmfkc.kernel.beta.cv} selects the optimal beta parameter of the kernel function by applying cross-validation over a set of candidate values.
#'
#' @param Y Observation matrix \eqn{Y(P,N)}.
#' @param Q Rank of the basis matrix.
#' @param U Covariate matrix \eqn{U(K,N) = (u_1, \dots, u_N)}. Each row may be normalized in advance.
#' @param V Covariate matrix \eqn{V(K,M) = (v_1, \dots, v_M)}, typically used for prediction. If \code{NULL}, the default is \code{U}.
#' @param beta A numeric vector of candidate kernel parameters to evaluate via cross-validation.
#' @param plot Logical. If TRUE (default), plots the objective function values for each candidate \code{beta}.
#' @param ... Additional arguments passed to \code{nmfkc.cv}.
#'
#' @return A list with components:
#' \item{beta}{The beta value that minimizes the cross-validation objective function.}
#' \item{objfunc}{Objective function values for each candidate \code{beta}.}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' nmfkc.kernel.beta.cv(Y,Q=1,U,V,beta=25:30/1000)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel.beta.cv <- function(Y,Q=2,U,V=NULL,beta=NULL,plot=TRUE,...){
  extra_args <- list(...)
  kernel_arg_names <- names(formals(nmfkc.kernel))
  cv_arg_names <- names(formals(nmfkc.cv))
  kernel_args_for_call <- extra_args[names(extra_args) %in% kernel_arg_names]
  cv_args_for_call <- extra_args[names(extra_args) %in% cv_arg_names]

  if(is.null(beta)){
    if(is.null(V)) V <- U
    med_args <- c(list(U = V), extra_args[names(extra_args) %in% names(formals(nmfkc.kernel.beta.nearest.med))])
    result.beta <- do.call("nmfkc.kernel.beta.nearest.med", med_args)
    beta <- result.beta$beta_candidates
    if (is.null(beta) || length(beta) == 0) stop("Failed to determine beta candidates from nearest-neighbor median.")
  }
  objfuncs <- numeric(length(beta))
  for(i in seq_along(beta)){
    start.time <- Sys.time()
    message(paste0("beta=",beta[i],"..."),appendLF=FALSE)

    kernel_args <- c(
      list(U = U, V = V, beta = beta[i]),
      kernel_args_for_call
    )
    A <- do.call("nmfkc.kernel", kernel_args)

    cv_args <- c(
      list(Y = Y, A = A, Q = Q),
      cv_args_for_call
    )
    result <- do.call("nmfkc.cv", cv_args)

    objfuncs[i] <- result$objfunc
    end.time <- Sys.time()
    diff.time <- difftime(end.time,start.time,units="sec")
    diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                           paste0(round(diff.time/60,1),"min"))
    message(diff.time.st)
  }
  i0 <- which.min(objfuncs)
  beta.best <- beta[i0]
  if(plot){
    plot(beta,objfuncs,type="l",col=2,xlab="beta",ylab="objfunc",log="x")
    graphics::points(beta[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(beta,objfuncs,format(beta,scientific=TRUE,digits=5))
  }
  names(objfuncs) <- beta
  result <- list(beta=beta.best,objfunc=objfuncs)
  return(result)
}






#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------




#' @title Initialize W (X) matrix using NNDSVDar (Vectorized)
#' @description
#'   Internal function to compute the NNDSVDar (Nonnegative Double SVD and Random)
#'   initialization for the basis matrix X. This method fills the zero entries
#'   generated by basic NNDSVD with small random values scaled by the matrix average,
#'   improving stability and convergence.
#' @param Y Input matrix (P x N)
#' @param Q Rank (number of components)
#' @return X (P x Q) non-negative initial basis matrix
#' @keywords internal
#' @noRd
.nndsvdar <- function(Y, Q) {
  P <- nrow(Y)
  N <- ncol(Y)
  s <- svd(Y, nu = Q, nv = Q)
  W <- matrix(0, P, Q)
  W[, 1] <- sqrt(s$d[1]) * abs(s$u[, 1])
  if (Q > 1) {
    idx <- 2:Q
    U <- s$u[, idx, drop = FALSE]
    V <- s$v[, idx, drop = FALSE]
    D_sqrt <- sqrt(s$d[idx])
    U_p <- pmax(U, 0)
    U_n <- pmax(-U, 0)
    V_p <- pmax(V, 0)
    V_n <- pmax(-V, 0)
    norm_Up <- sqrt(colSums(U_p^2))
    norm_Un <- sqrt(colSums(U_n^2))
    norm_Vp <- sqrt(colSums(V_p^2))
    norm_Vn <- sqrt(colSums(V_n^2))
    Mp <- norm_Up * norm_Vp
    Mn <- norm_Un * norm_Vn
    eps <- .Machine$double.eps
    norm_Up_safe <- norm_Up; norm_Up_safe[norm_Up == 0] <- eps
    norm_Un_safe <- norm_Un; norm_Un_safe[norm_Un == 0] <- eps
    W_pos <- sweep(U_p, 2, (D_sqrt * Mp / norm_Up_safe), FUN = "*")
    W_neg <- sweep(U_n, 2, (D_sqrt * Mn / norm_Un_safe), FUN = "*")
    use_pos <- (Mp > Mn)
    W_combined <- W_pos
    W_combined[, !use_pos] <- W_neg[, !use_pos]
    W[, idx] <- W_combined
  }
  W[is.nan(W)] <- 0
  W[W < 0] <- 0
  avg_Y <- mean(Y)
  idx_zero <- which(W == 0)
  n_zero <- length(idx_zero)
  fill_values <- stats::runif(n_zero) * avg_Y / 100
  W[idx_zero] <- fill_values
  return(W)
}





#' @title Compute a simplified/approximate silhouette coefficient (Internal)
#' @description
#' This internal function computes an approximate version of the silhouette coefficient
#' for clustering results. Unlike the standard definition (e.g., \code{cluster::silhouette}),
#' it determines the nearest neighboring cluster (\code{qn}) by comparing the
#' **Euclidean distances between each sample and the cluster mean vectors (centroids)**,
#' rather than comparing the average distance to all members of other clusters.
#'
#' However, the calculation of \code{a(i)} (average distance to members of the
#' assigned cluster) and \code{b(i)} (average distance to members of the nearest
#' neighboring cluster) uses the pairwise distances of all samples, consistent
#' with the standard definition. This combined approach reduces the computational
#' complexity of the \code{qn} determination step while maintaining consistency
#' for \code{a(i)} and \code{b(i)} using pre-computed distance matrix \code{D}.
#'
#' @param B.prob The coefficient probability matrix (Q x N).
#' @param B.cluster The hard clustering vector derived from \code{B.prob}.
#'
#' @return A list containing:
#' \item{cluster}{Sorted cluster labels.}
#' \item{silhouette}{Sorted silhouette values.}
#' \item{silhouette.means}{Mean of the cluster indices.}
#' \item{silhouette.mean}{Overall mean silhouette score.}
#' @keywords internal
#' @noRd
.silhouette.simple <- function(B.prob,B.cluster){
  if(is.matrix(B.prob)){Q <- nrow(B.prob)}else{Q <- 1}
  if(Q==1){
    return(list(cluster=NA,silhouette=NA,
                silhouette.means=NA,silhouette.mean=NA))
  }else{
    index <-!is.na(B.cluster)
    B.prob <- B.prob[,index, drop=FALSE]
    B.cluster <- B.cluster[index]
    N_samples <- ncol(B.prob)
    D <- as.matrix(stats::dist(t(B.prob)))
    cluster.means <- matrix(0,nrow=Q,ncol=Q)
    ns <- NULL
    cluster.list <- NULL
    for(q in 1:Q){
      ns <- c(ns,sum(B.cluster==q))
      cluster.list <- c(cluster.list,list(which(B.cluster==q)))
      cluster.means[,q] <- rowMeans(B.prob[,B.cluster==q,drop=FALSE])
    }
    B_prob_sq_sum <- colSums(B.prob^2)
    cluster_means_sq_sum <- colSums(cluster.means^2)
    d2_matrix <- outer(B_prob_sq_sum, rep(1, Q)) +
      outer(rep(1, N_samples), cluster_means_sq_sum) -
      2 * t(B.prob) %*% cluster.means
    d2_matrix[d2_matrix < 0] <- 0
    si <- 0*B.cluster
    neighbor.cluster <- 0*B.cluster
    for(q in 1:Q){
      q_indices <- cluster.list[[q]]
      for(i_idx in seq_along(q_indices)){
        i <- q_indices[i_idx]
        di <- d2_matrix[i, ]
        qn <- order(di)[ifelse(order(di)[1]==q, 2, 1)]
        neighbor.cluster[i] <- qn
        if(ns[q]==1){
          si[i] <- 0
        }else{
          ai_distances <- D[i, q_indices]
          ai <- sum(ai_distances) / (ns[q]-1)
          qn_indices <- cluster.list[[qn]]
          bi_distances <- D[i, qn_indices]
          bi <- sum(bi_distances) / ns[qn]
          si[i] <- (bi-ai)/max(ai,bi)
        }
      }
    }
    si.mean <- mean(si)
    si.sort.cluster.means <- 0*ns
    for(q in 1:Q){
      si.sort.cluster.means[q] <- mean(cluster.list[[q]])
    }
    si.sort <- NULL
    si.sort.cluster <- NULL
    for(q in 1:Q){
      si.sort <- c(si.sort,sort(si[cluster.list[[q]]],decreasing=TRUE))
      si.sort.cluster <- c(si.sort.cluster,rep(q,length(cluster.list[[q]])))
    }
    return(list(cluster=si.sort.cluster,silhouette=si.sort,
                silhouette.means=si.sort.cluster.means,silhouette.mean=si.mean))
  }
}





#' @title Parse formula and prepare Y and A matrices
#' @description
#' Internal function to handle formula input, parse variables, and generate
#' the transposed Y and A matrices (features/covariates in rows, samples in columns)
#' as required by the core \code{nmfkc} function. Supports 'data' mode
#' (variable names and dot notation) and direct expression evaluation mode.
#' @param formula A formula object.
#' @param data Optional data frame or environment.
#' @return A list containing the prepared matrices: \code{Y} (P x N) and \code{A} (R x N or NULL).
#' @keywords internal
#' @noRd
.nmfkc_parse_formula <- function(formula, data) {

  # Helper function to abbreviate column messages (shows max_show columns)
  .abbreviate_msg <- function(cols, prefix, max_show = 5) {
    total_len <- base::length(cols)
    if (total_len > max_show) {
      msg <- base::paste0(base::paste(cols[1:max_show], collapse = ", "),
                          base::sprintf("... (Total: %d)", total_len))
    } else {
      msg <- base::paste(cols, collapse = ", ")
    }
    base::message(base::paste(prefix, "created using columns:", msg))
  }

  # Capture the formula object and set the environment
  f <- formula
  f_env <- if (base::missing(data)) base::environment(f) else base::environment()

  # Extract the expressions for Y and A
  Y_expr <- f[[2]]
  A_expr <- f[[3]]

  # --- Input Mode Branching ---

  if (!base::missing(data)) {
    # MODE 1: data provided (Variable names and dot notation)
    data <- base::as.data.frame(data)
    all_cols <- base::names(data)

    # Check if A_expr is structurally present in the formula object.
    A_is_structurally_present <- base::length(f) >= 3

    # Extract Y and A expressions as clean strings using deparse
    Y_expr_str <- base::paste(base::deparse(Y_expr), collapse = " ")

    A_expr_str <- if (A_is_structurally_present) {
      base::paste(base::deparse(A_expr), collapse = " ")
    } else {
      NA # Use simple NA for missing structural part
    }

    # Determine dot notation and missing status based on the safe string A_expr_str
    Y_is_dot <- (Y_expr_str == ".")
    A_is_missing <- base::is.na(A_expr_str)
    A_is_dot <- (!A_is_missing && A_expr_str == ".")

    # --- Check for explicit A omission symbols (0 or -1) ---
    A_is_explicitly_omitted <- (!A_is_missing &&
                                  (A_expr_str == "0" || A_expr_str == "-1" ||
                                     A_expr_str == "0 + ." || A_expr_str == "-1 + ." ||
                                     A_expr_str == ". + 0" || A_expr_str == ". + -1"))

    # Update A_is_missing status if explicitly omitted
    A_is_missing <- A_is_missing || A_is_explicitly_omitted

    Y_cols <- NULL
    A_cols <- NULL

    # Helper to clean and split variable names from expression string
    .extract_cols <- function(expr_str) {
      # Split by '+', '~', ' ' and remove empty elements. This handles "A1 + A2" structure.
      cols <- base::unlist(base::strsplit(expr_str, "[ +~]", fixed = FALSE))
      cols <- cols[cols != ""]
      # Remove explicit omission symbols if they were part of the split result
      cols <- cols[!(cols %in% base::c("0", "-1"))]
      return(cols)
    }

    if (!Y_is_dot) { Y_cols <- .extract_cols(Y_expr_str) }

    # A_cols is only extracted if A is NOT missing and NOT explicitly omitted (and not dot)
    if (!A_is_missing && !A_is_dot) {
      A_cols <- .extract_cols(A_expr_str)
    }

    # Error Check: Both sides cannot be '.'
    if (Y_is_dot && A_is_dot) {
      base::stop("Formula error: '.' cannot be used on both the left (Y) and right (A) sides simultaneously when 'data' is provided.")
    }

    used_cols <- base::unique(base::c(Y_cols, A_cols))
    remaining_cols <- base::setdiff(all_cols, used_cols)

    # Assign dot notation variables
    if (Y_is_dot) {
      Y_cols <- remaining_cols
      if (base::length(Y_cols) == 0) { base::stop("Formula error: '.' for Y resulted in no remaining variables.") }
      .abbreviate_msg(Y_cols, "Y")
    } else if (A_is_dot) {
      A_cols <- remaining_cols
      if (base::length(A_cols) == 0) { base::stop("Formula error: '.' for A resulted in no remaining variables.") }
      # NOTE: Message output for A is suppressed here to avoid duplication
      #       and is handled in the A_mat creation block below.
    }

    # Validate column existence
    missing_Y <- base::setdiff(Y_cols, all_cols)
    if (base::length(missing_Y) > 0) {
      base::stop(base::paste0("Formula error: Y columns not found in data: ",
                               base::paste(missing_Y, collapse = ", ")))
    }
    if (!A_is_missing) {
      missing_A <- base::setdiff(A_cols, all_cols)
      if (base::length(missing_A) > 0) {
        base::stop(base::paste0("Formula error: A columns not found in data: ",
                                 base::paste(missing_A, collapse = ", ")))
      }
    }

    # Matrix Creation
    Y_mat <- data[, Y_cols, drop = FALSE]

    # --- A Matrix Finalization and Message Output ---
    if (A_is_missing) { # Catches structural omission, NA, or explicit 0/-1
      A_mat <- NULL
      base::message("A (covariate matrix) is omitted. Performing standard NMF (Y ~ X B).")
    } else {
      # This block now runs if A is explicitly defined with variables or is dot notation
      A_mat <- data[, A_cols, drop = FALSE]
      .abbreviate_msg(A_cols, "A") # Message output for A is performed here ONCE
    }

  } else {
    # MODE 2: data omitted (Direct matrix expression evaluation)
    Y_cols <- NULL
    A_cols <- NULL

    if (base::is.symbol(Y_expr) && base::as.character(Y_expr) == ".") {
      base::stop("Formula error: '.' is not supported for direct matrix evaluation mode (without 'data' argument).")
    }

    Y_mat <- base::tryCatch({ base::as.matrix(base::eval(Y_expr, envir = f_env)) },
                            error = function(e) { base::stop(base::paste("Error evaluating Y expression:", base::conditionMessage(e))) })

    A_mat <- base::tryCatch({
      # Check if A_expr exists (length >= 3) or is explicitly 0/-1
      A_is_structurally_present <- base::length(f) >= 3
      if (!A_is_structurally_present || (base::is.numeric(A_expr) && A_expr == 0) || (base::is.numeric(A_expr) && A_expr == -1)) { NULL
      } else { base::as.matrix(base::eval(A_expr, envir = f_env)) }
    }, error = function(e) { base::stop(base::paste("Error evaluating A expression:", base::conditionMessage(e))) })

    if (base::is.null(A_mat)) {
      base::message("A (covariate matrix) is omitted. Performing standard NMF (Y ~ X B).")
    }
    # NOTE: No column abbreviation message needed for Mode 2 as matrices are evaluated directly.
  }

  # Transpose: R rows=samples -> nmfkc columns=samples (P x N or R x N)
  # NOTE: R data frames/matrices are typically samples x variables. nmfkc requires variables x samples.
  Y <- base::t(Y_mat)

  if (!base::is.null(A_mat)) {
    A <- base::t(A_mat)
    if (base::ncol(Y) != base::ncol(A)) {
      base::stop(base::paste0("Dimension error: Number of columns (samples) in Y (", base::ncol(Y), ") must match number of columns (samples) in A (", base::ncol(A), ")."))
    }
  } else {
    A <- NULL
  }

  return(base::list(Y = Y, A = A, Y_cols = Y_cols, A_cols = A_cols))
}


#' Resolve formula input to Y/A matrices with metadata
#'
#' Internal helper that detects formula input, parses it, and returns
#' the Y/A matrices along with formula metadata for downstream use.
#' The \code{data_missing} flag solves the \code{missing(data)} propagation
#' problem in R (cannot check \code{missing()} through nested calls).
#'
#' @param Y A formula or matrix.
#' @param A A matrix or NULL.
#' @param data_missing Logical; TRUE if \code{data} was not supplied by the caller.
#' @param data_value The value of \code{data} (or NULL if missing).
#' @return A list with \code{Y}, \code{A}, and \code{formula.meta} (NULL if not formula mode).
#' @keywords internal
#' @noRd
.nmfkc_resolve_formula <- function(Y, A, data_missing, data_value) {
  if (!base::inherits(Y, "formula")) {
    return(base::list(Y = Y, A = A, formula.meta = NULL))
  }
  if (data_missing) {
    data_list <- .nmfkc_parse_formula(formula = Y)
  } else {
    data_list <- .nmfkc_parse_formula(formula = Y, data = data_value)
  }
  formula.meta <- base::list(
    formula = Y,
    Y_cols  = data_list$Y_cols,
    A_cols  = data_list$A_cols
  )
  base::list(Y = data_list$Y, A = data_list$A, formula.meta = formula.meta)
}










#' @title Optimize NMF with kernel covariates (Full Support for Missing Values)
#' @description
#' \code{nmfkc} fits a nonnegative matrix factorization with kernel covariates
#' under the tri-factorization model \eqn{Y \approx X C A = X B}.
#'
#' This function supports two major input modes:
#' 1. **Matrix Mode (Existing)**: \code{nmfkc(Y=matrix, A=matrix, ...)}
#' 2. **Formula Mode (New)**: \code{nmfkc(formula=Y_vars ~ A_vars, data=df, rank=Q, ...)}
#'
#' The rank of the basis matrix can be specified using either the \code{rank} argument
#' (preferred for formula mode) or the hidden \code{Q} argument (for backward compatibility).
#'
#' @param Y Observation matrix (P x N), OR a formula object for Formula Mode.
#'   In Formula Mode, use \code{Y1 + Y2 ~ A1 + A2} with \code{data}, or
#'   \code{Y_matrix ~ A_matrix} for direct matrix evaluation.
#'   Supports dot notation (\code{. ~ A1 + A2}) when \code{data} is supplied.
#' @param A Covariate matrix. Default is \code{NULL} (no covariates).
#'   Ignored when \code{Y} is a formula.
#' @param rank Integer. The rank of the basis matrix \eqn{X} (Q). Preferred over \code{Q}.
#' @param data Optional. A data frame from which variables in the formula should be taken.
#' @param epsilon Positive convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @param ... Additional arguments passed for fine-tuning regularization, initialization, constraints,
#'   and output control. This includes the backward-compatible arguments \code{Q} and \code{method}.
#'   \itemize{
#'     \item \code{Y.weights}: Optional numeric matrix (P x N) or vector (length N).
#'       0 indicates missing/ignored values. If NULL (default), weights are automatically
#'       set to 0 for NAs in Y, and 1 otherwise.
#'     \item \code{X.L2.ortho}: Nonnegative penalty parameter for the orthogonality of \eqn{X} (default: 0).
#'       It minimizes the off-diagonal elements of the Gram matrix \eqn{X^\top X}, reducing the correlation
#'       between basis vectors (conceptually minimizing \eqn{\| X^\top X - \mathrm{diag}(X^\top X) \|_F^2}).
#'       (Formerly \code{lambda.ortho}).
#'     \item \code{B.L1}: Nonnegative penalty parameter for L1 regularization on \eqn{B = C A} (default: 0).
#'       Promotes **sparsity in the coefficients**. (Formerly \code{gamma}).
#'     \item \code{C.L1}: Nonnegative penalty parameter for L1 regularization on \eqn{C} (default: 0).
#'       Promotes **sparsity in the parameter matrix**. (Formerly \code{lambda}).
#'     \item \code{Q}: Backward-compatible name for the rank of the basis matrix (Q).
#'     \item \code{method}: Objective function: Euclidean distance \code{"EU"} (default) or Kullback–Leibler divergence \code{"KL"}.
#'     \item \code{X.restriction}: Constraint for columns of \eqn{X}. Options: \code{"colSums"} (default), \code{"colSqSums"}, \code{"totalSum"}, or \code{"fixed"}.
#'     \item \code{X.init}: Method for initializing the basis matrix \eqn{X}. Options: \code{"kmeans"} (default), \code{"runif"}, \code{"nndsvd"}, or a user-specified matrix.
#'     \item \code{nstart}: Number of random starts for \code{kmeans} when initializing \eqn{X} (default: 1).
#'     \item \code{seed}: Integer seed for reproducibility (default: 123).
#'     \item \code{C.init}: Optional numeric matrix giving the initial value of the parameter matrix \eqn{C}
#'       (i.e., \eqn{\Theta}). If \code{A} is \code{NULL}, \code{C} has dimension \eqn{Q \times N} (equivalently \eqn{B});
#'       otherwise, \code{C} has dimension \eqn{Q \times K} where \eqn{K = nrow(A)}. Default initializes all entries to 1.
#'     \item \code{prefix}: Prefix for column names of \eqn{X} and row names of \eqn{B} (default: "Basis").
#'     \item \code{print.trace}: Logical. If \code{TRUE}, prints progress every 10 iterations (default: \code{FALSE}).
#'     \item \code{print.dims}: Logical. If \code{TRUE} (default), prints matrix dimensions and elapsed time.
#'     \item \code{save.time}: Logical. If \code{TRUE} (default), skips some post-computations (e.g., CPCC, silhouette) to save time.
#'     \item \code{save.memory}: Logical. If \code{TRUE}, performs only essential computations (implies \code{save.time = TRUE}) to reduce memory usage (default: \code{FALSE}).
#'   }
#' @return A list with components:
#' \item{call}{The matched call, as captured by `match.call()`.}
#' \item{dims}{A character string summarizing the matrix dimensions of the model.}
#' \item{runtime}{A character string summarizing the computation time.}
#' \item{X}{Basis matrix. Column normalization depends on \code{X.restriction}.}
#' \item{B}{Coefficient matrix \eqn{B = C A}.}
#' \item{XB}{Fitted values for \eqn{Y}.}
#' \item{C}{Parameter matrix.}
#' \item{B.prob}{Soft-clustering probabilities derived from columns of \eqn{B}.}
#' \item{B.cluster}{Hard-clustering labels (argmax over \eqn{B.prob} for each column).}
#' \item{X.prob}{Row-wise soft-clustering probabilities derived from \eqn{X}.}
#' \item{X.cluster}{Hard-clustering labels (argmax over \eqn{X.prob} for each row).}
#' \item{A.attr}{List of attributes of the input covariate matrix \code{A}, containing metadata like lag order and intercept status if created by \code{nmfkc.ar} or \code{nmfkc.kernel}.}
#' \item{formula.meta}{If fitted via Formula Mode, a list with \code{formula}, \code{Y_cols}, and \code{A_cols}; otherwise \code{NULL}.}
#' \item{objfunc}{Final objective value.}
#' \item{objfunc.iter}{Objective values by iteration.}
#' \item{r.squared}{Coefficient of determination \eqn{R^2} between \eqn{Y} and \eqn{X B}.}
#' \item{method}{Character string indicating the optimization method used (\code{"EU"} or \code{"KL"}).}
#' \item{n.missing}{Number of missing (or zero-weighted) elements in \eqn{Y}.}
#' \item{n.total}{Total number of elements in \eqn{Y}.}
#' \item{rank}{The rank \eqn{Q} used in the factorization.}
#' \item{sigma}{The residual standard error, representing the typical deviation of the observed values \eqn{Y} from the fitted values \eqn{X B}.}
#' \item{mae}{Mean Absolute Error between \eqn{Y} and \eqn{X B}.}
#' \item{criterion}{A list of selection criteria, including \code{ICp} (\code{ICp1}, \code{ICp2}, \code{ICp3}), \code{CPCC}, \code{silhouette}, \code{AIC}, \code{BIC}, \code{dist.cor}, \code{B.prob.sd.min}, \code{B.prob.max.mean}, and \code{B.prob.entropy.mean}.}
#' @seealso \code{\link{nmfkc.cv}}, \code{\link{nmfkc.rank}}, \code{\link{nmfkc.kernel}}, \code{\link{nmfkc.ar}}, \code{\link{predict.nmfkc}}
#' @export
#' @source Satoh, K. (2024). Applying Non-negative Matrix Factorization with Covariates
#'   to the Longitudinal Data as Growth Curve Model. arXiv:2403.05359.
#'   \url{https://arxiv.org/abs/2403.05359}
#' @references
#' Ding, C., Li, T., Peng, W., & Park, H. (2006). Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering.
#'   In \emph{Proceedings of the 12th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining} (pp. 126–135).
#'   \doi{10.1145/1150402.1150420}
#' Potthoff, R. F., & Roy, S. N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems.
#'   \emph{Biometrika}, 51, 313–326. \doi{10.2307/2334137}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1. Matrix Mode (Existing)
#' library(nmfkc)
#' X <- cbind(c(1,0,1),c(0,1,0))
#' B <- cbind(c(1,0),c(0,1),c(1,1))
#' Y <- X %*% B
#' rownames(Y) <- paste0("P",1:nrow(Y))
#' colnames(Y) <- paste0("N",1:ncol(Y))
#' print(X); print(B); print(Y)
#' library(nmfkc)
#' res <- nmfkc(Y,Q=2,epsilon=1e-6)
#' res$X
#' res$B
#'
#' # Example 2. Formula Mode
#' set.seed(1)
#' dummy_data <- data.frame(Y1=rpois(10,5), Y2=rpois(10,10),
#'                          A1=abs(rnorm(10,5)), A2=abs(rnorm(10,3)))
#' res_f <- nmfkc(Y1 + Y2 ~ A1 + A2, data=dummy_data, rank=2)
#'
nmfkc <- function(Y, A=NULL, rank=NULL, data, epsilon=1e-4, maxit=5000, ...){
  # A small constant for numerical stability to prevent division by zero and log(0).
  .eps <- 1e-10

  extra_args <- base::list(...)

  # --- 1. Parameter Extraction ---
  Q_hidden <- if (!base::is.null(extra_args$Q)) extra_args$Q else NULL
  Q_val <- if (!base::is.null(rank)) rank else if (!base::is.null(Q_hidden)) Q_hidden else 2
  Q <- Q_val

  C.L1 <- if (!base::is.null(extra_args$C.L1)) extra_args$C.L1 else 0
  B.L1 <- if (!base::is.null(extra_args$B.L1)) extra_args$B.L1 else 0
  X.L2.ortho <- if (!base::is.null(extra_args$X.L2.ortho)) extra_args$X.L2.ortho else 0

  if (C.L1 == 0 && !base::is.null(extra_args$lambda)) C.L1 <- extra_args$lambda
  if (B.L1 == 0 && !base::is.null(extra_args$gamma)) B.L1 <- extra_args$gamma
  if (X.L2.ortho == 0 && !base::is.null(extra_args$lambda.ortho)) X.L2.ortho <- extra_args$lambda.ortho

  method <- if (!base::is.null(extra_args$method)) extra_args$method else "EU"
  X.restriction <- if (!base::is.null(extra_args$X.restriction)) extra_args$X.restriction else "colSums"
  X.init <- if (!base::is.null(extra_args$X.init)) extra_args$X.init else "kmeans"
  nstart <- if (!base::is.null(extra_args$nstart)) extra_args$nstart else 1
  seed <- if (!base::is.null(extra_args$seed)) extra_args$seed else 123
  C.init <- if (!is.null(extra_args$C.init)) extra_args$C.init else NULL

  prefix <- if (!base::is.null(extra_args$prefix)) extra_args$prefix else "Basis"
  print.trace <- if (!base::is.null(extra_args$print.trace)) extra_args$print.trace else FALSE
  print.dims <- if (!base::is.null(extra_args$print.dims)) extra_args$print.dims else TRUE
  save.time <- if (!base::is.null(extra_args$save.time)) extra_args$save.time else TRUE
  save.memory <- if (!base::is.null(extra_args$save.memory)) extra_args$save.memory else FALSE

  Y.weights <- if (!base::is.null(extra_args$Y.weights)) extra_args$Y.weights else NULL

  # --- 2. Input Data Preparation ---
  formula.meta <- NULL
  if (base::inherits(Y, "formula")) {
    data_missing <- base::missing(data)
    resolved <- .nmfkc_resolve_formula(Y, A, data_missing, if (!data_missing) data else NULL)
    Y <- resolved$Y
    A <- resolved$A
    formula.meta <- resolved$formula.meta
  } else {
    if(base::is.vector(Y)) Y <- base::matrix(Y,nrow=1)
    if(!base::is.matrix(Y)) Y <- base::as.matrix(Y)
  }

  # --- Input Validation (after formula dispatch) ---
  if(!base::is.null(A)) {
    if(any(is.na(A))) base::stop("Covariate matrix A contains NAs. Please impute or remove them.")
    if(base::min(A, na.rm=TRUE)<0) base::stop("The matrix A should be non-negative.")
  }
  if(base::min(Y, na.rm=TRUE)<0) base::stop("The matrix Y should be non-negative.")

  # === Weights Handling ===
  # 1. Vector Expansion (Column-wise weights)
  if (!is.null(Y.weights) && is.vector(Y.weights)) {
    if (length(Y.weights) == ncol(Y)) {
      # Expand vector to matrix: each row gets the same weight vector
      Y.weights <- matrix(Y.weights, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
    } else if (length(Y.weights) == 1) {
      Y.weights <- matrix(Y.weights, nrow = nrow(Y), ncol = ncol(Y))
    } else {
      stop("Length of Y.weights vector must match ncol(Y) (or be 1).")
    }
  }

  # 2. Check Dimensions & Handle NAs
  if (is.null(Y.weights)) {
    if (any(is.na(Y))) {
      Y.weights <- matrix(1, nrow=nrow(Y), ncol=ncol(Y))
      Y.weights[is.na(Y)] <- 0
      Y[is.na(Y)] <- 0
      if(print.dims) message("Notice: Missing values (NA) in Y were treated as weights=0.")
    } else {
      Y.weights <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))
    }
  } else {
    if (!is.matrix(Y.weights)) Y.weights <- as.matrix(Y.weights)
    if (!all(dim(Y.weights) == dim(Y))) stop("Dimension mismatch between Y and Y.weights.")
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y)] <- 0
    Y[Y.weights == 0] <- 0
  }

  # --- 3. Algorithm Setup ---
  X.restriction <- base::match.arg(X.restriction, base::c("colSums", "colSqSums", "totalSum","fixed"))
  xnorm <- base::switch(X.restriction,
                        colSums   = function(X) base::sweep(X, 2, base::colSums(X), "/"),
                        colSqSums = function(X) base::sweep(X, 2, base::sqrt(base::colSums(X^2)), "/"),
                        totalSum  = function(X) X / base::sum(X),
                        fixed = function(X) X
  )

  if(base::is.null(A)){
    dims <- base::sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)",
                          base::nrow(Y),base::ncol(Y),base::nrow(Y),Q,Q,base::ncol(Y))
  }else{
    dims <- base::sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)",
                          base::nrow(Y),base::ncol(Y),base::nrow(Y),Q,Q,base::nrow(A),base::nrow(A),base::ncol(Y),Q,base::ncol(Y))
  }
  if(print.dims) base::message(base::paste0(dims,"..."),appendLF=FALSE)
  start.time <- base::Sys.time()

  # Initialize X
  is.X.scalar <- FALSE
  if(nrow(Y)>=2){
    # [Fix 2] Create Y_init for initialization (impute NAs with row means)
    Y_init <- Y
    if (is.matrix(Y.weights) && any(Y.weights == 0)) {
      # Calculate weighted mean (Since Y is zero-filled, sum(Y) is equivalent to sum(Y*W))
      row_means <- rowSums(Y) / (rowSums(Y.weights) + .eps)
      mask_binary <- (Y.weights > 0)
      idx_missing <- which(!mask_binary, arr.ind = TRUE)
      if (nrow(idx_missing) > 0) {
        Y_init[idx_missing] <- row_means[idx_missing[,1]]
      }
    }

    if (is.matrix(X.init)) {
      X <- X.init
    } else if (is.character(X.init)) {
      if (!is.null(seed)) set.seed(seed)
      if (X.init == "kmeans") {
        # kmeans requires ncol(Y) >= Q (N >= Q)
        res.kmeans <- if (ncol(Y) >= Q) {
          tryCatch(stats::kmeans(t(Y_init), centers = Q, iter.max = maxit, nstart = nstart),
                   error = function(e) NULL)
        } else { NULL }
        if(!is.null(res.kmeans)) X <- t(res.kmeans$centers) else X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
      } else if (X.init == "runif") {
        X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
      } else {
        # nndsvd requires Q <= min(P, N)
        if (Q > min(nrow(Y), ncol(Y))) {
          X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
        } else {
          X <- .nndsvdar(Y_init, Q)
        }
      }
    } else if (ncol(Y) == Q) {
      X <- Y_init
    } else {
      X <- matrix(stats::runif(nrow(Y) * Q), nrow = nrow(Y), ncol = Q)
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  X <- xnorm(X)

  # [FIX: Initialization of tX]
  # Initialize tX here so it exists even if the X update loop is skipped (e.g., scalar X)
  tX <- t(X)

  if(is.null(A)){
    if(is.null(C.init)) C <- matrix(1, nrow=Q, ncol=ncol(Y)) else C <- C.init
  }else{
    if(is.null(C.init)) C <- matrix(1, nrow=Q, ncol=nrow(A)) else C <- C.init
  }
  hasA <- !is.null(A)

  ones_QN <- matrix(1, nrow=Q, ncol=ncol(Y))
  if(hasA) {
    At <- t(A)
    ones_QN_At <- ones_QN %*% At
  }

  epsilon.iter <- Inf
  objfunc.iter <- 0*(1:maxit)
  i_end <- NULL

  # --- 4. Main Loop (Weighted) ---
  for(i in 1:maxit){
    if(is.null(A)) B <- C else B <- C %*% A
    XB <- X %*% B
    if(print.trace && i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))

    if(method=="EU"){
      if(!is.X.scalar && X.restriction!="fixed"){
        num_X <- (Y.weights * Y) %*% t(B)
        den_X <- (Y.weights * XB) %*% t(B)
        if (X.L2.ortho > 0) {
          XtX <- crossprod(X); diag(XtX) <- 0
          den_X <- den_X + X.L2.ortho * (X %*% XtX)
        }
        X <- X * (num_X / (den_X + .eps))
        X <- xnorm(X)
        tX <- t(X)
      }
      if(is.null(A)) {
        num_C <- tX %*% (Y.weights * Y)
        den_C <- tX %*% (Y.weights * XB)
        if (C.L1 != 0) den_C <- den_C + (C.L1/2) * ones_QN
        if (B.L1 != 0) den_C <- den_C + (B.L1/2) * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        num_C <- tX %*% (Y.weights * Y) %*% At
        den_C <- tX %*% (Y.weights * XB) %*% At
        if (C.L1 != 0) den_C <- den_C + (C.L1/2) * matrix(1, nrow=Q, ncol=nrow(A))
        if (B.L1 != 0) den_C <- den_C + (B.L1/2) * ones_QN_At
        C <- C * (num_C / (den_C + .eps))
      }
      resid <- Y.weights * (Y - XB)
      obj <- sum(resid^2)

    }else{ # KL
      if(!is.X.scalar && X.restriction!="fixed"){
        ratio <- Y.weights * (Y / (XB + .eps))
        num_X <- ratio %*% t(B)
        den_X <- Y.weights %*% t(B)
        if (X.L2.ortho > 0) {
          XtX <- crossprod(X); diag(XtX) <- 0
          den_X <- den_X + X.L2.ortho * (X %*% XtX)
        }
        X <- X * (num_X / (den_X + .eps))
        X <- xnorm(X)
        tX <- t(X)
      }
      if(is.null(A)) {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- tX %*% ratio
        den_C <- tX %*% Y.weights
        if (C.L1 != 0) den_C <- den_C + C.L1 * ones_QN
        if (B.L1 != 0) den_C <- den_C + B.L1 * ones_QN
        C <- C * (num_C / (den_C + .eps))
      } else {
        ratio <- Y.weights * (Y / (XB + .eps))
        num_C <- tX %*% ratio %*% At
        den_C <- tX %*% Y.weights %*% At
        if (C.L1 != 0) den_C <- den_C + C.L1 * matrix(1, nrow=Q, ncol=nrow(A))
        if (B.L1 != 0) den_C <- den_C + B.L1 * ones_QN_At
        C <- C * (num_C / (den_C + .eps))
      }
      term1 <- - (Y.weights * Y) * log(XB + .eps)
      term2 <- Y.weights * XB
      obj <- sum(term1 + term2)
    }

    if (C.L1 != 0) obj <- obj + C.L1 * sum(C)
    if (B.L1 != 0) obj <- obj + if (hasA) B.L1 * sum(C %*% A) else B.L1 * sum(C)
    if (X.L2.ortho != 0) {
      XtX <- crossprod(X); diag(XtX) <- 0
      obj <- obj + (X.L2.ortho / 2) * sum(XtX^2)
    }
    objfunc.iter[i] <- obj

    if(i>=10){
      #epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1]) / pmax(abs(objfunc.iter[i]), 1)
      if(epsilon.iter <= abs(epsilon)){ i_end <- i; break }
    }
  }

  if(is.null(A)) B <- C else B <- C %*% A
  XB <- X %*% B

  if(method=="EU"){
    resid <- Y.weights * (Y - XB)
    objfunc <- sum(resid^2)
  } else {
    term1 <- - (Y.weights * Y) * log(XB + .eps)
    term2 <- Y.weights * XB
    objfunc <- sum(term1 + term2)
  }

  if(!is.null(i_end)){ objfunc.iter <- objfunc.iter[10:i_end]
  } else if (i >= 10){ objfunc.iter <- objfunc.iter[10:i]
  } else { objfunc.iter <- objfunc.iter[1:i] }

  if(ncol(X) > 1 && X.restriction != "fixed"){
    index <- order(matrix(1:nrow(X)/nrow(X),nrow=1) %*% X)
    X <- X[,index,drop=FALSE]; B <- B[index,,drop=FALSE]; C <- C[index,,drop=FALSE]
  }
  rownames(C) <- paste0(prefix,1:nrow(C))
  rownames(X) <- rownames(Y); colnames(X) <- paste0(prefix,1:ncol(X))
  rownames(B) <- paste0(prefix,1:nrow(B)); colnames(B) <- colnames(Y)

  N_obs <- sum(Y.weights > 0)

  if(X.restriction=="fixed") nparam.X <- 0 else if(X.restriction=="totalSum") nparam.X <- prod(dim(X)) - 1 else nparam.X <- (nrow(X)-1)*ncol(X)
  if(is.null(A)) nparam <- nparam.X + prod(dim(B)) else nparam <- nparam.X + prod(dim(C))

  # Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
  if (method == "EU") {
    sigma2 <- objfunc / N_obs
    P_dim <- nrow(Y)
    T_dim <- ncol(Y)
    g_icp1 <- (P_dim + T_dim) / (P_dim * T_dim) * log(P_dim * T_dim)
    g_icp2 <- (P_dim + T_dim) / (P_dim * T_dim) * log(min(P_dim, T_dim))
    g_icp3 <- log(min(P_dim, T_dim)) / min(P_dim, T_dim)
    ICp1 <- log(sigma2) + nparam * g_icp1
    ICp2 <- log(sigma2) + nparam * g_icp2
    ICp3 <- log(sigma2) + nparam * g_icp3
    ICp <- min(ICp1, ICp2, ICp3)
    AIC <- N_obs * log(sigma2) + 2 * nparam
    BIC <- N_obs * log(sigma2) + nparam * log(N_obs)
  } else {
    ICp <- NA;ICp1 <- NA;ICp2 <-NA;ICp3 <- NA;
    AIC <- NA
    BIC <- NA
  }

  if(save.memory==FALSE){
    # [Fix 3] Calculate statistics using valid data only
    valid_idx <- (Y.weights > 0)
    if(any(valid_idx)) {
      r2 <- stats::cor(XB[valid_idx], Y[valid_idx])^2
      sigma <- stats::sd(Y[valid_idx] - XB[valid_idx])
      mae <- mean(abs(Y[valid_idx] - XB[valid_idx]))
    } else { r2 <- NA; sigma <- NA; mae <- NA }

    B.prob <- t( t(B) / (colSums(B) + .eps) )
    if(Q > 1){
      B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
      B.prob.max.mean <- mean(apply(B.prob,2,max))
      p <- B.prob + .eps # avoid log(0)
      B.prob.entropy.mean <- -mean(colSums(p * log(p))) / log(Q)
    } else {
      B.prob.sd.min <- 0
      B.prob.entropy.mean <- 0
      B.prob.max.mean <- 1
    }
    B.cluster <- apply(B.prob,2,which.max)
    B.cluster[colSums(B.prob)==0] <- NA
    X.prob <- X / (rowSums(X) + .eps)
    X.cluster <- apply(X.prob,1,which.max)

    # [Fix 4] Skip dist.cor if there are missing values
    if(save.time || any(Y.weights == 0)){
      silhouette <- NA; CPCC <- NA; dist.cor <- NA
    }else{
      silhouette <- .silhouette.simple(B.prob,B.cluster)
      dist.cor <- stats::cor(as.vector(stats::dist(t(Y))),as.vector(stats::dist(t(B))))
      if(Q>=2){
        M <- t(B.prob) %*% B.prob
        h.dist <- as.matrix(stats::cophenetic(stats::hclust(stats::as.dist(1-M))))
        up <- upper.tri(M)
        CPCC <- stats::cor(h.dist[up],(1-M)[up])
      }else{ CPCC <- NA }
    }
  }else{
    r2 <- NA; sigma <- NA; mae <- NA
    B.prob <- NA; B.cluster <- NA
    B.prob.sd.min <- NA; B.prob.max.mean <- NA; B.prob.entropy.mean <- NA
    XB <- NA; X.prob <- NA; X.cluster <- NA;
    silhouette <- NA; CPCC <- NA; dist.cor <- NA
  }
  if(epsilon.iter > abs(epsilon)) warning(paste0("maximum iterations (",maxit,") reached..."))
  end.time <- Sys.time()
  diff.time.st <- paste0(round(difftime(end.time,start.time,units="sec"),1),"sec")
  if(print.dims) message(diff.time.st)

  n.missing <- sum(Y.weights == 0)
  n.total <- prod(dim(Y))
  A.attr <- NULL
  if (!is.null(A)) A.attr <- attributes(A)

  result <- list(
    call      = match.call(),
    dims      = dims,
    runtime   = diff.time.st,
    method    = method,
    X         = X,
    B         = B,
    XB        = XB,
    C         = C,
    B.prob    = B.prob,
    B.cluster = B.cluster,
    X.prob    = X.prob,
    X.cluster = X.cluster,
    A.attr    = A.attr,
    formula.meta = formula.meta,
    n.missing = n.missing,
    n.total   = n.total,
    rank      = Q,
    objfunc   = objfunc,
    objfunc.iter = objfunc.iter,
    r.squared = r2,
    sigma     = sigma,
    mae =mae,
    criterion = list(
      B.prob.sd.min = B.prob.sd.min,
      B.prob.max.mean = B.prob.max.mean,
      B.prob.entropy.mean= B.prob.entropy.mean,
      ICp1 = ICp1, ICp2 = ICp2, ICp3 = ICp3, ICp = ICp,
      AIC  = AIC,  BIC  = BIC,
      silhouette = silhouette,
      CPCC       = CPCC,
      dist.cor   = dist.cor
    )
  )
  class(result) <- "nmfkc"
  return(result)
}















#' @title Plot method for objects of class \code{nmfkc}
#' @description
#' \code{plot.nmfkc} produces a diagnostic plot for the return value of
#' \code{nmfkc}, showing the objective function across iterations.
#'
#' @param x An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param ... Additional arguments passed to the base \code{\link{plot}} function.
#' @return Called for its side effect (a plot). Returns \code{NULL} invisibly.
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' plot(result)
#'
#' @export
plot.nmfkc <- function(x,...){
  extra_args <- list(...)
  args <- list(x = x$objfunc.iter)
  if(is.null(extra_args$main)) args$main <- paste0("r.squared=",round(x$r.squared,3))
  if(is.null(extra_args$xlab)) args$xlab <- "iter"
  if(is.null(extra_args$ylab)) args$ylab <- "objfunc"
  all_args <- c(args, extra_args)
  do.call("plot",all_args)
}




#' @title Summary method for objects of class \code{nmfkc}
#' @description
#' Produces a summary of an \code{nmfkc} object, including matrix dimensions,
#' runtime, fit statistics, and diagnostics.
#'
#' @param object An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{summary.nmfkc}, containing summary statistics.
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' summary(result)
#'
#' @export
summary.nmfkc <- function(object, ...) {
  ans <- list()
  ans$call <- object$call
  ans$dims <- object$dims
  ans$rank <- object$rank
  ans$runtime <- object$runtime

  # Missing values
  if(!is.null(object$n.missing) && !is.null(object$n.total)){
    ans$n.missing <- object$n.missing
    ans$prop.missing <- object$n.missing / object$n.total * 100
  } else {
    ans$n.missing <- NULL
  }

  ans$formula.meta <- object$formula.meta
  ans$method <- object$method
  ans$iter <- length(object$objfunc.iter)
  ans$objfunc <- object$objfunc
  ans$r.squared <- object$r.squared
  ans$sigma <- object$sigma
  ans$mae <- object$mae
  if(!is.null(object$criterion)){
    ans$ICp <- object$criterion$ICp
  } else {
    ans$ICp <- NULL
  }

  # --- Diagnostics (Sparsity & Clustering Quality) ---

  # 1. Basis (X)
  if (!is.null(object$X) && is.matrix(object$X)) {
    # Sparsity: Proportion of elements close to zero (< 1e-4)
    ans$X.sparsity <- mean(object$X < 1e-4)
  }

  # 2. Probabilities (B.prob)
  if (!is.null(object$B.prob)){
    # Sparsity
    ans$B.prob.sparsity <- mean(object$B.prob < 1e-4)
    ans$B.prob.entropy.mean <- object$criterion$B.prob.entropy.mean
    ans$B.prob.max.mean <- object$criterion$B.prob.max.mean
  }

  class(ans) <- "summary.nmfkc"
  return(ans)
}

#' @title Print method for \code{summary.nmfkc} objects
#' @description
#' Prints a formatted summary of an \code{nmfkc} model fit.
#'
#' @param x An object of class \code{summary.nmfkc}.
#' @param digits Minimum number of significant digits to be used.
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @examples
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' print(summary(result))
#'
#' @export
print.summary.nmfkc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if (!is.null(x$formula.meta)) {
    cat("Formula:    ", deparse(x$formula.meta$formula), "\n")
  }
  cat("Dimensions:", x$dims, "\n")
  if(!is.null(x$rank)) cat("Rank (Q):   ", x$rank, "\n")
  cat("Runtime:    ", x$runtime, "\n")
  if (!is.null(x$method)) cat("Method:     ", x$method, "\n")
  cat("Iterations: ", x$iter, "\n")

  if (!is.null(x$n.missing)) {
    cat("Missing:    ", x$n.missing,
        sprintf("(%.1f%%)", x$prop.missing), "\n")
  }

  cat("\nStatistics:\n")
  cat("  Objective function:  ", format(x$objfunc, digits = digits), "\n")
  cat("  Multiple R-squared:  ", format(x$r.squared, digits = digits), "\n")
  cat("  Residual Std Error:  ", format(x$sigma, digits = digits), "\n")
  cat("  Mean Absolute Error: ", format(x$mae, digits = digits), "\n")

  if (!is.null(x$ICp)) {
    cat("  ICp:                 ", format(x$ICp, digits = digits), "\n")
  }

  cat("\nStructure Diagnostics:\n")
  if (!is.null(x$X.sparsity)) {
    cat("  Basis (X) Sparsity:   ", sprintf("%.1f%%", x$X.sparsity * 100), "(< 1e-4)\n")
  }
  if (!is.null(x$B.prob.sparsity)) {
    cat("  Coef (B) Sparsity:    ", sprintf("%.1f%%", x$B.prob.sparsity * 100), "(< 1e-4)\n")
  }
  if(!is.null(x$B.prob.entropy.mean)){
    cat("  Clustering Entropy:   ", format(x$B.prob.entropy.mean, digits = digits),
        "(range: 0-1, closer to 0 is better)\n")
    cat("  Clustering Crispness: ", format(x$B.prob.max.mean, digits = digits),
        "(range: 0-1, closer to 1 is better)\n")
  }
  cat("\n")
  invisible(x)
}














#' @title Normalize a matrix to the range \eqn{[0,1]}
#' @description
#' \code{nmfkc.normalize} rescales the values of a matrix to lie between 0 and 1
#' using the column-wise minimum and maximum values of a reference matrix.
#'
#' @param x A numeric matrix (or vector) to be normalized.
#' @param ref A reference matrix from which the column-wise minima and maxima are taken.
#'   Default is \code{x}.
#'
#' @return A matrix of the same dimensions as \code{x}, with each column rescaled to the \eqn{[0,1]} range.
#' @seealso \code{\link{nmfkc.denormalize}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' x <- nmfkc.normalize(iris[,-5])
#' apply(x,2,range)
nmfkc.normalize <- function(x,ref=x){
  if(is.vector(x)){
    x <- matrix(x,ncol=1)
    ref <- matrix(ref,ncol=1)
  }
  r <- apply(ref,2,range)
  denom <- r[2, ] - r[1, ]
  denom[denom == 0] <- 1   # leave zero-width columns unchanged
  y <- sweep(x, 2, r[1, ], FUN = "-")
  y <- sweep(y, 2, denom,   FUN = "/")
  return(y)
}




#' @title Denormalize a matrix from \eqn{[0,1]} back to its original scale
#' @description
#' \code{nmfkc.denormalize} rescales a matrix with values in \eqn{[0,1]} back to its
#' original scale using the column-wise minima and maxima of a reference matrix.
#'
#' @param x A numeric matrix (or vector) with values in \eqn{[0,1]} to be denormalized.
#' @param ref A reference matrix used to obtain the original column-wise minima
#'   and maxima. Must have the same number of columns as \code{x}.
#'
#' @return A numeric matrix with values transformed back to the original scale.
#' @seealso \code{\link{nmfkc.normalize}}
#' @export
#' @examples
#' x <- nmfkc.normalize(iris[, -5])
#' x_recovered <- nmfkc.denormalize(x, iris[, -5])
#' apply(x_recovered - iris[, -5], 2, max)
nmfkc.denormalize <- function(x, ref=x) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    ref <- matrix(ref, ncol = 1)
  }
  r <- apply(ref, 2, range)
  y <- sweep(x, 2, r[2, ] - r[1, ], FUN = "*")
  y <- sweep(y, 2, r[1, ], FUN = "+")
  return(y)
}




#' @title Create a class (one-hot) matrix from a categorical vector
#' @description
#' \code{nmfkc.class} converts a categorical or factor vector into a class matrix
#' (one-hot encoded representation), where each row corresponds to a category
#' and each column corresponds to an observation.
#'
#' @param x A categorical vector or a factor.
#'
#' @return A binary matrix with one row per unique category and one column per observation. Each column has exactly one entry equal to 1, indicating the category of the observation.
#' @seealso \code{\link{nmfkc}}
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- nmfkc.class(iris$Species)
#' Y[,1:6]
nmfkc.class <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  lev <- levels(x)
  X <- outer(lev, x, "==")
  mode(X) <- "numeric"
  rownames(X) <- lev
  if(!is.null(names(x))) colnames(X) <- names(x) else colnames(X) <- seq_along(x)
  X
}





#' @title Prediction method for objects of class \code{nmfkc}
#' @description
#' \code{predict.nmfkc} generates predictions from an object of class \code{nmfkc},
#' either using the fitted covariates or a new covariate matrix.
#'
#' When the model was fitted using a formula (Formula Mode), a \code{newdata}
#' data frame can be supplied instead of \code{newA}; the covariate matrix is
#' then constructed automatically from the stored formula metadata.
#'
#' @param object An object of class \code{nmfkc}, i.e., the return value of \code{nmfkc}.
#' @param newA Optional. A new covariate matrix to be used for prediction.
#' @param newdata Optional data frame. Only available when the model was fitted
#'   using a formula. Covariate columns are extracted automatically using the
#'   stored formula metadata. If both \code{newdata} and \code{newA} are
#'   supplied, \code{newdata} takes precedence (with a warning).
#' @param type Type of prediction to return. Options are "response" (fitted values matrix),
#'   "prob" (soft-clustering probabilities), or "class" (hard-clustering labels based on row names of X).
#' @param ... Further arguments passed to or from other methods.
#' @return Depending on \code{type}: a numeric matrix (\code{"response"} or \code{"prob"})
#'   or a character vector of class labels (\code{"class"}).
#' @examples
#' # Prediction with newA
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' result <- nmfkc(Y, A, Q = 1)
#' newA <- rbind(1, c(10, 20, 30))
#' predict(result, newA = newA)
#'
#' @export
predict.nmfkc <- function(object, newA = NULL, newdata = NULL, type = "response", ...) {
  x <- object
  .eps <- 1e-10

  # --- newdata handling (Formula Mode) ---
  if (!is.null(newdata)) {
    if (is.null(x$formula.meta)) {
      stop("'newdata' can only be used when the model was fitted with a formula.")
    }
    if (!is.null(newA)) {
      warning("Both 'newdata' and 'newA' supplied; 'newdata' takes precedence.")
    }
    A_cols <- x$formula.meta$A_cols
    if (is.null(A_cols)) {
      stop("Formula metadata has no A columns. Model was fitted without covariates (standard NMF).")
    }
    missing_cols <- setdiff(A_cols, names(newdata))
    if (length(missing_cols) > 0) {
      stop(paste0("'newdata' is missing required columns: ", paste(missing_cols, collapse = ", ")))
    }
    newdata <- as.data.frame(newdata)
    newA <- t(as.matrix(newdata[, A_cols, drop = FALSE]))
  }

  if(is.null(newA)){
    if(type=="response"){
      result <- x$X %*% x$B
    }else{
      XB.prob <- x$X %*% x$B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }else{
    B <- x$C %*% newA
    if(type=="response"){
      result <- x$X %*% B
    }else{
      B.prob <- sweep(B, 2, colSums(B) + .eps, "/")
      XB.prob <- x$X %*% B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }
  return(result)
}








#' @title Perform k-fold cross-validation for NMF with kernel covariates
#'
#' @description
#' \code{nmfkc.cv} performs k-fold cross-validation for the tri-factorization model
#' \eqn{Y \approx X C A = X B}, where
#' \itemize{
#'   \item \eqn{Y(P,N)} is the observation matrix,
#'   \item \eqn{A(R,N)} is the covariate (or kernel) matrix,
#'   \item \eqn{X(P,Q)} is the basis matrix,
#'   \item \eqn{C(Q,R)} is the parameter matrix, and
#'   \item \eqn{B(Q,N)} is the coefficient matrix (\eqn{B = C A}).
#' }
#' Given \eqn{Y} (and optionally \eqn{A}), \eqn{X} and \eqn{C} are fitted on each
#' training split and predictive performance is evaluated on the held-out split.
#'
#' @param Y Observation matrix, or a formula (see \code{\link{nmfkc}} for Formula Mode).
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#'   Ignored when \code{Y} is a formula.
#' @param Q Rank of the basis matrix \eqn{X}.
#' @param data A data frame (required when \code{Y} is a formula with column names).
#' @param ... Additional arguments controlling CV and the internal \code{\link{nmfkc}} call:
#'   \describe{
#'     \item{\code{Y.weights}}{Optional numeric matrix or vector; 0 indicates missing/ignored values.}
#'     \item{\code{div}}{Number of folds (\eqn{k}); default: \code{5}.}
#'     \item{\code{seed}}{Integer seed for reproducible partitioning; default: \code{123}.}
#'     \item{\code{shuffle}}{Logical. If \code{TRUE} (default), randomly shuffles samples (standard CV);
#'       if \code{FALSE}, splits sequentially (block CV; recommended for time series).}
#'     \item{\emph{Arguments passed to} \code{\link{nmfkc}}}{e.g., \code{gamma} (\code{B.L1}), \code{epsilon},
#'       \code{maxit}, \code{method} (\code{"EU"} or \code{"KL"}), \code{X.restriction}, \code{X.init}, etc.}
#'   }
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{objfunc}}{Mean loss per valid entry over all folds (MSE for \code{method="EU"}).}
#'     \item{\code{sigma}}{Residual standard error (RMSE). Available only if \code{method="EU"}; on the same scale as \code{Y}.}
#'     \item{\code{objfunc.block}}{Loss for each fold.}
#'     \item{\code{block}}{Vector of fold indices (1, …, \code{div}) assigned to each column of \eqn{Y}.}
#'   }
#'
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.kernel.beta.cv}}, \code{\link{nmfkc.ar.degree.cv}}
#'
#' @examples
#' # Example 1 (with explicit covariates):
#' Y <- matrix(cars$dist, nrow = 1)
#' A <- rbind(1, cars$speed)
#' res <- nmfkc.cv(Y, A, Q = 1)
#' res$objfunc
#'
#' # Example 2 (kernel A and beta sweep):
#' Y <- matrix(cars$dist, nrow = 1)
#' U <- matrix(c(5, 10, 15, 20, 25), nrow = 1)
#' V <- matrix(cars$speed, nrow = 1)
#' betas <- 25:35/1000
#' obj <- numeric(length(betas))
#' for (i in seq_along(betas)) {
#'   A <- nmfkc.kernel(U, V, beta = betas[i])
#'   obj[i] <- nmfkc.cv(Y, A, Q = 1, div = 10)$objfunc
#' }
#' betas[which.min(obj)]
#'
#' @export

nmfkc.cv <- function(Y, A=NULL, Q=2, data, ...){
  # --- Formula Mode ---
  if (base::inherits(Y, "formula")) {
    resolved <- .nmfkc_resolve_formula(Y, A, base::missing(data), if (!base::missing(data)) data else NULL)
    Y <- resolved$Y
    A <- resolved$A
  }

  # A small constant for numerical stability to prevent division by zero and log(0).
  .eps <- 1e-10

  extra_args <- list(...)

  div <- if (!is.null(extra_args$div)) extra_args$div else 5
  seed <- if (!is.null(extra_args$seed)) extra_args$seed else 123
  shuffle <- if (!is.null(extra_args$shuffle)) extra_args$shuffle else TRUE

  gamma <- if (!is.null(extra_args$B.L1)) extra_args$B.L1 else if (!is.null(extra_args$gamma)) extra_args$gamma else 0
  epsilon <- if (!is.null(extra_args$epsilon)) extra_args$epsilon else 1e-4
  maxit   <- if (!is.null(extra_args$maxit))   extra_args$maxit   else 5000
  method  <- if (!is.null(extra_args$method))  extra_args$method  else "EU"
  Y.weights <- if (!is.null(extra_args$Y.weights)) extra_args$Y.weights else NULL

  if(!is.matrix(Y)) Y <- as.matrix(Y)
  P <- nrow(Y)
  N <- ncol(Y)

  # === Weights Preparation (Same as nmfkc) ===
  # Expand vector to matrix for easier slicing
  if (!is.null(Y.weights) && is.vector(Y.weights)) {
    if (length(Y.weights) == N) {
      Y.weights <- matrix(Y.weights, nrow = P, ncol = N, byrow = TRUE)
    } else if (length(Y.weights) == 1) {
      Y.weights <- matrix(Y.weights, nrow = P, ncol = N)
    } else {
      stop("Length of Y.weights vector must match ncol(Y) (or be 1).")
    }
  }
  # Handle NAs
  if (is.null(Y.weights)) {
    if (any(is.na(Y))) {
      Y.weights <- matrix(1, nrow=P, ncol=N)
      Y.weights[is.na(Y)] <- 0
      Y[is.na(Y)] <- 0
    } else {
      Y.weights <- matrix(1, nrow=P, ncol=N) # Full 1s for easier slicing logic
    }
  } else {
    if (!is.matrix(Y.weights)) Y.weights <- as.matrix(Y.weights)
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y)] <- 0
    Y <- Y * Y.weights
  }

  # --- Helper: Weighted Optimization of B given fixed X and Y_test ---
  optimize.B.from.Y <- function(result, Y_test, W_test, gamma, epsilon, maxit, method){
    X <- result$X
    # Initialize C (which acts as B here)
    C <- matrix(1, nrow=ncol(X), ncol=ncol(Y_test))

    # Precompute ones for penalty (must match test dimension)
    ones_QN <- matrix(1, nrow=ncol(X), ncol=ncol(Y_test))

    oldSum <- 0
    epsilon.iter <- Inf

    for(l in 1:maxit){
      B <- C
      XB <- X %*% B

      if(method=="EU"){
        # Weighted Update for B (EU)
        # Num: X^T (W * Y)
        num <- t(X) %*% (W_test * Y_test)
        # Den: X^T (W * XB)
        den <- t(X) %*% (W_test * XB)

        if(gamma != 0) den <- den + (gamma/2) * ones_QN
        C <- C * ( num / (den + .eps) )

      }else{ # KL
        # Weighted Update for B (KL)
        # Num: X^T (W * (Y/XB))
        ratio <- W_test * (Y_test / (XB + .eps))
        num <- t(X) %*% ratio
        # Den: X^T W
        den <- t(X) %*% W_test

        if(gamma != 0) den <- den + gamma * ones_QN
        C <- C * ( num / (den + .eps) )
      }

      newSum <- sum(C)
      if(l>=10){
        #epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
        epsilon.iter <- abs(newSum-oldSum) / pmax(abs(newSum), 1)
        if(epsilon.iter <= abs(epsilon)) break
      }
      oldSum <- sum(C)
    }

    B <- C
    XB <- X %*% B
    # Note: colnames handling omitted for speed in CV
    return(list(B=B, XB=XB))
  }

  is.identity.matrix <- function(A, tol = .Machine$double.eps) {
    if (nrow(A) != ncol(A)) return(FALSE)
    isTRUE(all.equal(A, diag(nrow(A)), tolerance = tol))
  }

  if(is.null(A)){
    is_identity <- TRUE
    is_symmetric.matrix <- FALSE
  }else{
    is_identity <- is.identity.matrix(A)
    is_symmetric.matrix <- isSymmetric(A, tol=.Machine$double.eps)
    A.function <- attr(A, "function.name")
    is_kernel_matrix <- !is.null(A.function) && A.function == "nmfkc.kernel"
  }

  # Create Folds
  remainder <- N %% div
  division <- N %/% div
  block <- 0*(1:N)

  if(shuffle){
    set.seed(seed)
    perm_index <- sample(1:N, N, replace=FALSE)
  } else {
    perm_index <- 1:N
  }

  processed_count <- 0
  for(i in 1:(div-1)){
    plus <- ifelse(i <= remainder, 1, 0)
    chunk_size <- division + plus
    target_indices <- perm_index[(processed_count + 1):(processed_count + chunk_size)]
    block[target_indices] <- i
    processed_count <- processed_count + chunk_size
  }
  target_indices <- perm_index[(processed_count + 1):N]
  block[target_indices] <- div

  objfunc.block <- numeric(div)
  total_valid_obs <- 0 # Denominator for weighted mean error

  for(j in 1:div){
    # Slice Y and Weights
    # Train
    Y_train <- Y[,block!=j, drop=FALSE]
    W_train <- Y.weights[,block!=j, drop=FALSE]

    # Test
    Y_test <- Y[,block==j, drop=FALSE]
    W_test <- Y.weights[,block==j, drop=FALSE]

    # Slice A
    if(is_identity){
      A_train <- NULL
    }else{
      if(is_symmetric.matrix && is_kernel_matrix){
        A_train <- A[block!=j,block!=j, drop=FALSE]
        A_test  <- A[block!=j,block==j, drop=FALSE]
      }else{
        A_train <- A[,block!=j, drop=FALSE]
        A_test  <- A[,block==j, drop=FALSE]
      }
    }

    # Run NMF on Training set (passing W_train)
    nmfkc_args <- c(
      list(...),
      list(Y = Y_train, A = A_train, Q = Q, Y.weights = W_train, # Pass weights!
           seed=NULL, print.trace = FALSE, print.dims = FALSE,
           save.time = TRUE, save.memory = TRUE)
    )
    nmfkc_args$shuffle <- NULL

    # Suppress messages from inner nmfkc calls
    res_j <- suppressMessages(do.call("nmfkc", nmfkc_args))

    # Predict on Test set
    if(is_identity){
      # Standard NMF: Optimize B for test set using weights
      resj <- optimize.B.from.Y(res_j, Y_test, W_test, gamma, epsilon, maxit, method)
      XB_test <- resj$XB
    }else{
      # Covariate NMF: Predict using A_test
      # XB = X * C * A_test
      XB_test <- res_j$X %*% res_j$C %*% A_test
    }

    # Evaluate Error (Weighted)
    if(method=="EU"){
      resid <- W_test * (Y_test - XB_test)
      objfunc.block[j] <- sum(resid^2)
    }else{
      term1 <- - (W_test * Y_test) * log(XB_test + .eps)
      term2 <- W_test * XB_test
      objfunc.block[j] <- sum(term1 + term2)
    }

    total_valid_obs <- total_valid_obs + sum(W_test > 0)
  }

  # Mean error per valid observation
  # (Avoid division by zero if all weights are 0, though unlikely)
  objfunc <- sum(objfunc.block) / max(total_valid_obs, 1)

  # Calculate RMSE for EU (This corresponds to 'sigma')
  sigma <- if(method == "EU") sqrt(objfunc) else NA

  return(list(objfunc=objfunc, sigma=sigma, objfunc.block=objfunc.block, block=block))
}





#' @title Perform Element-wise Cross-Validation (Wold's CV)
#' @description
#' \code{nmfkc.ecv} performs k-fold cross-validation by randomly holding out
#' individual elements of the data matrix (element-wise), assigning them a
#' weight of 0 via `Y.weights`, and evaluating the reconstruction error on
#' those held-out elements.
#'
#' This method (also known as Wold's CV) is theoretically robust for determining
#' the optimal rank (Q) in NMF. This function supports vector input for `Q`,
#' allowing simultaneous evaluation of multiple ranks on the same folds.
#'
#' @param Y Observation matrix, or a formula (see \code{\link{nmfkc}} for Formula Mode).
#' @param A Covariate matrix. Ignored when \code{Y} is a formula.
#' @param Q Vector of ranks to evaluate (e.g., 1:5).
#' @param div Number of folds (default: 5).
#' @param seed Integer seed for reproducibility.
#' @param data A data frame (required when \code{Y} is a formula with column names).
#' @param ... Additional arguments passed to \code{\link{nmfkc}} (e.g., \code{method="EU"}).
#'
#' @return A list with components:
#' \item{objfunc}{Numeric vector containing the Mean Squared Error (MSE) for each Q.}
#' \item{sigma}{Numeric vector containing the Residual Standard Error (RMSE) for each Q. Only available if method="EU".}
#' \item{objfunc.fold}{List of length equal to Q vector. Each element contains the MSE values for the k folds.}
#' \item{folds}{A list of length \code{div}, containing the linear indices of held-out elements for each fold (shared across all Q).}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.cv}}
#' @examples
#' # Element-wise CV to select rank
#' Y <- t(iris[1:30, 1:4])
#' res <- nmfkc.ecv(Y, Q = 1:2, div = 3)
#' res$objfunc
#'
#' @export
nmfkc.ecv <- function(Y, A=NULL, Q=1:3, div=5, seed=123, data, ...){
  # --- Formula Mode ---
  if (base::inherits(Y, "formula")) {
    resolved <- .nmfkc_resolve_formula(Y, A, base::missing(data), if (!base::missing(data)) data else NULL)
    Y <- resolved$Y
    A <- resolved$A
  }

  if(!is.matrix(Y)) Y <- as.matrix(Y)
  P <- nrow(Y)
  N <- ncol(Y)

  # --- Argument Handling ---
  extra_args <- list(...)

  # If user mistakenly passed 'rank' in ..., treat it as Q
  if (!is.null(extra_args$rank)) {
    Q <- extra_args$rank
  }

  # 1. Create Folds
  if (!is.null(seed)) set.seed(seed)
  valid_indices <- which(!is.na(Y))
  n_valid <- length(valid_indices)

  perm_indices <- sample(valid_indices)
  folds <- vector("list", div)
  chunk_size <- n_valid %/% div
  remainder <- n_valid %% div

  start_idx <- 1
  for(k in 1:div){
    current_size <- chunk_size + ifelse(k <= remainder, 1, 0)
    end_idx <- start_idx + current_size - 1
    folds[[k]] <- perm_indices[start_idx:end_idx]
    start_idx <- end_idx + 1
  }

  method <- if(!is.null(extra_args$method)) extra_args$method else "EU"

  # 2. Loop over Q vectors
  num_Q <- length(Q)
  result_objfunc <- numeric(num_Q)
  result_sigma   <- numeric(num_Q)
  result_fold    <- vector("list", num_Q)

  names(result_objfunc) <- paste0("Q=", Q)
  names(result_sigma)   <- paste0("Q=", Q)
  names(result_fold)    <- paste0("Q=", Q)

  message(paste0("Performing Element-wise CV for Q = ", paste(Q, collapse=","), " (", div, "-fold)..."))

  for(i in 1:num_Q){
    q_curr <- Q[i]
    objfunc.fold <- numeric(div)

    for(k in 1:div){
      test_idx <- folds[[k]]
      weights_train <- matrix(1, nrow=P, ncol=N)
      if(any(is.na(Y))) weights_train[is.na(Y)] <- 0
      weights_train[test_idx] <- 0

      # Prepare arguments for nmfkc
      # Remove 'Q' and 'rank' from extra_args to avoid conflicts
      nmfkc_clean_args <- extra_args
      nmfkc_clean_args$Q <- NULL
      nmfkc_clean_args$rank <- NULL

      nmfkc_clean_args$print.trace <- NULL
      nmfkc_clean_args$print.dims <- NULL
      nmfkc_clean_args$save.time <- NULL
      nmfkc_clean_args$save.memory <- NULL
      nmfkc_args <- c(list(Y=Y, A=A, Q=q_curr, Y.weights=weights_train,
                           print.trace=FALSE, print.dims=FALSE,
                           save.time=TRUE), nmfkc_clean_args)

      fit <- suppressMessages(do.call("nmfkc", nmfkc_args))

      pred <- fit$XB
      if(method == "KL"){
        .eps <- 1e-10
        term1 <- -Y[test_idx] * log(pred[test_idx] + .eps)
        term2 <- pred[test_idx]
        objfunc.fold[k] <- mean(term1 + term2)
      }else{
        residuals <- Y[test_idx] - pred[test_idx]
        objfunc.fold[k] <- mean(residuals^2)
      }
    }

    result_fold[[i]]  <- objfunc.fold
    result_objfunc[i] <- mean(objfunc.fold)
    result_sigma[i]   <- if(method == "EU") sqrt(result_objfunc[i]) else NA
  }

  return(list(objfunc = result_objfunc,
              sigma = result_sigma,
              objfunc.fold = result_fold,
              folds = folds))
}



#' @title Rank selection diagnostics with graphical output
#' @description
#' \code{nmfkc.rank} provides diagnostic criteria for selecting the rank (\eqn{Q})
#' in NMF with kernel covariates. Several model selection measures are computed
#' (e.g., R-squared, silhouette, CPCC, ARI), and results can be visualized in a plot.
#'
#' By default (\code{save.time = FALSE}), this function also computes the
#' Element-wise Cross-Validation error (Wold's CV Sigma) using \code{\link{nmfkc.ecv}}.
#'
#' The plot explicitly marks the "BEST" rank based on two criteria:
#' \enumerate{
#'   \item **Elbow Method (Red)**: Based on the curvature of the R-squared values (always computed if Q > 2).
#'   \item **Min RMSE (Blue)**: Based on the minimum Element-wise CV Sigma (only if \code{save.time=FALSE}).
#' }
#'
#' @param Y Observation matrix, or a formula (see \code{\link{nmfkc}} for Formula Mode).
#' @param A Covariate matrix. If \code{NULL}, the identity matrix is used.
#'   Ignored when \code{Y} is a formula.
#' @param rank A vector of candidate ranks to be evaluated.
#' @param save.time Logical. If \code{TRUE}, skips heavy computations like Element-wise CV.
#'   Default is \code{FALSE} (computes everything).
#' @param plot Logical. If \code{TRUE} (default), draws a plot of the diagnostic criteria.
#' @param data A data frame (required when \code{Y} is a formula with column names).
#' @param ... Additional arguments passed to \code{\link{nmfkc}} and \code{\link{nmfkc.ecv}}.
#'   \itemize{
#'     \item \code{Q}: (Deprecated) Alias for \code{rank}.
#'   }
#'
#' @return A list containing:
#' \item{rank.best}{The estimated optimal rank. Prioritizes ECV minimum if available, otherwise R-squared Elbow.}
#' \item{criteria}{A data frame containing diagnostic metrics for each rank.}
#' @seealso \code{\link{nmfkc}}, \code{\link{nmfkc.ecv}}
#' @export
#' @references
#' Brunet, J.P., Tamayo, P., Golub, T.R., Mesirov, J.P. (2004).
#' Metagenes and molecular pattern discovery using matrix factorization.
#' \emph{Proc. Natl. Acad. Sci. USA}, 101, 4164–4169.
#' \doi{10.1073/pnas.0308531101}
#' Punera, K., & Ghosh, J. (2008).
#' Consensus-based ensembles of soft clusterings.
#' \emph{Applied Artificial Intelligence}, 22(7–8), 780–810.
#' \doi{10.1080/08839510802170546}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' # Full run (default)
#' nmfkc.rank(Y, rank=1:4)
#' # Fast run (skip ECV)
#' nmfkc.rank(Y, rank=1:4, save.time=TRUE)

nmfkc.rank <- function(Y, A=NULL, rank=1:2, save.time=FALSE, plot=TRUE, data, ...){
  # --- Formula Mode ---
  if (base::inherits(Y, "formula")) {
    resolved <- .nmfkc_resolve_formula(Y, A, base::missing(data), if (!base::missing(data)) data else NULL)
    Y <- resolved$Y
    A <- resolved$A
  }

  extra_args <- list(...)

  # Backward Compatibility for Q
  if (!is.null(extra_args$Q)) rank <- extra_args$Q
  Q <- rank
  # ---------------------------------------------
  AdjustedRandIndex <- function(x){
    choose2 <- function(n) choose(n,2)
    a <- sum(apply(x,c(1,2),choose2))
    ab <- sum(sapply(rowSums(x),choose2))
    b <- ab-a
    ac <- sum(sapply(colSums(x),choose2))
    c <- ac-a
    total <- choose2(sum(x))
    d <- total-a-b-c
    (ri <- (a+d)/total)
    e <- ab*ac/total+(total-ab)*(total-ac)/total
    (ari <- (a+d-e)/(total-e))
    return(list(RI=ri,ARI=ari))
  }

  num_q <- length(Q)
  results_df <- data.frame(
    rank = Q,
    r.squared = numeric(num_q),
    ICp = numeric(num_q),
    AIC = numeric(num_q),
    BIC = numeric(num_q),
    B.prob.sd.min = numeric(num_q),
    B.prob.entropy.mean = numeric(num_q),
    B.prob.max.mean = numeric(num_q),
    ARI = numeric(num_q),
    silhouette = numeric(num_q),
    CPCC = numeric(num_q),
    dist.cor = numeric(num_q),
    sigma.ecv = numeric(num_q)
  )

  cluster.old <- NULL

  # --- Main Loop for Standard Metrics ---
  for(q_idx in 1:num_q){
    current_Q <- Q[q_idx]

    extra_args_nmfkc <- extra_args
    extra_args_nmfkc$save.memory <- FALSE
    extra_args_nmfkc$save.time <- save.time
    extra_args_nmfkc$Q <- NULL

    nmfkc_args <- c(list(Y = Y, A = A, rank = current_Q), extra_args_nmfkc)
    result <- do.call("nmfkc", nmfkc_args)

    results_df$r.squared[q_idx] <- result$r.squared
    results_df$ICp[q_idx] <- result$criterion$ICp
    results_df$AIC[q_idx] <- result$criterion$AIC
    results_df$BIC[q_idx] <- result$criterion$BIC
    results_df$B.prob.sd.min[q_idx] <- result$criterion$B.prob.sd.min
    results_df$B.prob.max.mean[q_idx] = result$criterion$B.prob.max.mean
    results_df$B.prob.entropy.mean[q_idx] = result$criterion$B.prob.entropy.mean

    if(save.time){
      results_df$CPCC[q_idx] <- NA
      results_df$dist.cor[q_idx] <- NA
      results_df$silhouette[q_idx] <- NA
    } else {
      results_df$CPCC[q_idx] <- result$criterion$CPCC
      results_df$dist.cor[q_idx] <- result$criterion$dist.cor
      results_df$silhouette[q_idx] <- result$criterion$silhouette$silhouette.mean
    }

    if(is.null(cluster.old)){
      results_df$ARI[q_idx] <- NA
    } else {
      df <- data.frame(old=cluster.old, new=result$B.cluster)
      df <- df[stats::complete.cases(df),]
      if (nrow(df) > 0 && length(unique(df$old)) > 1 && length(unique(df$new)) > 1) {
        f <- table(df$old, df$new)
        results_df$ARI[q_idx] <- AdjustedRandIndex(f)$ARI
      } else {
        results_df$ARI[q_idx] <- NA
      }
    }
    cluster.old <- result$B.cluster
  }

  # --- Element-wise CV (Wold's CV) ---
  rank.best.ecv <- NA
  if(!save.time){
    ecv_args <- list(Y = Y, A = A, Q = Q)
    extra_args_ecv <- extra_args
    extra_args_ecv$save.time <- NULL
    extra_args_ecv$save.memory <- NULL
    extra_args_ecv$Q <- NULL

    ecv_full_args <- c(ecv_args, extra_args_ecv)

    message("Running Element-wise CV (this may take time)...")
    ecv_res <- do.call("nmfkc.ecv", ecv_full_args)
    results_df$sigma.ecv <- ecv_res$sigma

    # Determine best rank by ECV
    idx_best_ecv <- which.min(results_df$sigma.ecv)
    rank.best.ecv <- if(length(idx_best_ecv) > 0) results_df$rank[idx_best_ecv] else NA
  } else {
    results_df$sigma.ecv <- NA
  }

  # --- Determine R-squared Best Rank (Elbow) ---
  rank.best.r2 <- NA
  if(num_q > 2){
    x <- 1:num_q
    y <- results_df$r.squared
    y_range <- max(y) - min(y)
    x_range <- max(x) - min(x)
    y_norm <- if(y_range > 0) (y - min(y)) / y_range else rep(0.5, length(y))
    x_norm <- if(x_range > 0) (x - min(x)) / x_range else rep(0.5, length(x))
    x1 <- x_norm[1]; y1 <- y_norm[1]
    x2 <- x_norm[num_q]; y2 <- y_norm[num_q]

    distances <- numeric(num_q)
    denom <- sqrt((y2-y1)^2 + (x2-x1)^2)
    if (denom > 0) {
      for(i in 1:num_q){
        distances[i] <- abs((y2-y1)*x_norm[i] - (x2-x1)*y_norm[i] + x2*y1 - y2*x1) / denom
      }
    }
    idx_best_r2 <- which.max(distances)
    rank.best.r2 <- results_df$rank[idx_best_r2]
  }

  # Decide final recommended rank (Prioritize ECV)
  rank.final <- if(!is.na(rank.best.ecv)) rank.best.ecv else rank.best.r2

  # --- Plotting ---
  if(plot){
    old_par <- graphics::par(mar = c(5, 4, 4, 5) + 0.1)
    on.exit(graphics::par(old_par))

    # 1. Left Axis Plot: 0-1 Metrics
    plot(results_df$rank, results_df$r.squared, type="l", col=2, lwd=3,
         xlab="Rank (Q)", ylab="Fit / Stability (0-1)", ylim=c(0,1),
         main="Rank Selection Diagnostics")
    for(q in results_df$rank) graphics::abline(v=q,col="gray90",lwd=0.5)
    graphics::points(results_df$rank, results_df$r.squared, pch=16, col=2, cex=0.8)
    graphics::text(results_df$rank, results_df$r.squared, results_df$rank, pos=3, col=2, cex=0.8)

    # Always mark R2 Elbow BEST if found
    if(!is.na(rank.best.r2)){
      idx_r2 <- which(results_df$rank == rank.best.r2)
      graphics::points(rank.best.r2, results_df$r.squared[idx_r2], pch=16, col="red", cex=1.5)
      graphics::text(rank.best.r2, results_df$r.squared[idx_r2], "Best (Elbow)", pos=1, col="red", cex=0.8)
    }

    legend_txt <- c("r.squared")
    legend_col <- c(2)
    legend_lty <- c(1)

    graphics::lines(results_df$rank, results_df$B.prob.sd.min, col="green1", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.sd.min"); legend_col <- c(legend_col, "green1"); legend_lty <- c(legend_lty, 1)

    graphics::lines(results_df$rank, results_df$B.prob.max.mean, col="green3", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.max.mean"); legend_col <- c(legend_col, "green3"); legend_lty <- c(legend_lty, 1)

        graphics::lines(results_df$rank, results_df$B.prob.entropy.mean, col="green4", lwd=2)
    legend_txt <- c(legend_txt, "B.prob.entropy.mean"); legend_col <- c(legend_col, "green4"); legend_lty <- c(legend_lty, 1)


    if (any(!is.na(results_df$ARI))) {
      graphics::lines(results_df$rank, results_df$ARI, col=4, lwd=2)
      legend_txt <- c(legend_txt, "ARI"); legend_col <- c(legend_col, 4); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$silhouette))) {
      graphics::lines(results_df$rank, results_df$silhouette, col=7, lwd=2)
      legend_txt <- c(legend_txt, "Silhouette"); legend_col <- c(legend_col, 7); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$CPCC))) {
      graphics::lines(results_df$rank, results_df$CPCC, col=6, lwd=2)
      legend_txt <- c(legend_txt, "CPCC"); legend_col <- c(legend_col, 6); legend_lty <- c(legend_lty, 1)
    }
    if (any(!is.na(results_df$dist.cor))) {
      graphics::lines(results_df$rank, results_df$dist.cor, col=5, lwd=2)
      legend_txt <- c(legend_txt, "dist.cor"); legend_col <- c(legend_col, 5); legend_lty <- c(legend_lty, 1)
    }

    # 2. Right Axis Plot: Sigma ECV (RMSE)
    if (any(!is.na(results_df$sigma.ecv))) {
      graphics::par(new = TRUE)
      graphics::plot(results_df$rank, results_df$sigma.ecv, type="l", col="blue", lwd=3,
                     axes=FALSE, xlab="", ylab="")

      graphics::points(results_df$rank, results_df$sigma.ecv, pch=16, col="blue", cex=0.8)
      graphics::text(results_df$rank, results_df$sigma.ecv, results_df$rank, pos=3, col="blue", cex=0.8)

      # Always mark ECV Min BEST
      if(!is.na(rank.best.ecv)){
        idx_ecv <- which(results_df$rank == rank.best.ecv)
        graphics::points(results_df$rank, results_df$sigma.ecv, pch=1, col="blue")
        graphics::points(rank.best.ecv, results_df$sigma.ecv[idx_ecv], pch=16, col="blue", cex=1.5)
        graphics::text(rank.best.ecv, results_df$sigma.ecv[idx_ecv], "Best (Min)", pos=1, col="blue", cex=0.8)
      }

      graphics::axis(side=4, col="blue", col.axis="blue")
      graphics::mtext("ECV Sigma (RMSE)", side=4, line=3, col="blue")

      legend_txt <- c(legend_txt, "ECV Sigma")
      legend_col <- c(legend_col, "blue")
      legend_lty <- c(legend_lty, 1)
    }

    graphics::legend("right", legend=legend_txt, col=legend_col, lty=legend_lty, lwd=2, bg="white", cex=0.7)
  }

  return(list(rank.best = rank.final, criteria = results_df))
}







#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



#' @title Plot Diagnostics: Original, Fitted, and Residual Matrices as Heatmaps
#' @description
#' This function generates a side-by-side plot of three heatmaps: the original
#' observation matrix Y, the fitted matrix XB (from NMF), and the residual matrix E (Y - XB).
#' This visualization aids in diagnosing whether the chosen rank Q is adequate
#' by assessing if the residual matrix E appears to be random noise.
#'
#' The axis labels (X-axis: Samples, Y-axis: Features) are integrated into the main title of each plot
#' to maximize the plot area, reflecting the compact layout settings.
#'
#' @param Y The original observation matrix (P x N).
#' @param result The result object returned by the nmfkc function.
#' @param Y_XB_palette A vector of colors used for Y and XB heatmaps. Defaults to a white-orange-red gradient.
#' @param E_palette A vector of colors used for the residuals (E) heatmap. Defaults to a blue-white-red gradient.
#' @param ... Additional graphical parameters passed to the internal image calls.
#'
#' @return NULL. The function generates a plot.
#' @examples
#' Y <- t(iris[1:30, 1:4])
#' result <- nmfkc(Y, Q = 2)
#' nmfkc.residual.plot(Y, result)
#'
#' @export
nmfkc.residual.plot <- function(Y, result,
                                Y_XB_palette = grDevices::colorRampPalette(c("white", "orange", "red"))(256),
                                E_palette = grDevices::colorRampPalette(c("blue", "white", "red"))(256), ...){
  if (!inherits(result, "nmfkc")) {
    stop("The 'result' argument must be an object of class 'nmfkc'.")
  }
  XB <- result$XB
  if (is.null(XB) || (length(XB) == 1 && is.na(XB))) {
    stop("'result$XB' is not available. The model may have been fitted with save.memory=TRUE.")
  }
  E <- Y - XB
  if(nrow(Y) != nrow(E) || ncol(Y) != ncol(E)){
    stop("Dimension mismatch between Y and result$XB. Check input matrices.")
  }
  old_par <- graphics::par(mfrow = c(1, 3), mar = c(1,0.5,5,0.5) + 0.1)
  on.exit(graphics::par(old_par))
  min_YX <- min(Y, XB, na.rm = TRUE)
  max_YX <- max(Y, XB, na.rm = TRUE)
  max_abs_E <- max(abs(E), na.rm = TRUE)
  min_E <- -max_abs_E
  max_E <- max_abs_E
  graphics::image(t(Y)[, nrow(Y):1],
                  col = Y_XB_palette,
                  zlim = c(min_YX, max_YX),
                  main = "1. Original Matrix Y\n X-axis: Samples (N), Y-axis: Features (P)",
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
  graphics::image(t(XB)[, nrow(XB):1],
                  col = Y_XB_palette,
                  zlim = c(min_YX, max_YX),
                  main = paste0("2. Fitted Matrix XB (Q=",ncol(result$X),")\n X-axis: Samples (N), Y-axis: Features (P)"),
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
  graphics::image(t(E)[, nrow(E):1],
                  col = E_palette,
                  zlim = c(min_E, max_E),
                  main = "3. Residual Matrix E (Y - XB)\n X-axis: Samples (N), Y-axis: Features (P)",
                  xlab = "",
                  ylab = "",
                  axes = FALSE, ...)
  graphics::box()
}



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
.nmfkc_dot_digits_from_threshold <- function(threshold) {
  if (!is.finite(threshold) || threshold <= 0) {
    return(2L)
  }
  if (threshold >= 1) {
    return(0L)
  }
  s <- format(threshold, scientific = FALSE, trim = TRUE)
  # Examples: "0.01" -> "01", "0.005" -> "005"
  s_dec <- sub("^0\\.", "", s)
  nchar(s_dec)
}


#' Format coefficient values for DOT edge labels
#'
#' If \code{digits} is \code{NULL}, the function uses the traditional
#' magnitude-based formatting rules.
#' If \code{digits} is provided, the coefficient is formatted with the
#' specified number of decimal places (typically derived from the threshold).
#'
#' @param x Numeric scalar.
#' @param digits Integer or \code{NULL}; number of decimal places.
#' @return Character string suitable for use inside DOT \code{label="..."}.
#' @keywords internal
#' @noRd
.nmfkc_dot_format_coef <- function(x, digits = NULL) {
  if (!is.finite(x)) return("")

  # Default behavior: magnitude-dependent auto formatting
  if (is.null(digits)) {
    ax <- abs(x)
    if (ax == 0) return("0")

    if (ax >= 1) {
      sprintf("%.2f", x)
    } else if (ax >= 1e-1) {
      sprintf("%.3f", x)
    } else if (ax >= 1e-2) {
      sprintf("%.4f", x)
    } else {
      sprintf("%.3e", x)
    }

    # Threshold-based fixed-digit formatting
  } else {
    if (digits <= 0) {
      sprintf("%.0f", x)
    } else {
      fmt <- paste0("%.", digits, "f")
      sprintf(fmt, x)
    }
  }
}


#' Sanitize character vectors for DOT node / cluster IDs
#'
#' Replace characters that are not alphanumeric, underscore, or dot
#' with underscores to generate safe DOT identifiers.
#'
#' @param x Character vector.
#' @return Character vector with sanitized identifiers.
#' @keywords internal
#' @noRd
.nmfkc_dot_sanitize_id <- function(x) {
  gsub("[^[:alnum:]_.]", "_", x, perl = TRUE)
}


#' Generate a standard DOT graph header for NMF-related diagrams
#'
#' Provides a consistent DOT header used by NMF visualization functions.
#'
#' @param graph_name Character; name of the DOT graph.
#' @param rankdir Character; Graphviz rank direction (e.g., "LR", "RL", "TB", "BT").
#' @param fontname Character; default font name used in the graph.
#'
#' @return Character scalar containing DOT header lines.
#' @keywords internal
#' @noRd
.nmfkc_dot_header <- function(graph_name = "NMF_GRAPH",
                              rankdir    = "LR",
                              fontname   = "Arial") {
  paste0(
    "digraph ", graph_name, " {\n",
    "  graph [rankdir=", rankdir, " compound=true];\n",
    "  splines=true; nodesep=0.4; ranksep=0.7; fontname=\"", fontname, "\";\n"
  )
}


#' Generate a DOT subgraph cluster with nodes
#'
#' Internal helper to define a cluster (subgraph) of nodes with shared
#' visual properties.
#' If \code{cluster_style = "none"}, nodes are listed individually without
#' creating a DOT subgraph.
#'
#' @param cluster_id Character; cluster identifier suffix (e.g. "Y", "X", "F").
#' @param title Character; cluster label (ignored if \code{cluster_style = "none"}).
#' @param node_ids Character vector; internal DOT node IDs.
#' @param node_labels Character vector; display labels for nodes.
#' @param shape Character; Graphviz node shape.
#' @param fill Logical; whether to use filled node shapes.
#' @param fillcolor Character; fill color (if \code{fill = TRUE}).
#' @param line_width Numeric; node border width.
#' @param indent Character; indentation prefix for formatting.
#' @param cluster_style Character; boundary style ("rounded", "dashed", "invis", "none").
#' @param cluster_color Character; cluster boundary color.
#' @param cluster_penwidth Numeric; cluster boundary width.
#'
#' @return Character scalar containing DOT code for the cluster.
#' @keywords internal
#' @noRd
.nmfkc_dot_cluster_nodes <- function(cluster_id,
                                     title,
                                     node_ids,
                                     node_labels,
                                     shape      = "box",
                                     fill       = TRUE,
                                     fillcolor  = "white",
                                     line_width = 1.5,
                                     indent     = "  ",
                                     cluster_style    = "rounded",
                                     cluster_color    = "black",
                                     cluster_penwidth = 1.0) {

  if (length(node_ids) == 0L) return("")
  if (length(node_ids) != length(node_labels)) {
    stop("node_ids and node_labels must have the same length.")
  }

  # Node style definition
  if (fill) {
    node_style <- sprintf(
      '%snode [shape=%s, style="filled,rounded", fillcolor="%s", color=black, penwidth=%.1f];\n',
      indent, shape, fillcolor, line_width
    )
  } else {
    node_style <- sprintf(
      '%snode [shape=%s, style="rounded", color=black, penwidth=%.1f];\n',
      indent, shape, line_width
    )
  }

  # Case 1: no cluster → output nodes directly
  if (identical(cluster_style, "none")) {
    scr <- node_style
    for (i in seq_along(node_ids)) {
      scr <- paste0(
        scr,
        sprintf('%s%s [label="%s"];\n', indent, node_ids[i], node_labels[i])
      )
    }
    return(scr)
  }

  # Case 2: cluster with boundary
  scr <- paste0(
    indent, "subgraph cluster_", cluster_id,
    '{label="', title, '"',
    ' style="', cluster_style, '"',
    ' color="', cluster_color, '"',
    ' penwidth=', sprintf("%.1f", cluster_penwidth), ";\n",
    node_style
  )

  for (i in seq_along(node_ids)) {
    scr <- paste0(
      scr,
      sprintf('%s  %s [label="%s"];\n', indent, node_ids[i], node_labels[i])
    )
  }
  paste0(scr, indent, "}\n")
}


#' Scale DOT edge width based on coefficient magnitude
#'
#' Computes an appropriate \code{penwidth} value for Graphviz edges
#' by scaling the coefficient relative to the maximum coefficient in
#' its category.
#'
#' @param value Numeric; coefficient value.
#' @param max_value Numeric; maximum coefficient within its group.
#' @param weight_scale Numeric; global scale factor for edge width.
#' @param min_pw Numeric; minimum penwidth.
#'
#' @return Numeric penwidth value.
#' @keywords internal
#' @noRd
.nmfkc_dot_penwidth <- function(value,
                                max_value,
                                weight_scale = 5,
                                min_pw       = 0.5) {
  if (!is.finite(max_value) || max_value <= 0) return(min_pw)
  if (!is.finite(value)     || value <= 0)     return(min_pw)
  max(min_pw, value * weight_scale / max_value)
}



