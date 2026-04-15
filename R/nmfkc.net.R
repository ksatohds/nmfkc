# =====================================================
# nmfkc.net.R — Network analysis for symmetric NMF
# =====================================================


############################################################
## nmfkc.net.inference
############################################################

#' Statistical Inference for Symmetric NMF Parameters
#'
#' @description
#' Performs statistical inference on the parameter matrix \eqn{C} of a
#' symmetric NMF model (\eqn{Y \approx XCX^\top}).
#'
#' This is a wrapper around \code{\link{nmfkc.inference}} that automatically
#' sets the covariate matrix to \eqn{A = X^\top}, which is the defining
#' property of symmetric NMF.
#'
#' @param object A fitted model returned by \code{nmfkc()} with
#'   \code{Y.symmetric = "bi"} or \code{"tri"}.
#' @param Y The original symmetric matrix (P x P).
#' @param wild.bootstrap Logical; if TRUE (default), perform wild bootstrap
#'   inference for the C matrix.
#' @param ... Additional arguments passed to \code{nmfkc.inference()}
#'   (e.g., \code{wild.B}, \code{wild.seed}, \code{C.p.side}).
#'
#' @return The object augmented with inference fields
#'   (same as \code{nmfkc.inference}), with class
#'   \code{c("nmfkc.net.inference", "nmfkc.inference", "nmf.inference", "nmfkc", "nmf")}.
#'
#' @seealso \code{\link{nmfkc.inference}}, \code{\link{nmfkc.net.DOT}},
#'   \code{\link{summary.nmfkc.net.inference}}
#'
#' @examples
#' \donttest{
#' library(nmfkc)
#' Y <- matrix(c(0,1,1,0,0,0,
#'               1,0,1,0,0,0,
#'               1,1,0,1,0,0,
#'               0,0,1,0,1,1,
#'               0,0,0,1,0,1,
#'               0,0,0,1,1,0), 6, 6)
#' res <- nmfkc(Y, rank = 2, Y.symmetric = "tri", nstart = 20)
#' res_inf <- nmfkc.net.inference(res, Y)
#' summary(res_inf)
#' }
#'
#' @export
nmfkc.net.inference <- function(object, Y, wild.bootstrap = TRUE, ...) {
  if (!inherits(object, "nmfkc"))
    stop("object must be of class 'nmfkc'.")
  if (is.null(object$X))
    stop("object must contain X matrix.")

  ## Key: A = t(X) for symmetric NMF (Y = X C X')
  A <- t(object$X)
  result <- nmfkc.inference(object, Y, A = A, wild.bootstrap = wild.bootstrap, ...)

  ## Relabel coefficients for symmetric context: Covariate -> Basis.col
  if (!is.null(result$coefficients)) {
    cf <- result$coefficients
    names(cf)[names(cf) == "Basis"]     <- "Basis.row"
    names(cf)[names(cf) == "Covariate"] <- "Basis.col"

    ## Keep only strict upper triangle (off-diagonal) since C is symmetric
    ## Diagonal (within-group) is always large and trivially significant;
    ## only between-group interactions are meaningful to test.
    Q <- ncol(object$X)
    rlabs <- if (!is.null(rownames(object$C))) rownames(object$C)
             else paste0("Basis", 1:Q)
    clabs <- if (!is.null(colnames(object$C))) colnames(object$C)
             else paste0("Basis", 1:Q)
    row_idx <- match(cf$Basis.row, rlabs)
    col_idx <- match(cf$Basis.col, clabs)
    keep <- !is.na(row_idx) & !is.na(col_idx) & (row_idx < col_idx)
    cf <- cf[keep, , drop = FALSE]
    ## Unify Basis.col labels to match Basis.row naming
    cf$Basis.col <- rlabs[col_idx[keep]]
    rownames(cf) <- NULL
    result$coefficients <- cf
  }

  class(result) <- c("nmfkc.net.inference", "nmfkc.inference",
                      "nmf.inference", "nmfkc", "nmf")
  result
}


############################################################
## summary / print for nmfkc.net.inference
############################################################

#' @title Summary method for nmfkc.net.inference objects
#' @description
#' Produces a summary of a symmetric NMF model with inference results,
#' including the coefficients table for \eqn{C}.
#'
#' @param object An object of class \code{"nmfkc.net.inference"}.
#' @param ... Additional arguments passed to \code{summary.nmfkc}.
#' @return An object of class \code{"summary.nmfkc.net.inference"}.
#' @seealso \code{\link{nmfkc.net.inference}}
#' @export
summary.nmfkc.net.inference <- function(object, ...) {
  ans <- summary.nmfkc(object, ...)
  ans$sigma2.used  <- object$sigma2.used
  ans$coefficients <- object$coefficients
  ans$C.p.side     <- object$C.p.side
  ans$C            <- object$C
  ans$Y.symmetric  <- tryCatch(eval(object$call$Y.symmetric),
                               error = function(e) NULL)
  class(ans) <- "summary.nmfkc.net.inference"
  ans
}


#' @title Print method for summary.nmfkc.net.inference objects
#' @description
#' Prints a formatted summary of a symmetric NMF model with inference results.
#' The coefficients table uses \code{Basis.row} and \code{Basis.col} labels
#' reflecting the symmetric structure \eqn{C_{qr}}.
#'
#' @param x An object of class \code{"summary.nmfkc.net.inference"}.
#' @param digits Minimum number of significant digits.
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfkc.net.inference}}
#' @export
print.summary.nmfkc.net.inference <- function(x,
    digits = max(3L, getOption("digits") - 3L), ...) {

  # Print base summary
  print.summary.nmfkc(x, digits = digits, ...)

  # Model type
  sym_type <- if (!is.null(x$Y.symmetric)) x$Y.symmetric else "unknown"
  cat(sprintf("Symmetric NMF type: %s\n", sym_type))

  # Inference section
  cat("Inference (conditional on X):\n")
  cat("  sigma^2:  ", format(x$sigma2.used, digits = digits), "\n")

  # C matrix
  if (!is.null(x$C)) {
    Q <- nrow(x$C)
    cat(sprintf("  C matrix:  %d x %d", Q, Q))
    if (sym_type == "tri") {
      sym_check <- max(abs(x$C - t(x$C)))
      cat(sprintf(" (max asymmetry: %.1e)", sym_check))
    }
    cat("\n")
  }

  # Coefficients table
  if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
    cf <- x$coefficients

    p_side <- if (!is.null(x$C.p.side)) x$C.p.side else "one.sided"
    p_header <- if (p_side == "one.sided") "Pr(>z)" else "Pr(>|z|)"

    sig_stars <- function(p) {
      ifelse(!is.finite(p), " ",
        ifelse(p < 0.001, "***",
          ifelse(p < 0.01, "**",
            ifelse(p < 0.05, "*",
              ifelse(p < 0.1, ".", " ")))))
    }
    format_pval <- function(p) {
      ifelse(!is.finite(p), "      NA",
        ifelse(p < 2.2e-16, "  <2e-16",
          formatC(p, format = "g", digits = 4, width = 8)))
    }

    n_total <- nrow(cf)
    n_sig <- sum(cf$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("\nC coefficients: %d total, %d significant\n",
                n_total, n_sig))

    rnames <- paste0(cf$Basis.row, ":", cf$Basis.col)

    est <- formatC(cf$Estimate, format = "f", digits = 3, width = 9)
    se  <- formatC(cf$SE, format = "f", digits = 3, width = 10)
    bse <- formatC(cf$BSE, format = "f", digits = 3, width = 6)
    zv  <- formatC(cf$z_value, format = "f", digits = 2, width = 7)
    pv_str <- format_pval(cf$p_value)
    stars <- sig_stars(cf$p_value)

    max_lw <- max(nchar(rnames))
    hdr <- sprintf("%s %s %s %s %s %s",
                   formatC("Estimate", width = 9),
                   formatC("Std. Error", width = 10),
                   formatC("(Boot)", width = 6),
                   formatC("z value", width = 7),
                   formatC(p_header, width = 8), "")
    cat(sprintf("%s %s\n", formatC("Row:Col", width = max_lw), hdr))
    for (i in seq_len(n_total)) {
      cat(sprintf("%s %s %s %s %s %s %s\n",
                  formatC(rnames[i], width = max_lw),
                  est[i], se[i], bse[i], zv[i], pv_str[i], stars[i]))
    }
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }

  cat("\n")
  invisible(x)
}


############################################################
## nmfkc.net.DOT
############################################################

#' Generate a Graphviz DOT Diagram for a Symmetric NMF Network
#'
#' @description
#' Creates a Graphviz DOT script that visualizes the two-layer structure
#' of a symmetric NMF model (\eqn{Y \approx X C X^\top}).
#'
#' The resulting diagram displays:
#' \itemize{
#'   \item outer nodes: the original network nodes (rows/columns of Y),
#'   \item inner nodes: the latent basis/group nodes (columns of X),
#'   \item directed edges from basis to node: membership weights from X,
#'   \item undirected edges between basis nodes: inter-group interactions
#'         from C matrix (for tri-symmetric models only).
#' }
#'
#' @param result A list returned by \code{nmfkc()} with
#'   \code{Y.symmetric = "bi"} or \code{"tri"}.
#'   If inference results are present (from \code{\link{nmfkc.net.inference}}),
#'   C edges are decorated with significance stars.
#' @param threshold Minimum coefficient value to display an edge.
#' @param sig.level Significance level for filtering C edges (if inference
#'   results are present). Set to \code{NULL} to show all edges above threshold.
#' @param weight_scale Base scaling factor for edge widths.
#' @param weight_scale_xy Scaling factor for X edges (basis -> node).
#' @param weight_scale_xx Scaling factor for C edges (basis <-> basis).
#' @param rankdir Graphviz rank direction (\code{"TB"} or \code{"LR"}).
#' @param fill Logical; whether nodes are filled with color.
#' @param hide.isolated Logical; if TRUE, omit outer nodes with no X edge
#'   above threshold.
#' @param Y.label Character vector of labels for outer nodes.
#' @param X.label Character vector of labels for basis nodes.
#' @param Y.title Cluster title for outer nodes.
#' @param X.title Cluster title for basis nodes.
#' @param show.theta Logical or NULL. Whether to draw C edges
#'   between basis nodes. NULL = auto-detect (TRUE for tri, FALSE for bi).
#' @param cluster.box Style of cluster box: \code{"none"}, \code{"normal"},
#'   \code{"faint"}, \code{"invisible"}.
#' @param layout Graphviz layout engine: \code{"fdp"}, \code{"dot"},
#'   \code{"neato"}, \code{"circo"}, \code{"twopi"}.
#' @param X.color Color palette for basis nodes (length Q).
#' @param Y.cluster Coloring mode for outer nodes: \code{"soft"} (weighted mix)
#'   or \code{"hard"} (most probable basis color).
#'
#' @return A character string of class \code{c("nmfkc.net.DOT", "nmfkc.DOT")}
#'   representing a valid Graphviz DOT script.
#'   Use \code{plot()} to render (requires DiagrammeR).
#'
#' @seealso \code{\link{nmfkc.net.inference}}, \code{\link{nmfkc.DOT}},
#'   \code{\link{plot.nmfkc.DOT}}
#'
#' @examples
#' \donttest{
#' library(nmfkc)
#' Y <- matrix(c(0,1,1,0,0,0,
#'               1,0,1,0,0,0,
#'               1,1,0,1,0,0,
#'               0,0,1,0,1,1,
#'               0,0,0,1,0,1,
#'               0,0,0,1,1,0), 6, 6)
#' res <- nmfkc(Y, rank = 2, Y.symmetric = "tri", nstart = 20)
#' dot <- nmfkc.net.DOT(res)
#' plot(dot)
#' }
#'
#' @export
nmfkc.net.DOT <- function(
    result,
    threshold       = 0.01,
    sig.level       = 0.1,
    weight_scale    = 5,
    weight_scale_xy = 1,
    weight_scale_xx = weight_scale,
    rankdir         = "TB",
    fill            = TRUE,
    hide.isolated   = TRUE,
    Y.label         = NULL,
    X.label         = NULL,
    Y.title         = "Nodes (Y)",
    X.title         = "Basis (X)",
    show.theta      = NULL,
    cluster.box     = c("none", "normal", "faint", "invisible"),
    layout          = c("fdp", "dot", "neato", "circo", "twopi"),
    X.color         = NULL,
    Y.cluster       = c("soft", "hard")
) {

  ## ---------------------------------------------------------
  ## Extract matrices
  ## ---------------------------------------------------------
  X <- result$X          # P x Q
  C_mat <- result$C      # Q x Q

  if (is.null(X)) stop("result must contain X (basis matrix).")
  if (is.null(C_mat)) stop("result must contain C (parameter matrix).")

  P <- nrow(X)
  Q <- ncol(X)

  ## ---------------------------------------------------------
  ## Auto-detect bi vs tri model
  ## ---------------------------------------------------------
  if (is.null(show.theta)) {
    sym_arg <- tryCatch(eval(result$call$Y.symmetric), error = function(e) NULL)
    if (!is.null(sym_arg)) {
      show.theta <- (sym_arg == "tri")
    } else {
      is_identity <- (all(dim(C_mat) == c(Q, Q)) &&
                      isTRUE(all.equal(C_mat, diag(Q), tolerance = 1e-6)))
      show.theta <- !is_identity
    }
  }

  ## Symmetrize C (numerical optimization may leave tiny asymmetry)
  if (show.theta) {
    C_mat <- (C_mat + t(C_mat)) / 2
  }

  ## ---------------------------------------------------------
  ## Argument matching
  ## ---------------------------------------------------------
  layout      <- match.arg(layout)
  Y.cluster   <- match.arg(Y.cluster)
  cluster.box <- match.arg(cluster.box)

  cluster_style <- switch(cluster.box,
                          normal    = "rounded",
                          faint     = "rounded,dashed",
                          invisible = "rounded",
                          none      = "none")
  cluster_color <- switch(cluster.box,
                          normal    = "black",
                          faint     = "gray80",
                          invisible = "none",
                          none      = "none")

  ## ---------------------------------------------------------
  ## Labels
  ## ---------------------------------------------------------
  Y_labels <- if (!is.null(Y.label)) Y.label else rownames(X)
  X_labels <- if (!is.null(X.label)) X.label else colnames(X)
  if (is.null(Y_labels)) Y_labels <- as.character(seq_len(P))
  if (is.null(X_labels)) X_labels <- paste0("Basis", seq_len(Q))

  Y_ids <- .nmfkc_dot_sanitize_id(paste0("Y_", seq_len(P)))
  X_ids <- .nmfkc_dot_sanitize_id(paste0("X_", seq_len(Q)))

  ## ---------------------------------------------------------
  ## Filter isolated outer nodes
  ## ---------------------------------------------------------
  idx_N <- seq_len(P)
  if (isTRUE(hide.isolated)) {
    used <- apply(X, 1L, function(row) any(row >= threshold, na.rm = TRUE))
    idx_N <- which(used)
  }

  ## ---------------------------------------------------------
  ## Membership probability: B.prob (Q x P), columns sum to 1
  ## ---------------------------------------------------------
  B_prob <- result$B.prob
  if (is.null(B_prob)) {
    B_raw <- result$B
    if (!is.null(B_raw)) {
      cs <- colSums(B_raw)
      cs[cs == 0] <- 1
      B_prob <- sweep(B_raw, 2, cs, "/")
    } else {
      stop("result must contain B.prob or B.")
    }
  }
  H <- t(B_prob)   # P x Q, each row sums to 1

  ## ---------------------------------------------------------
  ## Colors: basis nodes and outer nodes
  ## ---------------------------------------------------------
  if (is.null(X.color)) X.color <- grDevices::hcl.colors(Q, "Set 2")
  if (length(X.color) < Q) X.color <- rep_len(X.color, Q)

  basis_rgb <- t(grDevices::col2rgb(X.color)) / 255  # Q x 3

  node_hex <- character(P)
  if (Y.cluster == "soft") {
    for (i in seq_len(P)) {
      w <- H[i, ]
      mixed <- as.vector(w %*% basis_rgb)
      mixed <- pmin(pmax(mixed, 0), 1)
      node_hex[i] <- grDevices::rgb(mixed[1], mixed[2], mixed[3])
    }
  } else {
    for (i in seq_len(P)) {
      node_hex[i] <- X.color[which.max(H[i, ])]
    }
  }

  ## ---------------------------------------------------------
  ## DOT header
  ## ---------------------------------------------------------
  scr <- .nmfkc_dot_header(graph_name = "NMF_NET", rankdir = rankdir)
  scr <- paste0(scr, '  layout=', layout, ';\n',
                '  start=123;\n',
                '  node [fontsize=8];\n')

  ## ---------------------------------------------------------
  ## Basis node cluster
  ## ---------------------------------------------------------
  scr <- paste0(scr, '\n  // Basis nodes (latent groups)\n')
  if (!identical(cluster_style, "none")) {
    scr <- paste0(scr,
      sprintf('  subgraph cluster_Basis{label="%s" style="%s" color="%s" penwidth=1.0;\n',
              X.title, cluster_style, cluster_color))
  }
  for (q in seq_len(Q)) {
    scr <- paste0(scr,
      sprintf('  %s [label="%s", shape=ellipse, style="filled,rounded", fillcolor="%s"];\n',
              X_ids[q], X_labels[q], X.color[q]))
  }
  if (!identical(cluster_style, "none")) {
    scr <- paste0(scr, '  }\n')
  }

  ## ---------------------------------------------------------
  ## Outer node cluster
  ## ---------------------------------------------------------
  scr <- paste0(scr, '\n  // Original network nodes\n')
  if (!identical(cluster_style, "none")) {
    scr <- paste0(scr,
      sprintf('  subgraph cluster_Node{label="%s" style="%s" color="%s" penwidth=1.0;\n',
              Y.title, cluster_style, cluster_color))
  }
  for (i in idx_N) {
    if (fill) {
      scr <- paste0(scr,
        sprintf('  %s [label="%s", shape=circle, style="filled,rounded", fillcolor="%s", color=black, penwidth=0.5];\n',
                Y_ids[i], Y_labels[i], node_hex[i]))
    } else {
      scr <- paste0(scr,
        sprintf('  %s [label="%s", shape=circle, style="rounded", color=black, penwidth=0.5];\n',
                Y_ids[i], Y_labels[i]))
    }
  }
  if (!identical(cluster_style, "none")) {
    scr <- paste0(scr, '  }\n')
  }

  ## ---------------------------------------------------------
  ## Edge defaults
  ## ---------------------------------------------------------
  scr <- paste0(scr,
    '\n  edge [fontname="Arial", fontsize=8, arrowhead=open];\n')

  pw     <- .nmfkc_dot_penwidth
  digits <- .nmfkc_dot_digits_from_threshold(threshold)
  fmtc   <- function(x) .nmfkc_dot_format_coef(x, digits)

  ## ---------------------------------------------------------
  ## X edges: Basis -> Node (membership probability)
  ## ---------------------------------------------------------
  scr <- paste0(scr,
    '\n  // Membership edges: Basis -> Node (B.prob)\n',
    '  edge [color="gray0", fontcolor="gray0", style=solid];\n')

  if (max(H, na.rm = TRUE) > 0) {
    for (i in idx_N) {
      for (q in seq_len(Q)) {
        val <- H[i, q]
        if (is.finite(val) && val >= threshold) {
          pen <- pw(val, 1, weight_scale_xy)
          lab <- sprintf("%.2f", val)
          scr <- paste0(scr,
            sprintf('  %s -> %s [label="%s", penwidth=%.2f];\n',
                    X_ids[q], Y_ids[i], lab, pen))
        }
      }
    }
  }

  ## ---------------------------------------------------------
  ## Significance stars for C edges (from nmfkc.net.inference)
  ## ---------------------------------------------------------
  C_stars <- NULL
  C_show  <- NULL
  if (!is.null(result$coefficients) && show.theta) {
    C_stars <- matrix("", nrow = Q, ncol = Q)
    C_pval  <- matrix(NA_real_, nrow = Q, ncol = Q)
    cf <- result$coefficients

    ## Column names depend on whether inference was via nmfkc.net.inference
    row_col <- if ("Basis.row" %in% names(cf)) "Basis.row" else "Basis"
    col_col <- if ("Basis.col" %in% names(cf)) "Basis.col" else "Covariate"

    rlabs <- if (!is.null(rownames(result$C))) rownames(result$C)
             else paste0("Basis", 1:Q)
    clabs <- if (!is.null(colnames(result$C))) colnames(result$C)
             else paste0("Basis", 1:Q)

    for (k in seq_len(nrow(cf))) {
      q <- match(cf[[row_col]][k], rlabs)
      r <- match(cf[[col_col]][k], clabs)
      if (!is.na(q) && !is.na(r) && !is.na(cf$p_value[k])) {
        p <- cf$p_value[k]
        C_pval[q, r] <- p
        if      (p < 0.001) C_stars[q, r] <- "***"
        else if (p < 0.01)  C_stars[q, r] <- "**"
        else if (p < 0.05)  C_stars[q, r] <- "*"
      }
    }
    if (!is.null(sig.level)) {
      C_show <- !is.na(C_pval) & C_pval < sig.level
    }
  }

  ## ---------------------------------------------------------
  ## C edges: Basis <-> Basis (inter-group interaction, tri only)
  ## ---------------------------------------------------------
  if (isTRUE(show.theta) && Q >= 2) {
    scr <- paste0(scr,
      '\n  // Inter-basis interaction edges (C, off-diagonal)\n',
      '  edge [color=black, fontcolor=black, style=dashed, dir=none];\n')

    max_C_off <- suppressWarnings(
      max(C_mat[upper.tri(C_mat)], na.rm = TRUE))
    if (is.finite(max_C_off) && max_C_off > 0) {
      for (q in seq_len(Q - 1)) {
        for (r in (q + 1):Q) {
          val <- C_mat[q, r]
          show_edge <- if (!is.null(C_show)) C_show[q, r]
                       else is.finite(val) && val >= threshold
          if (show_edge) {
            pen <- pw(val, max_C_off, weight_scale_xx)
            lab <- fmtc(val)
            if (!is.null(C_stars)) lab <- paste0(lab, C_stars[q, r])
            scr <- paste0(scr,
              sprintf('  %s -> %s [label="%s", penwidth=%.2f, dir=none];\n',
                      X_ids[q], X_ids[r], lab, pen))
          }
        }
      }
    }
  }

  ## ---------------------------------------------------------
  ## Close and return
  ## ---------------------------------------------------------
  result <- paste0(scr, "}\n")
  class(result) <- c("nmfkc.net.DOT", "nmfkc.DOT")
  result
}
