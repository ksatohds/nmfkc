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
#'   \code{Y.symmetric = "bi"} or \code{"tri"}, or the newer
#'   \code{\link{nmfkc.net}} with \code{type = "bi"} or \code{"tri"}.
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
#' res <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 20)
#' res_inf <- nmfkc.net.inference(res, Y)
#' summary(res_inf)
#' }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
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

  class(result) <- c("nmfkc.net.inference", "nmfkc.net",
                      "nmfkc.inference", "nmf.inference", "nmfkc", "nmf")
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
#' @param by Character; grouping order of the coefficients table. Because the
#'   symmetric model sets \eqn{A = X^\top}, the column factor \code{Basis.col}
#'   plays the covariate role: \code{"covariate"} (default) lists all
#'   \code{Basis.row} within each \code{Basis.col} (1-1, 2-1, ...), while
#'   \code{"basis"} lists all \code{Basis.col} within each \code{Basis.row}
#'   (1-1, 1-2, ...).
#' @param ... Additional arguments (currently unused).
#' @return Called for its side effect (printing). Returns \code{x} invisibly.
#' @seealso \code{\link{summary.nmfkc.net.inference}}
#' @export
print.summary.nmfkc.net.inference <- function(x,
    digits = max(3L, getOption("digits") - 3L),
    by = c("covariate", "basis"), ...) {
  by <- match.arg(by)

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
    cf <- cf[.coef.order.by(cf, by), , drop = FALSE]   # grouping order (by)

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
#'   \code{Y.symmetric = "bi"} or \code{"tri"}, or the newer
#'   \code{\link{nmfkc.net}} with \code{type = "bi"} or \code{"tri"}.
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
#' @param signed Logical. If \code{TRUE}, \eqn{C} is treated as a signed
#'   matrix: positive entries rendered as solid edges and negative as
#'   dashed, with edge visibility threshold on \eqn{|C|}.  Default is
#'   \code{inherits(result, "nmfkc.net.signed")} so that results from
#'   \code{\link{nmfkc.net}} are auto-detected.
#' @param cluster.box Style of cluster box: \code{"none"}, \code{"normal"},
#'   \code{"faint"}, \code{"invisible"}.
#' @param layout Graphviz layout engine, in recommended order:
#'   \code{"neato"} (default; spring model, clearest for small/medium community
#'   graphs), \code{"fdp"} (force-directed, scales to larger graphs),
#'   \code{"twopi"} (radial), \code{"circo"} (circular), \code{"dot"}
#'   (hierarchical).  For community networks, \code{"neato"} or \code{"fdp"} with
#'   a raised \code{threshold} (e.g.\ 0.2--0.3) separate the groups best.
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
#' res <- nmfkc.net(Y, rank = 2, type = "tri", nstart = 20)
#' dot <- nmfkc.net.DOT(res)
#' plot(dot)
#' }
#'
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
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
    signed          = inherits(result, "nmfkc.net.signed"),
    cluster.box     = c("none", "normal", "faint", "invisible"),
    layout          = c("neato", "fdp", "twopi", "circo", "dot"),
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
    if (!is.null(result$type)) {
      ## New nmfkc.net API: $type is "bi" / "tri" / "signed".  Only "bi"
      ## has C = I (no inter-class interaction), so it draws no C edges.
      show.theta <- (result$type != "bi")
    } else {
      sym_arg <- tryCatch(eval(result$call$Y.symmetric), error = function(e) NULL)
      if (!is.null(sym_arg)) {
        show.theta <- (sym_arg == "tri")
      } else {
        ## Legacy fallback: detect C = I.  Drop dimnames first, otherwise
        ## all.equal() reports a names mismatch and mis-detects bi as tri.
        is_identity <- (all(dim(C_mat) == c(Q, Q)) &&
                        isTRUE(all.equal(unname(C_mat), diag(Q),
                                         tolerance = 1e-6)))
        show.theta <- !is_identity
      }
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
  ## Membership probability: H (P x Q, each row sums to 1)
  ## New nmfkc.net objects carry $X.prob; legacy nmfkc(Y.symmetric=...)
  ## objects carry $B.prob (Q x P) which we transpose.
  ## ---------------------------------------------------------
  if (!is.null(result$X.prob)) {
    H <- result$X.prob
  } else if (!is.null(result$B.prob)) {
    H <- t(result$B.prob)
  } else if (!is.null(result$B)) {
    B_raw <- result$B
    cs <- colSums(B_raw); cs[cs == 0] <- 1
    H <- t(sweep(B_raw, 2, cs, "/"))
  } else {
    stop("result must contain X.prob, B.prob, or B.")
  }

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
  ## When `signed = TRUE`: positive C rendered solid, negative C dashed.
  ## When `signed = FALSE` (default for non-signed): all edges dashed.
  ## ---------------------------------------------------------
  if (isTRUE(show.theta) && Q >= 2) {
    default_style <- if (isTRUE(signed)) "solid" else "dashed"
    scr <- paste0(scr,
      '\n  // Inter-basis interaction edges (C, off-diagonal)\n',
      sprintf('  edge [color=black, fontcolor=black, style=%s, dir=none];\n',
              default_style))

    ## Visibility threshold: abs(val) >= threshold for signed; val >= threshold otherwise
    visible_val <- function(v) {
      if (isTRUE(signed)) abs(v) else v
    }
    max_C_off <- suppressWarnings(
      max(visible_val(C_mat[upper.tri(C_mat)]), na.rm = TRUE))
    if (is.finite(max_C_off) && max_C_off > 0) {
      for (q in seq_len(Q - 1)) {
        for (r in (q + 1):Q) {
          val <- C_mat[q, r]
          v_show <- visible_val(val)
          show_edge <- if (!is.null(C_show)) C_show[q, r]
                       else is.finite(v_show) && v_show >= threshold
          if (show_edge) {
            pen <- pw(v_show, max_C_off, weight_scale_xx)
            lab <- fmtc(val)
            if (!is.null(C_stars)) lab <- paste0(lab, C_stars[q, r])
            if (isTRUE(signed) && val < 0) {
              ## Negative C: dashed style, overriding the solid default
              scr <- paste0(scr,
                sprintf('  %s -> %s [label="%s", penwidth=%.2f, style=dashed, dir=none];\n',
                        X_ids[q], X_ids[r], lab, pen))
            } else {
              scr <- paste0(scr,
                sprintf('  %s -> %s [label="%s", penwidth=%.2f, dir=none];\n',
                        X_ids[q], X_ids[r], lab, pen))
            }
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


.nmfkc.net.run_once <- function(Y, X, C, maxit, epsilon,
                                X.restriction, C.L1, X.L2.ortho,
                                Wmat, has.weights) {
  N <- nrow(Y); Q <- ncol(X)
  small <- 1e-16
  norm_XC <- function(X, C) {
    if (X.restriction %in% c("fixed", "none")) return(list(X = X, C = C))
    cs <- colSums(X) + small
    if (X.restriction == "colSums") {
      d <- cs
    } else if (X.restriction == "colSqSums") {
      d <- sqrt(colSums(X * X)) + small
    } else d <- cs
    X <- sweep(X, 2, d, "/")
    C <- sweep(C, 1, d, "*"); C <- sweep(C, 2, d, "*")
    list(X = X, C = C)
  }
  cc <- norm_XC(X, C); X <- cc$X; C <- cc$C
  compute_obj <- function(X, C) {
    Yhat <- X %*% C %*% t(X)
    ## lm()-style weighted least squares: L = sum(W * (Y - Yhat)^2).
    ## W = W^2 for binary {0,1} masks, so this is unchanged for ECV.
    if (has.weights) sum(Wmat * (Y - Yhat)^2) else sum((Y - Yhat)^2)
  }
  obj_prev <- compute_obj(X, C)
  objfunc.iter <- numeric(maxit)
  iter <- 0L
  if (has.weights) WmatY <- Wmat * Y
  for (iter in seq_len(maxit)) {
    tX <- t(X)

    ## ---- C update (symmetric tri) ----
    if (has.weights) {
      Yhat <- X %*% C %*% tX
      num_C <- tX %*% WmatY %*% X
      den_C <- tX %*% (Wmat * Yhat) %*% X
    } else {
      num_C <- tX %*% Y %*% X
      P <- tX %*% X
      den_C <- P %*% C %*% P
    }
    if (C.L1 > 0) den_C <- den_C + (C.L1 / 2)
    C <- C * (num_C / (den_C + small))

    ## ---- X update (bilateral gradient for Y = X C X^T, symmetric C) ----
    ##   num = (W*Y) X (C + C^T)
    ##   den = (W*Yhat) X (C + C^T)
    Yhat <- X %*% C %*% tX
    Ssym <- C + t(C)
    if (has.weights) {
      WY <- WmatY; WYhat <- Wmat * Yhat
    } else {
      WY <- Y; WYhat <- Yhat
    }
    num_X <- WY    %*% X %*% Ssym
    den_X <- WYhat %*% X %*% Ssym
    if (X.L2.ortho > 0) {
      XtX <- crossprod(X); diag(XtX) <- 0
      den_X <- den_X + X.L2.ortho * (X %*% XtX)
    }
    X <- X * (num_X / (den_X + small))

    cc <- norm_XC(X, C); X <- cc$X; C <- cc$C

    obj_cur <- compute_obj(X, C)
    objfunc.iter[iter] <- obj_cur
    if (iter > 1) {
      rel <- abs(obj_prev - obj_cur) / (abs(obj_prev) + small)
      if (rel < epsilon) break
    }
    obj_prev <- obj_cur
  }
  if (iter == maxit && exists("rel") && rel >= epsilon)
    warning(paste0("maximum iterations (", maxit, ") reached..."))
  list(X = X, C = C, iter = iter,
       objfunc.iter = objfunc.iter[seq_len(iter)],
       objfunc = obj_prev)
}


## ==============================================================
## Internal runner for bi-symmetric: Y ~ X X^T (C = I_Q fixed)
## Cube-root damping (He et al. 2011) for monotone descent.
## ==============================================================

.nmfkc.net.run_once_bi <- function(Y, X, maxit, epsilon,
                                    X.restriction, X.L2.ortho,
                                    Wmat, has.weights) {
  N <- nrow(Y); Q <- ncol(X)
  small <- 1e-16
  compute_obj <- function(X) {
    Yhat <- tcrossprod(X)
    ## lm()-style weighted least squares (linear W).
    if (has.weights) sum(Wmat * (Y - Yhat)^2) else sum((Y - Yhat)^2)
  }
  obj_prev <- compute_obj(X)
  objfunc.iter <- numeric(maxit)
  iter <- 0L
  if (has.weights) WmatY <- Wmat * Y
  for (iter in seq_len(maxit)) {
    Yhat <- tcrossprod(X)
    if (has.weights) {
      num_X <- WmatY    %*% X
      den_X <- (Wmat * Yhat) %*% X
    } else {
      num_X <- Y %*% X
      den_X <- X %*% crossprod(X)
    }
    if (X.L2.ortho > 0) {
      XtX <- crossprod(X); diag(XtX) <- 0
      den_X <- den_X + X.L2.ortho * (X %*% XtX)
    }
    update_ratio <- num_X / (den_X + small)
    X <- X * update_ratio^(1/3)

    if (X.restriction == "colSums") {
      cs <- colSums(X) + small
      X <- sweep(X, 2, cs, "/")
    } else if (X.restriction == "colSqSums") {
      d <- sqrt(colSums(X * X)) + small
      X <- sweep(X, 2, d, "/")
    }
    obj_cur <- compute_obj(X)
    objfunc.iter[iter] <- obj_cur
    if (iter > 1) {
      rel <- abs(obj_prev - obj_cur) / (abs(obj_prev) + small)
      if (rel < epsilon) break
    }
    obj_prev <- obj_cur
  }
  if (iter == maxit && exists("rel") && rel >= epsilon)
    warning(paste0("maximum iterations (", maxit, ") reached..."))
  list(X = X, C = diag(Q), iter = iter,
       objfunc.iter = objfunc.iter[seq_len(iter)],
       objfunc = obj_prev)
}


## ==============================================================
#' @title Symmetric NMF for networks (tri / bi / signed)
#' @description
#' Single entry point for symmetric NMF of network data with correct
#' multiplicative updates.  Three model types are supported via \code{type}:
#' \itemize{
#'   \item \strong{tri} (\code{type="tri"}, default):
#'     \eqn{Y \approx X C X^\top} with \eqn{X, C \ge 0} (both non-negative;
#'     \eqn{C} symmetric by design).  Uses Frobenius-full bilateral gradient.
#'   \item \strong{bi} (\code{type="bi"}):
#'     \eqn{Y \approx X X^\top} (\eqn{C} fixed to \eqn{I_Q}),
#'     cube-root damping (He et al. 2011).
#'   \item \strong{signed} (\code{type="signed"}):
#'     \eqn{Y \approx X (C_{+} - C_{-}) X^\top} with \eqn{X \ge 0} and signed
#'     \eqn{C = C_{+} - C_{-}}.  Preserves the soft-clustering interpretation
#'     of \eqn{X} while allowing negative off-diagonals of \eqn{C}
#'     (inter-cluster \emph{repulsion}).
#' }
#'
#' \strong{Non-negative adjacency matrix assumption.} All three types
#' assume \eqn{Y \ge 0} (a non-negative adjacency/affinity matrix).
#' The qualifier \dQuote{signed} in \code{type = "signed"} refers to the
#' middle coefficient \eqn{C}, not to \eqn{Y} itself.  The underlying
#' Ding, Li & Jordan (2010) sign-splitting updates require \eqn{Y \ge 0}
#' to guarantee monotone descent; supplying a signed \eqn{Y} triggers an
#' error.  For a signed data matrix, see \code{\link{nmfkc.signed}}.
#'
#' @param Y Symmetric (N x N) \strong{non-negative} adjacency matrix.
#'   \code{NA} entries are automatically treated as masked edges
#'   (equivalent to supplying \code{Y.weights} with 0 at those positions);
#'   see the note on \code{Y.weights} below.
#' @param rank Integer Q.
#' @param type \code{"tri"} (default), \code{"bi"}, or \code{"signed"}.
#' @param epsilon,maxit,verbose Standard.
#' @param ... Hidden options: \code{nstart} (default 1; see note below),
#'   \code{seed} (default 123), \code{X.restriction}, \code{X.init},
#'   \code{C.init} (tri only) or \code{Cp.init}/\code{Cn.init} (signed only),
#'   \code{Y.weights}, \code{C.L1} (tri only), \code{X.L2.ortho}, \code{prefix}.
#'
#' \strong{\code{Y.weights}} is an optional non-negative N x N weight
#' matrix (symmetric, same shape as \code{Y}).  When supplied, the loss
#' becomes \eqn{\sum W_{ij} \, (Y_{ij} - \hat Y_{ij})^2}
#' (\code{\link[stats]{lm}}()-style, \strong{linear} in \eqn{W}).
#' Logical matrices (\code{TRUE} / \code{FALSE}) are also accepted.
#' Typical usage by \code{\link{nmfkc.net.ecv}} is a binary mask
#' (\eqn{W \in \{0,1\}}) holding out test edges on the upper triangle;
#' real-valued weights for edge-level importance weighting are also
#' supported.  If \code{Y.weights} is \code{NULL} (default) and \code{Y}
#' contains \code{NA}, a binary mask is auto-generated (0 at \code{NA}
#' positions, 1 elsewhere), and the \code{NA} entries in \code{Y} are
#' replaced by 0 so the multiplicative updates can proceed.
#'
#' \strong{\code{X.init}} controls the initialization of the N x Q basis
#' matrix \eqn{X}.  Accepted values:
#' \itemize{
#'   \item \code{"kmeans"} (default): k-means on the rows of \eqn{Y}
#'     (equivalently columns, since \eqn{Y} is symmetric); the Q
#'     cluster centers become the columns of \eqn{X}.  Each node is
#'     treated as an N-dimensional connectivity profile, so clusters
#'     correspond to nodes with similar neighborhood structure --
#'     essentially a fast proxy for spectral clustering (Kuang, Yun &
#'     Park 2015, \emph{SymNMF}).  Scales well and is the recommended
#'     default for network data.
#'   \item \code{"kmeansar"}: \code{"kmeans"} followed by filling zero
#'     entries of \eqn{X} with \eqn{\mathrm{Uniform}(0, \bar Y / 100)}
#'     to escape trivial stationary points.
#'   \item \code{"nndsvd"}: Non-negative Double SVD with additive
#'     randomness (NNDSVDar).  Requires a full SVD of \eqn{Y}, so for
#'     very large networks (N > a few thousand) \code{"kmeans"} is
#'     preferable.
#'   \item \code{"runif"}: Uniform random entries in \eqn{[0, 1]}.
#'   \item \code{"random"}: Legacy default (pre-v0.6.8), equivalent to
#'     \code{abs(rnorm(N * Q)) * 0.1}.  Kept for backward compatibility.
#'   \item A numeric N x Q matrix supplied by the user (used as-is).
#' }
#' When \code{nstart > 1}, each restart uses a distinct seed so that
#' k-means / runif / NNDSVDar produce different candidate initial
#' values across the multi-start loop.
#'
#' \strong{Multi-start recommendation.} For \code{type = "signed"} the
#' \eqn{C = C_{+} - C_{-}} bottleneck can take both positive and negative
#' values, so the objective has more local minima than for \code{"tri"}
#' or \code{"bi"}.  A larger \code{nstart} (e.g., 10-50) is recommended
#' during exploration to reduce the chance of being trapped at a
#' suboptimal stationary point.  The default 1 is intended for fast
#' development; raise for publication-grade runs.
#' @return Object of class \code{c("nmfkc.net.<type>", "nmfkc.net", "nmfkc")}.
#'   For \code{type = "signed"} the return also carries \code{$Cp, $Cn}.
#' @seealso \code{\link{nmfkc.net.ecv}}, \code{\link{nmfkc.net.DOT}}
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @export
nmfkc.net <- function(Y, rank = 2, type = c("tri", "bi", "signed"),
                      epsilon = 1e-4, maxit = 5000,
                      verbose = FALSE, ...) {
  cl <- match.call()
  type <- match.arg(type)
  ex <- list(...)

  nstart        <- if (!is.null(ex$nstart))        ex$nstart        else 1L
  seed          <- if (!is.null(ex$seed))          ex$seed          else 123L
  X.restriction <- if (!is.null(ex$X.restriction)) ex$X.restriction
                   else if (type == "bi")          "none"
                   else                            "colSums"
  X.restriction <- match.arg(X.restriction,
                             c("colSums", "colSqSums", "none", "fixed"))
  X.init        <- if (!is.null(ex$X.init))        ex$X.init        else "kmeans"
  C.init        <- ex$C.init       # tri only
  Cp.init       <- ex$Cp.init      # signed only
  Cn.init       <- ex$Cn.init      # signed only
  Y.weights     <- ex$Y.weights
  C.L1          <- if (!is.null(ex$C.L1))          ex$C.L1          else 0
  X.L2.ortho    <- if (!is.null(ex$X.L2.ortho))    ex$X.L2.ortho    else 0
  prefix        <- if (!is.null(ex$prefix))        ex$prefix        else "Basis"

  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  if (nrow(Y) != ncol(Y)) stop("Y must be square for nmfkc.net.")
  if (!isSymmetric(Y, tol = 1e-10)) {
    Y <- (Y + t(Y)) / 2
    if (verbose) message("nmfkc.net: symmetrized Y <- (Y + Y^T) / 2")
  }
  N <- nrow(Y); Q <- as.integer(rank)
  small <- 1e-16; t0 <- proc.time()

  ## ---- Y.weights handling (supports auto-NA mask, user mask, logical) ----
  ## Matches the convention of nmfkc() / nmfkc.signed(): if Y has NA and
  ## Y.weights is NULL, a binary mask is auto-generated with 0 at the
  ## NA positions and 1 elsewhere.  NA is then replaced by 0 in Y so the
  ## multiplicative updates can proceed (masked entries contribute
  ## nothing to the objective since W == 0 there).
  has.weights <- !is.null(Y.weights)
  if (is.null(Y.weights) && anyNA(Y)) {
    Y.weights <- matrix(1, nrow = N, ncol = N)
    Y.weights[is.na(Y)] <- 0
    ## Also zero the symmetric mirror so W stays symmetric
    Y.weights[is.na(t(Y))] <- 0
    Y[is.na(Y)] <- 0
    has.weights <- TRUE
    if (verbose) message("nmfkc.net: auto-generated Y.weights mask from ",
                         sum(Y.weights == 0), " NA entries in Y.")
  } else if (has.weights) {
    Y.weights <- as.matrix(Y.weights)
    if (!all(dim(Y.weights) == c(N, N))) stop("Y.weights must be N x N.")
    storage.mode(Y.weights) <- "double"  # coerce logical T/F to 1/0
    if (!isSymmetric(Y.weights, tol = 1e-10))
      Y.weights <- (Y.weights + t(Y.weights)) / 2
    Y.weights[is.na(Y.weights)] <- 0
    Y[is.na(Y) | Y.weights == 0] <- 0
  } else if (anyNA(Y)) {
    stop("Y contains NA; please impute, remove, or supply Y.weights.")
  }
  Wmat <- if (has.weights) Y.weights else NULL

  if (min(Y) < 0) {
    ## nmfkc.net is designed for non-negative adjacency matrices in network
    ## analysis.  `type = "signed"` means the middle coefficient matrix C is
    ## allowed to be signed (enabling inhibitory/excitatory latent factor
    ## interactions), NOT that Y itself can be signed.  The Ding-style
    ## multiplicative updates require Y >= 0 in all three types; allowing a
    ## signed Y here would silently violate the monotone descent guarantee
    ## (the C+/C- numerators could go negative).
    stop("nmfkc.net(type='", type, "') requires Y >= 0 ",
         "(Y is expected to be a non-negative adjacency matrix; ",
         "`type = \"signed\"` refers to the middle coefficient C, not Y).")
  }

  ## ---- run_once dispatches to the appropriate MU kernel ----
  run_once <- function(X, C_or_Cp, Cn = NULL) {
    if (type == "bi") {
      .nmfkc.net.run_once_bi(Y, X, maxit, epsilon,
                              X.restriction, X.L2.ortho, Wmat, has.weights)
    } else if (type == "signed") {
      .nmfkc.net.signed.run_once(Y, X, C_or_Cp, Cn, maxit, epsilon,
                                  X.restriction, Wmat, has.weights)
    } else {
      .nmfkc.net.run_once(Y, X, C_or_Cp, maxit, epsilon,
                          X.restriction, C.L1, X.L2.ortho, Wmat, has.weights)
    }
  }

  ## ---- multi-start ----
  best <- NULL
  explicit_X <- is.matrix(X.init) || (is.numeric(X.init) && length(X.init) > 1)
  ## "random" is the legacy default -- kept as a backward-compat alias for
  ## abs(rnorm(N*Q)) * 0.1, since pre-v0.6.8 nmfkc.net() used this form.
  ## Other strings ("kmeans", "kmeansar", "nndsvd", "runif") delegate to
  ## the shared .init_X_method() helper for consistency with nmfkc() /
  ## nmf.sem().
  for (s in seq_len(nstart)) {
    s_seed <- seed + 7919L * (s - 1L)
    if (explicit_X) {
      X0 <- as.matrix(X.init)
      if (!identical(dim(X0), c(N, Q)))
        stop("X.init must have dimensions (nrow(Y), rank).")
    } else if (identical(X.init, "random")) {
      set.seed(s_seed)
      X0 <- matrix(abs(stats::rnorm(N * Q)) * 0.1, N, Q)
    } else if (is.character(X.init)) {
      ## "kmeans" / "kmeansar" / "nndsvd" / "runif" via the shared helper.
      ## Each restart uses a distinct seed so k-means and runif produce
      ## different candidates across nstart; nndsvdar's additive randomness
      ## also varies.
      X0 <- .init_X_method(X.init, Y, Q, seed = s_seed)
    } else {
      stop("X.init must be one of \"kmeans\", \"kmeansar\", \"nndsvd\", ",
           "\"runif\", \"random\", or a numeric (N x Q) matrix; got ",
           class(X.init)[1], ".")
    }

    if (type == "bi") {
      out <- run_once(X0, diag(Q))
    } else if (type == "signed") {
      if (!is.null(Cp.init) && !is.null(Cn.init)) {
        Cp0 <- as.matrix(Cp.init); Cn0 <- as.matrix(Cn.init)
        if (!isSymmetric(Cp0, tol = 1e-10) || !isSymmetric(Cn0, tol = 1e-10)) {
          warning("nmfkc.net(type='signed'): Cp.init/Cn.init symmetrized via (M + M^T)/2.")
          Cp0 <- (Cp0 + t(Cp0)) / 2
          Cn0 <- (Cn0 + t(Cn0)) / 2
        }
      } else {
        set.seed(s_seed + 1L)
        C0 <- matrix(stats::runif(Q * Q, min = -1, max = 1), Q, Q)
        C0 <- (C0 + t(C0)) / 2       # symmetrize before sign split
        Cp0 <- pmax(C0, 0); Cn0 <- pmax(-C0, 0)
      }
      out <- run_once(X0, Cp0, Cn0)
    } else {  # tri
      if (!is.null(C.init)) {
        C0 <- as.matrix(C.init)
        if (!identical(dim(C0), c(Q, Q)))
          stop("C.init must have dimensions (rank, rank).")
        if (any(C0 < 0)) stop("C.init must be non-negative for type='tri'.")
        if (!isSymmetric(C0, tol = 1e-10)) {
          warning("nmfkc.net: C.init symmetrized via (C + C^T)/2.")
          C0 <- (C0 + t(C0)) / 2
        }
      } else {
        set.seed(s_seed + 1L)
        C0 <- matrix(abs(stats::rnorm(Q * Q)) * 0.1, Q, Q)
        C0 <- (C0 + t(C0)) / 2       # symmetrize for exact bilateral gradient
      }
      out <- run_once(X0, C0)
    }
    if (is.null(best) || out$objfunc < best$objfunc) best <- out
  }

  ## ---- extract best result ----
  X <- best$X
  objfunc.iter <- best$objfunc.iter
  iter <- best$iter
  if (type == "signed") {
    Cp <- best$Cp; Cn <- best$Cn; C <- Cp - Cn
  } else {
    C <- best$C; Cp <- NULL; Cn <- NULL
  }

  ## ---- centroid reorder ----
  if (Q > 1 && X.restriction != "fixed") {
    idx <- order(as.vector(matrix(seq_len(N) / N, nrow = 1) %*% X))
    X <- X[, idx, drop = FALSE]
    C <- C[idx, idx, drop = FALSE]
    if (!is.null(Cp)) Cp <- Cp[idx, idx, drop = FALSE]
    if (!is.null(Cn)) Cn <- Cn[idx, idx, drop = FALSE]
  }

  ## ---- names ----
  rownames(X) <- rownames(Y)
  colnames(X) <- paste0(prefix, seq_len(Q))
  rownames(C) <- colnames(C) <- paste0(prefix, seq_len(Q))
  if (!is.null(Cp)) { rownames(Cp) <- colnames(Cp) <- paste0(prefix, seq_len(Q)) }
  if (!is.null(Cn)) { rownames(Cn) <- colnames(Cn) <- paste0(prefix, seq_len(Q)) }

  ## ---- reconstruction and fit statistics ----
  Y1hat <- X %*% C %*% t(X)
  resid <- Y - Y1hat
  ## lm()-style weighted least squares (matches compute_obj inside the loop).
  objfunc <- if (has.weights) sum(Wmat * resid^2) else sum(resid^2)
  ## Unified three-variant R^2 (cor^2, uncentered,
  ## row-mean centered), all respecting Y.weights == 0 masking.
  r2_all <- .r.squared.all(Y, Y1hat,
                           Y.weights = if (has.weights) Wmat else NULL)
  r.squared          <- r2_all$r.squared
  r.squared.uncentered     <- r2_all$r.squared.uncentered
  r.squared.centered <- r2_all$r.squared.centered
  mae <- if (has.weights) sum(Wmat * abs(resid)) / max(sum(Wmat), small)
         else mean(abs(resid))
  ## sigma (RMSE) for parity with nmfkc / nmfkc.signed fit objects
  sigma <- if (has.weights) sqrt(sum(Wmat * resid^2) / max(sum(Wmat), small))
           else sqrt(mean(resid^2))

  X.prob <- X / (rowSums(X) + small)
  X.cluster <- apply(X.prob, 1, which.max)

  result <- list(
    call = cl,
    X = X, C = C, Cp = Cp, Cn = Cn,
    XB = Y1hat, Y1hat = Y1hat,          # XB alias for fitted.nmf compatibility
    X.prob = X.prob, X.cluster = X.cluster,
    rank = Q, dims = c(N = N),
    type = type,
    objfunc = objfunc, objfunc.iter = objfunc.iter,
    r.squared          = r.squared,
    r.squared.uncentered     = r.squared.uncentered,
    r.squared.centered = r.squared.centered,
    mae = mae,
    sigma = sigma,
    iter = iter,
    runtime = as.numeric((proc.time() - t0)[3]),
    X.restriction = X.restriction
  )
  class(result) <- c(paste0("nmfkc.net.", type), "nmfkc.net", "nmfkc", "nmf")
  result
}


## ==============================================================
## Signed: symmetric tri-NMF with signed C = Cp - Cn
## (X >= 0 for soft clustering; C signed for inter-cluster repulsion)
## ==============================================================

.nmfkc.net.signed.run_once <- function(Y, X, Cp, Cn, maxit, epsilon,
                                       X.restriction, Wmat, has.weights) {
  N <- nrow(Y); Q <- ncol(X)
  small <- 1e-16
  norm_XC <- function(X, Cp, Cn) {
    if (X.restriction %in% c("fixed", "none"))
      return(list(X = X, Cp = Cp, Cn = Cn))
    cs <- colSums(X) + small
    d <- switch(X.restriction,
                colSums   = cs,
                colSqSums = sqrt(colSums(X * X)) + small,
                cs)
    X <- sweep(X, 2, d, "/")
    Cp <- sweep(Cp, 1, d, "*"); Cp <- sweep(Cp, 2, d, "*")
    Cn <- sweep(Cn, 1, d, "*"); Cn <- sweep(Cn, 2, d, "*")
    list(X = X, Cp = Cp, Cn = Cn)
  }
  cc <- norm_XC(X, Cp, Cn); X <- cc$X; Cp <- cc$Cp; Cn <- cc$Cn
  compute_obj <- function(X, Cp, Cn) {
    Yhat <- X %*% (Cp - Cn) %*% t(X)
    ## lm()-style weighted least squares (linear W).
    if (has.weights) sum(Wmat * (Y - Yhat)^2) else sum((Y - Yhat)^2)
  }
  obj_prev <- compute_obj(X, Cp, Cn)
  objfunc.iter <- numeric(maxit)
  iter <- 0L
  if (has.weights) WmatY <- Wmat * Y

  for (iter in seq_len(maxit)) {
    tX <- t(X)
    P <- tX %*% X

    ## ---- Cp, Cn update (Ding sign split) ----
    if (has.weights) {
      Yh_p <- X %*% Cp %*% tX
      Yh_n <- X %*% Cn %*% tX
      G    <- tX %*% WmatY %*% X
      Hp   <- tX %*% (Wmat * Yh_p) %*% X
      Hn   <- tX %*% (Wmat * Yh_n) %*% X
    } else {
      G  <- tX %*% Y %*% X
      Hp <- P %*% Cp %*% P
      Hn <- P %*% Cn %*% P
    }
    Cp <- Cp * (G + Hn) / (Hp + small)
    if (has.weights) {
      Yh_p <- X %*% Cp %*% tX
      Hp   <- tX %*% (Wmat * Yh_p) %*% X
    } else Hp <- P %*% Cp %*% P
    Cn <- Cn * (Hp) / (G + Hn + small)

    ## ---- X update (bilateral + Ding split, exact for symmetric Cp, Cn) ----
    Sp <- Cp + t(Cp); Sn <- Cn + t(Cn)
    Yh_p <- X %*% Cp %*% tX
    Yh_n <- X %*% Cn %*% tX
    if (has.weights) {
      WY   <- WmatY;    WYhp <- Wmat * Yh_p;    WYhn <- Wmat * Yh_n
    } else {
      WY <- Y; WYhp <- Yh_p; WYhn <- Yh_n
    }
    num_X <- WY %*% X %*% Sp + WYhp %*% X %*% Sn + WYhn %*% X %*% Sp
    den_X <- WY %*% X %*% Sn + WYhp %*% X %*% Sp + WYhn %*% X %*% Sn
    X <- X * (num_X / (den_X + small))

    cc <- norm_XC(X, Cp, Cn); X <- cc$X; Cp <- cc$Cp; Cn <- cc$Cn

    obj_cur <- compute_obj(X, Cp, Cn)
    objfunc.iter[iter] <- obj_cur
    if (iter > 1) {
      rel <- abs(obj_prev - obj_cur) / (abs(obj_prev) + small)
      if (rel < epsilon) break
    }
    obj_prev <- obj_cur
  }
  if (iter == maxit && exists("rel") && rel >= epsilon)
    warning(paste0("maximum iterations (", maxit, ") reached..."))
  list(X = X, Cp = Cp, Cn = Cn, iter = iter,
       objfunc.iter = objfunc.iter[seq_len(iter)],
       objfunc = obj_prev)
}


## ECV: upper-triangle element-wise cross-validation
## ==============================================================

.nmfkc.net.make_uppertri_folds <- function(N, div = 5, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  ut_idx <- which(upper.tri(matrix(0, N, N), diag = TRUE))
  perm <- sample(ut_idx)
  n <- length(perm)
  chunk <- n %/% div; rem <- n %% div
  folds <- vector("list", div)
  s <- 1L
  for (k in 1:div) {
    sz <- chunk + ifelse(k <= rem, 1L, 0L)
    folds[[k]] <- perm[s:(s + sz - 1L)]
    s <- s + sz
  }
  folds
}

## Mirror a 0/1 mask to be symmetric via AND (upper & lower both retained).
## For general (non-binary) weights use (W + t(W))/2 instead.
.nmfkc.net.mirror_mask <- function(W) W * t(W)

#' Element-wise cross-validation for nmfkc.net (upper-triangle folds)
#'
#' @description k-fold CV with folds taken over the upper triangle of the
#' symmetric \eqn{Y} (mirrored to the lower triangle) to prevent information
#' leakage through symmetry.  A single entry point covers all three symmetric
#' NMF variants; \code{type} selects the fitting function:
#' \itemize{
#'   \item \code{"tri"} (default): \code{\link{nmfkc.net}} with \eqn{C \ge 0}
#'   \item \code{"bi"}: \code{\link{nmfkc.net}} with \eqn{C = I_Q}
#'   \item \code{"signed"}: \code{\link{nmfkc.net}} with signed \eqn{C = C_{+} - C_{-}}
#' }
#'
#' @param Y Symmetric N x N non-negative matrix.
#' @param rank Integer vector of ranks to evaluate. Default \code{1:3}.
#' @param type Model type: \code{"tri"} (default), \code{"bi"}, or \code{"signed"}.
#' @param ... Passed to the underlying fitter; also accepts \code{nfolds}
#'   (default 5; \code{div} alias), \code{seed} (default 123).
#'
#' @return A list with \code{objfunc}, \code{sigma}, \code{objfunc.fold},
#'   \code{folds}, \code{Q.grid}, \code{type}.
#' @seealso \code{\link{nmfkc.net}}
#' @section Lifecycle:
#' This function is \strong{experimental}. The interface may change in
#' future versions; details are to be described in an upcoming paper.
#'
#' @export
nmfkc.net.ecv <- function(Y, rank = 1:3,
                          type = c("tri", "bi", "signed"), ...) {
  type <- match.arg(type)
  ex <- list(...)
  nfolds <- if (!is.null(ex$nfolds)) ex$nfolds
            else if (!is.null(ex$div)) ex$div else 5
  seed <- if (!is.null(ex$seed)) ex$seed else 123
  Y <- as.matrix(Y); N <- nrow(Y)
  folds <- .nmfkc.net.make_uppertri_folds(N, div = nfolds, seed = seed)

  fit_args <- ex; fit_args$nfolds <- NULL; fit_args$div <- NULL
  fit_args$rank <- NULL; fit_args$type <- NULL

  run_one <- function(q, k) {
    test_ut <- folds[[k]]
    W <- matrix(1, N, N); W[test_ut] <- 0
    W <- .nmfkc.net.mirror_mask(W)
    fit <- suppressMessages(do.call(nmfkc.net,
             c(list(Y = Y, rank = q, type = type,
                    Y.weights = W, verbose = FALSE), fit_args)))
    Yhat <- fit$X %*% fit$C %*% t(fit$X)
    mean((Y[test_ut] - Yhat[test_ut])^2)
  }

  message(sprintf("nmfkc.net ECV (type=%s): %d ranks, %d-fold, upper-triangle.",
                  type, length(rank), nfolds))
  cv <- .ecv.run(sprintf("Q=%d", rank), nfolds,
                 run_one = function(i, k) run_one(rank[i], k),
                 progress = function(i, o, s)
                   message(sprintf("  Q=%d: MSE=%.6f, sigma=%.4f", rank[i], o, s)))
  cls <- if (type == "signed")
           c("nmfkc.net.signed.ecv", "nmfkc.net.ecv", "nmfkc.ecv")
         else
           c("nmfkc.net.ecv", "nmfkc.ecv")
  structure(list(objfunc = cv$objfunc, sigma = cv$sigma,
                 rank = rank, nfolds = nfolds,
                 objfunc.fold = cv$objfunc.fold, folds = folds,
                 Q.grid = rank, type = type),
            class = cls)
}


#' @title Rank selection for nmfkc.net (concise diagnostics)
#' @description
#' Fits \code{\link{nmfkc.net}} across a range of ranks and reports the
#' three rank-selection criteria -- \code{r.squared}, the effective rank
#' (utilization), and the element-wise CV error \code{sigma.ecv} -- with
#' the same concise diagnostics plot as \code{\link{nmfkc.rank}}.
#' @param Y Symmetric (network) observation matrix.
#' @param rank Integer vector of ranks to evaluate.
#' @param type One of \code{"tri"}, \code{"bi"}, \code{"signed"}.
#' @param detail \code{"full"} (default) also runs element-wise CV
#'   (\code{sigma.ecv}); \code{"fast"} skips it (plots r.squared and
#'   eff.rank only, and recommends the R-squared elbow).
#' @param plot Logical; draw the diagnostics plot (default \code{TRUE}).
#' @param ... Passed on to \code{\link{nmfkc.net}} and
#'   \code{\link{nmfkc.net.ecv}} (e.g.\ \code{nstart}, \code{maxit},
#'   \code{nfolds}, \code{seed}).
#' @return A list with \code{rank.best} (ECV minimum, or the R-squared
#'   elbow under \code{detail = "fast"}) and \code{criteria} (data
#'   frame: \code{rank}, \code{effective.rank}, \code{effective.rank.ratio},
#'   \code{r.squared}, \code{sigma.ecv}).
#' @seealso \code{\link{nmfkc.net}}, \code{\link{nmfkc.net.ecv}},
#'   \code{\link{nmfkc.rank}}
#' @references
#' Roy, O., & Vetterli, M. (2007).  The effective rank: A measure of
#' effective dimensionality.  \emph{Proc. EUSIPCO}, 606--610.
#' (\code{effective.rank})
#' Wold, S. (1978).  Cross-validatory estimation of the number of
#' components in factor and principal components models.
#' \emph{Technometrics}, 20(4), 397--405. (\code{sigma.ecv})
#' @export
#' @examples
#' \donttest{
#' Y <- matrix(c(0,1,1,0,0,0, 1,0,1,0,0,0, 1,1,0,1,0,0,
#'               0,0,1,0,1,1, 0,0,0,1,0,1, 0,0,0,1,1,0), 6, 6)
#' nmfkc.net.rank(Y, rank = 1:3, type = "tri", nstart = 5, nfolds = 3)
#' }
nmfkc.net.rank <- function(Y, rank = 1:5, type = c("tri", "bi", "signed"),
                           detail = c("full", "fast"), plot = TRUE, ...) {
  type <- match.arg(type); detail <- match.arg(detail)
  Y <- as.matrix(Y)
  rs <- numeric(length(rank)); er <- numeric(length(rank))
  for (i in seq_along(rank)) {
    f <- suppressMessages(nmfkc.net(Y, rank = rank[i], type = type,
                                    verbose = FALSE, ...))
    rs[i] <- f$r.squared
    er[i] <- .effective.rank(t(f$X))
  }
  ecv <- if (detail == "full")
    suppressMessages(nmfkc.net.ecv(Y, rank = rank, type = type, ...))$sigma
    else rep(NA_real_, length(rank))
  criteria <- data.frame(rank = rank, effective.rank = er,
                         effective.rank.ratio = er / rank,
                         r.squared = rs, sigma.ecv = as.numeric(ecv))
  .rank.finish(criteria, plot = plot,
               main = sprintf("nmfkc.net rank selection (type=%s)", type))
}


## ==============================================================
## Summary methods for nmfkc.net and nmfkc.net.signed
## (Format unified with print.summary.nmfkc)
## ==============================================================

#' Summary method for nmfkc.net objects
#' @param object An \code{nmfkc.net} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfkc.net"}.
#' @seealso \code{\link{nmfkc.net}}
#' @export
summary.nmfkc.net <- function(object, ...) {
  N <- object$dims["N"]; Q <- object$rank
  small <- 1e-16
  dims_str <- sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)X'", N, N, N, Q, Q, Q)

  ans <- list(
    call = object$call,
    dims = dims_str,
    rank = Q,
    type = if (!is.null(object$type)) object$type else "tri",
    iter = object$iter,
    runtime = sprintf("%.2fsec", object$runtime),

    objfunc   = object$objfunc,
    r.squared          = object$r.squared,
    r.squared.uncentered     = object$r.squared.uncentered,
    r.squared.centered = object$r.squared.centered,
    mae       = object$mae,
    ## Effective rank over the Q factors' membership across the N nodes
    ## (X is N x Q, so transpose to factors x nodes).
    effective.rank = .effective.rank(t(object$X)),

    X.sparsity      = mean(object$X  < 1e-4),
    C.sparsity      = mean(abs(object$C) < 1e-4),

    X.range = range(object$X),
    C.range = range(object$C),
    X.cluster.counts = table(object$X.cluster)
  )
  class(ans) <- "summary.nmfkc.net"
  ans
}

#' Print method for summary.nmfkc.net objects
#' @param x An object of class \code{"summary.nmfkc.net"}.
#' @param digits Minimum number of significant digits.
#' @param ... Unused.
#' @return Invisible \code{x}.
#' @seealso \code{\link{summary.nmfkc.net}}
#' @export
print.summary.nmfkc.net <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Dimensions:", x$dims, "\n")
  cat("Rank (Q):   ", x$rank, "\n")
  cat("Type:       ", x$type, "\n")
  cat("Runtime:    ", x$runtime, "\n")
  cat("Iterations: ", x$iter, "\n")

  .print.fit.statistics(x, digits = digits)

  .print.structure.diagnostics(
    sparsity = c("Basis (X)" = x$X.sparsity, "C" = x$C.sparsity))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "X range:", x$X.range[1], x$X.range[2]))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "C range:", x$C.range[1], x$C.range[2]))

  cat("\nCluster sizes (hard assignment by argmax of X):\n  ")
  print(x$X.cluster.counts)
  cat("\n")
  invisible(x)
}

#' Summary method for nmfkc.net.signed objects
#' @param object An \code{nmfkc.net.signed} object.
#' @param ... Unused.
#' @return An object of class \code{"summary.nmfkc.net.signed"}.
#' @seealso \code{\link{nmfkc.net}}
#' @export
summary.nmfkc.net.signed <- function(object, ...) {
  ans <- summary.nmfkc.net(object, ...)
  ans$Cp.range <- range(object$Cp)
  ans$Cn.range <- range(object$Cn)
  ans$neg.frac <- sum(pmax(-object$C, 0)) / (sum(abs(object$C)) + 1e-16)
  class(ans) <- c("summary.nmfkc.net.signed", "summary.nmfkc.net")
  ans
}

#' Print method for summary.nmfkc.net.signed objects
#' @param x An object of class \code{"summary.nmfkc.net.signed"}.
#' @param digits Minimum number of significant digits.
#' @param ... Unused.
#' @return Invisible \code{x}.
#' @seealso \code{\link{summary.nmfkc.net.signed}}
#' @export
print.summary.nmfkc.net.signed <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Dimensions:", x$dims, "\n")
  cat("Rank (Q):   ", x$rank, "\n")
  cat("Type:       ", "signed (C = Cp - Cn)\n")
  cat("Runtime:    ", x$runtime, "\n")
  cat("Iterations: ", x$iter, "\n")

  .print.fit.statistics(x, digits = digits)

  .print.structure.diagnostics(
    sparsity = c("Basis (X)" = x$X.sparsity, "C" = x$C.sparsity))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "X range:",  x$X.range[1],  x$X.range[2]))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "Cp range:", x$Cp.range[1], x$Cp.range[2]))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "Cn range:", x$Cn.range[1], x$Cn.range[2]))
  cat(sprintf("  %-24s [%.3g, %.3g]\n", "C range:",  x$C.range[1],  x$C.range[2]))
  cat(sprintf("  %-24s %.1f%% (of sum |C|)\n", "Negative mass:", 100 * x$neg.frac))

  cat("\nCluster sizes (hard assignment by argmax of X):\n  ")
  print(x$X.cluster.counts)
  cat("\n")
  invisible(x)
}
