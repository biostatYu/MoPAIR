#' Fit module-level paired-path association inference
#'
#' See the package manual and README for input definitions and examples.
fit_mopair <- function(X, M, Y, C = NULL,
                       family_y = c("gaussian", "binomial", "cox"),
                       alpha = 0.10, modules = NULL, corr_cut = 0.45,
                       min_module_size = 2L, hclust_method = "average",
                       composite_method = c("dact", "maxp"),
                       combination = c("acat", "simes"), leading_edge_n = 3L,
                       binomial_test = c("score", "wald"),
                       cox_ties = c("efron", "breslow", "exact"),
                       calibration = c("none", "dual_permutation"),
                       n_perm = 200L, strata = NULL, seed = 1L) {
  family_y <- match.arg(family_y)
  composite_method <- match.arg(composite_method)
  combination <- match.arg(combination)
  calibration <- match.arg(calibration)
  binomial_test <- match.arg(binomial_test)
  cox_ties <- match.arg(cox_ties)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).", call. = FALSE)
  if (length(leading_edge_n) != 1L || !is.finite(leading_edge_n) || leading_edge_n < 0) stop("leading_edge_n must be nonnegative.", call. = FALSE)
  if (length(seed) != 1L || !is.finite(seed) || abs(seed) > .Machine$integer.max - 1000000L) stop("seed must be a finite integer in range.", call. = FALSE)
  if (family_y == "cox" && calibration != "none") stop("Cox models support BH selection only; censoring-aware permutation calibration is not implemented.", call. = FALSE)
  dat <- .mopair_validate(X, M, Y, C, family_y)
  X <- dat$X; M <- dat$M; Y <- dat$Y; C <- dat$C
  set.seed(seed)
  info <- .mopair_prepare_modules(X, M, C, modules, corr_cut, as.integer(min_module_size), hclust_method)
  paths <- .mopair_fit_paths(X, M, Y, C, family_y, binomial_test, cox_ties)
  paths$module <- unname(info$module[paths$mediator])
  paths$composite_p <- .mopair_composite(paths$a_p, paths$b_p, composite_method)
  combine <- if (combination == "acat") .mopair_acat else .mopair_simes
  groups <- split(paths$composite_p, paths$module)
  module_results <- data.frame(module = names(groups), n_mediators = lengths(groups),
    n_tested = vapply(groups, function(z) sum(is.finite(z)), integer(1)),
    module_p = vapply(groups, combine, numeric(1)), stringsAsFactors = FALSE, row.names = NULL)
  module_results$module_q <- stats::p.adjust(module_results$module_p, method = "BH")
  module_results$selected_bh <- is.finite(module_results$module_q) & module_results$module_q <= alpha
  paths <- merge(paths, module_results, by = "module", all.x = TRUE, sort = FALSE)
  paths <- paths[order(paths$module_q, paths$composite_p, paths$mediator), , drop = FALSE]
  paths$leading_edge_rank <- stats::ave(paths$composite_p, paths$module, FUN = function(z) rank(z, ties.method = "first", na.last = "keep"))
  paths$leading_edge <- paths$selected_bh & is.finite(paths$leading_edge_rank) & paths$leading_edge_rank <= leading_edge_n
  rownames(paths) <- NULL
  module_results <- module_results[order(module_results$module_q, module_results$module), , drop = FALSE]
  rownames(module_results) <- NULL
  fit <- list(call = match.call(), module_results = module_results, mediator_results = paths,
    modules = info$module, selected_modules = module_results[module_results$selected_bh, , drop = FALSE],
    leading_edge = paths[paths$leading_edge, , drop = FALSE], calibration = NULL,
    metadata = list(method = "MoPAIR", version = "0.2.0", family_y = family_y, alpha = alpha,
      binomial_test = if (family_y == "binomial") binomial_test else NA_character_,
      cox_ties = if (family_y == "cox") cox_ties else NA_character_,
      module_source = info$source, corr_cut = if (is.null(modules)) corr_cut else NA_real_,
      hclust_method = if (is.null(modules)) hclust_method else NA_character_,
      composite_method = composite_method, dact_correction = if (composite_method == "dact") "none" else NA_character_,
      combination = combination, inference_target = "module", individual_confirmed = FALSE,
      n_input = dat$n_input, n_complete = length(X), n_dropped = dat$n_input - length(X),
      row_index = dat$row_index, n_unavailable_members = sum(!is.finite(paths$composite_p)), seed = as.integer(seed)),
    input = list(X = X, M = M, Y = Y, C = C))
  class(fit) <- "mopair_fit"
  if (calibration == "dual_permutation") {
    if (!is.null(strata)) {
      if (length(strata) != dat$n_input || anyNA(strata)) stop("strata must be complete and have the original sample length.", call. = FALSE)
      strata <- strata[dat$row_index]
    }
    fit <- .mopair_calibrate(fit, n_perm, strata, seed)
  }
  fit$input <- NULL
  fit
}

print.mopair_fit <- function(x, ...) {
  cat("MoPAIR module fit\n")
  cat("  outcome:", x$metadata$family_y, "\n")
  if (x$metadata$family_y == "binomial") cat("  outcome-path test:", x$metadata$binomial_test, "\n")
  cat("  complete observations:", x$metadata$n_complete, "of", x$metadata$n_input, "\n")
  cat("  modules:", nrow(x$module_results), "\n")
  cat("  BH-selected:", sum(x$module_results$selected_bh), "at q <=", x$metadata$alpha, "\n")
  if (!is.null(x$calibration)) cat("  sensitivity selections:", sum(x$module_results$selected_calibrated), "\n")
  invisible(x)
}
