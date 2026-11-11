fit_mopair_grid <- function(X, M, Y, ...) {
  if (!is.numeric(X)) stop("X must be a numeric sample-by-exposure matrix.", call. = FALSE)
  X <- as.matrix(X)
  if (!ncol(X)) stop("X needs at least one exposure.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (anyDuplicated(colnames(X))) stop("Exposure names must be unique.", call. = FALSE)
  args <- list(...)
  if (is.null(args$modules)) stop("fit_mopair_grid requires common fixed modules.", call. = FALSE)
  if (!is.null(args$calibration) && args$calibration != "none") stop("Grid analysis supports family-wide BH only; a complete-family calibration driver is required for sensitivity calibration.", call. = FALSE)
  fits <- lapply(seq_len(ncol(X)), function(j) do.call(fit_mopair, c(list(X = X[, j], M = M, Y = Y), args)))
  names(fits) <- colnames(X)
  family <- do.call(rbind, lapply(seq_along(fits), function(j) {
    out <- fits[[j]]$module_results
    out$exposure <- names(fits)[j]
    names(out)[names(out) == "selected_bh"] <- "selected_within_exposure"
    out
  }))
  family$family_q <- stats::p.adjust(family$module_p, method = "BH")
  alpha <- fits[[1L]]$metadata$alpha
  family$selected_family <- is.finite(family$family_q) & family$family_q <= alpha
  family <- family[order(family$family_q, family$exposure, family$module), , drop = FALSE]
  rownames(family) <- NULL
  list(fits = fits, family_results = family, selected_hypotheses = family[family$selected_family, , drop = FALSE],
    metadata = list(method = "MoPAIR", alpha = alpha, n_exposures = ncol(X), n_hypotheses = nrow(family)))
}
