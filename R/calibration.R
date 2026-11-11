.mopair_refit_p <- function(fit, X, Y) {
  input <- fit$input
  path <- .mopair_fit_paths(X, input$M, Y, input$C, fit$metadata$family_y,
    fit$metadata$binomial_test, fit$metadata$cox_ties)
  path$module <- unname(fit$modules[path$mediator])
  path$composite_p <- .mopair_composite(path$a_p, path$b_p, fit$metadata$composite_method)
  combine <- if (fit$metadata$combination == "acat") .mopair_acat else .mopair_simes
  vapply(split(path$composite_p, path$module), combine, numeric(1))[fit$module_results$module]
}

.mopair_calibrate <- function(fit, n_perm, strata, seed) {
  if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 1L || n_perm != as.integer(n_perm)) stop("n_perm must be a positive integer.", call. = FALSE)
  if (n_perm < 20L) warning("n_perm < 20 is suitable only for software checks.", call. = FALSE)
  input <- fit$input
  if (is.null(strata)) strata <- rep("all", length(input$X))
  strata <- split(seq_along(input$X), as.factor(strata))
  permute <- function() {
    index <- seq_along(input$X)
    for (ii in strata) index[ii] <- ii[sample.int(length(ii))]
    index
  }
  h10 <- matrix(NA_real_, n_perm, nrow(fit$module_results)); h01 <- h10
  for (b in seq_len(n_perm)) {
    set.seed(seed + b); index <- permute()
    h10[b, ] <- .mopair_refit_p(fit, input$X, input$Y[index])
    set.seed(seed + 1000000L + b); index <- permute()
    h01[b, ] <- .mopair_refit_p(fit, input$X[index], input$Y)
  }
  if (any(!is.finite(h10)) || any(!is.finite(h01))) stop("Unavailable permutation outputs prevent complete sensitivity calibration; no missing values are counted as zero discoveries.", call. = FALSE)
  thresholds <- sort(unique(fit$module_results$module_p[is.finite(fit$module_results$module_p)]))
  observed <- vapply(thresholds, function(t) sum(fit$module_results$module_p <= t, na.rm = TRUE), integer(1))
  v10 <- vapply(thresholds, function(t) mean(rowSums(h10 <= t)), numeric(1))
  v01 <- vapply(thresholds, function(t) mean(rowSums(h01 <= t)), numeric(1))
  curve <- data.frame(threshold = thresholds, observed_discoveries = observed, expected_false_h10 = v10,
    expected_false_h01 = v01, expected_false = pmax(v10, v01))
  curve$empirical_fdr <- pmin(1, curve$expected_false / pmax(observed, 1))
  eligible <- which(observed > 0 & curve$empirical_fdr <= fit$metadata$alpha)
  choice <- if (length(eligible)) eligible[which.max(observed[eligible])] else NA_integer_
  threshold <- if (is.na(choice)) NA_real_ else thresholds[choice]
  fit$module_results$selected_calibrated <- is.finite(threshold) & is.finite(fit$module_results$module_p) & fit$module_results$module_p <= threshold
  fit$calibration <- list(type = "dual_permutation_max_envelope", curve = curve, threshold = threshold,
    empirical_fdr = if (is.na(choice)) NA_real_ else curve$empirical_fdr[choice], n_perm = as.integer(n_perm),
    null_h10_module_p = h10, null_h01_module_p = h01, sensitivity_only = TRUE)
  fit
}
