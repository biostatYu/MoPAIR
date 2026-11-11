.mopair_residualize <- function(X, M, C) {
  design <- .mopair_design(X, C)
  residual <- matrix(NA_real_, nrow(M), ncol(M), dimnames = dimnames(M))
  for (j in seq_len(ncol(M))) {
    ok <- stats::complete.cases(design, M[, j])
    if (sum(ok) <= ncol(design)) next
    residual[ok, j] <- stats::lm.fit(design[ok, , drop = FALSE], M[ok, j])$residuals
  }
  residual
}

.mopair_prepare_modules <- function(X, M, C, modules, corr_cut, min_module_size, hclust_method) {
  if (!is.null(modules)) {
    ids <- names(modules)
    modules <- as.character(modules)
    if (is.null(ids)) {
      if (length(modules) != ncol(M)) stop("Unnamed modules must have length ncol(M).", call. = FALSE)
      ids <- colnames(M)
    }
    if (anyNA(ids) || anyDuplicated(ids)) stop("Named modules must have unique mediator identifiers.", call. = FALSE)
    names(modules) <- ids
    missing <- setdiff(colnames(M), ids)
    if (length(missing)) stop("modules is missing mediators: ", paste(missing, collapse = ", "), call. = FALSE)
    modules <- modules[colnames(M)]
    if (anyNA(modules) || any(!nzchar(modules))) stop("Every mediator needs a nonempty module label.", call. = FALSE)
    return(list(module = modules, source = "provided_fixed"))
  }
  if (length(corr_cut) != 1L || !is.finite(corr_cut) || corr_cut <= 0 || corr_cut >= 1) stop("corr_cut must be in (0, 1).", call. = FALSE)
  if (length(min_module_size) != 1L || is.na(min_module_size) || min_module_size < 1L) stop("min_module_size must be positive.", call. = FALSE)
  if (ncol(M) == 1L) return(list(module = stats::setNames("module_1", colnames(M)), source = "same_sample_outcome_free"))
  residual <- .mopair_residualize(X, M, C)
  correlation <- suppressWarnings(stats::cor(residual, use = "pairwise.complete.obs"))
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  tree <- stats::hclust(stats::as.dist(1 - abs(correlation)), method = hclust_method)
  groups <- stats::cutree(tree, h = 1 - corr_cut)
  sizes <- table(groups)
  labels <- paste0("module_", groups)
  small <- as.integer(sizes[as.character(groups)]) < min_module_size
  labels[small] <- paste0("singleton_", which(small))
  list(module = stats::setNames(labels, colnames(M)), source = "same_sample_outcome_free")
}
