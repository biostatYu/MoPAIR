.mopair_clip_p <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  pmin(pmax(p, .Machine$double.xmin), 1)
}

.mopair_validate <- function(X, M, Y, C, family_y) {
  if (!is.numeric(X) || !is.numeric(M)) stop("X and M must be numeric.", call. = FALSE)
  M <- as.matrix(M)
  X <- as.numeric(X)
  n <- length(X)
  if (family_y == "cox") {
    if (!inherits(Y, "Surv")) stop("For Cox models, Y must be survival::Surv(time, status).", call. = FALSE)
    if (!identical(attr(Y, "type"), "right")) stop("Only right-censored survival outcomes are supported.", call. = FALSE)
    ny <- nrow(Y)
    y_ok <- apply(is.finite(as.matrix(Y)), 1L, all)
    if (any(Y[, 1L] < 0, na.rm = TRUE)) stop("Survival times must be nonnegative.", call. = FALSE)
  } else {
    if (!is.numeric(Y) || is.matrix(Y)) stop("Y must be a numeric outcome vector.", call. = FALSE)
    Y <- as.numeric(Y)
    ny <- length(Y)
    y_ok <- is.finite(Y)
  }
  if (!n || nrow(M) != n || ny != n) stop("X, Y, and rows of M must have the same nonzero sample size.", call. = FALSE)
  if (ncol(M) < 1L) stop("M must contain at least one mediator.", call. = FALSE)
  if (is.null(colnames(M))) colnames(M) <- paste0("M", seq_len(ncol(M)))
  if (anyNA(colnames(M)) || any(!nzchar(colnames(M))) || anyDuplicated(colnames(M))) stop("Mediator names must be unique and nonempty.", call. = FALSE)
  C <- if (is.null(C)) NULL else as.data.frame(C)
  if (!is.null(C) && nrow(C) != n) stop("C must have one row per sample.", call. = FALSE)
  c_ok <- rep(TRUE, n)
  if (!is.null(C)) {
    if (!all(vapply(C, function(z) is.numeric(z) || is.factor(z) || is.logical(z), logical(1)))) stop("Covariates must be numeric, logical, or factors.", call. = FALSE)
    names(C) <- paste0("C", seq_len(ncol(C)))
    c_ok <- stats::complete.cases(C)
    for (z in C) if (is.numeric(z)) c_ok <- c_ok & is.finite(z)
  }
  keep <- is.finite(X) & y_ok & c_ok
  if (sum(keep) < 4L) stop("Too few complete exposure/outcome/covariate observations.", call. = FALSE)
  X <- X[keep]; M <- M[keep, , drop = FALSE]
  Y <- if (family_y == "cox") Y[keep, , drop = FALSE] else Y[keep]
  if (!is.null(C)) C <- droplevels(C[keep, , drop = FALSE])
  M[!is.finite(M)] <- NA_real_
  if (length(unique(X)) < 2L) stop("X must vary among complete observations.", call. = FALSE)
  if (family_y == "binomial" && !setequal(unique(Y), c(0, 1))) stop("Binary Y must contain both 0 and 1, with no other values.", call. = FALSE)
  if (family_y == "cox" && sum(Y[, 2L]) < 1L) stop("At least one observed survival event is required.", call. = FALSE)
  list(X = X, M = M, Y = Y, C = C, n_input = n, row_index = which(keep))
}

.mopair_design <- function(X, C = NULL) {
  data <- data.frame(X = X)
  if (!is.null(C)) data <- cbind(data, C)
  stats::model.matrix(~ ., data = data)
}

.mopair_safe_term <- function(fit, term) {
  empty <- c(estimate = NA_real_, se = NA_real_, p = NA_real_)
  if (is.null(fit)) return(empty)
  sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(sm) || !term %in% rownames(sm)) return(empty)
  pc <- grep("^Pr\\(", colnames(sm), value = TRUE)[1L]
  se <- if ("Std. Error" %in% colnames(sm)) "Std. Error" else "se(coef)"
  ec <- if ("Estimate" %in% colnames(sm)) "Estimate" else "coef"
  out <- c(estimate = unname(sm[term, ec]), se = unname(sm[term, se]), p = if (is.na(pc)) NA_real_ else unname(sm[term, pc]))
  if (!all(is.finite(out)) || out["se"] <= 0) return(empty)
  out
}

.mopair_acat <- function(p) {
  p <- .mopair_clip_p(p)
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  .mopair_clip_p(0.5 - atan(mean(tan((0.5 - p) * pi))) / pi)
}

.mopair_simes <- function(p) {
  p <- sort(.mopair_clip_p(p))
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  min(1, min(length(p) * p / seq_along(p)))
}

.mopair_composite <- function(a_p, b_p, method) {
  a_p <- .mopair_clip_p(a_p); b_p <- .mopair_clip_p(b_p)
  ok <- is.finite(a_p) & is.finite(b_p)
  out <- rep(NA_real_, length(a_p))
  if (!any(ok)) return(out)
  if (method == "maxp") {
    out[ok] <- pmax(a_p[ok], b_p[ok])
  } else {
    if (!requireNamespace("DACT", quietly = TRUE)) stop("Install DACT to use composite_method = 'dact'; no alternative is substituted automatically.", call. = FALSE)
    value <- tryCatch(DACT::DACT(a_p[ok], b_p[ok], correction = "none"), error = function(e) stop("DACT failed: ", conditionMessage(e), call. = FALSE))
    if (length(value) != sum(ok)) stop("DACT returned an unexpected number of values.", call. = FALSE)
    out[ok] <- as.numeric(value)
  }
  .mopair_clip_p(out)
}
