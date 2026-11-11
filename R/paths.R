.mopair_score_null <- function(X, Y, C) {
  if (!setequal(unique(Y), c(0, 1))) return(NULL)
  design <- .mopair_design(X, C)
  if (nrow(design) <= ncol(design) || qr(design)$rank != ncol(design)) return(NULL)
  fit <- suppressWarnings(stats::glm.fit(x = design, y = Y, family = stats::binomial(), intercept = FALSE, control = stats::glm.control(maxit = 100L)))
  if (!isTRUE(fit$converged)) return(NULL)
  probability <- fit$fitted.values
  weight <- probability * (1 - probability)
  information <- crossprod(design, design * weight)
  inverse <- tryCatch(solve(information), error = function(e) NULL)
  if (is.null(inverse) || any(!is.finite(inverse))) return(NULL)
  direct <- unname(fit$coefficients["X"])
  direct_se <- sqrt(inverse["X", "X"])
  list(design = design, probability = probability, weight = weight, inverse = inverse,
    direct = c(estimate = direct, se = direct_se, p = 2 * stats::pnorm(-abs(direct / direct_se))))
}

.mopair_score_member <- function(mediator, Y, null) {
  if (is.null(null)) return(NULL)
  cross <- crossprod(null$design, mediator * null$weight)
  information <- sum(mediator^2 * null$weight) - as.numeric(crossprod(cross, null$inverse %*% cross))
  if (!is.finite(information) || information <= 1e-12) return(NULL)
  score <- sum((Y - null$probability) * mediator)
  stat <- score / sqrt(information)
  c(estimate = score / information, se = 1 / sqrt(information), p = 2 * stats::pnorm(-abs(stat)))
}

.mopair_fit_paths <- function(X, M, Y, C, family_y, binomial_test, cox_ties) {
  rows <- vector("list", ncol(M))
  score_mode <- family_y == "binomial" && binomial_test == "score"
  common_null <- if (score_mode) .mopair_score_null(X, Y, C) else NULL
  empty <- c(estimate = NA_real_, se = NA_real_, p = NA_real_)
  for (j in seq_len(ncol(M))) {
    ok <- is.finite(M[, j])
    xx <- X[ok]; mm <- M[ok, j]
    yy <- if (family_y == "cox") Y[ok, , drop = FALSE] else Y[ok]
    cc <- if (is.null(C)) NULL else droplevels(C[ok, , drop = FALSE])
    a <- b <- direct <- empty
    state <- "ok"
    events <- if (family_y == "cox") sum(yy[, 2L]) else NA_integer_
    n_used <- sum(ok)
    design <- tryCatch(.mopair_design(xx, cc), error = function(e) NULL)
    if (is.null(design) || n_used <= ncol(design) + 1L) {
      state <- "insufficient_complete_cases"
    } else if (length(unique(mm)) < 2L) {
      state <- "constant_mediator"
    } else {
      ad <- data.frame(M = mm, X = xx)
      bd <- data.frame(X = xx, M = mm)
      if (!is.null(cc)) { ad <- cbind(ad, cc); bd <- cbind(bd, cc) }
      bd$Y <- yy
      afit <- tryCatch(stats::lm(M ~ ., data = ad), error = function(e) NULL)
      a <- .mopair_safe_term(afit, "X")
      if (score_mode) {
        null <- if (all(ok)) common_null else .mopair_score_null(xx, yy, cc)
        value <- .mopair_score_member(mm, yy, null)
        if (is.null(value)) {
          state <- "score_information_unavailable"
        } else {
          b <- value; direct <- null$direct
        }
      } else {
        bad_fit <- FALSE
        bfit <- tryCatch(withCallingHandlers({
          if (family_y == "gaussian") stats::lm(Y ~ ., data = bd)
          else if (family_y == "binomial") stats::glm(Y ~ ., data = bd, family = stats::binomial(), control = stats::glm.control(maxit = 100L))
          else survival::coxph(Y ~ ., data = bd, ties = cox_ties, singular.ok = FALSE)
        }, warning = function(w) {
          if (grepl("converg|infinite|0 or 1", conditionMessage(w), ignore.case = TRUE)) bad_fit <<- TRUE
          invokeRestart("muffleWarning")
        }), error = function(e) NULL)
        if (family_y == "binomial" && (!setequal(unique(yy), c(0, 1)) || (!is.null(bfit) && !isTRUE(bfit$converged)))) bad_fit <- TRUE
        if (!bad_fit) { b <- .mopair_safe_term(bfit, "M"); direct <- .mopair_safe_term(bfit, "X") }
      }
      if (!is.finite(a["p"])) state <- "a_fit_unavailable"
      else if (!is.finite(b["p"]) && state == "ok") state <- "b_fit_unavailable"
    }
    b_scale <- switch(family_y, gaussian = "outcome", binomial = "log_odds", cox = "log_hazard")
    b_estimator <- if (score_mode) "score_one_step" else if (family_y == "cox") "partial_likelihood" else "regression"
    b_test <- if (score_mode) "score" else if (family_y == "gaussian") "t" else "wald"
    hr <- if (family_y == "cox" && is.finite(b["estimate"])) exp(b["estimate"]) else NA_real_
    rows[[j]] <- data.frame(mediator = colnames(M)[j], a = unname(a["estimate"]), a_se = unname(a["se"]), a_p = unname(a["p"]),
      b = unname(b["estimate"]), b_se = unname(b["se"]), b_p = unname(b["p"]),
      indirect_proxy = unname(a["estimate"] * b["estimate"]), direct = unname(direct["estimate"]), direct_p = unname(direct["p"]),
      b_test = b_test, b_estimator = b_estimator, b_scale = b_scale, b_hr = unname(hr),
      b_hr_lower = if (is.finite(hr)) unname(exp(b["estimate"] - stats::qnorm(0.975) * b["se"])) else NA_real_,
      b_hr_upper = if (is.finite(hr)) unname(exp(b["estimate"] + stats::qnorm(0.975) * b["se"])) else NA_real_,
      n_used = n_used, events = events, path_status = state, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
