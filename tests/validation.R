source_dir <- Sys.getenv("MOPAIR_TEST_SOURCE")
if (nzchar(source_dir)) {
  for (f in list.files(file.path(source_dir, "R"), "[.]R$", full.names = TRUE)) sys.source(f, globalenv())
} else {
  library(MoPAIR)
}
stopifnot(exists("fit_mopair", mode = "function"), exists("fit_mopair_grid", mode = "function"))

check <- function(label, code) {
  force(code)
  cat("PASS:", label, "\n")
}
near <- function(x, y, tolerance = 1e-9) stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = tolerance)))
fails <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}

set.seed(731)
n <- 180L
p <- 20L
X <- rbinom(n, 1, 0.45)
C <- data.frame(age = rnorm(n), group = factor(rep(c("A", "B", "C"), length.out = n)))
M <- matrix(rnorm(n * p), n, p)
M[, 1:3] <- M[, 1:3] + 0.5 * X
M[, 5] <- 0.7 * M[, 4] + 0.4 * rnorm(n)
colnames(M) <- paste0("feature_", seq_len(p))
modules <- setNames(rep(c("G1", "G2"), each = 10L), colnames(M))
Y <- 0.3 * X + 0.7 * M[, 1] + 0.2 * C$age + rnorm(n)
Y_binary <- rbinom(n, 1, plogis(-0.4 + 0.3 * X + 0.9 * M[, 1] + 0.3 * C$age))
args <- list(X = X, M = M, C = C, modules = modules, composite_method = "maxp", alpha = 0.10)
fit <- do.call(fit_mopair, c(args, list(Y = Y)))

check("Gaussian path coefficients and t tests agree with lm", {
  ref <- summary(lm(Y ~ X + M[, 1] + C$age + C$group))$coefficients[3, ]
  row <- fit$mediator_results[match("feature_1", fit$mediator_results$mediator), ]
  near(row$b, ref[1]); near(row$b_se, ref[2]); near(row$b_p, ref[4])
  near(fit$module_results$module_q, p.adjust(fit$module_results$module_p, "BH"))
  stopifnot(inherits(fit, "mopair_fit"), fit$metadata$method == "MoPAIR")
})

binary <- do.call(fit_mopair, c(args, list(Y = Y_binary, family_y = "binomial")))
check("default binary test agrees with the independent glm Rao score test", {
  null <- glm(Y_binary ~ X + C$age + C$group, family = binomial())
  full <- glm(Y_binary ~ X + C$age + C$group + M[, 1], family = binomial())
  ref_p <- anova(null, full, test = "Rao")$`Pr(>Chi)`[2]
  row <- binary$mediator_results[match("feature_1", binary$mediator_results$mediator), ]
  near(row$b_p, ref_p, 1e-7)
  stopifnot(binary$metadata$binomial_test == "score", row$b_estimator == "score_one_step", row$b_test == "score")
})

wald <- do.call(fit_mopair, c(args, list(Y = Y_binary, family_y = "binomial", binomial_test = "wald")))
check("explicit logistic Wald option remains available", {
  ref <- summary(glm(Y_binary ~ X + M[, 1] + C$age + C$group, family = binomial()))$coefficients[3, ]
  row <- wald$mediator_results[match("feature_1", wald$mediator_results$mediator), ]
  near(row$b, ref[1]); near(row$b_se, ref[2]); near(row$b_p, ref[4])
  stopifnot(wald$metadata$binomial_test == "wald")
})

set.seed(732)
event_time <- rexp(n, exp(0.3 * X + 0.5 * M[, 1] + 0.2 * C$age))
censor_time <- rexp(n, 0.4)
time <- pmax(0.05, round(pmin(event_time, censor_time), 1))
status <- as.integer(event_time <= censor_time)
Y_surv <- survival::Surv(time, status)
cox <- do.call(fit_mopair, c(args, list(Y = Y_surv, family_y = "cox")))
check("Cox coefficient, Wald p, HR and confidence limits agree with survival", {
  ref <- summary(survival::coxph(Y_surv ~ X + M[, 1] + C$age + C$group, ties = "efron"))
  row <- cox$mediator_results[match("feature_1", cox$mediator_results$mediator), ]
  near(row$b, ref$coefficients[2, "coef"])
  near(row$b_se, ref$coefficients[2, "se(coef)"])
  near(row$b_p, ref$coefficients[2, "Pr(>|z|)"])
  near(row$b_hr, ref$conf.int[2, "exp(coef)"])
  near(row$b_hr_lower, ref$conf.int[2, "lower .95"])
  near(row$b_hr_upper, ref$conf.int[2, "upper .95"])
  stopifnot(row$b_scale == "log_hazard", cox$metadata$cox_ties == "efron")
})

check("Cox tie method is configurable", {
  alt <- do.call(fit_mopair, c(args, list(Y = Y_surv, family_y = "cox", cox_ties = "breslow")))
  ref <- summary(survival::coxph(Y_surv ~ X + M[, 1] + C$age + C$group, ties = "breslow"))
  row <- alt$mediator_results[match("feature_1", alt$mediator_results$mediator), ]
  near(row$b_p, ref$coefficients[2, "Pr(>|z|)"])
})

check("named fixed modules align by mediator identifier", {
  reordered <- do.call(fit_mopair, c(within(args, modules <- rev(modules)), list(Y = Y)))
  stopifnot(identical(fit$modules, reordered$modules))
  near(fit$module_results$module_p, reordered$module_results$module_p)
})

check("multiple exposures use one family-wide BH adjustment", {
  g <- fit_mopair_grid(cbind(X1 = X, X2 = C$age), M, Y_binary, C = C["group"],
    modules = modules, family_y = "binomial", composite_method = "maxp", alpha = 0.05)
  near(g$family_results$family_q, p.adjust(g$family_results$module_p, "BH"))
  stopifnot(nrow(g$family_results) == 4, all(g$family_results$selected_family == (g$family_results$family_q <= 0.05)))
  fails(fit_mopair_grid(cbind(X), M, Y, composite_method = "maxp"), "fixed modules")
})

check("missing mediator observations use matched a and b samples", {
  mm <- M; mm[1:12, 1] <- NA_real_
  cc <- C; cc$age[15] <- NA_real_
  yy <- Y_binary; yy[16] <- NA_real_
  z <- fit_mopair(X, mm, yy, cc, family_y = "binomial", modules = modules, composite_method = "maxp")
  ok <- complete.cases(X, mm[, 1], yy, cc)
  null <- glm(yy[ok] ~ X[ok] + cc$age[ok] + cc$group[ok], family = binomial())
  full <- glm(yy[ok] ~ X[ok] + cc$age[ok] + cc$group[ok] + mm[ok, 1], family = binomial())
  row <- z$mediator_results[match("feature_1", z$mediator_results$mediator), ]
  near(row$b_p, anova(null, full, test = "Rao")$`Pr(>Chi)`[2], 1e-7)
  stopifnot(row$n_used == sum(ok), z$metadata$n_complete == n - 2L)
})

check("invalid inputs and unsupported survival schemes fail explicitly", {
  fails(fit_mopair(X[-1], M, Y, composite_method = "maxp"), "same")
  fails(fit_mopair(X, M, rep(1, n), family_y = "binomial", composite_method = "maxp"), "0 and 1")
  fails(fit_mopair(X, M, cbind(time, status), family_y = "cox", composite_method = "maxp"), "Surv")
  fails(fit_mopair(X, M, survival::Surv(time, rep(0, n)), family_y = "cox", composite_method = "maxp"), "event")
  fails(fit_mopair(X, M, survival::Surv(rep(0, n), time, status), family_y = "cox", composite_method = "maxp"), "right")
  fails(fit_mopair(X, M, Y_surv, family_y = "cox", calibration = "dual_permutation", composite_method = "maxp"), "censor")
  fails(fit_mopair(X, M, Y, alpha = 1, composite_method = "maxp"), "alpha")
})

check("unavailable member evidence is retained as missing", {
  mm <- M; mm[, 1] <- 0; mm[, 11:20] <- NA_real_
  z <- fit_mopair(X, mm, Y, C, modules = modules, composite_method = "maxp")
  stopifnot(is.na(z$mediator_results$b_p[z$mediator_results$mediator == "feature_1"]),
    is.na(z$module_results$module_p[z$module_results$module == "G2"]),
    !z$module_results$selected_bh[z$module_results$module == "G2"])
})

check("a feature observed in only one binary outcome class has unavailable outcome evidence", {
  mm <- M
  mm[Y_binary == 1, 1] <- NA_real_
  z <- fit_mopair(X, mm, Y_binary, C, modules = modules,
    family_y = "binomial", composite_method = "maxp")
  stopifnot(is.na(z$mediator_results$b_p[z$mediator_results$mediator == "feature_1"]))
})

check("Gaussian, binary and Cox examples are available", {
  for (family in c("gaussian", "binomial", "cox")) {
    dat <- mopair_example_data(n = 100, p = 20, family_y = family, seed = 24)
    z <- fit_mopair(dat$X, dat$M, dat$Y, dat$C, family_y = family, modules = dat$modules, composite_method = "maxp")
    stopifnot(inherits(z, "mopair_fit"), nrow(z$mediator_results) == 20)
  }
})

if (requireNamespace("DACT", quietly = TRUE)) check("DACT uses the uncorrected external implementation", {
  z <- fit_mopair(X, M, Y, C, modules = modules, composite_method = "dact")
  members <- z$mediator_results[match(colnames(M), z$mediator_results$mediator), ]
  expected <- DACT::DACT(members$a_p, members$b_p, correction = "none")
  near(members$composite_p, expected, 1e-10)
})

check("plots and print methods run without altering a fit", {
  outfile <- tempfile(fileext = ".pdf")
  grDevices::pdf(outfile)
  plot_mopair(fit, "modules"); plot_mopair(fit, "leading_edge")
  plot(fit)
  grDevices::dev.off()
  stopifnot(file.info(outfile)$size > 0, inherits(print(fit), "mopair_fit"))
})
