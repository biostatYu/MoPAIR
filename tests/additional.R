source_dir <- Sys.getenv("MOPAIR_TEST_SOURCE")
if (nzchar(source_dir)) {
  for (f in list.files(file.path(source_dir, "R"), "[.]R$", full.names = TRUE)) sys.source(f, globalenv())
} else library(MoPAIR)

dat <- mopair_example_data(n = 90, p = 20, seed = 510)
base <- list(X = dat$X, M = dat$M, Y = dat$Y, C = dat$C, modules = dat$modules, composite_method = "maxp")
fit <- do.call(fit_mopair, base)
alt <- do.call(fit_mopair, c(base, list(combination = "simes")))
for (g in alt$module_results$module) {
  pp <- sort(alt$mediator_results$composite_p[alt$mediator_results$module == g])
  expected <- min(1, min(length(pp) * pp / seq_along(pp)))
  stopifnot(abs(expected - alt$module_results$module_p[alt$module_results$module == g]) < 1e-12)
}
cat("PASS: Simes comparison matches its definition\n")

learned <- fit_mopair(dat$X, dat$M, dat$Y, dat$C, composite_method = "maxp")
changed_y <- fit_mopair(dat$X, dat$M, rev(dat$Y), dat$C, composite_method = "maxp")
stopifnot(identical(learned$modules, changed_y$modules))
singleton <- fit_mopair(dat$X, dat$M[, 1, drop = FALSE], dat$Y, dat$C, composite_method = "maxp")
stopifnot(nrow(singleton$module_results) == 1)
cat("PASS: module learning does not use outcomes and supports one feature\n")

cal <- suppressWarnings(do.call(fit_mopair, c(base, list(calibration = "dual_permutation", n_perm = 2L, seed = 20))))
stopifnot(identical(fit$module_results$selected_bh, cal$module_results$selected_bh),
  all(is.finite(cal$calibration$curve$empirical_fdr)), isTRUE(cal$calibration$sensitivity_only))
single_strata <- suppressWarnings(do.call(fit_mopair, c(base,
  list(calibration = "dual_permutation", n_perm = 2L, strata = seq_along(dat$X), seed = 20))))
expected <- single_strata$module_results$module_p
stopifnot(max(abs(single_strata$calibration$null_h10_module_p[1, ] - expected)) < 1e-12,
  max(abs(single_strata$calibration$null_h01_module_p[1, ] - expected)) < 1e-12)
cat("PASS: sensitivity calibration preserves primary BH and singleton strata\n")

outfile <- tempfile(fileext = ".pdf")
grDevices::pdf(outfile)
plot_mopair(cal, "calibration")
missing <- dat$M; missing[,] <- NA_real_
empty <- fit_mopair(dat$X, missing, dat$Y, dat$C, modules = dat$modules, composite_method = "maxp")
plot_mopair(empty, "modules"); plot_mopair(empty, "leading_edge")
grDevices::dev.off()
stopifnot(file.info(outfile)$size > 0)
cat("PASS: sensitivity and unavailable-evidence plotting\n")

expect_error <- function(expr) stopifnot(inherits(tryCatch({ force(expr); NULL }, error = identity), "error"))
expect_error(fit_mopair(as.character(dat$X), dat$M, dat$Y))
expect_error(fit_mopair(dat$X, matrix("x", 90, 20), dat$Y))
expect_error(fit_mopair(dat$X, dat$M, letters[1:90]))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, dat$C[-1, , drop = FALSE]))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, C = data.frame(text = rep("x", 90))))
expect_error(fit_mopair(rep(1, 90), dat$M, dat$Y))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, alpha = NA_real_))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, leading_edge_n = -1))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, seed = Inf))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, corr_cut = 1))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, min_module_size = 0))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, modules = rep("G", 2)))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, modules = setNames(rep("G", 20), rep("same", 20))))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, modules = setNames(rep("G", 20), paste0("wrong", 1:20))))
expect_error(fit_mopair(dat$X, dat$M, dat$Y, modules = rep("", 20)))
expect_error(do.call(fit_mopair, c(base, list(calibration = "dual_permutation", n_perm = 0))))
expect_error(do.call(fit_mopair, c(base, list(calibration = "dual_permutation", strata = 1:2))))
expect_error(plot_mopair(list()))
expect_error(plot_mopair(fit, top = 0))
expect_error(plot_mopair(fit, "calibration"))
expect_error(mopair_example_data(n = 2))
expect_error(mopair_example_data(p = 5))
expect_error(fit_mopair_grid(cbind(dat$X), dat$M, dat$Y, modules = dat$modules, calibration = "dual_permutation"))
cat("PASS: public argument validation\n")
