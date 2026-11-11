library(MoPAIR)

# Primary DACT analysis; install DACT before running.
for (family in c("gaussian", "binomial", "cox")) {
  dat <- mopair_example_data(n = 240, p = 60, family_y = family, seed = 2026)
  fit <- fit_mopair(dat$X, dat$M, dat$Y, dat$C,
    family_y = family, modules = dat$modules,
    composite_method = "dact", combination = "acat",
    calibration = "none", alpha = 0.10,
    binomial_test = "score", cox_ties = "efron")
  print(fit)
  print(fit$selected_modules)
}
