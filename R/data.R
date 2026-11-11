mopair_example_data <- function(n = 240L, p = 60L, family_y = c("gaussian", "binomial", "cox"), seed = 2026L) {
  family_y <- match.arg(family_y)
  n <- as.integer(n); p <- as.integer(p)
  if (length(n) != 1L || is.na(n) || n < 20L) stop("n must be at least 20.", call. = FALSE)
  if (length(p) != 1L || is.na(p) || p < 20L || p %% 10L != 0L) stop("p must be at least 20 and divisible by 10.", call. = FALSE)
  set.seed(seed)
  X <- stats::rbinom(n, 1, 0.4)
  C <- data.frame(age_z = stats::rnorm(n))
  module <- rep(paste0("block_", seq_len(p / 10L)), each = 10L)
  latent <- matrix(stats::rnorm(n * (p / 10L)), n, p / 10L)
  M <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) M[, j] <- 0.75 * latent[, ceiling(j / 10L)] + 0.45 * stats::rnorm(n)
  active <- seq_len(min(5L, p))
  M[, active] <- M[, active] + 0.45 * X
  linear <- 0.10 * X + 0.35 * rowMeans(M[, active, drop = FALSE]) + 0.15 * C$age_z
  Y <- if (family_y == "gaussian") linear + stats::rnorm(n) else if (family_y == "binomial") {
    stats::rbinom(n, 1, stats::plogis(linear - 0.5))
  } else {
    event <- stats::rexp(n, rate = exp(linear))
    censor <- stats::rexp(n, rate = 0.4)
    survival::Surv(pmin(event, censor), as.integer(event <= censor))
  }
  colnames(M) <- paste0("metabolite_", seq_len(p))
  list(X = X, M = scale(M), Y = Y, C = C, modules = stats::setNames(module, colnames(M)))
}
