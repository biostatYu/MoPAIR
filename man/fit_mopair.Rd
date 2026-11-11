\name{fit_mopair}
\alias{fit_mopair}
\title{Fit Module-Level Paired-Path Association Inference}
\description{Evaluate member-aligned exposure and outcome paths, aggregate evidence within outcome-free or fixed modules, and apply module-level BH selection.}
\usage{
fit_mopair(X, M, Y, C = NULL,
  family_y = c("gaussian", "binomial", "cox"),
  alpha = 0.10, modules = NULL, corr_cut = 0.45,
  min_module_size = 2L, hclust_method = "average",
  composite_method = c("dact", "maxp"),
  combination = c("acat", "simes"), leading_edge_n = 3L,
  binomial_test = c("score", "wald"),
  cox_ties = c("efron", "breslow", "exact"),
  calibration = c("none", "dual_permutation"),
  n_perm = 200L, strata = NULL, seed = 1L)
}
\arguments{
\item{X}{Numeric exposure vector.}
\item{M}{Numeric sample-by-mediator matrix with unique column names.}
\item{Y}{Numeric continuous or binary outcome, or a right-censored \code{survival::Surv(time, status)} object for Cox analysis. Binary outcomes use 0 and 1.}
\item{C}{Optional covariate matrix or data frame containing numeric, logical, or factor covariates.}
\item{family_y}{Outcome family. The default is Gaussian.}
\item{alpha}{Primary module BH q-value cutoff.}
\item{modules}{Optional fixed module labels. Named labels are matched to mediator identifiers; unnamed labels follow the columns of M.}
\item{corr_cut}{Absolute residual-correlation threshold; the clustering tree is cut at height \code{1 - corr_cut}.}
\item{min_module_size}{Clusters below this size become singleton modules.}
\item{hclust_method}{Linkage passed to \code{stats::hclust}.}
\item{composite_method}{Uncorrected DACT for the primary analysis, or an explicit MaxP comparison. DACT requires the DACT package.}
\item{combination}{Equal-weight ACAT, or the Simes comparison rule.}
\item{leading_edge_n}{Maximum number of descriptive representatives per selected module.}
\item{binomial_test}{Efficient logistic score test (default), or explicit full-model Wald test.}
\item{cox_ties}{Cox tie handling, passed to \code{survival::coxph}.}
\item{calibration}{Primary BH analysis without empirical calibration, or an optional Gaussian/binary dual-permutation sensitivity analysis. Cox permutation calibration is not supported.}
\item{n_perm}{Permutations per null generator in sensitivity calibration.}
\item{strata}{Optional complete permutation-stratum labels, of the original sample length.}
\item{seed}{Random seed.}
}
\details{
The main procedure is uncorrected DACT, equal-weight ACAT, and BH. Gaussian paths use linear regression t tests. Binary score tests use the reduced model containing exposure and covariates and nuisance-adjusted efficient information. In score mode, \code{b} is the one-step proxy \eqn{U/I}, \code{b_se} is \eqn{1/\sqrt{I}}, and \code{direct} is the exposure coefficient from the reduced model. The optional Wald mode uses a separate full logistic model for each mediator.

Cox paths use separate exposure- and covariate-adjusted proportional hazards models with two-sided Wald tests. The b coefficient is on the log-hazard scale; HRs and nominal 95 percent confidence limits are also returned. Survival inputs support right censoring only. Prepare cause-specific event and censoring times before fitting.

Exposure, outcome, and covariate missingness is excluded jointly. Mediator missingness uses feature-specific matched a/b samples. Unavailable fits remain missing. Measurements are not standardized or imputed automatically. The path product and leading members are descriptive and do not establish individual mediator significance or a causal indirect effect. FDR calibration for fitted marginal paths can differ from evaluation against structural mediation truth.
}
\value{An object of class \code{mopair_fit} containing \code{module_results}, \code{mediator_results}, \code{modules}, \code{selected_modules}, \code{leading_edge}, \code{calibration}, and settings/availability \code{metadata}. Primary BH discoveries remain in \code{selected_modules}; optional sensitivity flags are in \code{module_results$selected_calibrated}.}
\seealso{\code{\link{fit_mopair_grid}}, \code{\link{mopair_example_data}}, \code{\link{plot_mopair}}}
\examples{
dat <- mopair_example_data(n = 100, p = 20, seed = 7)
fit <- fit_mopair(dat$X, dat$M, dat$Y, dat$C,
  composite_method = "maxp", modules = dat$modules)
fit$module_results

binary <- mopair_example_data(n = 100, p = 20, family_y = "binomial")
fit_binary <- fit_mopair(binary$X, binary$M, binary$Y, binary$C,
  family_y = "binomial", composite_method = "maxp")

surv <- mopair_example_data(n = 100, p = 20, family_y = "cox")
fit_survival <- fit_mopair(surv$X, surv$M, surv$Y, surv$C,
  family_y = "cox", composite_method = "maxp")
}
