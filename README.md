# MoPAIR

**Module-level Paired-path Association Inference for Related mediators**

MoPAIR organizes correlated molecular features into outcome-free modules,
combines exposure and outcome path evidence within each member, aggregates
member-level evidence using ACAT, and selects module hypotheses using BH.

The main analysis is **uncorrected DACT → ACAT → BH**.

## Installation

```r
install.packages("remotes")
remotes::install_github("zhonghualiu/DACT", upgrade = "never")
remotes::install_github("biostatYu/MoPAIR", upgrade = "never")
library(MoPAIR)
```

DACT is required for the main analysis. An explicit `composite_method = "maxp"`
option provides a different comparison rule; it is never substituted for DACT
automatically. The package is developed for R 4.2 or later.

## Continuous outcomes

```r
dat <- mopair_example_data(family_y = "gaussian", seed = 2026)
fit <- fit_mopair(
  X = dat$X, M = dat$M, Y = dat$Y, C = dat$C,
  family_y = "gaussian",
  composite_method = "dact", combination = "acat",
  calibration = "none", alpha = 0.10,
  corr_cut = 0.45, hclust_method = "average", seed = 2026
)
fit$selected_modules
head(fit$mediator_results)
plot_mopair(fit, type = "modules")
```

`X` is a numeric exposure vector; `M` is a sample-by-mediator numeric matrix;
`Y` is the outcome; and optional `C` contains numeric, logical, or factor
covariates. Inputs must share sample order. Apply study-specific transformations
and quality control before fitting; the fitting function does not standardize
or impute measurements automatically.

## Binary outcomes: efficient score test

```r
dat <- mopair_example_data(family_y = "binomial", seed = 2026)
fit_binary <- fit_mopair(
  dat$X, dat$M, dat$Y, dat$C,
  family_y = "binomial", binomial_test = "score",
  modules = dat$modules,
  composite_method = "dact", combination = "acat", alpha = 0.10
)
```

Binary `Y` must contain 0 and 1. The default `binomial_test = "score"` fits the
reduced exposure- and covariate-adjusted logistic model, then uses each
mediator's score and nuisance-adjusted efficient information. It matches the
outcome-path calculation used for the PRISM score analysis.

In this mode, `b` is the one-step score proxy `U / I`, `b_se` is `1 / sqrt(I)`,
and `b_p` is the two-sided score-test p value. The output identifies this
estimate as `b_estimator = "score_one_step"`; it is not the full-model logistic
maximum-likelihood coefficient. The `direct` field in score mode is the reduced
model's exposure coefficient.

Use `binomial_test = "wald"` explicitly for earlier analyses based on
single-mediator logistic Wald tests. The change of test can change discoveries.

## Right-censored survival outcomes

```r
dat <- mopair_example_data(family_y = "cox", seed = 2026)
fit_survival <- fit_mopair(
  dat$X, dat$M, dat$Y, dat$C,
  family_y = "cox", cox_ties = "efron",
  modules = dat$modules,
  composite_method = "dact", combination = "acat", alpha = 0.10
)
head(fit_survival$mediator_results[
  c("mediator", "b", "b_hr", "b_hr_lower", "b_hr_upper", "b_p")
])
```

For your data, construct `Y = survival::Surv(time, status)` with `status = 1`
for the event of interest and `status = 0` for censoring. The a path remains a
linear regression. Each b path uses a separate Cox model adjusting for exposure
and covariates, with a two-sided Wald test. `b` is on the log-hazard scale;
`b_hr` and its nominal 95% limits are member-level regression summaries.

The interface supports right-censored outcomes. Counting-process, interval-
censored, and multistate responses are rejected. For a cause-specific analysis,
prepare event and censoring times for the chosen endpoint before fitting.
Censoring-aware permutation calibration is not implemented; Cox analyses use BH.

## Fixed modules and multiple exposures

With `modules = NULL`, MoPAIR learns modules from mediator residuals adjusted
for `X` and `C`, without `Y`. It clusters distance `1 - abs(correlation)` using
average linkage. `corr_cut = 0.45` corresponds to tree height `0.55`. Small
clusters become singleton modules rather than being discarded.

Supply a named module vector to reuse a prespecified or independently learned
partition. Names are matched to `colnames(M)`; unnamed labels follow column
order. A multiple-exposure analysis requires one common fixed partition.

```r
dat <- mopair_example_data(seed = 2026)
X_matrix <- cbind(exposure_1 = dat$X, exposure_2 = dat$C$age_z)
grid <- fit_mopair_grid(
  X_matrix, dat$M, dat$Y, C = NULL,
  modules = dat$modules, family_y = "gaussian",
  composite_method = "dact", combination = "acat", alpha = 0.10
)
grid$selected_hypotheses
```

`family_q` and `selected_family` apply BH across the entire exposure–module
family. The `module_q` and `selected_within_exposure` columns describe the
separate within-exposure calculations. PRISM-style analyses use a common
outcome-free metabolite partition, then family-wide BH across all genera.

## Results and settings

* `module_results`: module p values, BH q values, selection flags, total members,
  and members with available evidence.
* `mediator_results`: a/b coefficients, standard errors, p values, DACT statistic,
  test and scale labels, sample counts, fit status, and module membership.
* `selected_modules`: the primary BH discoveries.
* `leading_edge`: descriptive representatives ranked within selected modules.
* `metadata`: settings, complete-case counts, excluded rows, and availability.

Rows missing exposure, outcome, or required covariates are excluded jointly.
Mediator missingness is handled per feature using matched a/b samples. Failed
fits remain unavailable and are not converted to evidence or zero p values.

Set `alpha = 0.05` for BH q ≤ 0.05. `calibration = "none"` retains the primary
BH analysis. The optional `dual_permutation` sensitivity interface supports
Gaussian and binary outcomes with explicit strata and `n_perm`; its validity
depends on the study's permutation assumptions and it is not the primary rule.
The grid interface does not apply separate per-exposure empirical calibrations.

MoPAIR reports module-level paired-path associations. Leading members and the
`a * b` product are descriptive summaries; they do not establish individual
mediator significance or a causal indirect effect. Calibration against
structural simulation truth can differ from calibration for fitted marginal
path hypotheses, particularly with correlated proxy signals.

## Migration from MOSAICCFDR 0.1.0

| Previous name | Current name |
|---|---|
| `MOSAICCFDR` | `MoPAIR` |
| `fit_mosaic()` | `fit_mopair()` |
| `fit_mosaic_grid()` | `fit_mopair_grid()` |
| `plot_mosaic()` | `plot_mopair()` |
| `mosaic_example_data()` | `mopair_example_data()` |
| `mosaic_fit` | `mopair_fit` |

Continuous path models and the uncorrected DACT–ACAT–BH sequence are retained.
The binary default is now score-based; select `"wald"` to recover the previous
test choice. Historical study outputs are not regenerated by this update.

## Tests

```sh
R CMD build .
R CMD check --no-manual MoPAIR_0.2.0.tar.gz
```

Tests compare Gaussian paths with `lm`, binary score p values with the Rao test
from `anova.glm`, and Cox estimates, p values, and confidence limits with
`survival::coxph`. They also cover fixed-module alignment, family-wide BH,
missingness, unsupported outcomes, and unavailable fits. All bundled examples
are synthetic; participant-level study data are not distributed in this package.
