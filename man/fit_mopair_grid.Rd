\name{fit_mopair_grid}
\alias{fit_mopair_grid}
\title{Fit a Common Family of Exposure-Module Hypotheses}
\description{Fit each exposure using a common fixed module partition and apply BH to the entire exposure-module family.}
\usage{fit_mopair_grid(X, M, Y, ...)}
\arguments{
\item{X}{Numeric sample-by-exposure matrix with unique column names.}
\item{M}{Numeric sample-by-mediator matrix.}
\item{Y}{Outcome vector or right-censored Surv object.}
\item{...}{Arguments passed to \code{fit_mopair}, including required common fixed \code{modules}, optional \code{C}, \code{family_y}, and the testing settings.}
}
\details{Final selection uses \code{family_q} and \code{selected_family}. The columns \code{module_q} and \code{selected_within_exposure} retain within-exposure summaries. The interface supports BH selection; empirical calibration of a complete multi-exposure family requires a separate study-specific driver.}
\value{A list with exposure-specific \code{fits}, \code{family_results}, \code{selected_hypotheses}, and \code{metadata}.}
\seealso{\code{\link{fit_mopair}}}
\examples{
dat <- mopair_example_data(n = 100, p = 20)
g <- fit_mopair_grid(cbind(X1 = dat$X, X2 = dat$C$age_z),
  dat$M, dat$Y, modules = dat$modules, composite_method = "maxp")
g$family_results
}
