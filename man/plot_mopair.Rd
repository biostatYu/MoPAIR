\name{plot_mopair}
\alias{plot_mopair}
\alias{plot.mopair_fit}
\alias{print.mopair_fit}
\title{Display a MoPAIR Analysis}
\description{Display module evidence, descriptive leading members, or an optional sensitivity calibration.}
\usage{
plot_mopair(x, type = c("modules", "calibration", "leading_edge"),
  top = 20L, ...)
\method{plot}{mopair_fit}(x, ...)
\method{print}{mopair_fit}(x, ...)
}
\arguments{
\item{x}{A \code{mopair_fit} object.}
\item{type}{Display type.}
\item{top}{Maximum number of displayed modules or members.}
\item{...}{Additional plotting arguments.}
}
\value{The fit, invisibly.}
\examples{
dat <- mopair_example_data(n = 100, p = 20)
fit <- fit_mopair(dat$X, dat$M, dat$Y, dat$C, composite_method = "maxp")
plot_mopair(fit)
}
