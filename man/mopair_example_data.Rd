\name{mopair_example_data}
\alias{mopair_example_data}
\title{Generate Synthetic MoPAIR Example Data}
\description{Generate correlated mediator blocks with continuous, binary, or right-censored survival outcomes for software examples.}
\usage{mopair_example_data(n = 240L, p = 60L,
  family_y = c("gaussian", "binomial", "cox"), seed = 2026L)}
\arguments{
\item{n}{Sample size, at least 20.}
\item{p}{Mediator count, at least 20 and divisible by 10.}
\item{family_y}{Outcome family.}
\item{seed}{Random seed.}
}
\value{A list with \code{X}, standardized \code{M}, \code{Y}, \code{C}, and known example \code{modules}. These are synthetic demonstrations rather than the manuscript's simulation benchmark.}
\examples{
dat <- mopair_example_data(n = 100, p = 20, family_y = "cox")
head(dat$Y)
}
