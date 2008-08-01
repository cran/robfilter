\name{critvals}
\docType{data}
\alias{critvals}
\title{Critical Values for the RM Goodness of Fit Test}
\description{
  This matrix contains critical values for the goodness of fit test
  for the last fitted value of a Repeated Median regression fit
  to a sample of size \code{n}. The critical values are based on the 0.95-quantiles
  of the distribution of a test statistic corresponding to the absolute value of
  the sum of a subset of residual signs. The critical value for a test based on the
  last \code{nI} out of \code{n} observations corresponds to \code{critvals[n,nI]}.
}
\usage{critvals}
\format{A (600x61)-matrix containing 30550 observations.}
\source{Simulation.}
\references{
Schettlinger, K., Fried, R., Gather, U. (2008) 
Real Time Signal Processing by Adaptive Repeated Median Filters,
\emph{International Journal of Adaptive Control and Signal Processing}, submitted.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians,
\emph{Biometrika} \bold{69} (1), 
242-244.\cr
}
\keyword{datasets}
