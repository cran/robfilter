\name{var.n}
\alias{var.n}
\docType{data}

\title{
Variance of the Repeated Median slope estimator.
}

\description{
Empirical variance of the Repeated Median slope estimator, computed on standard normal noise in moving windows of width n; required by the function \code{scarm.filter}.
}

\usage{data(var.n)}
\format{
  The format is:
 num [1:300] NA NA NA NA 0.138 ...
}

\details{
The variance of the Repeated Median slope estimator depends on the size of the window sample, i.e. the window width n. These are results from simulations, where the Repeated Median slope has been computed on standard normal noise in moving windows of width n. The value \code{var.n[n]} is the variance for the window width n.
}

\source{
The empirical variances have been obtained by simulations.
}

\references{
Borowski, M. and Fried, R. (2011)
Robust moving window regression for online signal extraction from
non-stationary time series: online window width adaption by testing for
signal changes, \emph{submitted}.
}

\keyword{datasets}
