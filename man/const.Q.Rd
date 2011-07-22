\name{const.Q}
\alias{const.Q}
\docType{data}

\title{
Correction factors to achieve unbiasedness of the regression-free Q scale estimator
}

\description{
Correction factors for the regression-free Q scale estimator (Rousseeuw and Hubert, 1996, Gelper et al., 2009) to achieve unbiasedness under Gaussian noise; required by the function \code{scarm.filter}.
}

\usage{data(const.Q)}

\format{
  The format is:
 num [1:151] NA NA NA NA 2.17 ...
}

\source{
The correction factors have been obtained by simulations.
}

\references{
Rousseeuw, P. and Hubert, M. (1996)
Regression-free and robust estimation of scale for bivariate data,
\emph{Computational Statistics and Data Analysis}, \bold{21}(1), 67-85.\cr

Gelper, S., Schettlinger, K., Croux, C., and Gather, U. (2009)
Robust online scale estimation in time series: A model-free approach,
\emph{Journal of Statistical Planning and Inference}, \bold{139}(2), 335-349.\cr
}

\keyword{datasets}
