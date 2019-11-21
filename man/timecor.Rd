\encoding{latin1}

\name{timecorrection}

\docType{data}

\alias{timecorrection}

\title{Correction factors for the scale estimation of the filtering procedure proposed by Fried (2004).}

\description{
Fried's (2004) signal extraction procedure includes optional rules for outlier replacement based on local scale estimation. Since detected outliers are treated as missing values, the finite sample correction for the scale estimation is adjusted for the reduced sample size, using the correction factors in the dataset 'timecorrection'.
}

\usage{timecorrection}

\format{
A (250x16)-matrix containing the correction factors for the scale estimators MAD, Qn, Sn, and LSH and for the outlier treatments 'trimming', 'downsizing large values', 'downsizing moderate values', and 'winsorization'.
}

\source{The correction factors have been obtained by simulations.}

\references{
Fried, R. (2004), Robust Filtering of Time Series with 
Trends, \emph{Journal of Nonparametric Statistics} \bold{16}, 
313-328.\cr
}

\seealso{\code{\link{robust.filter}}.}

\keyword{datasets}
