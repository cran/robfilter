\encoding{latin1}

\name{const}

\docType{data}

\alias{const}

\title{Correction factors to achieve unbiasedness of the Qn scale estimator}

\description{
This matrix contains correction factors for the univariate Qn scale estimator (Rousseeuw, Croux, 1993) to achieve unbiasedness under Gaussian noise. The \code{madore.filter} estimates the local error covariance matrix by the orthogonalized Gnanadesikan-Kettenring estimator (Gnanadesikan, Kettenring, 1972, Maronna, Zamar, 2002) which is based on the Qn scale estimator.
}

\usage{const}

\format{
A (96x2)-matrix containing the correction factors for the univariate Qn scale estimator for the samples sizes \eqn{n=10,11,...,100,200,300,400,500,1000}.
}

\source{The correction factors have been obtained by simulations.}

\references{
Gnanadesikan, R., Kettenring, J.R. (1972)
Robust Estimates, Residuals, and Outlier Detection with Multiresponse Data,
\emph{Biometrics} \bold{28}, 81-124.\cr

Maronna, R.A., Zamar, R.H. (2002)
Robust Estimates of Location and Dispersion for High-Dimensional Datasets,
\emph{Technometrics} \bold{44}, 307-317.\cr

Rousseeuw, P.J., Croux, C. (1993)
Alternatives to the Median Absolute Deviation,
\emph{Journal of the American Statistical Association} \bold{88}, 1273-1283.
}

\keyword{datasets}
