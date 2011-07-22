\name{dfs}
\alias{dfs}
\docType{data}

\title{
Degrees of freedom for the SCARM test statistic.
}

\description{
This matrix contains degrees of freedom for the t-distributed SCARM test statistic; required by the function \code{scarm.filter}.
}

\usage{data(dfs)}

\format{
  A data frame with 20 observations on the following 20 variables.
  \describe{
    \item{\code{X5}}{a numeric vector}
    \item{\code{X10}}{a numeric vector}
    \item{\code{X15}}{a numeric vector}
    \item{\code{X20}}{a numeric vector}
    \item{\code{X25}}{a numeric vector}
    \item{\code{X30}}{a numeric vector}
    \item{\code{X35}}{a numeric vector}
    \item{\code{X40}}{a numeric vector}
    \item{\code{X45}}{a numeric vector}
    \item{\code{X50}}{a numeric vector}
    \item{\code{X55}}{a numeric vector}
    \item{\code{X60}}{a numeric vector}
    \item{\code{X65}}{a numeric vector}
    \item{\code{X70}}{a numeric vector}
    \item{\code{X75}}{a numeric vector}
    \item{\code{X80}}{a numeric vector}
    \item{\code{X85}}{a numeric vector}
    \item{\code{X90}}{a numeric vector}
    \item{\code{X95}}{a numeric vector}
    \item{\code{X100}}{a numeric vector}
  }
}

\details{
The SCARM test from the function \code{scarm.filter} is based on the difference of Repeated Median slopes computed in a left-hand and right-hand window. The distribution of the SCARM test statistic is approximated by a t-distribution where the degrees of freedom depend on the width of the left- and right-hand window. This matrix delivers suitable degrees of freedom, obtained by simulations.
}

\source{
The degrees of freedom have been obtained by simulations.
}

\references{
Borowski, M. and Fried, R. (2011)
Robust moving window regression for online signal extraction from
non-stationary time series: online window width adaption by testing for
signal changes, \emph{submitted}.
}

\keyword{datasets}
