\encoding{latin1}

\name{madore.filter}

\alias{madore.filter}

\title{A Robust Adaptive Online-Filter for Multivariate Time Series}

\description{
Procedure for robust online extraction of low frequency components (the \emph{signal}) 
from a multivariate time series by a moving window technique with adaptive window width selection
(\emph{m}ultivariate \emph{ad}aptive \emph{o}nline \emph{re}peated median filter).
}

\usage{madore.filter(Y, byrow=FALSE, 
                 min.width=20, max.width=200, start.width=min.width, 
                 test.sample.size=min.width/2, width.search="geometric", 
                 rtr.size=10, extraction.lag=0, 
                 NA.sample.size=10, minNonNAs=5)
}

\arguments{
\item{Y}{a numeric matrix or (multivariate) time series object.}

\item{byrow}{logical. If \code{FALSE} (the default), the filtering is done by columns, otherwise the filtering is done by rows.}

\item{min.width}{a positive integer \eqn{\ge 10} specifying the minimal width of the moving time window.}

\item{max.width}{a positive integer \eqn{\ge} \code{min.width} specifying the maximal width of the moving time window. If \code{min.width} = \code{max.width}, the window width is fixed.}

\item{start.width}{a positive integer \eqn{\ge} \code{min.width}, \eqn{\le} \code{max.width} specifying the width of the first time window and thus the time point of the first signal estimation.}

\item{test.sample.size}{a positive integer \eqn{\ge 5}, \eqn{\le} \code{min.width} defining the number of the most recent \emph{Repeated Median} (RM) regression residuals within each time window. The RM residuals are used for the goodness of fit test from the univariate \code{adore.filter}. At each time point \eqn{t} this test serves for determining an adequate window width \eqn{n(t)}. For more details, see Schettlinger, Fried, Gather (2008).} 

\item{width.search}{a character specifying the way of determining an adequate window width \eqn{n(t)}; \code{width.search} can be either \code{"linear"}, \code{"binary"} or \code{"geometric"}. If the null hypothesis of a 'good fit' (see \code{test.sample.size}) is rejected for the window width \eqn{n*(t)=min[n(t-1)+1, \code{max.width}]}, where \eqn{n(t-1)} is the adapted window width at the previous time point, \eqn{n(t)} must be chosen smaller than \eqn{n*(t)} to obtain an adequate regression fit.\cr
If \code{width.search} = \code{"linear"}, \eqn{n(t)} is determined by decreasing \eqn{n*(t)} stepwise by one until the null hypothesis cannot be rejected or \eqn{n(t)} is equal to \code{min.width}.\cr
If \code{width.search} = \code{"binary"}, \eqn{n(t)} is determined by nested intervals within \eqn{[\code{min.width}, n*(t)]} such that \eqn{n(t)} leads to an adequate regression fit.\cr
If \code{width.search} = \code{"geometric"}, the smallest \eqn{i=1, 2, \ldots} is determined such that \eqn{n*(t)-2^i} is adequate or such that \eqn{n*(t)-2^i} is not bigger than \code{min.width}, respectively. Then \eqn{n(t)} is determined by nested intervals within \eqn{[n*(t)-2^i, n*(t)-2^{i-1}]} such that \eqn{n(t)} leads to an adequate regression fit.}

\item{rtr.size}{a non-negative integer specifying the size of a subset of the most recent observations within each window. The signal estimation is restricted to the range of the observations within this subset.}

\item{extraction.lag}{a non-negative integer \eqn{\le \code{min.width}/2}. At each time t the signal is estimated at the position \eqn{t-\code{extraction.lag}} within the time window \eqn{(t-n(t)+1,\ldots,t)}.}

\item{NA.sample.size}{a positive integer \eqn{\ge 10}, \eqn{\le \code{min.width}} specifying the size of a subset of the most recent observations within each window. See \code{minNonNAs}.}

\item{minNonNAs}{a positive integer \eqn{\ge 5}, \eqn{\le \code{NA.sample.size}}. If a variable does not offer at least \code{minNonNAs} non-missing observations within the subset specified by \code{NA.sample.size}, the signal is not estimated for this variable at this time point \eqn{t}.}
}

\details{
The \code{madore.filter} is based on \emph{Repeated Median} regression (Siegel, 1982) in moving time windows and serves for separating signals from noise and outliers in multivariate time series. At each time point \eqn{t} the test procedure of the \emph{adaptive online Repeated Median} filter (Schettlinger, Fried, Gather, 2008) is used to determine an appropriate window width \eqn{n(t) \in [\code{min.width}, \code{max.width}]}. Then the signal vector at time t is estimated within the time window \eqn{(t-n(t)+1,\ldots,t)} by a slight modification of the multivariate \emph{Trimmed Repeated Median-Least Squares} regression (Lanius, Gather, 2004). A more detailed description of the filter can be found in Borowski, Schettlinger, Gather (2008).
}

\value{
\code{madore.filter} returns an object of class \code{madore.filter}. An object of class \code{madore.filter} is a list containing the following components:
\item{signals}{a matrix containing the estimated signal vectors at each time point \eqn{t}.}
\item{widths}{a matrix containing the individual window widths of each variable at each time point \eqn{t}.}
\item{overall.width}{a vector containing the overall window widths at each time point \eqn{t}.}
In addition, the original input data is returned as list member \code{Y}, and the settings used for the analysis are returned as the list members \code{byrow}, \code{min.width}, \code{max.width}, \code{start.width}, \code{test.sample.size}, \code{width.search}, \code{rtr.size}, \code{extraction.lag}, \code{NA.sample.size}, and \code{minNonNAs}. Application of the function \code{plot} to an object of class \code{madore.filter} returns a plot showing the original multivariate time series with the filtered output.
}

\references{ 
Borowski, M., Schettlinger, K., Gather, U. (2008)
Multivariate Real Time Signal Extraction by a Robust Adaptive Regression Filter, submitted.\cr

Lanius, V., Gather, U. (2007)
Robust Online Signal Extraction from Multivariate Time Series, 
\emph{Technical Report 38/07, SFB 475, Universität Dortmund, Germany.}\cr

Schettlinger, K., Fried, R., Gather, U. (2008) 
Real Time Signal Processing by Adaptive Repeated Median Filters,
\emph{International Journal of Adaptive Control and Signal Processing}, submitted.

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians, 
\emph{Biometrika} \bold{69} (1), 242-244.\cr
}

\author{Matthias Borowski}

\seealso{\code{\link{robreg.filter}}, \code{\link{adore.filter}}.}

\examples{
# load multivariate time series sample
data(multi.ts)

# extract signals from 'multi.ts' by madore.filter
# this may take some time, depending on your system
extr <- madore.filter(multi.ts, min.width=30, max.width=100, extraction.lag=5)
plot(extr)
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}
\keyword{multivariate}
\keyword{regression}
