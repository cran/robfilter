\name{scarm.filter}
\alias{scarm.filter}       
\title{SCARM (Slope Comparing Adaptive Repeated Median)}

\description{A procedure for robust online signal extraction from univariate time series ("smoothing") by a moving window technique with adaptive window width selection based on Repeated Median regression}

\usage{
scarm.filter(time.series,
             right.width=30, min.left.width=right.width,
             min.width=floor(right.width/3), max.width=200,
             sign.level=0.001, bound.noise.sd=0.01, rtr=TRUE,
             autocorrelations="automatic")
}

\arguments{
\item{time.series}{a numeric vector or (univariate) time series object.}
\item{right.width}{a positive integer >=5 defining the fixed width of the right-hand window used for testing; the choice of \code{right.width} is crucial to distinguish between a patch of outliers and a signal change.}
\item{min.left.width}{a positive integer \eqn{\ge} \code{right.width} defining the minimum width of the left-hand window used for testing.}
\item{min.width}{a positive integer \eqn{\ge 5} specifying the minimum window width.}
\item{max.width}{a positive integer \eqn{\ge} \code{min.width} and \eqn{\ge} \code{right.width} + \code{min.left.width} specifying the maximum window width.}
\item{sign.level}{significance level of the test procedure; must be a value in \eqn{(0,0.5)}.}
\item{bound.noise.sd}{a lower bound for the estimate of the noise standard deviation; this bound ensures that the noise scale estimate cannot be zero due to ties in the data; must be a value \eqn{> 0}.}
\item{rtr}{if \code{rtr=TRUE}, the signal estimation is restricted to the range of the rightmost \code{min.width} observations.}
\item{autocorrelations}{the \code{scarm.filter} is developed for non-autocorrelated data, but can be adapted to work for AR(1) processes with parameter \eqn{\phi = -0.9,-0.6,...,0.9}; \code{autocorrelations} must be either "no" (\eqn{\phi=0}), "high.positive" (\eqn{\phi=0.9}), "moderate.positive" (\eqn{\phi=0.6}), "small.positive" (\eqn{\phi=0.3}), "small.negative" (\eqn{\phi=-0.3}), "moderate.negative (\eqn{\phi=-0.6})", "high.negative (\eqn{\phi=-0.9})" or "automatic"; if \code{autocorrelations="automatic"}, the true parameter \eqn{\phi} is estimated at each time point.}
}

\details{
The \code{scarm.filter} fits a Repeated Median (RM, Siegel, 1982) regression line to a moving window sample with length varying between \code{min.width} and \code{max.width}.

For each time point, the window width is adapted to the current data situation by a test comparing two RM slopes estimated in separated sub-windows, a right-hand and a left-hand window.

A more detailed description of the filter can be found in Borowski and Fried (2011).
}

\value{
\code{scarm.filter} returns an object of class \code{scarm.filter}.  
An object of class \code{scarm.filter} is a list containing the 
following components:

\item{signal.est}{a vector containing the signal estimations}
\item{slope.est}{a vector containing the slope (or trend) estimations}
\item{adapted.width}{a vector containing the adapted window widths}
\item{test.statistic}{a vector containing the SCARM test statistics}
\item{critvals}{a vector containing the critical values for test decision}
\item{noise.sd}{a vector containing the estimates of the noise standard deviation}
\item{slope.diff}{a vector containing the differences of the Repeated Median slopes estimated in the left-hand and right-hand window}
\item{acf.lag.one}{a vector containing the estimated autocorrelations at lag one for each time point; estimation is done on the recent \code{max.width} observations at each time point}
\item{time.series}{the original input data}

In addition, the input arguments used for the analysis are returned as list members.

Application of the function \code{plot} to an object of class
\code{scarm.filter} returns a plot showing the original time series 
with the filtered output. If \code{info==TRUE} (default), a plot of the adapted window widths is also given.
}


\references{ 
Borowski, M. and Fried, R. (2011)
Robust repeated median regression in moving windows with data-adaptive width selection, \emph{Discussion Paper 28/2011, SFB 823, TU Dortmund University}.\cr

Gelper, S., Schettlinger, K., Croux, C., and Gather, U. (2009)
Robust online scale estimation in time series: A model-free approach,
\emph{Journal of Statistical Planning and Inference}, \bold{139}(2), 335-349.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians,
\emph{Biometrika} \bold{69}(1), 242-244.\cr
}

\author{Matthias Borowski}

\seealso{\code{\link{robreg.filter}}, \code{\link{adore.filter}}, \code{\link{madore.filter}}, \code{\link{mscarm.filter}}.}

\examples{
# Time series
data(multi.ts)
x <- multi.ts[,1]

# apply SCARM Filter 
scarm.extr <- scarm.filter(x)
plot(scarm.extr)
}
\keyword{robust}
\keyword{smooth}
\keyword{ts}
