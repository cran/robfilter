\name{scarm.filter}
\alias{scarm.filter}       
\title{SCARM (Slope Comparing Adaptive Repeated Median)}

\description{A procedure for robust online signal extraction from univariate time series by a moving window technique with adaptive window width selection based on the Repeated Median regression}

\usage{scarm.filter(x, right.width=15, min.left.width=right.width, min.width=floor(right.width/3), max.width=200, sign.level=0.001, bound.slope.diff=0, bound.noise.sd=0, rtr=TRUE)}

\arguments{
\item{x}{a numeric vector or (univariate) time series object.}
\item{right.width}{a positive integer >=5 defining the fixed width of the right-hand window used for testing; the choice of \code{right.width} is crucial to distinguish between a patch of outliers and a signal change.}
\item{min.left.width}{a positive integer \eqn{\ge} \code{right.width} defining the minimum width of the left-hand window used for testing.}
\item{min.width}{a positive integer \eqn{\ge 5} specifying the minimum window width.}
\item{max.width}{a positive integer \eqn{\ge} \code{min.width} and \eqn{\ge} \code{right.width} + \code{min.left.width} specifying the maximum window width.}
\item{sign.level}{significance level of the test procedure; must be a value in \eqn{(0,0.5)}.}
\item{bound.slope.diff}{a lower bound for the absolute difference of the Repeated Median slopes in the left-hand and right-hand window; if this absolute difference is smaller than 'bound.slope.diff', the window width is not decreased; must be a value \eqn{> 0}.}
\item{bound.noise.sd}{a lower bound for the estimation of the noise standard deviation; this bound ensures that the noise estimation is zero due to ties in the data; must be a value \eqn{> 0}.}
\item{rtr}{if \code{rtr=TRUE}, the signal estimation is restricted to the range of the rightmost \code{min.width} observations.}
}

\details{
The scarm.filter works by applying Repeated Median (RM) regression (Siegel, 1982) to a moving time window with a length varying between min.width and max.width.

For each point in time, the window width is adapted to the current data situation by a test which is based on the comparison of two RM slopes estimated in separated sub-windows, a right-hand and a left-hand window.

A more detailed description of the filter can be found in Borowski and Fried (2012).
}


\value{
\code{scarm.filter} returns an object of class \code{scarm.filter}.  
An object of class \code{scarm.filter} is a list containing the 
following components:

\item{signal.est}{a vector containing the signal estimations}
\item{slope.est}{a vector containing the slope estimations}
\item{adapted.width}{a vector containing the adapted window widths}
\item{test.statistic}{a vector containing the SCARM test statistics}
\item{critvals}{a vector containing the critical values for test decision}
\item{noise.sd}{a vector containing the noise standard deviation estimated by a regression-free scale estimator (Gelper et al., 2009)}
\item{slope.diff}{a vector containing the differences of the Repeated Median slopes estimated in the left-hand and right-hand window}

In addition, the original input data and the input arguments used for the analysis are returned as list members.

Application of the function \code{plot} to an object of class
\code{scarm.filter} returns a plot showing the original time series 
with the filtered output. 
}


\references{ 
Borowski, M., Fried, R. (2012)
Robust moving window regression for online signal extraction from
non-stationary time series: online window width adaption by testing for
signal changes, submitted.\cr

Gelper, S., Schettlinger, K., Croux, C., and Gather, U. (2009)
Robust online scale estimation in time series: A model-free approach,
\emph{Journal of Statistical Planning and Inference}, \bold{139} (2), 335-349.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians,
\emph{Biometrika} \bold{69} (1), 242-244.\cr
}

\author{Matthias Borowski}

\seealso{\code{\link{robreg.filter}}, \code{\link{adore.filter}}, \code{\link{madore.filter}}.}

\examples{
# # # # # # # # # #
# Short and noisy time series
data(multi.ts)
x <- multi.ts[,1]

# SCARM Filter 
scarm.extr <- scarm.filter(x)
plot(scarm.extr)
}
\keyword{robust}
\keyword{smooth}
\keyword{ts}
