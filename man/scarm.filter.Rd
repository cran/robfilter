\name{scarm.filter}
\alias{scarm.filter}       
\title{SCARM (Slope Comparing Adaptive Repeated Median)}

\description{A procedure for robust online signal extraction from univariate time series by a moving window technique with adaptive window width selection based on the detection of signal changes}

\usage{scarm.filter(x, right.width=30, min.left.width=right.width, 
                    min.width=right.width/3, max.width=180, 
                    sign.level=0.001, bound.noise.sd=1)}

\arguments{
\item{x}{a numeric vector or (univariate) time series object.}
\item{right.width}{a positive integer >=5 defining the fixed width of the right-hand window.}
\item{min.left.width}{a positive integer >=5 defining the minimum width of the left-hand window.}
\item{min.width}{a positive integer specifying the minimal width of the whole window.}
\item{max.width}{a positive integer \code{>= min.width} and \code{>= right.width+min.left.width} specifying the maximum width of the whole window.}
\item{sign.level}{significance level of the test procedure; must be a value in (0,0.5)}
\item{bound.noise.sd}{a lower bound for the estimation of the noise standard deviation; this bound ensures that the noise estimation is zero due to ties in the data; must be a value > 0}
}

\details{
The scarm.filter works by applying Repeated Median (RM) regression (Siegel, 1982) to a moving time window with a length varying between min.width and max.width.

For each point in time, the window width is adapted to the current data situation by a test which is based on the comparison of two RM slopes estimated in separated sub-windows, a right-hand and a left-hand window. The choice of the fixed width of the right-hand window \code{right.width} is crucial to distinguish between a patch of outliers and a signal change. If a patch of < k aberrant observations is to be treated as an outlier-patch, it is recommended to set \code{right.width = 3*k}. 

A more detailed description of the filter can be found in Borowski and Fried (2011).
}


\value{
\code{scarm.filter} returns an object of class \code{scarm.filter}.  
An object of class \code{scarm.filter} is a list containing the 
following components:

\item{signal}{a vector containing the signal extractions}
\item{adapted.width}{a vector containing the adapted window widths}
\item{scarm.statistic}{a vector containing the SCARM test statistics}
\item{critvals}{a vector containing the critical values for test decision}
\item{noise.sd}{a vector containing the noise standard deviation estimated by a regression-free scale estimator (Gelper et al., 2009)}

In addition, the original input data and the input arguments used for the analysis are returned as list members.

Application of the function \code{plot} to an object of class
\code{scarm.filter} returns a plot showing the original time series 
with the filtered output. 
}


\references{ 
Borowski, M., Fried, R. (2011)
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
set.seed(12)
N <- 300
noise <- rnorm(N)
signal <- c(rep(0,200),rep(10,100))
x <- signal+noise

# SCARM Filter 
signal.extr <- scarm.filter(x)

par(mfrow=c(2,1))
ts.plot(x, col="darkgrey", main="time series data and signal extraction")
lines(signal.extr$signal)
ts.plot(signal.extr$adapted.width, main="adapted window width")
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}
