\name{mscarm.filter}
\alias{mscarm.filter}       
\title{MSCARM (Multivariate Slope Comparing Adaptive Repeated Median)}

\description{A multivariate version of the \code{scarm.filter} which also gives information about parallel running components of the multivariate time series}

\usage{
mscarm.filter(time.series,
              right.width=30, min.left.width=right.width,
              min.width=floor(right.width/3), max.width=200,
              sign.level=0.001, bound.noise.sd=0.01,
              rtr=TRUE, autocorrelations="automatic",
              c.bound=3, r.bound=0)
}

\arguments{
\item{time.series}{a numeric matrix or multivariate time series object.}
\item{right.width}{a positive integer \eqn{\ge 5} defining the fixed width of the right-hand window used for the SCARM test; the choice of \code{right.width} is crucial to distinguish between a patch of outliers and a signal change.}
\item{min.left.width}{a positive integer \eqn{\ge} \code{right.width} defining the minimum width of the left-hand window used for testing.}
\item{min.width}{a positive integer \eqn{\ge 5} specifying the minimum window width.}
\item{max.width}{a positive integer \eqn{\ge} \code{min.width} and \eqn{\ge} \code{right.width} + \code{min.left.width} specifying the maximum window width.}
\item{sign.level}{significance level of the SCARM test procedure; must be a value in \eqn{(0,0.5)}.}
\item{bound.noise.sd}{a lower bound for the estimate of the noise standard deviation; this bound ensures that the noise scale estimate cannot be zero due to ties in the data; must be a value \eqn{> 0}.}
\item{rtr}{if \code{rtr=TRUE}, the signal estimation is restricted to the range of the rightmost \code{min.width} observations.}
\item{autocorrelations}{the \code{mscarm.filter} is developed for non-autocorrelated data, but can be adapted to work for AR(1) processes with parameter \eqn{\phi = -0.9,-0.6,...,0.9}; \code{autocorrelations} must be either "no" (\eqn{\phi=0}), "high.positive" (\eqn{\phi=0.9}), "moderate.positive" (\eqn{\phi=0.6}), "small.positive" (\eqn{\phi=0.3}), "small.negative" (\eqn{\phi=-0.3}), "moderate.negative (\eqn{\phi=-0.6})", "high.negative (\eqn{\phi=-0.9})" or "automatic"; if \code{autocorrelations="automatic"}, the true parameter \eqn{\phi} is estimated at each time point.}
\item{c.bound}{the bound for the SSM (Similar Slope Monitoring) statistic which is the absolute difference of RM (Repeated Median) trend estimates (of two univariate time series) relative to the estimated standard deviation of this difference; if the SSM statistic is not larger than \code{c.bound}, a relationship between the two time series is assumed; \code{c.bound} must be a value >0.}
\item{r.bound}{the bound for the ratio of adapted window widths (of two univariate time series); a relationship between the two time series is only possible, if this ratio (smaller window width divided by larger window width) is smaller than \code{r.bound}; \code{r.bound} must be a value in \eqn{[0,1]}.}
}

\details{
The \code{mscarm.filter} is a procedure for real-time signal extraction from noisy and outlier-contaminated instationary multivariate time series. It is based on \emph{Repeated Median} regression (Siegel, 1982) in moving time windows. At each time point \eqn{t} the test procedure of the SCARM filter (Borowski and Fried, 2011) is used to determine an appropriate window width \eqn{n(t)} in [\code{min.width}, \code{max.width}]. Then the signal vector at time \eqn{t} is estimated within the time window \eqn{(t-n(t)+1,\ldots,t)} by a slight modification of the multivariate \emph{Trimmed Repeated Median-Least Squares} regression (Lanius, Gather, 2010). At each time point \eqn{t}, the \code{mscarm.filter} uses the \emph{Similar Slope Monitoring} (SSM) method to build blocks of currently interrelated univariate time series. This information is given to the user and is used to improve the signal estimations. A detailed description of the \code{mscarm.filter} can be found in Borowski (2012).
}

\value{
\code{mscarm.filter} returns an object of class \code{mscarm.filter}.  
An object of class \code{mscarm.filter} is a list containing the 
following components:

\item{signal.est}{a matrix containing the signal estimations}
\item{slope.est}{a matrix containing the slope (or trend) estimations}
\item{adapted.width}{a matrix containing the adapted window widths}
\item{noise.sd.est}{a matrix containing the estimated noise standard deviations}
\item{scarm.signal.est}{a matrix containing the signal estimates of the univariate SCARM}
\item{scarm.width}{a matrix containing the adapted window widths of the univariate SCARM}
\item{scarm.statistic}{a matrix containing the SCARM test statistics}
\item{scarm.critval}{a matrix containing the critical values of the SCARM test}
\item{ssm.statistic}{a matrix containing the SSM statistics}
\item{blocks}{a matrix of the blocks built by the SSM procedure}
\item{acf.lag.one}{a matrix containing the estimated autocorrelations at lag one for each time point; estimation is done on the recent \code{max.width} observations at each time point}
\item{time.series}{the original input data}

In addition, the input arguments used for the analysis are returned as list members.

Application of the function \code{plot} to an object of class
\code{mscarm.filter} returns a plot showing the original time series 
with the filtered output. If \code{info==TRUE} (default), a plot showing the results of the SSM procedure is given.
}

\references{ 
Borowski, M. (2012)
Echtzeit-Extraktion relevanter Information aus multivariaten Zeitreihen basierend auf robuster Regression, \emph{PhD thesis, TU Dortmund University (in German)}.\cr

Borowski, M. and Fried, R. (2011)
Robust repeated median regression in moving windows with data-adaptive width selection, \emph{Discussion Paper 28/2011, SFB 823, TU Dortmund University}.\cr

Lanius, V., Gather, U. (2010)
Robust Online Signal Extraction from Multivariate Time Series, 
\emph{Computational Statistics and Data Analysis} \bold{54}(4), 966-975.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians,
\emph{Biometrika} \bold{69}(1), 242-244.\cr
}

\author{Matthias Borowski}

\seealso{\code{\link{robreg.filter}}, \code{\link{adore.filter}}, \code{\link{madore.filter}}, \code{\link{scarm.filter}}.}

\examples{
\dontrun{
# Multivariate time series
data(multi.ts)

# apply MSCARM Filter 
mscarm.extr <- mscarm.filter(multi.ts)
plot(mscarm.extr)
}
}
\keyword{robust}
\keyword{smooth}
\keyword{ts}
