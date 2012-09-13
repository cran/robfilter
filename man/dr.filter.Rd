\name{dr.filter}
\alias{dr.filter}

\title{
Deepest Regression (DR) filter
}
\description{
This function extracts signals from time series by means of Deepest regression in a moving time window.
}
\usage{
dr.filter(y, width, online = FALSE, extrapolate = TRUE)
}

\arguments{
\item{y}{a numeric vector or (univariate) time series object.} 
\item{width}{a positive integer defining the window width used for fitting.\cr
    If \code{online=FALSE} (see below) this needs to be an odd integer.}
\item{online}{a logical indicating whether the current level estimate is 
    evaluated at the most recent time within each time window 
    (\code{TRUE}) or centred within each window (\code{FALSE}). 
    Setting \code{online=FALSE} requires the \code{width} to be odd. 
    Default is \code{online=FALSE}.}  
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series. \cr
    If \code{online=FALSE} the extrapolation consists of the 
    fitted values within the first half of the first window and the 
    last half of the last window; if \code{online=TRUE} the 
    extrapolation consists of the fitted values within the first 
    time window.}
}

\details{
\code{dr.filter} is suitable for extracting low 
frequency components (the \emph{signal}) from a time series which 
may be contaminated with outliers and can contain level shifts. 
For this, robust Deepest Regression is applied to a moving 
window, and the signal level is estimated by the fitted value 
either at the end of each time window for online signal 
extraction without time delay (\code{online=TRUE}) or in the 
centre of each time window (\code{online=FALSE}).
}

\value{
\code{dr.filter} returns an object of class \code{robreg.filter}. 
An object of class \code{robreg.filter} is a list containing the 
following components:

\item{level}{a data frame containing the extracted signal level.}
\item{slope}{a data frame containing the corresponding slope within each time window.}

In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{width}, \code{online} and \code{extrapolate}. 

Application of the function \code{plot} to an object of class
\code{robreg.filter} returns a plot showing the original time series 
with the filtered output. 
}

\references{
Davies, P.L., Fried, R., Gather, U. (2004) 
Robust Signal Extraction for On-Line Monitoring Data, 
\emph{Journal of Statistical Planning and Inference} \bold{122}, 
65-78.\cr

Gather, U., Schettlinger, K., Fried, R. (2006)
Online Signal Extraction by Robust Linear Regression,
\emph{Computational Statistics} \bold{21}(1), 
33-51. \cr

Schettlinger, K., Fried, R., Gather, U. (2006) 
Robust Filters for Intensive Care Monitoring: Beyond the Running Median,
\emph{Biomedizinische Technik} \bold{51}(2), 
49-56.
}

\author{
Roland Fried, Karen Schettlinger and Matthias Borowski
}

\seealso{
\code{\link{robreg.filter}}
}

\examples{
# Generate random time series:
y <- cumsum(runif(500)) - .5*(1:500)
# Add jumps:
y[200:500] <- y[200:500] + 5
y[400:500] <- y[400:500] - 7
# Add noise:
n <- sample(1:500, 30)
y[n] <- y[n] + rnorm(30)

# Online filtering with DR filter:
y.rr <- dr.filter(y,width=41,online=TRUE)
plot(y.rr)
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}