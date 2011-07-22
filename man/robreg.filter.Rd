\name{robreg.filter} 
\alias{robreg.filter} 

\title{Robust Regression Filters for Univariate Time Series} 

\description{ Procedures for robust (online) extraction of low 
frequency components (the \emph{signal}) from a univariate time 
series by applying robust regression techniques to moving time 
windows. 
}

\usage{
robreg.filter(y, width, method = "all", h = floor(width/2)+1,   
                 minNonNAs = 5, online = FALSE, extrapolate = TRUE)

%lts.filter(y, width, h = floor(width/2)+1, 
%              online = FALSE, extrapolate = TRUE)
%              
%med.filter(y, width, minNonNAs = 5, online = FALSE, extrapolate = TRUE)
%rm.filter (y, width, minNonNAs = 5, online = FALSE, extrapolate = TRUE)
%
%dr.filter (y, width, online = FALSE, extrapolate = TRUE)
%lms.filter(y, width, online = FALSE, extrapolate = TRUE)
%lqd.filter(y, width, online = FALSE, extrapolate = TRUE)
}

\arguments{ 
\item{y}{a numeric vector or (univariate) time series object.} 
\item{width}{a positive integer defining the window width used for fitting.\cr
    If \code{online=FALSE} (see below) this needs to be an odd integer.}
\item{method}{a (vector of) character string(s) containing the method(s) to be used for robust 
    approximation of the signal within one time window. 
    It is possible to specify any combination of the values:
    \describe{\item{\code{"DR"}}{Deepest Regression}
              \item{\code{"LMS"}}{Least Median of Squares regression} 
              \item{\code{"LQD"}}{Least Quartile Difference regression}
              \item{\code{"LTS"}}{Least Trimmed Squares regression}
              \item{\code{"MED"}}{Median}
              \item{\code{"RM"}}{Repeated Median regression}
              \item{\code{"all"}}{all of the above (default)} 
             } 
    Using \code{dr.filter}, \code{lms.filter}, \code{lqd.filter}, \code{lts.filter}, \code{med.filter} or \code{rm.filter}
    forces \code{"DR"}, \code{"LMS"}, \code{"LQD"}, \code{"LTS"}, \code{"MED"} or \code{"RM"} respectively. \cr
    Currently, only \code{method="MED"} and \code{method="RM"} (\code{med.filter} / \code{rm.filter}) 
    can handle missing values in the input time series.
    For the other regression filters missing values have to be replaced before the analysis.
    }
\item{h}{a positive integer defining the trimming quantile for 
    LTS regression.}
\item{minNonNAs}{a positive integer defining the minimum number of
    non-missing observations within one window which is required 
    for a \sQuote{sensible} estimation. 
    Currently, this option only has an effect for the two methods \code{"MED"} and /or \code{"RM"} (see \code{method}).}
\item{online}{a logical indicating whether the current level estimate is 
    evaluated at the most recent time within each time window 
    (\code{TRUE}) or centred within each window (\code{FALSE}). 
    Setting \code{online=FALSE} requires the \code{width} to be odd. 
    Default is \code{online=FALSE}.
}
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series. \cr
    If \code{online=FALSE} the extrapolation consists of the 
    fitted values within the first half of the first window and the 
    last half of the last window; if \code{online=TRUE} the 
    extrapolation consists of the fitted values within the first 
    time window. }
}

\details{\code{robreg.filter} is suitable for extracting low 
frequency components (the \emph{signal}) from a time series which 
may be contaminated with outliers and can contain level shifts. 
For this, robust regression methods are applied to a moving 
window, and the signal level is estimated by the fitted value 
either at the end of each time window for online signal 
extraction without time delay (\code{online=TRUE}) or in the 
centre of each time window (\code{online=FALSE}).}

\value{
\code{robreg.filter} returns an object of class \code{robreg.filter}. 
An object of class \code{robreg.filter} is a list containing the 
following components:

\item{level}{a data frame containing the signal level extracted by the filter(s) specified in \code{method}.}
\item{slope}{a data frame containing the corresponding slope within each time window.}

In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{width}, \code{method}, 
\code{h}, \code{minNonNAs}, \code{online} and \code{extrapolate}. 

Application of the function \code{plot} to an object of class
\code{robreg.filter} returns a plot showing the original time series 
with the filtered output. 
}

\references{ 
Davies, P.L., Fried, R., Gather, U. (2004) 
Robust Signal Extraction for On-Line Monitoring Data, 
\emph{Journal of Statistical Planning and Inference} \bold{122}, 
65-78.\cr
(earlier version: \url{http://www.sfb475.uni-dortmund.de/berichte/tr02-02.ps})

Gather, U., Schettlinger, K., Fried, R. (2006)
Online Signal Extraction by Robust Linear Regression,
\emph{Computational Statistics} \bold{21}(1), 
33-51. \cr
(earlier version: \url{http://www.sfb475.uni-dortmund.de/berichte/tr53-04.ps})
  
Schettlinger, K., Fried, R., Gather, U. (2006) 
Robust Filters for Intensive Care Monitoring: Beyond the Running Median,
\emph{Biomedizinische Technik} \bold{51}(2), 
49-56.
}

\author{C++ code: Thorsten Bernholt and Robin Nunkesser\cr
Port to R: Roland Fried and Karen Schettlinger}

\note{
Missing values are treated by omitting them and thus by reducing 
the corresponding window width. The estimated signal level is 
only returned as \code{NA} if the window the estimation is based 
on contains less than \code{minNonNAs} non-missing values.
}

\seealso{\code{\link{wrm.filter}}, \code{\link{robust.filter}}, \code{\link{dw.filter}}, \code{\link{hybrid.filter}}.}

\examples{
# Generate random time series:
y <- cumsum(runif(500)) - .5*(1:500)
# Add jumps:
y[200:500] <- y[200:500] + 5
y[400:500] <- y[400:500] - 7
# Add noise:
n <- sample(1:500, 30)
y[n] <- y[n] + rnorm(30)

# Filtering with all methods:
y.rr <- robreg.filter(y, width=31, method=c("RM", "LMS", "LTS", "DR", "LQD"))
# Plot:
plot(y.rr)

# Delayed filtering with RM and LMS filter:
y2.rr <- robreg.filter(y,width=31,method=c("RM","LMS"))
plot(y2.rr)

# Online filtering with RM filter:
y3.rr <- rm.filter(y,width=41,online=TRUE)
plot(y3.rr)
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}
