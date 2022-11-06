\name{hybrid.filter} 
\alias{hybrid.filter} 

\title{Robust Hybrid Filtering Methods for Univariate Time Series} 

\description{
Procedures for robust extraction of low frequency components (the 
\emph{signal}) from a univariate time series based on a moving 
window technique using the median of several one-sided 
half-window estimates (subfilters) in each step. 
}

\usage{
hybrid.filter(y, width, method = "all", minNonNAs=3, extrapolate = TRUE)
}

\arguments{
\item{y}{a numeric vector or (univariate) time series object.} 
\item{width}{an odd positive integer (\eqn{\geq 3}{>=3}) 
defining the window width used for fitting.}
\item{method}{a (vector of) character string(s) containing the method(s) to be used for the estimation
    of the signal level.\cr
    It is possible to specify any combination of 
    \code{"MED"}, \code{"RM"}, \code{"MEAN"}, 
    \code{FMH}, \code{"PFMH"}, \code{"CFMH"}, 
    \code{"MH"}, \code{"PRMH"}, \code{"CRMH"}, 
    \code{"MMH"}, \code{"PRMMH"}, \code{"CRMMH"}, 
    and \code{"all"} (for all of the above). Default is \code{method="all"}. 
    For a detailed description see the section \sQuote{Methods} below.
    }
\item{minNonNAs}{a positive integer defining the minimum number 
    of non-missing observations within each window (half) which is required  
    for a \sQuote{sensible} estimation. Default: if a window (half) contains
    less than \code{minNonNAs = 3} observations an \code{NA} is returned (for that subfilter).}
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series.
    The extrapolation extends the first estimated value to the first 
    time in the first window and the last estimated value to 
    the last time in the last time window. Default is 
    \code{extrapolate=TRUE}.
    } 
}

\section{Methods}{ 
The following methods are available as \code{method} for signal extraction.

Filters applying only \emph{one} location or regression estimate 
to the whole window of length \code{width} and taking the 
location (in the centre of the time window) as final signal level
estimate:
\describe{
\item{\code{MED}}{ordinary running median filter.}
\item{\code{MEAN}}{ordinary moving average filter.}
\item{\code{RM}}{ordinary repeated median filter. \cr
    Applies repeated median regression to each time window.}
}

Filters applying several subfilters within one window, taking the 
median of the values listed below as the final signal level 
estimate:
\describe{
\item{\code{FMH}}{FIR median hybrid filter. \cr
Uses half-window averages and the central observation.} 
\item{\code{PFMH}}{predictive \code{FMH} filter. \cr
    Uses half-window least squares regression and the central 
    observation.} 
\item{\code{CFMH}}{combined \code{FMH} filter.\cr
    Uses half-window averages, half-window least squares 
    regression, and the central observation.} 
\item{\code{MH}}{median hybrid filter. \cr
    Uses half-window medians and the central observation.} 
\item{\code{PRMH}}{predictive repeated median hybrid filter. \cr 
    Uses half-window repeated median regression and the 
    central observation.} 
\item{\code{CRMH}}{combined repeated median hybrid filter. \cr 
    Uses half-window medians, half-window repeated median 
    regression, and the central observation.} 
\item{\code{MMH}}{median/median hybrid filter. \cr
    Uses half-window medians and the median of all 
    observations in the window.} 
\item{\code{PRMMH}}{predictive repeated median/median filter. \cr 
    Uses half-window repeated median regression and the 
    median of all observations in the window.} 
\item{\code{CRMMH}}{combined repeated median/median filter. \cr 
    Uses half-window medians, half-window repeated median 
    regression, and the median of all observations in the window.}
} % end itemize
} % end section 'Method'

\details{
\code{hybrid.filter} is suitable for extracting low frequency 
components (the \emph{signal}) from a time series which may be 
contaminated with outliers and can contain level shifts or local 
extremes. For this, moving window techniques are applied.

Within each time window several subfilters are applied to 
half-windows (left and right of the centre); the final signal 
level in the centre of the time window is then estimated by the 
median of the subfilter outputs.

For the subfilters, both, location-based and regression-based 
method are available, the former applying means or medians and 
the idea of a locally constant signal value, the latter using 
ordinary least squares (LS) regression or Siegel's (1982) 
repeated median (RM) and the idea of an underlying locally linear 
trend. 

The methods should be chosen based on an a-priori guess of the 
underlying signal and the data quality.  Location based methods 
(\code{MED}, \code{MEAN}, \code{FMH}, \code{MH}, \code{MMH})  are 
recommended in case of a locally (piecewise) constant signal. 
Regression based and predictive approaches (\code{RM}, 
\code{PFMH}, \code{PRMH}, \code{PRMMH}) in case of locally linear 
monotone trends.  The combined filters (\code{CFMH}, \code{CRMH}, 
\code{CRMMH}) can be seen as a compromise, but are 
computationally somewhat more expensive and may be inferior to 
the predictive filters during steep trends.

The approaches based on the median and RM are robust alternatives 
to the (in Gaussian samples) more efficient mean and least 
squares methods. The hybrid filters preserve shifts and local 
extremes much better than \code{MED}, \code{MEAN} or \code{RM} 
for the price of decreased robustness and / or Gaussian 
efficiency.
}

\value{\code{hybrid.filter} returns an object of class \code{hybrid.filter}.  
An object of class \code{hybrid.filter} is a list containing the 
following components:

\item{level}{a data frame containing the signal level extracted 
by the filter(s) specified in \code{method}.} 

\item{slope}{a data frame (possibly) containing \code{RM}, 
\code{RM.left}, \code{RM.right}, \code{LS.left} and 
\code{LS.right}: the slope estimated by Repeated Median regression 
in the whole window (for \code{method="RM"}) or in the left and right 
window half (for any \code{method} in \code{"PRMH"}, 
\code{"CRMH"}, \code{"PRMMH"} and \code{"CRMMH"}) or the least 
squares slope estimated from the left and right window half (for any \code{method} in 
\code{"PRFMH"} or \code{"CFMH"}).\cr Only those slopes are 
returned which are required by the filters specified in 
\code{method}. If only location-based filters are applied (i.e. 
\code{"MED"}, \code{"MEAN"}, \code{"FMH"}, \code{"MH"} and /or 
\code{"MMH"}) \code{NULL} is returned for the \code{slope}.}
In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{width}, \code{method} and 
\code{extrapolate}. 

Application of the function \code{plot} to an object of class
\code{hybrid.filter} returns a plot showing the original time series 
with the filtered output. 
}


\references{
Fried, R., Bernholt, T., Gather, U. (2006) 
Repeated Median and Hybrid Filters, \emph{Computational Statistics & Data Analysis} \bold{50}, 
2313-2338. \cr
(earlier version: \url{http://hdl.handle.net/2003/4866})

Schettlinger, K., Fried, R., Gather, U. (2006) 
Robust Filters for Intensive Care Monitoring: Beyond the Running Median, \emph{Biomedizinische Technik} \bold{51}(2), 
49-56.
} 

\author{Roland Fried and Karen Schettlinger}

\note{Missing values are treated by omitting them and thus by 
reducing the corresponding window width. \cr 
The \code{hybrid.filter} function only offers filters for signal 
extraction delayed by (\code{width}+1)/2 time units, in contrast 
to other filters available from the \code{robfilter} package 
which also offer online time series analysis without time delay.
}

\seealso{\code{\link{robreg.filter}}, \code{\link{robust.filter}}, \code{\link{dw.filter}}, \code{\link{wrm.filter}}.}

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
y.hy <- hybrid.filter(y, width=31)
# Plot:
plot(y.hy)

# Filtering with running median and PRMH only:
y2.hy <- hybrid.filter(y, width=31, method=c("MED","PRMH"))
plot(y2.hy)
}


\keyword{robust}
\keyword{smooth}
\keyword{ts}
