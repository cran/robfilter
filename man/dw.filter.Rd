\name{dw.filter} 
\alias{dw.filter}
\alias{dw.filter.online} 

\title{Robust Double Window Filtering Methods for Univariate Time Series} 

\description{
Procedures for robust (online) extraction of low frequency 
components (the \emph{signal}) from a univariate time series 
based on a moving window technique using two nested time windows 
in each step.
}

\usage{
dw.filter(y, outer.width, inner.width, method = "all", 
             scale = "MAD", d = 2, 
             minNonNAs = 5, online = FALSE, extrapolate = TRUE)
}

\arguments{
\item{y}{a numeric vector or (univariate) time series object.} 
\item{outer.width}{a positive integer specifying the window width 
    of the outer window used for determining the final estimate. \cr
    If \code{online=FALSE} (see below) this needs to be an odd integer.}
\item{inner.width}{a positive integer (not larger than \code{outer.width})
    specifying the window width of the inner window used for determining 
    the initial estimate and trimming features.\cr
    If \code{online=FALSE} (see below) this needs to be an odd integer.}
\item{method}{a (vector of) character string(s) containing the method(s) to be used for the estimation
    of the signal level.\cr
    It is possible to specify any combination of \code{"MED"}, \code{"RM"}, \code{"MTM"}, \code{"TRM"}, 
    \code{"MRM"}, \code{"DWRM"}, \code{"DWMTM"}, \code{"DWTRM"}, \code{"DWMRM"} and 
    \code{"all"} (for all of the above). Default is \code{method="all"}. 
    For a detailed description see the section \sQuote{Methods} below.}
\item{scale}{a character string specifying the method to be used for robust estimation of the local 
   variability (within one time window). Possible values are:
   \describe{\item{\code{"MAD"}}{Median absolute deviation about the median (default)}
             \item{\code{"QN"}}{Rousseeuw's and Croux' (1993) \eqn{Q_n}{Q_n} scale estimator}
             \item{\code{"SN"}}{Rousseeuw's and Croux' (1993) \eqn{S_n}{S_n} scale estimator}
            } }
\item{d}{a positive integer defining factor the current scale estimate is multiplied with for
    determining the trimming boundaries for outlier detection.\cr
    Observations deviating more than \eqn{d\cdot \hat{\sigma}_t}{d*\sigma_t} 
    from the current level approximation \eqn{\hat{\mu}_t}{\mu_t} 
    are replaced by \eqn{\hat{\mu}_t}{\mu_t} where 
    \eqn{\hat{\sigma}_t}{\sigma_t} denotes the current scale estimate.
    Default is \code{d = 2} meaning a \eqn{2\sigma}{2 \sigma} rule  for outlier 
    detection.}
\item{minNonNAs}{a positive integer defining the minimum number 
    of non-missing observations within each window which is required  
    for a \sQuote{sensible} estimation. Default: if windows contain
    less than \code{minNonNAs = 5} observations \code{NA}s are returned.}
\item{online}{a logical indicating whether the current level and 
    scale estimates are evaluated at the most recent time 
    within each (inner and outer) window (\code{TRUE}) or centred within 
    the windows (\code{FALSE}). Setting \code{online=FALSE} requires odd 
    \code{inner.width} and \code{outer.width}. Default is \code{online=FALSE}.}
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series. \cr
    If \code{online=FALSE} the extrapolation consists of the 
    fitted values within the first half of the first window and the 
    last half of the last window; if \code{online=TRUE} the 
    extrapolation consists of the all fitted values within the first 
    time window. }
}

\section{Methods}{ 
The following methods are available as \code{method} for signal extraction, 
whereby the prefix \code{DW} denotes the fact that different 
window widths are used in the first and second step of the 
calculations within one window (i.e. 
\code{inner.width}<\code{outer.width}) while for the 
methods \code{MED}, \code{RM}, \code{MTM}, \code{TRM} and \code{MRM}
the first and second step take place in a window of fixed 
length \code{outer.width}.

\describe{
\item{\code{MED}}{ordinary running median filter. \cr
     The simple median is applied to the observations within a moving 
     time window of length \code{outer.width}.
}
\item{\code{RM}}{ordinary repeated median filter. \cr
    Repeated median regression is applied to the observations within 
    a moving time window of length \code{outer.width}.
}
\item{\code{MTM}, \code{DWMTM}}{modified trimmed mean filters. \cr
    In a first step the median is applied to (\code{MTM}): the whole 
    window with \code{outer.width} or (\code{DWMTM}): the inner window 
    with \code{inner.width}; in a second step the mean is applied to 
    the (trimmed) observations in the whole window (with 
    \code{outer.width}).
}
\item{\code{TRM}, \code{DWTRM}}{trimmed repeated median filters. \cr
    In a first step repeated median regression is applied to 
    (\code{TRM}): the whole window with \code{outer.width} or 
    (\code{DWTRM}): the inner window with \code{inner.width}; in a 
    second step least squares regression is applied to the (trimmed) 
    observations in the whole window (with \code{outer.width}).
}
\item{\code{MRM}, \code{DWMRM}}{modified repeated median filters. \cr
     In a first step repeated median regression is applied to 
    (\code{MRM}): the whole window with \code{outer.width} or 
    (\code{DWMRM}): the inner window with \code{inner.width}; in a 
    second step another repeated median regression is applied to the 
    (trimmed) observations in the whole window (with 
    \code{outer.width}).
} 
\item{\code{DWRM}}{double window repeated median filter. \cr 
    In a first step repeated median regression is applied to the 
    inner window with \code{inner.width} to determine the trend 
    (slope); in a second step the median is applied to the trend 
    corrected observations in the whole window with 
    \code{outer.width} (without trimming). 
}
} % end itemize
} % end section 'Methods'

\details{\code{dw.filter} is suitable for extracting low 
frequency components (the \emph{signal}) from a time series 
which may be contaminated with outliers and can contain level 
shifts. For this, moving window techniques are applied.

A short inner window of length \code{inner.width} is used in each 
step for calculating an initial level estimate (by using either 
the median or a robust regression fit) and a robust estimate of 
the local standard deviation. Observations deviating strongly from this 
initial fit are trimmed from an outer time window of length 
\code{outer.width}, and the signal level is estimated from the 
remaining observations (by using either a location or regression 
estimator). Values specified in \code{method} determine which 
combination of estimation methods should be applied to the inner 
and outer window (see section \sQuote{Methods} below).

The applied \code{method} should be chosen based on an a-priori 
guess of the underlying signal and the data quality: Location 
based method (\code{MED} / \code{MTM}) are recommended in case 
of a locally (piecewise) constant signal, regression based 
approaches (\code{RM} / \code{DWRM} / \code{TRM} / \code{MRM}) in 
case of locally linear, monotone trends. 

Since no big differences have been reported between \code{TRM} 
and \code{MRM}, the quicker and somewhat more efficient 
\code{TRM} option might be preferred. \code{DWRM} is the quickest 
of all regression based methods and performs better than the 
ordinary \code{RM} at shifts, but it is the least robust and 
least efficient method.

If location based methods are used, the \code{inner.width} should 
be chosen at least twice the length of expected patches of 
subsequent outliers in the time series; if regression based 
methods are used, the \code{inner.width} should be at least three 
times that length, otherwise outlier patches can influence the 
estimations strongly. To increase the efficiency of the final 
estimates, \code{outer.width} can then be chosen rather large - 
provided that it is smaller than the time between subsequent 
level shifts.

For robust scale estimation, \code{MAD} is the classical choice; 
\code{SN} is a somewhat more efficient and almost equally robust 
alternative, while \code{QN} is much more efficient if the 
window widths are not too small, and it performs very well at the 
occurrence of level shifts.

The factor \code{d}, specifying the trimming boundaries as a 
multiple of the estimated scale, can be chosen similarly to 
classical rules for detecting unusual observations in a Gaussian 
sample. Choosing \code{d=3} instead of \code{d=2} increases 
efficiency, but decreases robustness; \code{d=2.5} might be seen 
as a compromise.}

\value{
\code{dw.filter} returns an object of class \code{dw.filter}.  
An object of class \code{dw.filter} is a list containing the 
following components:

\item{level}{a data frame containing the corresponding signal level extracted by the filter(s) specified in \code{method}.}
\item{slope}{a data frame containing the corresponding slope within each time window.}
\item{sigma}{a data frame containing 
             \code{inner.loc.sigma}, \code{inner.reg.sigma}, \code{outer.loc.sigma} and \code{outer.reg.sigma},
             the scale estimated from the observations (\code{loc}) or the residuals from the Repeated Median regression (\code{reg}) 
             within the inner window of length \code{inner.width} or the outer window of length \code{outer.width}, respectively.\cr
             \code{MTM} uses \code{outer.loc.sigma} for trimming outliers, 
             \code{MRM} and \code{TRM} use \code{outer.reg.sigma} for trimming outliers,\cr
             \code{DWMTM} uses \code{inner.loc.sigma} for trimming outliers, 
             \code{DWMRM} and \code{DWTRM} use \code{inner.reg.sigma} for trimming outliers;\cr
             \code{MED}, \code{RM} and \code{RM} require no scale estimation.\cr 
             The function only returns values for \code{inner.loc.sigma}, \code{inner.reg.sigma},
             \code{outer.loc.sigma} or \code{outer.reg.sigma} if any specified \code{method}
             requires their estimation; otherwise \code{NA}s are returned.}

In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{outer.width}, 
\code{inner.width}, \code{method}, \code{scale}, \code{d}, 
\code{minNonNAs}, \code{online} and \code{extrapolate}. 

Application of the function \code{plot} to an object of class
\code{dw.filter} returns a plot showing the original time series 
with the filtered output. 
}

\references{ 
Bernholt, T., Fried, R., Gather, U., Wegener, I. (2006) 
Modified Repeated Median Filters, 
\emph{Statistics and Computing} \bold{16}, 
177-192.\cr 
(earlier version: \url{http://www.sfb475.uni-dortmund.de/berichte/tr46-04.ps})
 
  
Schettlinger, K., Fried, R., Gather, U. (2006) 
Robust Filters for Intensive Care Monitoring: Beyond the Running Median,
\emph{Biomedizinische Technik} \bold{51}(2), 
49-56.
}

\author{Roland Fried and Karen Schettlinger}

\note{
Missing values are treated by omitting them and thus by reducing 
the corresponding window width. \cr
\code{MED}, \code{RM}, \code{MTM}, \code{TRM} and \code{MRM} 
require at least \code{minNonNAs} non-missing observations in each 
outer window; \code{DWRM}, \code{DWMTM}, \code{DWTRM} and 
\code{DWMRM} require at least \code{minNonNAs} non-missing 
observations in each inner window. Otherwise \code{NA}s are 
returned for \code{level}, \code{slope} and \code{sigma}.
}

\seealso{\code{\link{robreg.filter}}, \code{\link{robust.filter}}, \code{\link{hybrid.filter}}, \code{\link{wrm.filter}}.}



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
y.dw <- dw.filter(y, outer.width=31, inner.width=11, method="all")
# Plot:
plot(y.dw)

# Filtering with trimmed RM and double window TRM only:
y2.dw <- dw.filter(y, outer.width=31, inner.width=11, method=c("TRM","DWTRM"))
plot(y2.dw)
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}
