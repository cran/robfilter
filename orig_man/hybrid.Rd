\name{hybrid}
\alias{hybrid}

\title{function for edge preserving signal extraction from time series
  (smoothing)}

\description{Different filtering procedures using moving centered time
windows for robust delayed extraction of low frequency components ('the
signal') in the presence of outliers, level shifts and local
extremes. Several one-sided half-window estimates ('subfilters',
specified in brackets below) are calculated in a first step and the
median of these and a further estimate of the signal value in the window
center is used for the final estimate.  Both location based and
regression based methods are available, the former applying means or
medians and the idea of a locally constant signal value, the latter
using ordinary least squares (LS) regression or Siegel's (1982) repeated
median (RM) and the idea of an underlying locally linear trend.}

\usage{ 
  hybrid(y, n, methods="all", plot=FALSE, plot.methods)
}

\arguments{
  \item{y}{One dimensional data vector}
  \item{n}{Window width. Needs to be an impair positive integer.}
  \item{methods}{The method(s) used for the analyses. Possible
  values are "MED", "RM", "MEAN", "FMH", "PFMH", "CFMH", "RMH", "PRMH",
  "CRMH", "RMMH", "PRMMH" and "CRMMH". For a detailed description see
  the Methods section.

  Example: An analysis by the ordinary RM and the PRMH filter would be
  achieved by setting methods=c"RM", "PRMH")

  Default: methods="all"}
  \item{plot}{Plots some or all signal estimates obtained by the
    different methods in a time series plot if set to TRUE.}
  \item{plot.methods}{The estimates to be shown in the time series
    plot. Possible values are the same as in 'methods'.}
}

\section{Methods}{
  \itemize{
    \item{MED}{Ordinary median filter (one step only)}
    \item{MEAN}{Ordinary moving average (one step only)}
    \item{RM}{Ordinary repeated median filter (one step only)}
    \item{FMH}{FIR Median Hybrid (half-window averages and central
      observation)} 
    \item{PFMH}{Predictive FMH (half-window LS regression and central
      observation)} 
    \item{CFMH}{Combined FMH (half-window averages and LS regression and
      central observation)} 
    \item{RMH}{RM Hybrid (half-window medians and central observation)}
    \item{PRMH}{Predictive RMH (half-window RM regression and central
      observation)} 
    \item{CRMH}{Combined FMH (half-window medians and RM regression and
      central observation)} 
    \item{RMMH}{RM Hybrid (half-window medians and full-window median)} 
    \item{PRMMH}{Predictive RMH (half-window RM regression and
      full-window median)} 
    \item{CRMMH}{Combined FMH (half-window medians and RM regression and
      full-window median)}
  }
}

\details{The method should be chosen based on an a-priori guess of the
underlying signal and the data quality.  Location based (MED / MEAN /
FMH / RMH / RMMH) methods are recommended in case of a locally
(piecewise) constant signal, regression based (MEAN / RM / PFMH / PRMH /
PRMMH) in case of locally linear monotone trends.  The combined filters
(CFMH / CRMH / CRMMH) can be seen as a compromise, but are
computationally somewhat more expensive and may be inferior to the
predictive filters during steep trends.

The median and RM based approaches are robust alternatives to the more
efficient (in Gaussian samples) mean and least squares methods.  The
hybrid filters preserve shifts and local extremes much better than the
others (MED / MEAN / RM) for the price of decreased robustness and / or
Gaussian efficiency.}

\value{A matrix with one column per selected signal estimation method
  containing the smoothed data.}


\references{
  Fried, R., Bernholt, T., Gather, U. (2006), \dQuote{Repeated Median
  and Hybrid Filters}, \emph{Computational Statistics \& Data Analysis}, \bold{50}, 2313-2338;
  preliminary version available as technical report from
  \url{http://www.statistik.uni-dortmund.de/fixme}. 
}

\author{Roland Fried}

\seealso{ \code{\link{double.window}} for more efficient and more robust
  delayed filtering procedures also providing some edge preservation}

\examples{
data(serie1t)
hybrid(serie1t$y, n=31, methods="all", plot=TRUE, plot.methods=c("RM", "PFMH", "PRMH"))
}

\concept{FMH}
\concept{CFMH}
\concept{PFMH}
\concept{RMH}
\concept{CRMH}
\concept{PRMH}
\concept{FMH}
\concept{CFMH}
\concept{running median}
\concept{mean}
\concept{moving average}
\concept{RM}
\concept{robust} 
\concept{level shifts}
\concept{local extremes}
\concept{edges}
\keyword{ts}
