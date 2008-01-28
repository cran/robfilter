\name{double.window}
\alias{double.window}
\alias{double.window.online}
\title{Online / Offline signal extraction from time series}

%% FIXME: Zu lang. Lieber nur 1 bis 2 Saetze und den Rest in die
%% \details{} packen.
\description{ Different filtering procedures using moving centered time
  windows for robust delayed extraction of low frequency components
  (\sQuote{the signal}) in the presence of outliers and shifts. In a
  first step, a short inner window is used for calculating an initial
  estimate of the local signal value and a robust estimate of the local
  variance.  Based on these initial estimates, strongly deviating
  observations (outliers) are trimmed from a possibly longer outer
  window, and the final signal estimate is calculated from the remaining
  observations in a second step. Both location based and regression
  based methods are available, the former applying the median and the
  idea of a locally constant signal value, the latter using Siegel's
  (1982) repeated median (RM) and the idea of an underlying locally
  linear trend.}

\usage{
double.window(y, outer.n, inner.n, d = 2, scale = "MAD",
              mad.corr.file = NULL, methods = "all", plot = FALSE,
              plot.methods = methods, online = FALSE)

double.window.online(...)
}
\arguments{
  \item{y}{One dimensional data vector}

  \item{outer.n}{Width of the outer window used for the final
    estimate. Needs to be an impair positive integer.}

  \item{inner.n}{Width of the inner window used for the initial
    estimate. Needs to be an impair positive integer no larger than
    \code{outer.n}.}

  \item{d}{Factor the scale estimate is multiplied with in the trimming
    boundaries.  Default: \code{d = 2} meaning a \eqn{2\sigma}{2 sigma} rule.}
  
  \item{scale}{Scale estimator to be used for the estimation of the
    local variance.  Possible values: "MAD" (standard), Rousseeuw and
    Croux's (1993) "SN" and "QN".}

  \item{mad.corr.file}{File with correction factors to make the MAD
    scale estimate unbiased in finite samples.  Default: Internal
    simulation of these correction factors for inner.n and outer.n}

  \item{methods}{A vector of the method(s) used for signal estimattion.
    Possible values are "MED", "RM", "MTM", "TRM", "MRM", "DWRM",
    "DWMTM", "DWTRM" and "DWMRM". For a thorough description see the
    Methods section.}

  \item{plot}{Plots some or all signal estimates obtained by the
    different methods in a time series plot if set to TRUE.}

  \item{plot.methods}{The estimates to be shown in the time series
    plot. Possible values are the same as in 'methods'. Only those
    selected in the \code{methods} parameter may be named.}

  \item{online}{Use a window that only uses observations from the past.} 
}

\section{Methods}{
  The following methods for signal estimation are currently implemented:
  %% FIXME: first step = outer window und second step = inner window?
  \itemize{
    \item{MED}{Ordinary median filter (one step only)}
    \item{RM}{Ordinary repeated median filter (one step only)}
    \item{MTM}{Modified trimmed mean (median in the first / mean in the
      second step)}
    \item{TRM}{Trimming with repeated median (repeated median in the
      first and trimmed least squares in the second step)} 
    \item{MRM}{Modified repeated median (RM in the first and the second
      step)}
    \item{DWRM}{Double window RM (only RM slope in the first step, and
      median of the trend corrected observations in the second step
      without trimming)}
    \item{DWMTM}{Double window MTM (different window widths in the two
      steps)} 
    \item{DWTRM}{Double window TRM (different window widths in the two
      steps)} 
    \item{DWMRM}{Double window MRM (different window widths in the two
      steps)}
  }
}

\details{The method should be chosen based on an a-priori guess of the
underlying signal and the data quality: Location based (MED / MTM) is
recommended in case of a locally (piecewise) constant signal, regression
based (RM / DWRM / TRM / MRM) in case of locally linear, monotone
trends.  No big differences have been reported between TRM and MRM, so
preference might be given to the quicker and somewhat more efficient TRM
option. DWRM is the quickest of all regression based methods and
performs better than the ordinary RM at shifts, but it is the least
robust and least efficient method.

The MAD is the classical highly robust choice for the estimation of the
variance. The SN is a somewhat more efficient and almost equally robust
alternative, while the QN is much more efficient in case of not very
small window widths and performs very well at the occurrence of shifts.

The inner window width should be chosen at least twice the length of
outlier patches in the series to be ignored in case of the location
based, and at least three times this length in case of the regression
based methods. Otherwise the methods can be severely influenced by
outlier patches.  The outer window width can then be chosen rather large
to increase the efficiency of the final estimate, provided that it is
smaller than the time in between subsequent level shifts.

The factor d with which the scale is multiplied for fixing the trimming
boundaries can be chosen similar to classical rules for detecting
unusual observations in a Gaussian sample. Choosing d=3 instead of d=2
increases efficiency, but decreases robustness; d=2.5 might be seen as a
compromise.}

\value{
  \item{signal}{Matrix (?) with the different signal estimates obtained
    from the choosen methods in its columns.} 
  \item{slope}{Matrix (?) with the different estimates of the signal
    slope in its columns.} 
}

\references{
  Bernholt, T., Fried, R., Gather, U., Wegener, I. (2006), \dQuote{Modified
    Repeated Median Filters}, \emph{Statistics and Computing}, \bold{16},
  177-192; Preliminary version available as technical report from
  \url{http://www.statistik.uni-dortmund.de/fixme} 
}

\author{Roland Fried}

\seealso{\code{\link{hybrid}}} 

\examples{
data(serie1t)
double.window(serie1t$y,
              outer.n=25, inner.n=15,
              d=2.5, scale="QN",
              methods = c("RM", "TRM"))
}

\concept{time delay}
\concept{DWRM}
\concept{MAD} 
\concept{running median}
\concept{MTM} 
\concept{Qn}   
\concept{robust} 
\concept{Repeated median (RM)} 
\concept{Sn} 
\concept{level shifts}
\concept{robust smoothing}
\keyword{ts}

