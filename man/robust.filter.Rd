\name{robust.filter}
\alias{robust.filter}

\title{Robust Filtering Methods for Univariate Time Series}

\description{ 
Procedure for robust (online) extraction of low frequency
components (the \emph{signal}) from a univariate time series with
optional rules for outlier replacement and shift detection.
}

\usage{robust.filter(y, width, trend = "RM", scale = "QN", outlier = "T", 
                        shiftd = 2, wshift = floor(width/2), lbound = 0.1, p = 0.9,
                        adapt = 0, max.width = width, 
                        online = FALSE, extrapolate = TRUE)}

\arguments{
\item{y}{a numeric vector or (univariate) time series object.}
\item{width}{a positive integer defining the window width used for fitting.
             If \code{online=FALSE} (default) this needs to be an odd number.}
\item{trend}{a character string defining the method to be used for robust approximation of the signal 
             within one time window. Possible values are:\cr
    \describe{
    \item{\code{"MED"}:}{Median}
    \item{\code{"RM"}:}{Repeated Median regression (default)}
    \item{\code{"LTS"}:}{Least Trimmed Squares regression}
    \item{\code{"LMS"}:}{Least Median of Squares regression} 
     }
    }
\item{scale}{a character string defining the method to be used for robust estimation of the local 
             variability (within one time window). 
             Possible values are:\cr
    \describe{
    \item{\code{"MAD"}:}{Median absolute deviation about the median}
    \item{\code{"QN"}:}{Rousseeuw's and Croux' (1993) \eqn{Q_n} scale estimator (default)}
    \item{\code{"SN"}:}{Rousseeuw's and Croux' (1993) \eqn{S_n} scale estimator}
    \item{\code{"LSH"}:}{Length of the shortest half}   
     }
}
\item{outlier}{a single character defining the rule to be used for outlier detection and outlier treatment. 
    Observations deviating more than \eqn{d\cdot \hat{\sigma}_t}{d* \sigma_t} 
    from the current level approximation \eqn{\hat{\mu}_t}{\mu_t} 
    are replaced by \eqn{\hat{\mu}_t\pm k\hat{\sigma}_t}{\mu_t +/- k \sigma_t}
    where \eqn{\hat{\sigma}_t}{\sigma_t} denotes the current scale estimate. \cr
    Possible values are:\cr
    \describe{
    \item{\code{"T"}:}{Replace ('trim') large outliers detected by a \eqn{3\sigma}-rule 
        (\eqn{d=3}{d=3}) by the current level estimate (\eqn{k=0}{k=0}). (default)}
    \item{\code{"L"}:}{Shrink large outliers (\eqn{d=3}{d=3}) strongly 
        towards the current level estimate (\eqn{k=1}{k=1}).}
    \item{\code{"M"}:}{Shrink large and moderatly sized outliers (\eqn{d=2}{d=2}) strongly 
        towards the current level estimate (\eqn{k=1}{k=1}).}
    \item{\code{"W"}:}{Shrink large and moderatly sized outliers (\eqn{d=2}{d=2}) 
        towards the current level estimate (\eqn{k=2}{k=2}).}
    }
    \code{W} is the most efficient, \code{T} the most robust method (which should ideally 
        be combined with a suitable value of \code{lbound}).  }
\item{shiftd}{a positive numeric value defining the factor the current scale estimate is multiplied 
    with for shift detection. Default is \code{shiftd}=2 
    corresponding to a \eqn{2\sigma}{2\sigma} rule for shift detection.}
\item{wshift}{a positive integer specifying the number of the most recent observations used for shift detection
   (regulates therefore also the delay of shift detection). Only used
   in the \code{online} mode; should be less than half the (minimal)
   window width then. In the offline mode (\code{online=FALSE}, default), shift
   detection is based on the right half of the time window, i.e. \code{wshift=floor(width/2)} (default).}
\item{lbound}{a positive real value specifying an optional lower bound for the scale to prevent 
    the scale estimate from reaching zero (implosion). } 
\item{p}{a fraction \eqn{\in [2/3,1]}{in [2/3,1]} of observations 
    for additional rules in case of only two or three different values 
    within one window.\cr
    If 100 percent of the observations within one window take on 
    only two different values, the current level is estimated by the 
    mean of these values regardless of the \code{trend} 
    specification. In case of three differing values the median is 
    taken as the current level estimate.}
\item{adapt}{a numeric value  defining the fraction which regulates the adaption of the
   moving window width. \code{adapt} can be either 0 or a value \eqn{\in [0.6,1]}{in [0.6,1]} . 
   \code{adapt = 0} means that a fixed window width is used. 
   Otherwise, the window width is reduced whenever more than a fraction of 
   \code{adapt} \eqn{\in [0.6,1]}{in [0.6,1]} of the residuals in a
   certain part of the current time window are all positive or all
   negative.}
\item{max.width}{a positive integer (\code{>= width}) specifying the maximal width of the time window.\cr
   \code{width} specifies the minimal (and also the initial) width.}
\item{online}{a logical indicating whether the current level and 
    scale estimates are evaluated at the most recent time 
    within each window (\code{TRUE}) or centered within the window 
    (\code{FALSE}). \code{online=FALSE} (default) requires an odd 
    \code{width} for the window and means a time delay of 
    (\code{width}+1)/2 time units.}
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series. \cr
    If \code{online=FALSE} the extrapolation consists of the 
    fitted values within the first half of the first window and the 
    last half of the last window; if \code{online=TRUE} the 
    extrapolation consists of all fitted values within the first 
    time window. } 
}

\details{\code{robust.filter} works by applying the methods 
specified by \code{trend} and \code{scale} to a moving time 
window of length \code{width}.

Before moving the time window, it is checked whether the next 
(incoming) observation is considered an 'outlier' by applying the 
rule specified by \code{outlier}. Therefore, the trend in the 
current time window is extrapolated to the next point in time and 
the residual of the incoming observation is standardised by the 
current scale estimate.

After moving the time window, it can be tested whether a level 
shift has occurred within the window: If more than half of the 
residuals in the right part of the window are larger than 
\code{shiftd}\eqn{\cdot\sigma_t}{*\sigma_t}, a shift is detected and 
appropriate actions are taken. In
the \code{online} mode, the number of the rightmost residuals can be
chosen by \code{wshift} to regulate the resistance of the detection
rule against outliers, its power and the time delay of detection.

A more detailed description of the filter can be found in Fried 
(2004). The adaption of the window width is described by Gather and
Fried (2004). For more explanations on shift detection, see Fried
and Gather (2007).}


\value{ 
\code{robust.filter} returns an object of class \code{robust.filter}.  
An object of class \code{robust.filter} is a list containing the 
following components:

\item{level}{a numeric vector containing the signal level extracted by the (regression) filter specified by \code{trend}, \code{scale} and \code{outlier}.}
\item{slope}{a numeric vector containing the corresponding slope within each time window.}
\item{sigma}{a numeric vector containing the corresponding scale within each time window.}
\item{ol}{an outlier indicator. 
      0: no outlier, +1: positive outlier, -1: negative outlier}
\item{level.shift}{a level shift indicator.
                   0: no level shift, t: positive level shift detected at processing
                   time t, -t: negative level shift detected at processing time t (the
                   position in the vector gives an estimate of the point in time before
                   which the shift has occurred).
                   }

In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{width}, \code{trend}, 
\code{scale}, \code{outlier}, \code{shiftd}, 
\code{wshift}, \code{lbound}, 
\code{p}, \code{adapt}, \code{max.width}, \code{online} and \code{extrapolate}. 

Application of the function \code{plot} to an object of class
\code{robust.filter} returns a plot showing the original time series 
with the filtered output. 
} 

\note{
Missing values have to be replaced or removed from the time series 
before applying \code{robust.filter}.
}

\references{ 
Fried, R. (2004), Robust Filtering of Time Series with 
Trends, \emph{Journal of Nonparametric Statistics} \bold{16}, 
313-328.\cr
(earlier version: \url{http://www.sfb475.uni-dortmund.de/berichte/tr30-03.ps})

Fried, R., Gather, U. (2007), On Rank Tests for Shift Detection in Time Series, 
\emph{Computational Statistics and Data Analysis, Special Issue on Machine Learning and Robust Data Mining} \bold{52}, 221-233.\cr
(earlier version: \url{http://www.sfb475.uni-dortmund.de/berichte/tr48-06.pdf}) 

Gather, U., Fried, R. (2004), Methods and Algorithms for Robust Filtering, 
\emph{COMPSTAT 2004: Proceedings in Computational Statistics}, J. Antoch (eds.), Physika-Verlag, Heidelberg, 159-170. \cr 

Schettlinger, K., Fried, R., Gather, U. (2006) Robust Filters
for Intensive Care Monitoring: Beyond the Running Median,
\emph{Biomedizinische Technik} \bold{51}(2), 49-56.
}

\author{Roland Fried and Karen Schettlinger}

\seealso{\code{\link{robreg.filter}}, \code{\link{hybrid.filter}}, \code{\link{dw.filter}}, \code{\link{wrm.filter}}.}

\examples{
# Generate random time series:
y <- cumsum(runif(500)) - .5*(1:500)
# Add jumps:
y[200:500] <- y[200:500] + 5
y[400:500] <- y[400:500] - 7
# Add noise:
n <- sample(1:500, 30)
y[n] <- y[n] + rnorm(30)

# Delayed Filtering of the time series with window width 23:
y.rf <- robust.filter(y, width=23)
# Plot:
plot(y.rf)

# Delayed Filtering with different settings and fixed window width 31:
y.rf2 <- robust.filter(y, width=31, trend="LMS", scale="QN", outlier="W")
plot(y.rf2)

# Online Filtering with fixed window width 24:
y.rf3 <- robust.filter(y, width=24, online=TRUE)
plot(y.rf3)

# Delayed Filtering with adaptive window width (minimal width 11, maximal width 51):
y.rf4 <- robust.filter(y, width=11, adapt=0.7, max.width=51)
plot(y.rf4)
}



\keyword{robust}
\keyword{smooth}
\keyword{ts}







%References concerning the regression and scale estimation methods:
%
%Gr?bel, R. (1988) The Length of the Shorth,
%\emph{The Annals of Statistics} \bold{16}, 619-628.
%
%Hampel, F.R. (1975) Beyond Location Parameters: Robust Concepts
%and Methods, \emph{Bulletin of the International Statistical
%Institute} \bold{46}, 375-382.
%
%Rousseeuw, P.J.(1983) Multivariate Estimation with High Breakdown
%Point, \emph{Proceedings of the 4th Pannonian Symposium on Mathematical
%Statistics and Probability} Vol. B, Dordrecht: Reidel.
%
%Rousseeuw, P.J. (1984) Least Median of Squares Regression,
%\emph{Journal of the American Statistical Association} \bold{79},
%871-880.
%Rousseeuw, P.J., Croux, C. (1993) Alternatives to the Median Absolute Deviation,
%\emph{Journal of the American Statistical Association} \bold{88},
%1273-1283.
%
%Siegel, A.F. (1982) Robust Regression Using Repeated
%Medians, \emph{Biometrika} \bold{69}, 242-244.
