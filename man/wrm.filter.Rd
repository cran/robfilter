\name{wrm.filter}
\alias{wrm.filter}

\title{Weighted Repeated Median Filters for Univariate Time Series}

\description{
Filtering procedure based on a weighted version of Siegel's (1982) repeated median (RM) and
a moving time window for robust extraction  of low frequency components (the signal)
in the presence of outliers and shifts.
One of several weight functions can be chosen to weight the observations in each time window.
        }

\usage{
wrm.filter(y, width, weight = 1, del = floor(width/2), extrapolate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{a numeric vector or (univariate) time series object. }
\item{width}{a positive integer defining the window width used for fitting. \cr
            If \code{del = floor(width/2)} (default) this needs to be an odd number.}
\item{weight}{Indicates the weight function used.
\describe{
\item{\code{weight=0}:}{equal weighting}
\item{\code{weight=1}:}{triangular weights (default)}
\item{\code{weight=2}:}{Epanechnikov weights}
}}
\item{del}{ a positve integer (smaller than width) specifying the delay of the signal extraction. \cr
 \code{del=0} means online signal extraction without delay. \cr
 Default is \code{del=floor(width/2)}}.
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the edges of the time series. \cr
    If \code{del = floor(width/2)} (default) the extrapolation consists of the 
    fitted values within the first half of the first window and the 
    last half of the last window; if \code{del=0} the 
    extrapolation consists of the all fitted values within the first 
    time window. }
}
\details{For online signal extraction without time delay, weighted repeated median filtering with triangular weights is recommendable 
 in the presence of isolated outliers and abrupt level shifts
since it reacts more quickly to shifts than unweighted repeated median filtering and provides higher
efficiencies. The window width should be chosen based on a guess of the minimal time period 
in which the signal can be approximated by a straight line without abrupt shifts.
Better results can be obtained by increasing the delay, but often minimization of the time delay
itself is one of the objectives so that one prefers del=0.
The procedure replaces missing values by simple extrapolations if these are not within the
first time window used for initialization.

For "offline" situations, it is intuitive to set \code{del} roughly
equal to  \code{width/2}. If the focus is rather on smoothing than on
signal extraction, the Epanechnikov kernel should be used rather than
the triangular kernel. In this case one can also use directly function \code{\link{wrm.smooth}}.   
}

\value{ 
\code{wrm.filter} returns an object of class 
\code{wrm.filter}. An 
object of class \code{wrm.filter} is a list containing the 
following components:

\item{y}{the original input time series.}
\item{level}{the corresponding signal level extracted by the filter.}
\item{slope}{the corresponding slope within each time window.}
\item{del}{the parameter specifying the delay of the signal extraction.}
\item{width}{width of the time window.}
\item{weight}{name of the weight function used for the fit.}
The function \code{plot} returns a plot 
showing the original time series with the filtered output. } 

\references{These filtering procedures are described and investigated in \cr
Fried, R., Einbeck, J., Gather, U. (2007), Weighted Repeated Median Smoothing and Filtering, 
\emph{Journal of the American Statistical Association} \bold{102}, 1300-1308. \cr
Preliminary version available as technical report from \url{www.sfb475.uni-dortmund.de/berichte/tr33-05.pdf} 
}
\author{Roland Fried and Jochen Einbeck}


%% Other sections like Warning with \section{Warning }{....} are possible

\seealso{ \code{\link{dw.filter}}, \code{\link{hybrid.filter}}, \code{\link{wrm.smooth}}
}

\examples{
data(Nile)
nile <- as.numeric(Nile)
obj <- wrm.filter(nile, width=11)
plot(obj)
}
% muss ausfuehrbar sein!

\keyword{robust}
\keyword{ts}
