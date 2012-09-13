\name{wrm.smooth}
\alias{wrm.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Weighted Repeated Median Smoothing}

\description{
A robust smoothing tool using a kernel weighted version of Siegel's (1982) repeated median.
It can be seen as an alternative to local linear L1 regression.
}

\usage{
wrm.smooth(x, y, h, xgrid, weight = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Vector of predictors.}
\item{y}{Vector of responses, needs to have the same length as \code{x}.}
\item{h}{Bandwidth, measured in the same units as the explanatory (independent) variable x: 
(x[0]-\code{h},x[0]+\code{h}) is the range of x-values to be included in the local smoothing at x[0]. Needs to be a positive number.}
\item{xgrid}{Grid on which fitted values are to be evaluated. The default is here to take the input
 values \code{x} for a sample size of at most 100, and \code{seq(min(x),max(x), l=100)} otherwise.}
\item{weight}{Indicates the weight function used.
    \describe{
        \item{\code{weight=1}}{triangular weights}
        \item{\code{weight=2}}{Epanechnikov weights (default)}
        \item{\code{weight=3}}{Gaussian weights}
        \item{\code{weight=4}}{Biweight}
        \item{\code{weight=5}}{Uniform weights}
    }
}
}

\details{Weighted repeated median (WRM) smoothing was suggested in a signal extraction framework by Fried, Einbeck & Gather (2007). 
It combines the advantages of weighted and repeated medians, i.e. the WRM smoother is robust to outliers and adapts to linear trends (through the slope parameter of the repeated median, which is calculated by applying two consecutive weighted medians onto the pairwise slopes).  
The theory and simulations provided by Fried, Einbeck & Gather focus on online signal extraction from time series. Warning: The case of a kernel weighted repeated median smoother for arbitraty non-equidistant design (as implemented here) is not fully investigated yet.

The procedure copes with missing values by omitting them.
}

\value{ 
\code{wrm.smooth} returns an object of class 
\code{wrm.smooth}. An 
object of class \code{wrm.smooth} is a list containing the 
following components:

\item{y}{the original input time series.}
\item{level}{the corresponding signal level extraceted by the weighted Repeated Median filter.}
\item{slope}{the corresponding WRM slope within each time window.}
\item{h}{bandwidth.}
\item{xgrid}{vector with grid values.}
\item{weight}{name of the weight function used for the fit.}
 The function \code{plot} returns a plot 
showing the original data with the smoothed output.} 


\references{
Fried, R., Einbeck, J., Gather, U. (2007), Weighted Repeated Median Smoothing and Filtering, 
\emph{Journal of the American Statistical Association} \bold{102}, 1300-1308. \cr
Preliminary version available as technical report from \url{www.sfb475.uni-dortmund.de/berichte/tr33-05.pdf} 

Siegel, A.F. (1982). Robust regression using repeated medians. \emph{Biometrika} \bold{68}, 242-244.
}
\author{Jochen Einbeck and Roland Fried}
%\note{ further notes}

%% Other sections like Warning with \section{Warning }{....} are possible

\seealso{\code{\link{wrm.filter}}}

\examples{
data(faithful) # Old Faithful Geyser data
faith.WRM <- wrm.smooth(faithful$w, faithful$e,h=4)
plot(faith.WRM)
}

\keyword{smooth}
\keyword{robust}
