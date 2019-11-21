\encoding{latin1}

\name{madore.filter}

\alias{madore.filter}

\title{A multivariate adaptive online repeated median filter}

\description{Procedure for robust signal extraction from a multivariate time series by a moving window technique with adaptive window width selection (\emph{m}ultivariate \emph{ad}aptive \emph{o}nline \emph{re}peated median filter). The window width adaption is based on the univariate \code{adore.filter}.}

\usage{
madore.filter(Y, byrow=FALSE,
              min.width=10, max.width=200,
              test.sample.size=min.width/2,
              width.search="geometric",
              rtr.size=min.width, sign.level=0.1,
              NA.sample.size=min.width, minNonNAs=min.width/2)
}

\arguments{\item{Y}{a numeric matrix or (multivariate) time series object.}

\item{byrow}{logical. If \code{FALSE} (the default), the filtering is done by columns, otherwise the filtering is done by rows.}

\item{min.width}{a positive integer \eqn{\ge 10} specifying the minimal width of the moving time window.}

\item{max.width}{a positive integer \eqn{\ge} \code{min.width} specifying the maximal width of the moving time window. If \code{min.width} = \code{max.width}, the window width is fixed.}

\item{test.sample.size}{a positive integer in [5, \code{min.width}] defining a test window of the rightmost \code{test.sample.size} time points within the current time window. The \emph{Repeated Median} (RM) regression residuals within the test window are used for a goodness of fit test (see \code{adore.filter}) for finding an adequate window width. For more details about the test, see Schettlinger, Fried, Gather (2010).}

\item{width.search}{a character string defining the search algorithm used for finding an adequate window width at each point in time.
    \describe{
    \item{\code{"linear"}}{The linear search always results in the largest window width possible and hence yields 
        the smoothest online signal. 
        However, if sudden changes (like level shifts) appear in the signal it requires a lot of computation time
        and thus, an increased variability of the extracted signal may be observed.
        }
    \item{\code{"binary"}}{The binary search is recommended if it can be expected that the window width needs 
        to be reduced drastically from a large to a very small value at certain times (for example at level shifts or trend changes).
        However, it may not always result in the largest possible window width.
        }
    \item{\code{"geometric"} (default)}{The geometric search is as fast as the binary search
        but it puts more weight on large window widths. 
        It offers a good compromise between the linear and the binary search (computation time vs. smooth output signal).
        }
    }
}

\item{rtr.size}{a non-negative integer specifying the size of a subset of the most recent observations within each window. The signal estimation is restricted to the range of the observations within this subset.}

\item{sign.level}{the level of significance for the goodness of fit test (see \code{adore.filter}) for finding an adequate window width. For more details about the test, see Schettlinger, Fried, Gather (2010).}

\item{NA.sample.size}{a positive integer in [10, \code{min.width}] specifying the size of a subset of the most recent observations within each window. See \code{minNonNAs}.}

\item{minNonNAs}{a positive integer in [5, \code{NA.sample.size}]. If a variable does not offer at least \code{minNonNAs} non-missing observations within the subset specified by \code{NA.sample.size}, the signal is not estimated for this variable at this time point \eqn{t}.}
}

\details{The \code{madore.filter} is based on \emph{Repeated Median} regression (Siegel, 1982) in moving time windows and serves for separating signals from noise and outliers in multivariate time series. At each time point \eqn{t} the test procedure of the \emph{adaptive online Repeated Median} filter (Schettlinger, Fried, Gather, 2010) is used to determine an appropriate window width \eqn{n(t)} in [\code{min.width}, \code{max.width}]. Then the signal vector at time \eqn{t} is estimated within the time window \eqn{(t-n(t)+1,\ldots,t)} by a slight modification of the multivariate \emph{Trimmed Repeated Median-Least Squares} regression (Lanius, Gather, 2010). A more detailed description of the \code{madore.filter} can be found in Borowski, Schettlinger, Gather (2009).}

\value{
\code{madore.filter} returns an object of class \code{madore.filter}. An object of class \code{madore.filter} is a list containing the following components:
\item{signals}{a matrix containing the estimated signal vectors at each time point \eqn{t}.}
\item{widths}{a matrix containing the individual window widths of each variable at each time point \eqn{t}.}
\item{overall.width}{a vector containing the overall window widths at each time point \eqn{t}.}
In addition, the original input data is returned as list member \code{Y}, and the settings used for the analysis are returned as the list members \code{byrow}, \code{min.width}, \code{max.width}, \code{start.width}, \code{test.sample.size}, \code{width.search}, \code{rtr.size}, \code{extr.delay}, \code{NA.sample.size}, and \code{minNonNAs}. Application of the function \code{plot} to an object of class \code{madore.filter} returns a plot showing the original multivariate time series with the filtered output.
}

\references{
Borowski, M., Schettlinger, K., Gather, U. (2009)
Multivariate Real Time Signal Extraction by a Robust Adaptive Regression Filter, \emph{Communications in Statistics - Simulation and Computation} \bold{38}, 426-440.\cr

Lanius, V., Gather, U. (2010)
Robust Online Signal Extraction from Multivariate Time Series, 
\emph{Computational Statistics and Data Analysis} \bold{54}(4), 966-975.\cr

Schettlinger, K., Fried, R., Gather, U. (2010) 
Real Time Signal Processing by Adaptive Repeated Median Filters,
\emph{International Journal of Adaptive Control and Signal Processing} \bold{24}(5), 346-362.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians, 
\emph{Biometrika} \bold{69}(1), 242-244.\cr
}

\author{Matthias Borowski}

\seealso{\code{\link{robreg.filter}}, \code{\link{adore.filter}}, \code{\link{scarm.filter}}, \code{\link{mscarm.filter}}.}

\examples{
\dontrun{
data(multi.ts)
extr <- madore.filter(multi.ts)
plot(extr)
}   
}

\keyword{robust}
\keyword{smooth}
\keyword{ts}
\keyword{multivariate}
\keyword{regression}
