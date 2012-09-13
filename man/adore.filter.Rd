\name{adore.filter}
\alias{adore.filter}
\title{A Robust Adaptive Online Repeated Median Filter for Univariate Time Series}

\description{ 
Procedure for robust online extraction of low frequency
components (the \emph{signal}) from a univariate time series by a moving window technique with
adaptive window width selection (ADaptive Online REpeated median FILTER).
}

\usage{
adore.filter(y,
             p.test=15, minNonNAs=5,
             min.width=10, max.width=200,
             width.search="geometric",
             rtr=2, extrapolate=FALSE,
             calc.qn=FALSE, sign.level=0.1)
}

\arguments{
\item{y}{a numeric vector or (univariate) time series object.}
\item{p.test}{defines the number of most recent Repeated Median residuals within each window used to test the goodness of fit of the online signal level. \cr
It can be either a value in (0.25, 0.3, 0.5), meaning that \code{floor(p.test*width)} residuals are considered for the goodness of fit test, where \code{width} is the currently used window width, or it can also be a positive integer \code{>= 5} specifying a fixed number of most recent residuals (default).\cr
If the number of residuals considered for the test exceeds \code{width/2}, the procedure sets it to \code{floor(width/2)}, if it is smaller than five, the number is set to five.}
\item{minNonNAs}{a positive integer \code{>= 5} defining the minimum number of non-missing observations within one window which is required for a \sQuote{sensible} estimation.}
\item{min.width}{a positive integer \code{>= 5} specifying the minimal window width.}
\item{max.width}{a positive integer \code{>= min.width} specifying the maximal window width.}
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
\item{rtr}{a value in 0, 1, 2 specifying whether a 'restrict to range' rule should be applied.
    \describe{
    \item{\code{rtr=0}}{The estimated signal level consists of the last fitted value of a Repeated Median regression fit 
        within a time window of adequate width.}
    \item{\code{rtr=1}}{The level estimation is restricted to the range of the observations within each time window.}
    \item{\code{rtr=2} (default)}{The level estimation is restricted to the range of the most recent observations 
        (specified by \code{p.test}) i.e., to the range of the observations which are used to evaluate the goodness of fit.}
    }
}
\item{extrapolate}{a logical indicating whether the level 
    estimations should be extrapolated to the beginning of the time series. The 
    extrapolation consists of all fitted values within the first 
    time window. }
\item{calc.qn}{a logical indicating whether the Qn scale (Rousseeuw, Croux, 1993) should also be calculated along with the signal level 
as an estimate of the standard deviation in each window.
Here, the \code{Qn} command from the \code{robustbase} library is applied with the built-in finite sample correction.
}
\item{sign.level}{significance level of the test procedure; must be a value in \eqn{(0,0.5)}.}
}


\details{The \code{adore.filter} works by applying Repeated Median (RM)
regression (Siegel, 1982) to a moving time window with a length varying between 
\code{min.width} and \code{max.width}. 

For each point in time, the window width is adapted to the current data situation by a goodness of fit test for
the most recent signal level estimation. 
The test uses the absolute value of the sum of the RM residuals in the subset specified by \code{p.test}.
The critical value for the test decision corresponds to a slightly modified 0.95-quantile of the 
distribution of the test statistic and is stored in the data set \code{critvals}.

A more detailed description of the filter can be found in 
Schettlinger, Fried, Gather (2010).} 


\value{ 
\code{adore.filter} returns an object of class \code{adore.filter}.  
An object of class \code{adore.filter} is a list containing the 
following components:

\item{level}{a numeric vector containing the signal level extracted by the RM filter 
with adaptive window width.}
\item{slope}{a numeric vector containing the corresponding slope within each time window.}
\item{width}{a numeric vector containing the corresponding window width used for the level and slope estimations.}

\item{level.list}{a list which contains with as many elements as the length of the input time series.
If at time \code{t}, the window width was not reduced, the entry \code{level.list[[t]]} simply corresponds to \code{level[t]}.
However, if more than one iteration took place, \code{level.list[[t]]} is a vector which contains all level estimations which were
evaluated until the final estimate \code{mu[t]} passed the goodness of fit test and was stored.
%Note: The last value in the vector \code{level.list[[t]]} corresponds to the signal estimation
%without the 'restrict to observational range' rule and hence might differ from \code{mu[t]}.
}
\item{slope.list}{a list containing the slope estimations corresponding to the values in \code{level.list}.
}
\item{width.list}{a list containing the window widths used for the estimations in \code{level.list} and \code{slope.list}.
}
\item{sigma}{a numeric vector containing the corresponding scale within each time window estimated by the robust Qn estimator
(only calculated if \code{calc.qn = TRUE}, else \code{sigma} does not exist).}

In addition, the original input time series is returned as list 
member \code{y}, and the settings used for the analysis are 
returned as the list members \code{min.width}, \code{max.width}, \code{width.search},
\code{p.test}, \code{minNonNAs}, \code{rtr}, \code{extrapolate}, and \code{calc.qn}. 

Application of the function \code{plot} to an object of class
\code{aoRM} returns a plot showing the original time series 
with the filtered output. 
} 

%\note{
%}

\references{ 
Rousseeuw, P. J., Croux, C. (1993)
Alternatives to the Median Absolute Deviation,
\emph{Journal of the American Statistical Association} \bold{88},
1273-1283.

Schettlinger, K., Fried, R., Gather, U. (2010) 
Real Time Signal Processing by Adaptive Repeated Median Filters,
\emph{International Journal of Adaptive Control and Signal Processing} \bold{24}(5), 346-362.\cr

Siegel, A.F. (1982)
Robust Regression Using Repeated Medians,
\emph{Biometrika} \bold{69} (1), 
242-244.\cr
}

\author{Karen Schettlinger}

\seealso{\code{\link{robreg.filter}}, \code{\link{wrm.filter}}, \code{\link{madore.filter}}, \code{\link{scarm.filter}}.}

\examples{
# # # # # # # # # #
# Short and noise-free time series
series <- c(rep(0,30),rep(10,30),seq(10,5,length=20),seq(5,15,length=20))

# Adaptive online signal extraction without & with 'restrict to range' rule
t.without.rtr <- adore.filter(series, rtr=0)
plot(t.without.rtr)
t.with.rtr1 <- adore.filter(series, rtr=1)
lines(t.with.rtr1$level, col="blue")
t.with.rtr2 <- adore.filter(series)
lines(t.with.rtr2$level, col="green3",lty=2)
legend("top",c("Signal with rtr=1","Signal with rtr=2"),col=c("blue","green3"),lty=c(1,2),bty="n")

# # # # # # # # # #
# Short and noise-free time series + 1 outlier
ol.series <- series
ol.series[63] <- 3

# Adaptive online signal extraction without & with 'restrict to range' rule
t.without.rtr <- adore.filter(ol.series, rtr=0)
plot(t.without.rtr)
t.with.rtr1 <- adore.filter(ol.series, rtr=1)
lines(t.with.rtr1$level, col="blue")
t.with.rtr2 <- adore.filter(ol.series)
lines(t.with.rtr2$level, col="green3",lty=2)
legend("top",c("Signal with rtr=1","Signal with rtr=2"),col=c("blue","green3"),lty=c(1,2),bty="n")

# # # # # # # # # #
# Noisy time series with level shifts, trend changes and shifts in the scale of the error term
true.signal  <- c(rep(0,150),rep(10,150),seq(10,5,length=100),seq(5,15,length=100))
series2      <- true.signal + c(rnorm(250,sd=1), rnorm(200,sd=3), rnorm(50,sd=1))

# Adaptive online signal extraction with additional Qn scale estimation
s2 <- adore.filter(series2, calc.qn=TRUE)
par(mfrow=c(3,1))
plot(s2)
plot(s2$sigma,type="l",main="Corresponding Qn Scale Estimation",ylab="sigma",xlab="time")
lines(c(rep(1,250),rep(3,200),rep(1,150)),col="grey")
legend("topleft",c("True scale","Qn"),lty=c(1,1),col=c("grey","black"),bty="n")
plot(s2$width,type="l",main="Corresponding Window Width",ylab="width",xlab="time")

}



\keyword{robust}
\keyword{smooth}
\keyword{ts}
