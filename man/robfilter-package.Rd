\name{robfilter-package}
\alias{robfilter-package}
\alias{robfilter}
\docType{package}

\title{robfilter - Robust time series filters}

\description{
A set of functions for signal extraction from univariate and multivariate time series based on concepts from robust statistics.
}

\details{
\tabular{ll}{
Package: \tab robfilter\cr
Version: \tab 3.0\cr
Date: \tab 2011-07-22\cr
Depends: \tab R (>= 2.5.0), robustbase, MASS\cr
License: \tab GPL (>= 2)\cr
URL: \tab 
http://www.statistik.tu-dortmund.de/fried.html\cr
LazyData: \tab yes\cr
Repository: \tab CRAN\cr
Packaged: \tab 2011-07-22 11:20:52 UTC; borowski\cr
Date/Publication: \tab 2011-07-22 11:40:35\cr
Built: \tab R 2.13.1; i386-pc-mingw32; 2011-07-22 07:44:45 UTC; windows\cr
}

Index:
\preformatted{
adore.filter            A Robust Adaptive Online Repeated Median Filter
                        for Univariate Time Series
const.Q                 Correction factors to achieve unbiasedness of 
                        the regression-free Q scale estimator                       
const                   Correction factors to achieve unbiasedness of
                        the Qn scale estimator
critvals                Critical Values for the RM Goodness of Fit Test
dfs                     Degrees of freedom for the SCARM test statistic
dr.filter               Deepest Regression (DR) filter
dw.filter               Robust Double Window Filtering Methods for
                        Univariate Time Series
hybrid.filter           Robust Hybrid Filtering Methods for Univariate
                        Time Series
lms.filter              Least Median of Squares (LMS) filter
lqd.filter              Least Quartile Difference (LQD) filter
lts.filter              Least Trimmed Squares (LTS) filter
madore.filter           A Robust Adaptive Online Filter for
                        Multivariate Time Series
med.filter              Median (MED) filter
multi.ts                Generated Multivariate Time Series
rm.filter               Repeated Median (RM) filter
robreg.filter           Robust Regression Filters for Univariate Time
                        Series
robust.filter           Robust Filtering Methods for Univariate Time
                        Series
scarm.filter            SCARM (Slope Comparing Adaptive Repeated Median)                        
var.n                   Variance of the Repeated Median slope estimator
wrm.filter              Weighted Repeated Median Filters for Univariate
                        Time Series
wrm.smooth              Weighted Repeated Median Smoothing
}
}

\author{
Roland Fried <fried@statistik.tu-dortmund.de>, Karen Schettlinger <schettlinger@statistik.tu-dortmund.de> and Matthias Borowski <borowski@statistik.tu-dortmund.de>
Maintainer: Matthias Borowski <borowski@statistik.tu-dortmund.de> and Roland Fried <fried@statistik.tu-dortmund.de>
}

\keyword{robfilter}
\keyword{filter}
\keyword{robust}
\keyword{time series}