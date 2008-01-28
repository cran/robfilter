\name{robust.filter}
\alias{robust.filter}
\title{Robust filter for time series}
\description{
  What do I do?
}
\usage{
robust.filter(y, m=15,
              trend.extraction=c("RM", "MED", "LTS", "LMS"),
              scale.estimation=c("QN", "SN", "MAD", "LSH"),
              outlier.treatment="T",
              shiftd=2, lbound=1, p=0.9)
}
\arguments{
  \item{y}{A numeric vector or time series object to filter}
  \item{m}{Half the window width}
  \item{trend.extraction}{Method used to estimate trend}
  \item{scale.estimation}{Method used to estimate scale}
  \item{outlier.treatment}{Outlier substitution}
  \item{shiftd}{??}
  \item{lbound}{??}
  \item{p}{??}
}
\details{
  Details of how the function works.

  Probably also detailed descriptions of each outlier treatment method.
}
\value{
  Currently undocumented since it is a moving target...
}
\references{Publications?!}
\author{
  Prof. Dr. Roland Fried \email{<fried@statistik.uni-dortmund.de>} and
  Karen Schettlinger \email{<schettlinger@statistik.uni-dortmund.de>}.
}
\seealso{For a detailed description of the trend estimatores see \link{trend.extraction}}
\examples{
data(series1t)
robust.filter(series1t)
}
\keyword{robust}
\keyword{ts}