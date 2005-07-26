\name{ckFGN0}
\alias{ckFGN0}
\title{Covariances of a Fractional Gaussian Process}
\description{
  Compute the Autocovariances of a fractional Gaussian process
}
\usage{
ckFGN0(n, H)
}
\arguments{
  \item{n}{sample size (length of time series).}
  \item{H}{self-similarity (`Hurst') parameter.}
}
\value{
  numeric vector of covariances upto lag n-1.
}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\seealso{\code{\link{ckARMA0}} which does the same for a fractional
  ARIMA process.}
\examples{
str(C.8 <- ckFGN0(50, H = 0.8))
plot(0:49, C.8, type = "h", ylim = 0:1)
plot(0:49, C.8, type = "h", log = "xy",
     main = "Log-Log  ACF for frac.GaussNoise(H = 0.8)")
}
\keyword{ts}
