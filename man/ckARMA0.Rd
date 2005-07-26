\name{ckARMA0}
\alias{ckARMA0}
\title{Covariances of a Fractional ARIMA(0,d,0) Process}
\description{
  Compute the Autocovariances of a fractional ARIMA(0,d,0) process (d =
  H - 1/2).
}
\usage{
ckARMA0(n, H)
}
\arguments{
  \item{n}{sample size (length of time series).}
  \item{H}{self-similarity (`Hurst') parameter.}
}
\details{
  The theoretical formula,
  \deqn{C(k) = (-1)^k  \Gamma(1-2d) / (\Gamma(k+1-d) \Gamma(1-k-d)) ,}
  where \eqn{d = H - 1/2},
  leads to over-/underflow for larger lags \eqn{k};
  hence use the asymptotical formula there.
}
\value{
  numeric vector of length \code{n} of covariances
  \eqn{C(0) \ldots C(n-1)}.
}
\references{Jan Beran (1994), p.63, (2.35) and (2.39).}
\author{Jan Beran (principal) and Martin Maechler (speedup, fine tuning)}
\seealso{\code{\link{ckFGN0}} which does the same for fractional
  Gaussian noise.}
\examples{
str(C.8 <- ckARMA0(50, H = 0.8))
yl <- c(0,max(C.8))
plot(0:49, C.8, type = "h", ylim = yl)
plot(0:49, C.8, type = "h", log = "xy",
     main = "Log-Log  ACF for ARIMA(0,d,0)")
}
\keyword{ts}
