\name{simGauss}
\alias{simGauss}
\alias{simARMA0}
\alias{simFGN0}
\alias{simFGN.fft}
\title{Simulate (Fractional) Gaussian Processes}
\description{
  Simulation of a Gaussian series \eqn{X(1), \dots, X(n)}.  Whereas
  \code{simGauss} works from autocovariances, where \code{simFGN0} and
  \code{simARMA0} call it,
  for simulating a fractional ARIMA(0,d,0) process (\eqn{d = H-1/2}),
  or fractional Gaussian noise, respectively.
}
\usage{
simARMA0  (n, H)
simFGN0   (n, H)
simFGN.fft(n, H, \dots)
simGauss(autocov)
}
\arguments{
  \item{n}{length of time series}
  \item{H}{self-similarity parameter}
  \item{\dots}{optional arguments passed to \code{\link{B.specFGN}()}.}
  \item{autocov}{numeric vector of auto covariances
    \eqn{\gamma(0), \ldots, \gamma(n-1)}.}
}
\details{
  \code{simGauss} implements the method by Davies and Harte which is
  relatively fast using the FFT (\code{\link{fft}}) twice.

  To simulate ARIMA(p, d, q), (for d in (-1/2, 1,2), you can use
% use   \link[stats]{...} for R >= 1.9.0 :
  \code{\link{arima.sim}(n,  model = list(ar= .., ma = ..),
    innov= simARMA0(n,H=d+1/2) , n.start = 0)}.

  \code{simFGN.fft()} is about twice as fast as \code{simFGN0()} and
  uses Paxson's proposal, by default via
  \code{\link{B.specFGN}(*, k.approx=3, adjust=TRUE)}.
}
\value{
  The simulated series \eqn{X(1), \dots, X(n)}, an \R object of class
  \code{"ts"}, constructed from \code{\link[stats]{ts}()}.
}
\references{Beran (1994), 11.3.3, p.216 f, referring to

  Davis, R.B. and Harte, D.S. (1987).
  Tests for Hurst effect, \emph{Biometrika} \bold{74}, 95--102.

  Vern Paxson (1997).
  Fast, Approximate Synthesis of Fractional Gaussian Noise for
  Generating Self-Similar Network Traffic;
  \emph{Computer Communications Review} \bold{27} 5, 5--18.
%% MM: ~/save/papers/fracdiff-longrange/Paxson-1997_FGN-CCR.ps.gz
}
\author{Jan Beran (original) and Martin Maechler (\code{simGauss},
  speedup, simplication).
  \code{simFGN.fft}: Vern Paxson.
}
\seealso{\code{\link{ckARMA0}} on which \code{simARMA0} relies, and
  \code{\link{ckFGN0}} on which \code{simFGN0} relies.
}
\examples{
  x1 <- simFGN0(100, 0.7)
  x2 <- simARMA0(100, 0.7)
  plot(simFGN0(1000, 0.8)) #- time series plot
}
\keyword{ts}
