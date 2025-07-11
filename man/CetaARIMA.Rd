\name{CetaARIMA}
\alias{CetaARIMA}
\title{Covariance for fractional ARIMA}
\description{
  Compute the covariance matrix of \eqn{\hat{eta}}{eta^} for a fractional
  ARIMA process.
}
\usage{
CetaARIMA(eta, p, q, m = 10000, delta = 1e-9)
}
\arguments{
  \item{eta}{parameter vector \code{eta = c(H, phi, psi)}.}
  \item{p,q}{integer scalars giving the AR and MA order respectively.}
  \item{m}{integer specifying the length of the Riemann sum, with step
    size \code{2 * pi/m}.}
  \item{delta}{step size for numerical derivative computation.}
}
\details{
  builds on calling \code{\link{specARIMA}(eta,p,q,m)}
}
\value{
  the (square) matrix containg covariances up to ...
}
\references{Beran(1984), listing on p.224--225.}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\examples{
 (C.7  <- CetaARIMA(0.7, m = 256, p = 0, q = 0))
  C.5  <- CetaARIMA(eta = c(H = 0.5, phi=c(-.06, 0.42, -0.36), psi=0.776),
                    m = 256, p = 3, q = 1)
 ## add row and col names to the cov.matrix (as it is "normal" in R):
 nmC <- c("H", paste0("phi", 1:3), paste0("psi", 1)); dimnames(C.5) <- list(nmC, nmC)
 C.5
}
\keyword{ts}
\keyword{models}
