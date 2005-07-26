\name{CetaFGN}
\alias{CetaFGN}
\title{Covariance Matrix of Eta for Fractional Gaussian Noise}
\description{
  Covariance matrix of \eqn{\hat{\eta}}{eta^} for fractional Gaussian
  noise (fGn).
}
\usage{
CetaFGN(eta, m = 10000, delta = 1e-9)
}
\arguments{
  \item{eta}{parameter vector \code{eta = c(H, *)}.}
  \item{m}{integer specifying the length of the Riemann sum, with step
    size \code{2 * pi/m}.  The default (10000) is realistic.}
  \item{delta}{step size for numerical derivative computation.}
}
\details{
  Currently, the step size for numerical derivative is the same in all
  coordinate directions of \code{eta}.  In principle, this can be far
  from optimal.
}
\value{
  Variance-covariance matrix of the estimated parameter vector
  \eqn{\hat{\eta}}{eta^}.
}
%\references{ .. }
\author{Jan Beran (principal) and Martin Maechler (speedup, fine tuning)}
\seealso{\code{\link{specFGN}}}%% FIXME
\examples{
 (C.7  <- CetaFGN(0.7, m = 256))
 (C.5  <- CetaFGN(eta = c(H = 0.5), m = 256))
 (C.5. <- CetaFGN(eta = c(H = 0.5), m = 1024))
}
\keyword{ts}


