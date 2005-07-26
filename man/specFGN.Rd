\name{specFGN}
\alias{specFGN}
\title{Spectral Density of Fractional Gaussian Noise}
\description{
  Calculation of the spectral density \eqn{f} of
  normalized fractional Gaussian noise with self-similarity parameter
  \eqn{H} at the Fourier frequencies 2*pi*j/m (j=1,...,(m-1)).
}
\usage{
specFGN(eta, m, nsum = 200)
}
\arguments{
  \item{eta}{parameter vector \code{eta = c(H, *)}.}
  \item{m}{sample size determining Fourier frequencies.}
  \item{nsum}{length of approximating Riemans sum.}
}
\details{
  Note that
  \enumerate{
    \item cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
    \item f=theta1*spec and integral[log(spec)]=0.
    }
}
\value{
  an object of class \code{"spec"} (see also \code{\link[ts]{spectrum}})
  with components
  \item{freq}{the Fourier frequencies (in \eqn{(0,\pi)}) at which the
    spectrum is computed.}
  \item{spec}{the \emph{scaled} values spectral density \eqn{f(\lambda)}
    values at the \code{freq} values of \eqn{\lambda}.\cr
    \eqn{f^*(\lambda) = f(\lambda) / \theta_1}{f*(lambda) = f(lambda) / theta1}
    adjusted such \eqn{\int \log(f^*(\lambda)) d\lambda = 0}.}
  \item{theta1}{the scale factor \eqn{\theta_1}.}
  \item{H}{the self-similarity parameter from input.}
  \item{method}{a character indicating the kind of model used.}
}
\references{
  Jan Beran (1994).
  \emph{Statistics for Long-Memory Processes};
  Chapman \& Hall, NY.
}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\examples{
 str(r.7  <- specFGN(0.7, m = 100))
 str(r.7f <- specFGN(0.7, m = 100, nsum = 10000))
 all.equal(r.7, r.7f)# different in about 5th digit only
 str(r.5  <- specFGN(0.5, m = 100))

 try(plot(r.7)) ## work around plot.spec() `bug' in R < 1.6.0
 plot(r.5, add = TRUE, col = "blue")
}
\keyword{ts}
\keyword{models}
