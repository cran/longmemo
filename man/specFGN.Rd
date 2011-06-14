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
  an object of class \code{"spec"} (see also \code{\link{spectrum}})
  with components
  \item{freq}{the Fourier frequencies \eqn{\omega_j \in (0,\pi)}{om_j in
      (0,pi)}) at which the spectrum is computed.  Note that
    \eqn{omega_j = 2\pi j/m}{om_j = 2*pi*j/m} for \eqn{j=1,..., m-1},
    and \eqn{m = \left\lfloor\frac{n-1}{2}\right\rfloor}{m = floor((n-1)/2)}.
  }
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
  Chapman & Hall, NY.
}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\seealso{The spectral estimate for fractional ARIMA,
  \code{\link{specARIMA}}; more generally, \code{\link{spectrum}}.
}
\examples{
 str(rg.7  <- specFGN(0.7, m = 100))
 str(rg.7f <- specFGN(0.7, m = 100, nsum = 10000))
 all.equal(rg.7, rg.7f)# different in about 5th digit only
 str(rg.5  <- specFGN(0.5, m = 100))# { H = 0.5 <--> white noise ! }

 plot(rg.7) ## work around plot.spec() `bug' in R < 1.6.0
 plot(rg.5, add = TRUE, col = "blue")
%% TODO: show how to do a  log-log spectrum [nicely!]
}
\keyword{ts}
\keyword{models}
