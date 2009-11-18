\name{specARIMA}
\alias{specARIMA}
\title{Spectral Density of Fractional ARMA Process}
\description{
  Calculate the spectral density of a fractional ARMA process
  with standard normal innovations and self-similarity parameter H.
}
\usage{
specARIMA(eta, p, q, m)
}
\arguments{
  \item{eta}{parameter vector \code{eta = c(H, phi, psi)}.}
  \item{p, q}{integers giving AR and MA order respectively.}
  \item{m}{sample size determining Fourier frequencies.}
}
\details{
  at the Fourier frequencies \eqn{2*\pi*j/n}, (\eqn{j=1,\dots,(n-1)}),
%% FIXME: which of these two is ok
  cov(X(t),X(t+k)) = (sigma/(2*pi))*integral(exp(iuk)g(u)du).

  --- or rather -- FIXME --

  1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]

  2. f() = theta1 * f*() ; spec = f*(),  and integral[log(f*())] = 0

}
\value{
  an object of class \code{"spec"} (see also \code{\link{spectrum}})
  with components
  \item{freq}{the Fourier frequencies (in \eqn{(0,\pi)}) at which the
    spectrum is computed.}
  \item{spec}{the \emph{scaled} values spectral density \eqn{f(\lambda)}
    values at the \code{freq} values of \eqn{\lambda}.\cr
    \eqn{f^*(\lambda) = f(\lambda) / \theta_1}{f*(lambda) = f(lambda) / theta1}
    adjusted such \eqn{\int \log(f^*(\lambda)) d\lambda = 0}.}
  \item{theta1}{the scale factor \eqn{\theta_1}.}
  \item{pq}{a vector of length two, \code{= c(p,q)}.}
  \item{eta}{a named vector \code{c(H=H, phi=phi, psi=psi)} from input.}
  \item{method}{a character indicating the kind of model used.}
}
\references{Beran (1994) and more, see .... }%% FIXME
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\seealso{The spectral estimate for fractional Gaussian noise,
  \code{\link{specFGN}}.
  In general, \code{\link{spectrum}} and \code{\link{spec.ar}}.
}
\examples{
 str(r.7  <- specARIMA(0.7, m = 256, p = 0, q = 0))
 str(r.5  <- specARIMA(eta = c(H = 0.5, phi=c(-.06, 0.42, -0.36), psi=0.776),
                       m = 256, p = 3, q = 1))
 plot(r.7)
 plot(r.5)
}
\keyword{ts}
\keyword{models}
