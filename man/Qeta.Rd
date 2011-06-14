\name{Qeta}
\alias{Qeta}
\title{Approximate Log Likelihood for Fractional Gaussian Noise / Fractional ARIMA}
\description{
  \code{Qeta()} (\eqn{= \tilde{Q}(\eta)}{Q~(eta)} of Beran(1994), p.117)
  is up to scaling the negative log likelihood function of the specified
  model, i.e., fractional Gaussian noise or fractional ARIMA.
}
\usage{
Qeta(eta, model = c("fGn","fARIMA"), n, yper, pq.ARIMA, give.B.only = FALSE)
}
\arguments{
  \item{eta}{parameter vector = (H, phi[1:p], psi[1:q]).}
  \item{model}{character specifying the kind model class.}
  \item{n}{data length}
  \item{yper}{numeric vector of length \code{(n-1)\%/\% 2}, the
    periodogram of the (scaled) data, see \code{\link{per}}.}
  \item{pq.ARIMA}{integer, = c(p,q) specifying models orders of AR and
    MA parts --- only used when \code{model = "fARIMA"}.}
  \item{give.B.only}{logical, indicating if only the \code{B} component
    (of the \code{Values} list below) should be returned.  Is set to
    \code{TRUE} for the Whittle estimator minimization.}
}
\details{
  Calculation of \eqn{A, B} and \eqn{T_n = A/B^2} where
  \eqn{A = 2\pi/n \sum_j 2*[I(\lambda_j)/f(\lambda_j)]},
  \eqn{B = 2\pi/n \sum_j 2*[I(\lambda_j)/f(\lambda_j)]^2}
  and the sum is taken over all Fourier frequencies
  \eqn{\lambda_j = 2\pi*j/n}, (\eqn{j=1,\dots,(n-1)/2}).
  \eqn{f} is the spectral density of fractional Gaussian
  noise or fractional ARIMA(p,d,q) with self-similarity parameter \eqn{H}.

  \deqn{cov(X(t),X(t+k)) = \int \exp(iuk) f(u) du}
}
\note{yper[1] must be the periodogram \eqn{I(\lambda_1)} at
  the frequency \eqn{2\pi/n}, i.e., not the frequency zero !
}
\value{
  a list with components
  \item{n}{= input}
  \item{H}{(\emph{input}) Hurst parameter, = \code{eta[1]}.}
  \item{eta}{= input}
  \item{A,B}{defined as above.}
  \item{Tn}{the goodness of fit test statistic
    \eqn{Tn= A/B^2} defined in Beran (1992)}
  \item{z}{the standardized test statistic}
  \item{pval}{the corresponding p-value P(W > z)}
  \item{theta1}{the scale parameter
    \deqn{\hat{\theta_1} = \frac{{\hat\sigma_\epsilon}^2}{2\pi}}{theta1^ = sigma_e^2 / (2 pi)}
    such that \eqn{f()= \theta_1 f_1()} and \eqn{integral(\log[f_1(.)]) = 0}.
  }
  \item{spec}{scaled spectral density \eqn{f_1} at the Fourier frequencies
    \eqn{\omega_j}, see \code{\link{FEXPest}}; a numeric vector.}
}
\seealso{
  \code{\link{WhittleEst}} computes an approximate MLE for fractional
  Gaussian noise / fractional ARIMA, by minimizing \code{Qeta}.
}
\references{
  Jan Beran (1992).
  A Goodness-of-Fit Test for Time Series with Long Range Dependence.
  \emph{JRSS B} \bold{54}, 749--760.

  Beran, Jan (1994).
  \emph{Statistics for Long-Memory Processes};
  Chapman & Hall.
  (Section 6.1, p.116--119; 12.1.3, p.223 ff)
}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\examples{
data(NileMin)
y <- NileMin
n <- length(y)
yper <- per(scale(y))[2:(1+ (n-1) \%/\% 2)]
eta <- c(H = 0.3)
q.res <- Qeta(eta, n=n, yper=yper)
str(q.res)
}
\keyword{ts}
\keyword{models}
