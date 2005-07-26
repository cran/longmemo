\name{Qeta}
\alias{Qeta}
\title{Function to be minimized for approx. MLE of frARIMA or frGn}
\description{
  Qmeta() is up to scaling the log likelihood function of the two models
  indicated.
}
\usage{
Qeta(eta, model = c("fGn","fARIMA"), n, yper, pq.ARIMA,
     verbose = getOption("verbose"), give.B.only = FALSE)
}
\arguments{
  \item{eta}{parameter vector = (H, phi[1:p], psi[1:q]).}
  \item{model}{character specifying the kind model class.}
  \item{n}{data length}
  \item{yper}{numeric vector of length \code{(n-1)\%/\% 2}, the
    periodogram of the (scaled) data, see \code{\link{per}}.}
  \item{pq.ARIMA}{integer, = c(p,q) specifying models orders of AR and
    MA parts --- only used when \code{model = "fARIMA"}.}
  \item{verbose}{logical indicating if diagnostic output should be
    produced during fitting.}
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
%% FIXME --- the concept of a  GLOBAL  vector is  THE HORROR !!
% \section{INPUT}{H;
%   (\code{n, nhalfm = trunc[(n-1)/2]} and the
%   nhalfm-dimensional  GLOBAL vector \code{yper} must be defined.)
% }
\note{yper[1] must be the periodogram \eqn{I(\lambda_1)} at
  the frequency \eqn{2\pi/n} (i.e., not the frequency zero !).
}
\value{
  a list with components
  \item{n}{= input}
  \item{H}{Hurst parameter, = \code{eta[1]}.}
  \item{A,B}{defined as above.}
  \item{Tn}{the goodness of fit test statistic
    \eqn{Tn= A/B^2} defined in Beran (1992)}%<< Reference? [give below] FIXME!
  \item{z}{the standardized test statistic}
  \item{pval}{the corresponding p-value P(W > z)}
  \item{theta1}{the scale parameter such that \code{f=theta1*spec} and
    integral(log[spec])=0.}
  \item{spec}{}
}
% \references{Beran (1992), ???
%   __FIXME__
% }
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
