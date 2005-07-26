\name{FEXPest}
\alias{FEXPest}
\alias{print.FEXP}
\title{Fractional EXP (FEXP) Model Estimator}
\description{
  Computes Beran's Fractional EXP or \sQuote{FEXP} model estimator.
}
\usage{
FEXPest(x, order.poly, pvalmax, verbose = FALSE)
\method{print}{FEXP}(x, digits = getOption("digits"), \dots)
}
\arguments{
  \item{x}{numeric vector representing a time series.}
  \item{order.poly}{integer specifying the maximal polynomial order that
    should be taken into account.  \code{order.poly = 0} is equivalent
    to a FARIMA(0,d,0) model.}
  \item{pvalmax}{maximal P-value -- the other iteration stopping
    criterion.  Setting this to \code{1}, will use \code{order.poly} alone.}
  \item{verbose}{logical indicating if iteration output should be printed.}
  \item{digits,\dots}{optional arguments for \code{print} method, see
    \code{\link{print.default}}.}
}
\value{
  An object of class \code{FEXP} which is basically a list with components
  \item{call}{the function \code{\link{call}}.}
  \item{n}{time series length \code{length(x)}.}
  \item{H}{the \dQuote{Hurst} parameter which is simply \code{(1-theta[2])/2}.}
  \item{coefficients}{numeric 4-column matrix as returned from
    \code{\link{summary.glm}()}, with estimate of the full parameter
    vector \eqn{\theta}, its standard error estimates, t- and P-values,
    as from the \code{\link{glm}(*, family = Gamma)} fit.}
  \item{order.poly}{the effective polynomial order used.}
  \item{max.order.poly}{the original \code{order.poly} (argument).}
  \item{early.stop}{logical indicating if \code{order.poly} is less than
    \code{max.order.poly}, i.e., the highest order polynomial terms were
    dropped because of a non-significant P-value.}
  \item{freq, spec}{the Fourier frequencies \eqn{\omega_j} and spectral
    estimate \eqn{f(\omega_j)}.}
  \item{yper}{raw periodogram at Fourier frequencies \eqn{I(\omega_j)}.}

  There currently are a \code{\link{print}} and a \code{plot} method (see
  \code{\link{plot.FEXP}}) for objects of class \code{"FEXP"}.
}
\references{
  Beran, Jan (1993)
  Fitting long-memory models by generalized linear regression.
  \emph{Biometrika} \bold{80}, 817--822.

  Beran, Jan (1994). \emph{Statistics for Long-Memory Processes};
  Chapman \& Hall.
}
\author{Martin Maechler, using Beran's \dQuote{main program} in
  Beran(1994), p.234 ff
}
\seealso{\code{\link{plot.FEXP}}, ... }% < Whittle .. once we have it!
\examples{
data(videoVBR)
(fE <- FEXPest(videoVBR, order = 3, pvalmax = .5))
}
\keyword{ts}

