\name{FEXPest}
\title{Fractional EXP (FEXP) Model Estimator}
\alias{FEXPest}
\alias{print.FEXP}
\alias{.ffreq}
\description{
  \code{FEXPest(x, *)} computes Beran's Fractional EXP or \sQuote{FEXP}
  model estimator.

  \code{.ffreq(n)} returns the Fourier frequencies \eqn{\frac{2\pi j}{n}}{2pi j/n}
  (of a time series of length \code{n}).
}
\usage{
FEXPest(x, order.poly, pvalmax, verbose = FALSE)
\method{print}{FEXP}(x, digits = getOption("digits"), \dots)

.ffreq(n, full = FALSE)
}
\arguments{
  \item{x}{numeric vector representing a time series.}
  \item{order.poly}{integer specifying the maximal polynomial order that
    should be taken into account.  \code{order.poly = 0} is equivalent
    to a FARIMA(0,d,0) model.}
  \item{pvalmax}{maximal P-value -- the other iteration stopping
    criterion and \dQuote{model selection tuning parameter}.% << Martin Maechler
    Setting this to \code{1}, will use \code{order.poly} alone, and
    hence the final model order will be \code{= order.poly}.}
  \item{verbose}{logical indicating if iteration output should be printed.}
  \item{digits,\dots}{optional arguments for \code{print} method, see
    \code{\link{print.default}}.}
  %% freq
  \item{n}{a positive integer, typically the length of a time series.}
  \item{full}{logical indicating if \code{n \%/\% 2} or by default
    \dQuote{only} \code{(n-1) \%/\% 2} values should be returned.}
}
\value{
  \code{FEXPest(x,..)} returns an object of class \code{FEXP}, basically a list with components
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
  \item{spec}{the spectral estimate \eqn{f(\omega_j)}, at the Fourier
    frequencies \eqn{\omega_j}.  Note that \code{.ffreq(x$n)}
    recomputes the Fourier frequencies vector (from a fitted FEXP or
    WhittleEst model \code{x}).}
  \item{yper}{raw periodogram of (centered and scaled \code{x}) at
    Fourier frequencies \eqn{I(\omega_j)}.}

  There currently are methods for \code{\link{print}()},
  \code{\link{plot}} and \code{\link{lines}} (see
  \code{\link{plot.FEXP}}) for objects of class \code{"FEXP"}.
}
\references{
  Beran, Jan (1993)
  Fitting long-memory models by generalized linear regression.
  \emph{Biometrika} \bold{80}, 817--822.

  Beran, Jan (1994).
  \emph{Statistics for Long-Memory Processes};
  Chapman & Hall.
}
\author{Martin Maechler, using Beran's \dQuote{main program} in
  Beran(1994), p.234 ff
}
\seealso{
  \code{\link{WhittleEst}};
  the plot method, \code{\link{plot.FEXP}}.
}
\examples{
data(videoVBR)
(fE  <- FEXPest(videoVBR, order = 3, pvalmax = .5))
(fE3 <- FEXPest(videoVBR, order = 3, pvalmax = 1 ))

(fE7 <- FEXPest(videoVBR, order = 3, pvalmax = 0.10))
##--> this also choses order 2, as "FE" :
all.equal(fE $coef,
          fE7$coef) # -> TRUE

confint(fE)
confint(fE7, level = 0.99)

.ffreq(8)
.ffreq(8, TRUE)
stopifnot(all.equal((1:3)/4,
                    .ffreq(8) / pi))
}
\keyword{ts}

