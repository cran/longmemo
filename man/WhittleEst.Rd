\name{WhittleEst}
\alias{WhittleEst}
\alias{print.WhittleEst}
\alias{nobs}% for R < 2.13.0
\alias{nobs.WhittleEst}
\title{Whittle Estimator for Fractional Gaussian Noise / Fractional ARIMA}
\description{
  Computes Whittle's approximate MLE for fractional Gaussian noise or
  fractional ARIMA (=: fARIMA) models, according to Beran's prescript.
}
\usage{
WhittleEst(x,
           periodogr.x = per(if(scale) x/sd(x) else x)[2:((n+1) \%/\% 2)],
           n = length(x), scale = FALSE,
           model = c("fGn", "fARIMA"),
           p, q,
           start = list(H= 0.5, AR= numeric(), MA=numeric()),
           verbose = getOption("verbose"))

\method{print}{WhittleEst}(x, digits = getOption("digits"), \dots)
}
\arguments{
  \item{x}{numeric vector representing a time series.  Maybe omitted if
    \code{periodogr.x} and \code{n} are specified instead.}
  \item{periodogr.x}{the (raw) periodogram of \code{x}; the default, as by
    Beran, uses \code{\link{per}}, but tapering etc may be an
    alternative, see also \code{\link{spec.pgram}}.}
  \item{n}{length of the time series, \code{length(x)}.}
  \item{scale}{logical indicating if \code{x} should be standardized to
    (\code{\link{sd}}) scale 1; originally, \code{scale = TRUE} used to
    be built-in; for compatibility with other methods, notably plotting
    spectra, \code{scale = FALSE} seems a more natural default.}
  \item{model}{numeric vector representing a time series.}
  \item{p,q}{optional integers specifying the AR and MA orders of the
    fARIMA model, i.e., only applicable when \code{model} is \code{"fARIMA"}.}
  \item{start}{list of starting values; currently necessary for \code{model
      = "fARIMA"} and with a reasonable default for \code{model = "fGn"}.}
  \item{verbose}{logical indicating if iteration output should be printed.}
  \item{digits,\dots}{optional arguments for \code{print} method, see
    \code{\link{print.default}}.}
}
\value{
  An object of class \code{WhittleEst} which is basically a list with components
  \item{call}{the function \code{\link{call}}.}
  \item{model}{= input}
  \item{n}{time series length \code{length(x)}.}
  \item{p,q}{for "fARIMA": order of AR and MA parts, respectively.}
  \item{coefficients}{numeric 4-column matrix of coefficients
    with estimate of the full parameter vector \eqn{\eta}, its standard
    error estimates, z- and P-values.  This includes the Hurst parameter
    \eqn{H}.
  }
  \item{theta1}{the scale parameter \eqn{\hat{\theta_1}}{theta1^}, see
      \code{\link{Qeta}}.}
  \item{vcov}{the variance-covariance matrix for \eqn{\eta}.}
  \item{periodogr.x}{= input (with default).}
  \item{spec}{the spectral estimate \eqn{\hat{f}(\omega_j)}{f^(om_j)}.}

  There is a \code{\link{print}} method, and \code{\link{coef}},
  \code{\link{confint}} or \code{\link{vcov}} methods work as well for
  objects of class \code{"WhittleEst"}.
}
\references{
  Beran, Jan (1994).
  \emph{Statistics for Long-Memory Processes};
  Chapman & Hall.
  (Section 6.1, p.116--119; 12.1.3, p.223 ff)
}
\author{Martin Maechler, based on Beran's \dQuote{main program} in
  Beran(1994).
}
\seealso{
  \code{\link{Qeta}} is the function minimized by these Whittle
  estimators.

  \code{\link{FEXPest}} for an alternative model with Hurst parameter,
  also estimated by a \dQuote{Whittle} approximate MLE, i.e., a
  Whittle's estimator in the more general sense.

  The plot method, \code{\link{plot.WhittleEst}}.
}
\examples{
data(NileMin)
(f.Gn.N  <- WhittleEst(NileMin))                             # H = 0.837
(f.A00.N <- WhittleEst(NileMin, model = "fARIMA", p=0, q=0)) # H = 0.899
confint(f.Gn.N)
confint(f.A00.N)

data(videoVBR)
(f.GN    <- WhittleEst(videoVBR))

## similar {but faster !}:
(f.am00  <- WhittleEst(videoVBR, model = "fARIMA", p=0, q=0))
rbind(f.am00$coef,
      f.GN  $coef)# really similar

f.am11  <- WhittleEst(videoVBR, model = "fARIMA",
                      start= list(H= .5, AR = .5, MA= .5))
f.am11
vcov(f.am11)

op <- if(require("sfsmisc"))
  mult.fig(3, main = "Whittle Estimators for videoVBR data")$old.par  else
  par(mar = c(3,1), mgp = c(1.5, 0.6, 0), mar = c(4,4,2,1)+.1)
plot(f.GN)
plot(f.am00)
plot(f.am11)
par(op)

et <- as.list(coef(f.am11))
et$AR <- c(et$AR, 0, 0) # two more AR coefficients ..
f.am31  <- WhittleEst(videoVBR, model = "fARIMA", start = et)
lines(f.am31, col = "red3") ## drawing on top of ARMA(1,1) above - *small* diff

f.am31 # not all three are "significant"
round(cov2cor(vcov(f.am31)), 3) # and they are highly correlated

et <- as.list(coef(f.am31))
et$AR <- unname(unlist(et[c("AR1", "AR2")]))
f.am21  <- WhittleEst(videoVBR, model = "fARIMA",
                      start = c(et[c("H","AR", "MA")]))
f.am21
lines(f.am21, col = adjustcolor("gold", .4), lwd=4)

##--> ?plot.WhittleEst  for an example using  'periodogr.x'
%% ^^ "FIXME": do it here..
}
\keyword{ts}
