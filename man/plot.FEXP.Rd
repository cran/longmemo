\name{plot.FEXP}
\title{Plot Method for FEXP and WhittleEst Model Fits}
\alias{plot.FEXP}
\alias{plot.WhittleEst}
\alias{lines.FEXP}
\alias{lines.WhittleEst}
\description{
  (S3) methods for the generic functions \code{\link{plot}} (and
  \code{\link{lines}})   applied to \bold{f}ractional EXP (FEXP) and
  \code{"WhittleEst"} (\code{\link{FEXPest}},
  models.
  \code{plot()} plots the data periodogram and the \sQuote{FEXP'} model
  estimated spectrum, where \code{\link{lines}()} and does the latter.
}
\usage{
\method{plot}{FEXP}(x, log = "xy", type = "l",
      col.spec = 4, lwd.spec = 2, xlab = NULL, ylab = expression(hat(f)(nu)),
      main = paste(deparse(x$call)[1]), sub = NULL, \dots)

## (With identical argument list:)
\method{plot}{WhittleEst}(x, log = "xy", type = "l",
      col.spec = 4, lwd.spec = 2, xlab = NULL, ylab = expression(hat(f)(nu)),
      main = paste(deparse(x$call)[1]), sub = NULL, \dots)

\method{lines}{FEXP}      (x, type = "l", col = 4, lwd = 2, \dots)
\method{lines}{WhittleEst}(x, type = "l", col = 4, lwd = 2, \dots)
}
\arguments{
  \item{x}{an \R object of \code{class} \code{"FEXP"}, as from
    \code{\link{FEXPest}()}.}
  \item{log}{character specifying log scale should be used, see
    \code{\link{plot.default}}.  Note that the default log-log scale is
    particularly sensible for long-range dependence.}
  \item{type}{plot type for the periodogram, see \code{\link{plot.default}}.}
  \item{col.spec, lwd.spec, col, lwd}{graphical parameters used for drawing the
    estimated spectrum, see \code{\link{lines}}.}
  \item{xlab, ylab, main, sub}{labels for annotating the plot, see
    \code{\link{title}}, each with a sensible default.}
  \item{\dots}{further arguments passed to \code{\link{plot.default}}.}
}
\author{Martin Maechler}
\seealso{\code{\link{FEXPest}}, \code{\link{WhittleEst}};
  \code{\link{plot.default}} and \code{\link{spectrum}}.
}
\examples{
data(videoVBR)
fE <- FEXPest(videoVBR, order = 3, pvalmax = .5)
plot(fE)
fE3 <- FEXPest(videoVBR, order = 3, pvalmax = 1)#-> order 3
lines(fE3, col = "red3", lty=2)

f.GN    <- WhittleEst(videoVBR)
f.am21  <- WhittleEst(videoVBR, model = "fARIMA",
                      start= list(H= .5, AR = c(.5,0), MA= .5))
lines(f.GN,   col = "blue4")
lines(f.am21, col = "goldenrod")

##--- Using a tapered periodogram ---------
spvVBR <- spec.pgram(videoVBR, fast=FALSE, plot=FALSE)
fam21 <- WhittleEst(periodogr.x = head(spvVBR$spec, -1),
                    n = length(videoVBR), model = "fARIMA",
                    start= list(H= .5, AR = c(.5,0), MA= .5))
fam21
f.am21 # similar but slightly different

plot(fam21)

## Now, comparing to traditional ("log-X", not "log-log") spectral plot:
plot(fam21, log="y")

## compared to the standard R spectral plot : s
if(dev.interactive(TRUE)) getOption("device")()# a new graphics window
plot(spvVBR, log = "yes", col="gray50")
all.equal(.ffreq(fE$n) / (2*pi) -> ffr.,
          head(spvVBR$freq, -1))# TRUE
lines(ffr., fam21$spec, col=4, lwd=2)
## need to adjust for different 'theta1':
lines(ffr., f.am21$spec * fam21$theta1 / f.am21$theta1,
      col = adjustcolor("tomato", 0.6), lwd=2)
}
\keyword{hplot}
