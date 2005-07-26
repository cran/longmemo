\name{plot.FEXP}
\alias{plot.FEXP}
\title{Plot Method for FEXP Model Fits}
\description{
  This an (S3) method for the generic function \code{\link{plot}}
  applied to \bold{f}ractional EXP (FEXP) models.  It plots the data
  periodogram and the \sQuote{FEXP'} model estimated spectrum.
}
\usage{
\method{plot}{FEXP}(x, log = "xy", type = "l",
      col.spec = 4, lwd.spec = 2, xlab = NULL, ylab = expression(hat(f)(nu)),
      main = paste(deparse(x$call)[1]), sub = NULL, \dots)
}
\arguments{
  \item{x}{an \R object of \code{class} \code{"FEXP"}, as from
    \code{\link{FEXPest}()}.}
  \item{log}{character specifying log scale should be used, see
    \code{\link{plot.default}}.  Note that the default log-log scale is
    particularly sensible for long-range dependence.}
  \item{type}{plot type for the periodogram, see \code{\link{plot.default}}.}
  \item{col.spec, lwd.spec}{graphical parameters used for drawing the
    estimated spectrum, see \code{\link{lines}}.}
  \item{xlab, ylab, main, sub}{labels for annotating the plot, see
    \code{\link{title}}, each with a sensible default.}
  \item{\dots}{further arguments passed to \code{\link{plot.default}}.}
}
\author{Martin Maechler}
\seealso{\code{\link{FEXPest}}, \code{\link{plot.default}}
  \code{\link{spectrum}}.
}
\examples{
data(videoVBR)
fE <- FEXPest(videoVBR, order = 3, pvalmax = .5)
plot(fE)
}
\keyword{hplot}
