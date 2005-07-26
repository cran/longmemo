\name{per}
\alias{per}
\title{Simple Periodogram Estimate}
\description{
  Simply estimate the periodogram via the Fast Fourier Transform.
}
\usage{
per(z)
}
\arguments{
  \item{z}{numeric vector with the series to compute the periodogram
    from.}
}
\details{
  This is basically the same as\cr
  \code{spec.pgram(z, fast = FALSE, detrend = FALSE, taper = 0) $ spec},
  and not really recommended to use --- exactly for the reason that
  \code{\link[stats]{spec.pgram}} has the defaults differently,
  \code{fast = TRUE, detrend = TRUE, taper = 0.1}, see that help page.
}
\value{
  numeric vector of length \eqn{1 + floor(n/2)} where \code{n = length(z)}.
}
\author{Jan Beran (principal) and Martin Maechler (fine tuning)}
\seealso{a more versatile periodogram estimate by
  \code{\link[stats]{spec.pgram}}.
}
\examples{
 data(NileMin)
 plot(10*log10(per(NileMin)), type='l')
}
\keyword{ts}
