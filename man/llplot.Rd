\name{llplot}
\alias{llplot}
\alias{lxplot}
\title{Log-Log and Log-X Plot of Spectrum}
\description{
  Log-Log and \dQuote{Log-X} plot of spectrum.
  Very simple utilities, kept here mainly for back compatibility, as
  they appear in the book scripts.
}
\usage{
llplot(yper, spec)
lxplot(yper, spec)
}
\arguments{
  \item{yper}{periodogram values}
  \item{spec}{spectrum values}
}
%\value{
% }
\author{Jan Beran (principal) and Martin Maechler (speedup, fine tuning)}
\seealso{\code{\link[stats]{spectrum}()} from standard \R (package
  \pkg{stats}).
}
%\examples{
%}
\keyword{hplot}
\keyword{ts}
