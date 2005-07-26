\name{NileMin}
\alias{NileMin}
\docType{data}
\title{Nile River Minima, yearly 622--1284}
\description{
  Yearly minimal water levels of the Nile river for the years 622
  to 1281, measured at the Roda gauge near Cairo, (Tousson, p. 366--385).
}
\usage{data(NileMin)}
\format{
  Time-Series (\code{\link{ts}}) of length 663.
}
\source{
 The original Nile river data supplied by Beran only contained only
 500 observations (622 to 1121).  However, the book claimed to have
 660 observations (622 to 1281).  First added the remaining observations
 from the book by hand, and still came up short with only 653
 observations (622 to 1264).  Finally have 663 observations : years
 622--1284 (as in orig. source)
}
\references{
  Tousson, O. (1925)
  Mémoire sur l'Histoire du Nil;
  \emph{Mémoire de l'Institut d'Egypte}

  Jan Beran (1994). Dataset no.~1, p.20--22.
}
\examples{
data(NileMin)
plot(NileMin, main = "Nile River Minima 622 - 1284")
}
\keyword{datasets}
