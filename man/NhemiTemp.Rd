\name{NhemiTemp}
\alias{NhemiTemp}
\docType{data}
\title{Northern Hemisphere Temperature}
\description{
  Monthly temperature for the northern hemisphere for the years
  1854-1989, from the data base held at the Climate Research Unit of the
  University of East Anglia, Norwich, England.  The numbers consist of
  the temperature (degrees C) difference from the monthly average over
  the period 1950-1979.
}
\usage{data(NhemiTemp)}
\format{
  Time-Series (\code{\link{ts}}) of length 1632, frequency 12,
  starting 1854, ending 1990.
}
\source{Jan Beran and Brandon Whitcher by E-mail in fall 1995.
}
\references{
  Jones, P.D. and Briffa, K.R. (1992)
  Global surface air temperature variations during the twentieth
  century, part 1.
  \emph{The Holocene} \bold{2}, 165--179.

  Jan Beran (1994). Dataset no. 5, p.29--31.
}
\examples{
data(NhemiTemp)
plot(NhemiTemp)
mean(window(NhemiTemp, 1950,1979))# (about) 0 ``by definition''
}
\keyword{datasets}
