\name{NBSdiff1kg}
\alias{NBSdiff1kg}
\docType{data}
\title{NBS measurement deviations from 1 kg}
\encoding{UTF-8}
\description{
  NBS weight measurements - deviation from 1 kg in micrograms, see the
  references.  The data are listed in chronological sequence by row.
}
\usage{data(NBSdiff1kg)}
\format{
 A time series of length 289.
}
\source{Jan Beran and Brandon Whitcher by E-mail in fall 1995.
}
\references{
  H.P. Graf, F.R. Hampel, and J.Tacier (1984).
  The problem of unsuspected serial correlations.
  In J. Franke, W. Härdle, and R.D. Martin, editors, \emph{Robust
    and Nonlinear Time Series Analysis}, Lecture Notes in Statistics
  \bold{26}, 127--145; Springer.

  Pollak, M., Croakin, C., and Hagwood, C. (1993).
  \emph{Surveillance schemes with applications to mass calibration.}
  NIST report 5158; Gaithersburg, MD.
}
\examples{
data(NBSdiff1kg)
## maybe str(NBSdiff1kg) ; plot(NBSdiff1kg) ...
}
\keyword{datasets}
