\name{videoVBR}
\alias{videoVBR}
\docType{data}
\title{Video VBR data}
\description{
  Amount of coded information (\bold{v}ariable \bold{b}it \bold{r}ate)
  per frame for a certain video sequence. There were about 25 frames per
  second.
}
\usage{data(videoVBR)}
\format{
  a time-series of length 1000.
}
\references{
  Heeke, H. (1991)
  Statistical multiplexing gain for variable bit rate codecs in ATM
  networks. \emph{Int. J. Digit. Analog. Commun. Syst.} \bold{4},
  261--268.

  Heyman, D., Tabatabai, A., and Lakshman, T.V. (1991)
  Statistical analysis and simulation of video teleconferencing in ATM
  networks.
  \emph{IEEE Trans. Circuits. Syst. Video Technol.} \bold{2}, 49--59.

  Jan Beran (1994). Dataset no.~2, p.22--23.
}
\examples{
data(videoVBR)
plot(log(videoVBR), main="VBR Data (log)")
}
\keyword{datasets}
