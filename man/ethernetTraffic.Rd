\name{ethernetTraffic}
\alias{ethernetTraffic}
\docType{data}
\title{Ethernet Traffic Data Set}
\description{
  Ethernet traffic data from a LAN at Bellcore, Morristown (Leland et al.
  1993, Leland and Wilson 1991).  The data are listed in
  chronological sequence by row.
}
\usage{data(ethernetTraffic)}
\format{
 A times series of length 4000.
}
\source{Jan Beran and Brandon Whitcher by E-mail in fall 1995.
}
\references{%% FIXME ---
}
\examples{
data(ethernetTraffic)
str(ethernetTraffic)
plot(ethernetTraffic)## definitely special
}
\keyword{datasets}
