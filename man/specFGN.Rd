\name{specFGN}
\alias{specFGN}
\alias{B.specFGN}
\title{Spectral Density of Fractional Gaussian Noise}
\description{
  Calculation of the spectral density \eqn{f} of
  normalized fractional Gaussian noise with self-similarity parameter
  \eqn{H} at the Fourier frequencies 2*pi*j/m (j=1,...,(m-1)).

  \code{B.specFGN} computes (approximations of) the \eqn{B(\lambda, H)}
  component of the spectrum \eqn{f_H(\lambda)}.
}
\usage{
  specFGN(eta, m, \dots)
B.specFGN(lambd, H, k.approx=3, adjust = (k.approx == 3), nsum = 200)
}
\arguments{
  \item{eta}{parameter vector \code{eta = c(H, *)}.}
  \item{m}{sample size determining Fourier frequencies.}
  \item{\dots}{optional arguments for \code{B.specFGN()}:
    \code{k.approx}, etc}
  \item{lambd}{numeric vector of frequencies in [0, pi]}
  \item{H}{Hurst parameter in \eqn{(\frac 1 2, 1)}{(1/2, 1)}, (can be outside, here).}
  \item{k.approx}{either integer (the order of the Paxson approximation), or
    \code{NULL}, \code{NA} for choosing to use the slow direct sum (of
    \code{nsum} terms.)}
  \item{adjust}{logical indicating (only for \code{k.approx == 3}, the
    default) that Paxson's empirical adjustment should also be used.}
  \item{nsum}{if the slow sum is used (e.g. for k.approx = NA), the
    number of terms.}
}
\details{
  Note that
  \enumerate{
    \item cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
    \item f=theta1*spec and integral[log(spec)]=0.
  }

  Since \pkg{longmemo} version 1.1-0, a fast approximation is available
  (and default), using \code{k.approx} terms and an adjustment
  (\code{adjust=TRUE} in the default case of \code{k.approx=3}),
  which is due to the analysis and S code from Paxson (1997).
  %% MM: Found in fArma (Rmetrics) source, and then on Murad Taqqu's web page,
  %% -> "Statistical methods for long-range dependence"
  %%  -> "Whittle's Approximate MLE"
  %%  == http://math.bu.edu/people/murad/methods/whittle/
  % and its S source code \file{whittle.S} notes that the speedup is due
  % to Vern Paxson. MM found and confirmed the published reference Paxson (1997).

  %% see more comments in   ../R/WhittleEst.R
}
\value{
  \code{specFGN()} returns an object of class \code{"spec"} (see also
  \code{\link{spectrum}}) with components
  \item{freq}{the Fourier frequencies \eqn{\omega_j \in (0,\pi)}{om_j in
      (0,pi)}) at which the spectrum is computed.  Note that
    \eqn{omega_j = 2\pi j/m}{om_j = 2*pi*j/m} for \eqn{j=1,..., m-1},
    and \eqn{m = \left\lfloor\frac{n-1}{2}\right\rfloor}{m = floor((n-1)/2)}.
  }
  \item{spec}{the \emph{scaled} values spectral density \eqn{f(\lambda)}
    values at the \code{freq} values of \eqn{\lambda}.\cr
    \eqn{f^*(\lambda) = f(\lambda) / \theta_1}{f*(lambda) = f(lambda) / theta1}
    adjusted such \eqn{\int \log(f^*(\lambda)) d\lambda = 0}.}
  \item{theta1}{the scale factor \eqn{\theta_1}.}
  \item{H}{the self-similarity parameter from input.}
  \item{method}{a character indicating the kind of model used.}

  \code{B.specFGN()} returns a vector of (approximate) values
  \eqn{B(\lambda, H)}.
}

\references{
  Jan Beran (1994).
  \emph{Statistics for Long-Memory Processes};
  Chapman & Hall, NY.

  Vern Paxson (1997).
  Fast, Approximate Synthesis of Fractional Gaussian Noise for
  Generating Self-Similar Network Traffic;
  \emph{Computer Communications Review} \bold{27} 5, 5--18.
%% MM: ~/save/papers/fracdiff-longrange/Paxson-1997_FGN-CCR.pdf
}
\author{Jan Beran originally (using the slow sum);
  Martin Maechler, based on Vern Paxson (1997)'s code.
}
\seealso{The spectral estimate for fractional ARIMA,
  \code{\link{specARIMA}}; more generally, \code{\link{spectrum}}.
}
\examples{
 str(rg.7  <- specFGN(0.7, m = 100))
 str(rg.5  <- specFGN(0.5, m = 100))# { H = 0.5 <--> white noise ! }

 plot(rg.7) ## work around plot.spec() `bug' in R < 1.6.0
 plot(rg.5, add = TRUE, col = "blue")
 text(2, mean(rg.5$spec), "H = 0.5 [white noise]", col = "blue", adj = c(0,-1/4))
%% TODO: show how to do a  log-log spectrum [nicely!]
## This was the original method in longmemo, upto version 1.0-0 (incl):
 rg.7.o <- specFGN(0.7, m = 100, k.approx=NA, nsum = 200)
 ## quite accurate (but slightly slower):
 rg.7f  <- specFGN(0.7, m = 100, k.approx=NA, nsum = 10000)
 ## comparing old and new default :
 all.equal(rg.7, rg.7.o)# different in about 5th digit
 all.equal(rg.7, rg.7f )# ==> new default is *more* accurate: 1.42 e-6
\donttest{ ## takes about  7 sec {in 2011}:
 rg.7ff <- specFGN(0.7, m = 100, k.approx=NA, nsum = 500000)
 all.equal(rg.7f, rg.7ff)# ~ 10 ^ -7
 all.equal(rg.7  $spec, rg.7ff$spec)# ~ 1.33e-6 -- Paxson is accurate!
 all.equal(rg.7.o$spec, rg.7ff$spec)# ~ 2.40e-5 -- old default is less so
}
}
\keyword{ts}
\keyword{models}
