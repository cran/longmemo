####-*- mode: R; kept-old-versions: 12;  kept-new-versions: 20; -*-

### ------ =========================================
### 12.1.3 Whittle estimator for fractional Gaussian
### ------ noise and fractional ARIMA(p,d,q)		-- book, p. 223-233
###	   =========================================

## Splus-functions and program for the calculation of Whittle's
## estimator and the goodness of fit statistic defined in Beran (1992).
##
## The models are fractional Gaussian noise or fractional ARIMA.
##
## The data series may be divided
## into subseries for which the parameters are fitted separatly.
##____________________________________________________________________


CetaFGN <- function(eta, m = 10000, delta = 1e-9)
{
    ## Purpose: Covariance of hat{eta} for fGn = fractional Gaussian noise
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

    if(1 > (M <- length(eta)))
        stop("eta[] must have length at least 1") #FIXME
    if(2 > (m <- as.integer(m)))
        stop(sQuote("m")," must be at least 2")
    mhalfm <- (m-1) %/% 2:2

    ## partial derivatives of log f (at each Fourier frequency)

    lf <- matrix(ncol = M,nrow = mhalfm)
    f0 <- specFGN(eta,m)$spec

    for(j in (1:M)) { ## FIXME{MM}: better not using same delta for all [j]
        etaj <- eta
        etaj[j] <- etaj[j]+delta
        fj <- specFGN(etaj,m)$spec
        lf[,j] <- log(fj/f0)/delta
    }

    solve(crossprod(lf)) *m ## << improve: chol2inv(qr(lf)) !

}## CetaFGN()


CetaARIMA <- function(eta, p,q, m = 10000, delta = 1e-9)
{
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  if(1 > (M <- length(eta))) stop("eta[] must have length at least 1")#FIXME
  if(2 > (m <- as.integer(m)))stop(sQuote("m")," must be at least 2")#FIXME
  mhalfm <- (m-1) %/% 2:2

  ## partial derivatives of log f (at each Fourier frequency)
  lf <- matrix(ncol = M, nrow = mhalfm)

  f0 <- specARIMA(eta, p=p, q=q, m=m)$spec

  for(j in 1:M) {
    etaj <- eta
    etaj[j] <- etaj[j] + delta
    fj <- specARIMA(etaj, p=p, q=q, m=m)$spec
    lf[,j] <- log(fj/f0)/delta
  }

  ## Calculate D

  m* chol2inv(qr.R(qr(lf)))# == m* solve(crossprod(lf))

}## CetaARIMA()

specFGN <- function(eta, m, nsum = 200)
{
  ## Purpose: Calculation of the spectral density f of
  ## normalized fractional Gaussian noise with self-similarity parameter H=H
  ## at the Fourier frequencies 2*pi*j/m (j=1,...,(m-1)/2).
  ##
  ## Remarks:
  ## -------
  ## 1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
  ## 2. f=theta1*spec and integral[log(spec)]=0.
  ## -------------------------------------------------------------------------
  ## INPUT: m   = "sample size", 2 * #{Fourier frequencies}
  ##        eta = (H, *);  H = self-similarity parameter
  ##
  ## OUTPUT: list(spec=spec,theta1=theta1)
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  ##---------parameters for the calculation of f--------

  if(2 > (m <- as.integer(m)))
      stop(sQuote("m")," must be at least 2")

  if(length(eta) > 1) warning("eta[2..] are not used")
  H <- eta[1]
  hh <- -2*H-1

  mhalfm <- (m-1) %/% 2:2
  ##-- Fourier frequencies lambda_j = 2*pi*(j-1)/m (j=1,2,...,(m-1)/2)
  lambd <- 2*pi/m * (1:mhalfm)

  ##-----   calculation of f at Fourier frequencies   -------
  j <- 2*pi*(1:nsum)
  ## Using outer() instead of the loop is  *slower* (!)
  ##spec <- colSums(abs(outer(j, lambd, "+"))^hh + abs(outer(j, lambd, "-"))^hh)
  spec <- double(mhalfm)
  for(i in 1:mhalfm)
      spec[i] <- sum(abs(j + lambd[i])^hh + abs(j - lambd[i])^hh)

  spec <- sin(pi*H)/pi * gamma(-hh) * (1-cos(lambd)) * (lambd^hh + spec)

  ##----- adjust spectrum such that int{log(spec)} = 0 -------
  theta1 <- exp(2/m * sum(log(spec)))
  ## and specS = spec*() = spec / theta1 is adjusted :
  r <- list(freq = lambd, spec = spec/theta1, theta1 = theta1,
            H = H, method = paste("fGN( H=",format(H),")"))
  class(r) <- "spec"
  r
}## specFGN()

specARIMA <- function(eta, p,q, m)
{
  ## Purpose: Calculation of the spectral density of fractional ARMA
  ##	with standard normal innovations and self-similarity parameter H=H
  ##    at the Fourier frequencies 2*pi*j/m (j=1,..., (m-1)/2).
  ## cov(X(t),X(t+k)) = (sigma/(2*pi))*integral(exp(iuk)g(u)du).
  ##
  ## Remarks:
  ## -------
  ## 1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
  ## 2. f=theta1*spec and integral[log(spec)]=0.
  ## -------------------------------------------------------------------------
  ## INPUT: m = sample size
  ##        H = theta[1] = self-similarity parameter
  ##        phi = theta[    2:(p+1)]   = AR(p)-parameters
  ##        psi = theta[(p+2):(p+q+1)] = MA(q)-parameters
  ##
  ## OUTPUT: list(spec=spec,theta1=theta1)
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

  ##---------parameters for the calculation of f--------

  if(0 > (p <- as.integer(p)))
	stop("`p'  must be non-negative integer")
  if(0 > (q <- as.integer(q)))
        stop("`q'  must be non-negative integer")
  if(1+p+q != (M <- length(eta)))
        stop("eta[] must have length 1+p+q")
  if(2 > (m <- as.integer(m)))
        stop(sQuote("m")," must be at least 2")
  mhalfm <- (m-1) %/% 2:2

  H <- eta[1]

  ##------   Fourier frequencies   ----------------------
  ##------   x = 2*pi*(j-1)/m (j=1,2,...,(n-1)/2)   -----
  x <- 2*pi/m * (1:mhalfm)

  ##-----   calculation of f at Fourier frequencies   -------


  if(p > 0) {
    phi <- eta[2:(p+1)]
    px <- outer(x, 1:p)
    Rar <- cos(px) %*% phi
    Iar <- sin(px) %*% phi

    far <- (1-Rar)^2 + Iar^2
  } else {
      phi <- numeric(0)
      far <- 1
  }

  if(q > 0) {
    psi <- eta[(p+2):(p+q+1)]
    px <- outer(x, 1:q)
    Rma <- cos(px) %*% psi
    Ima <- sin(px) %*% psi

    fma <- (1+Rma)^2 + Ima^2
  } else {
      psi <- numeric(0)
      fma <- 1
  }

  spec <- fma/far * sqrt(2 - 2*cos(x))^(1-2*H)
  ##   			 2 - 2*cos(x) == (1-cos(x))^2 + sin(x)^2

  r <- list(freq = x, spec = spec, theta1 = 1/(2*pi),
            pq = c(p,q), eta = c(H = H, phi = phi, psi = psi),
            method = paste("ARIMA(",p,", ", H - 1/2,", ", q,")",sep = ""))
  class(r) <- "spec"
  r
}## specARIMA()


per <- function(z) {
  ## Purpose:  Computation of the periodogram via FFT
  ## -------------------------------------------------------------------------
  ## Arguments: z : time-series
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

### MM: per[-1] is the same as
### --    1/(2*pi) * spec.pgram(z, fast=FALSE, taper=0, detrend=FALSE)$spec
###    "[-1]" : spec.pgram() doesn't try to use the *wrong* value lambda=0

### Note the "fast = FALSE" needed to emulate this !
  n <- length(z)
  (Mod(fft(z))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
}



Qeta <- function(eta, model = c("fGn","fARIMA"), n, yper, pq.ARIMA,
                 verbose = getOption("verbose"), give.B.only = FALSE)
{
    ## Purpose: Calculation of A, B and Tn = A/B^2
    ## where A = 2pi/n sum 2*[I(lambda_j)/f(lambda_j)],
    ##       B = 2pi/n sum 2*[I(lambda_j)/f(lambda_j)]^2  and
    ## the sum is taken over all Fourier frequencies
    ## lambda_j = 2pi*j/n (j=1,...,(n-1)/2.
    ## f is the spectral density of fractional Gaussian
    ## noise or fractional ARIMA(p,d,q) with self-similarity parameter H=H.
    ## cov(X(t),X(t+k))=integral(exp(iuk)f(u)du)
    ##
    ## NOTE: yper[1] must be the periodogram I(lambda_1) at
    ## ----  the frequency 2pi/n (i.e. not the frequency zero !).
    ## -------------------------------------------------------------------------
    ## INPUT: H
    ##        (n,nhalfm = trunc[(n-1)/2] and the
    ##         nhalfm-dimensional  GLOBAL vector `yper' must be defined.)
    ##   (n,yper): now arguments {orig. w/ defaults 'n=n, yper=yper'}

    ## OUTPUT: list(n=n,H=H,A=A,B=B,Tn=Tn,z=z,pval=pval,
    ##              theta1=theta1,spec=spec)
    ##
    ##         Tn is the goodness of fit test statistic
    ##         Tn=A/B^2 defined in Beran (1992),
    ##         z is the standardized test statistic,
    ##         pval the corresponding p-value P(w>z).
    ##         theta1 is the scale parameter such that
    ##         f=theta1*spec and integral(log[spec])=0.
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

    if(verbose) cat("in function Qeta")
    model <- match.arg(model)

    stopifnot(is.numeric(yper),
              length(yper) == (n - 1) %/% 2)

    H <- eta[1]
    ##         spectrum at Fourier frequencies

### FIXME:  `n' is undefined    ; `p', `q' are undefined, too for ARIMA

### FIXME 2 : For MLE Minimization, only `B' is used --> do this via argument!

    sp <- switch(model,
                 "fGn"   = specFGN(eta,n),
                 "fARIMA" = specARIMA(eta,p,q,n)
                 )

    theta1 <- sp $theta1
    spec  <- sp $spec

    ## Tn := A / B^2

    yf <- yper/spec

    B <- 2*(2*pi/n)*sum(yf)

    if(give.B.only)
        return(B)
    ## else
    yfyf <- yf*yf

    A <- 2*(2*pi/n)*sum(yfyf)

    Tn <- A/(B^2)
    z <- sqrt(n)*(pi*Tn-1)/sqrt(2)

    pval <- 1-pnorm(z)
    theta1 <- B/(2*pi)

    list(n = n, H = H,
         eta = eta, A = A, B = B, Tn = Tn ,z = z, pval = pval,
         theta1 = theta1, spec = spec)
}## Qeta()



###MMMMMMMMMMMMMMMMMMMMMMM
###
### MAIN PROGRAM  (fARIMA)
###
###MMMMMMMMMMMMMMMMMMMMMMM


## ____ UNFINISHED ____  NOTA BENE: also need to "fix" Qeta() above!!
if(FALSE)
#########

## Idea:  define these *inside* fWhittleEst(.) below !
## ----   OTOH, have (somewhat) documented  ../man/Qeta.Rd
## (and Jan has it as public function too!) ==============

### FIXME:  Now copy the changes/technique from FEXPest() in ./polyFEXP.R
###         				        -------      ~~~~~~~~~~~~

## Martin's new main function  *INSTEAD* of any main program :

fWhittleEst <- function(x,
                        model = c("fGn","fARIMA"),
                        p, q,
                        start = list(H= numeric(), AR= numeric(), MA=numeric()),
                        verbose = getOption("verbose"))
{
    ## Author: Martin Maechler, 2004-2005,
    ## NOTA BENE: no "loop" over different data segments !!

    cl <- match.call()
    stopifnot(is.list(start),
              sapply(start, is.numeric),
              length(start$H) == 1)
    ## If 'fGn', don't need p nor q  -- maybe use default p = 0, q = 0
    model <- match.arg(model)
    if(model == "fGn") {
        stopifnot(length(start) == 1)   # namely only "$ H"
        if(!missing(p)) {
            if(p != 0) {
                p <- 0
                warning("'p' != 0 does not make senses in \"fGn\" model")
            }
        } else p <- 0
        if(!missing(q)) {
            if(q != 0) {
                q <- 0
                warning("'q' != 0 does not make senses in \"fGn\" model")
            }

        } else q <- 0
    }
    else { ## fARIMA
        if(0 > (p <- as.integer(p))) stop("must have integer p >= 0")
        if(0 > (q <- as.integer(q))) stop("must have integer q >= 0")
        stopifnot(length(start$AR) == p,
                  length(start$MA) == q)
    }

    eta <- unlist(start) # c(H, ar[1:p], ma[1:q]) # where p=0, q=0 is possible
    H0 <- eta[1]
    M <- length(eta)

    n <- length(x)
    stopifnot(is.numeric(x), n >= 4)
    nhalfm <- (n-1) %/% 2

    ## hmm: should give a warning when doing this:
    H <- max(0.2,min(H0,0.9))           # avoid extreme initial values
    eta[1] <- H

    ## standardize data
    y <- x
    y <- (y-mean(y))/sqrt(var(y))

    ## periodogram of data
    ffr <- 2*pi/n* (1:nhalfm)		# Fourier frequencies
    yper <- per(y)[2:(nhalfm+1)]


    ## find estimate

###################
## definition of function to be minimized
###################

    Qmin <- function(etatry) {
        ## Purpose: function to be minimized for MLE
        ## in R, we nicely have access to all the variables in "main" function
        ## ---------------------------------------------------------------------
        ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

        result <- Qeta(etatry, model = model, n = n, yper = yper,
                       pq.ARIMA = c(p,q), give.B.only = TRUE)
        if(getOption("verbose")) cat("Qmin(etatry =",etatry,") -> B =",result,"\n")
        drop(result)
    }


    s <- 2*(1 - H)

    etatry <- eta
    ##S+result <- nlmin(Qmin, etatry, xc.tol = 0.0000001, init.step = s)
    ##S+         ----- FIXME: use optim() instead
    ##S+eta <- result$x                   # = arg min _{eta} Qmin(eta, ..)
    if(M == 1) { # one-dimensional -- only 'H'
        res <- optimize(Qmin, lower = 0.1, upper = 0.99)
        eta <- res$minimum

    } else { ## M > 1
        res <- optim(Qmin, etatry)
        eta <- res$par
    }

    ## calculate goodness of fit statistic
    Qresult <- Qeta(eta, n=n, yper = yper)
    th1 <- Qresult$theta1
    theta <- c(th1, eta)

    ## output

    SD <- switch(model,
                 "fGn"   = CetaFGN  (eta),
                 "fARIMA" = CetaARIMA(eta,p,q)
                 )
    ## MM: Shouldn't SD be a symmetric matrix such that 'byrow' is not needed?
    SD <- matrix(SD, ncol = M ,nrow = M ,byrow = TRUE)/n

    cat("theta=",theta,fill = TRUE)
    cat("H=",eta[1],fill = TRUE)

    ssd <- sqrt(diag(SD))

    ##   Hlow <- eta[1]-1.96*sqrt(SD[1,1])
    ##   Hup  <- eta[1]+1.96*sqrt(SD[1,1])
    ##   if(verbose) cat("95%-C.I. for H: [",Hlow,",",Hup,"]\n")

    etalow <- eta - 1.96*ssd
    etaup <-  eta + 1.96*ssd

    ## if(verbose) cat("periodogram is in yper\n")

    fest <- th1 * Qresult$spec
    ## if(verbose) cat("spectral density is in fest\n")

    list(call = cl, model = model, n = n,
         theta = theta, eta = eta, SD = SD,
         ci.eta95 = cbind(etalow, etaup),
         fest = fest)

    ## FIXME: add a class !!

}# end fWhittleEst()
