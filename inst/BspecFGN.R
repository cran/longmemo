##' @title The function f(x) used in the approximation of  B.specFGN()
##' @param x
##' @param lambda
##' @param H
##' @return
##' @author Martin Maechler
fB <- function(x, lambda, H) {
    u <- 2*pi*x
    h <- -(2*H+1)
    (u + lambda)^h + (u - lambda)^h
}

## clearly has a pole "close to zero  --- actually we know the pole is where
##  2pi x = lambda  <==>  x = lambda / (2pi)
curve(fB(x, lambda=.001, H = 0.75), 0, 10)
curve(fB(x, lambda= 0.5, H = 0.75), 0, 10)
curve(fB(x, lambda= 1.5, H = 0.75), 0, 10)
curve(fB(x, lambda= 3.0, H = 0.75), 0, 10)

draw.f <- function(lambda, H, ..., log = "",
                   ylim= if(any("y" == strsplit(log,"")[[1]]))NULL else c(0, 2*fB1))
{
    fB1 <- fB(1, lambda=lambda, H = H)
    curve(fB(x, lambda=lambda, H = H), ylim=ylim, log=log, ...)
    mtext(bquote(list(lambda == .(lambda), H == .(H))))
    abline(v = lambda / (2*pi), lty=2, lwd=3, col = "blue3")
    xrng <- par("usr")[1:2]; if(par("xlog")) xrng <- 10^xrng
    abline(h = 0, lty=3, col="gray20")
    j <- floor(xrng[1]):ceiling(xrng[2])
    lines(j, fB(j, lambda=lambda, H = H), type = "h",
          lty=5, lwd= 0.75, col = "gray40")
    axis(1, at=j[j >= 1][1:2], col="green3")
}

draw.f(lambda=.001, H = 0.75, 0, 10)
draw.f(lambda= 0.5, H = 0.75, 0, 10)
draw.f(lambda= 1.5, H = 0.75, 0, 10)
draw.f(lambda= 3.0, H = 0.75, 0, 10)

draw.f(lambda=.001, H = 0.99, .001, 100, log="xy")
draw.f(lambda=.5,  H = 0.99, .001, 100, log="xy") ## !
draw.f(lambda= 3.0, H = 0.51, .001, 100, log="xy") ## !

##-->  Paxson's  formula using  \int_0^\infty f(x) dx ... is really not feasible
##  as f(x) is *not* defined (or \infty) for x \le \lambda/(2|pi)
