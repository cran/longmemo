library(longmemo)

(H0 <- c(seq(0.5, .999, by = 1/32), 0.98, 0.99))

mform <- function(x, digits=6, wid=9)
    formatC(x, flag="-", width=wid, digits=digits)

options(width = 99)
.proctime00 <- proc.time()

## CetaARIMA(*) does NOT depend on H !!
hiH <- H0 >= 0.9375
for(H in H0[!hiH]) {
    cat("H=", mform(H,wid=7),
        "; CetaFGN:",       mform(CetaFGN  (eta= c(H = H), m = 256)),
        "; C.ARIMA(*,0,0):",mform(CetaARIMA(eta= c(H = H), m = 256, p=0,q=0)),
        "\n")
}
## Ignoring differences for 0.9375 and 0.99 : -- only just for --no-long-double case
## IGNORE_RDIFF_BEGIN
for(H in H0[hiH]) {
    cat("H=", mform(H,wid=7),
        "; CetaFGN:",       mform(CetaFGN  (eta= c(H = H), m = 256)),
        "; C.ARIMA(*,0,0):",mform(CetaARIMA(eta= c(H = H), m = 256, p=0,q=0)),
        "\n")
}
## IGNORE_RDIFF_END

ems <- setNames(,6:18)
cFGN <- cT <- NA_real_ + ems
cARM <- array(NA, dim = c(3,3, length(ems)), dimnames = list(NULL, NULL, names(ems)))
for(em in ems) { ## oops -- becomes slow [specFGN() !] from about em = 13:
    m <- 2^em; ch.m <- as.character(em)
    cT[[ch.m]] <- system.time({
        cARM[ , , ch.m] <- CetaARIMA(eta = c(H = 0.7, phi=0.6, psi= -0.3), m = m, p=1,q=1)
        cFGN[     ch.m] <- CetaFGN  (eta = c(H = 0.7), m = m)
    })[1] # user.time
    cat("m= 2^",formatC(em, width = 2),": ",
        " C_{FGN}= ",      mform(cFGN[ch.m]),
        "; C_{AR; 1,1}= ", mform(cARM[1,1,ch.m]), "\n", sep="")
    if(interactive()) {
        print(cARM[,,ch.m], digits = 6) ## .. for --no-long-double
        cat("\n")
    }
}
## IGNORE_RDIFF_BEGIN
cat("Timings:\n"); print(cT)
(dcFGN <- diff(cFGN))
cA.18 <- matrix(c(6.6007,  -6.55447, 0.133132,
                 -6.55447, 11.2902, -5.43864,
                 0.133132, -5.43864, 6.80152), 3L)
all.equal(cA.18, cARM[,,"18"], tolerance = 0) # ~ 9.376e-7
(dcAdiffs <- apply(cARM[,,] - c(cA.18), 3, norm))
(mnD <- mean(abs(tail(cFGN, 3) - 0.42825))) # 0.0002317033
## IGNORE_RDIFF_END
stopifnot(exprs = {
    -0.075 < dcFGN ; dcFGN < 0
    mnD < 4e-4
    all.equal(cA.18, cARM[,,"18"], tolerance = 4e-6)
    tail(dcAdiffs, 5) < c(0.51, 0.28, 0.14, 0.05, 4e-5)
    diff(dcAdiffs) < 0
})

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
