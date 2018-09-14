==========================
Dependencies and Structure
==========================

R source in thematical parts --
according to the sections in Chapter 12 of the book.

File |		Section | Title | Pages
----	|	------- | ----- |
`R/simFracGauss.R`  12.1.1   Simulation of fractional Gaussian noise	218-220
`R/simFracARIMA.R`  12.1.2   Simulation of fractional ARIMA($0,d,0$)	220-223
`R/WhittleEst.R`	12.1.3   Whittle estimator for fractional Gaussian
	noise and fractional ARIMA(p,d,q)		223-233
`R/polyFEXP.R`  	12.1.4   Approximate MLE for polynomial FEXP-models	233-237

------------------------

simFracGauss.R:401:gkFGN0 <- function(n, H) {
simFracGauss.R:421:simFGN0 <- function(n,H) {

simFracARIMA.R:308:gkARMA0 <- function(n, H) {
simFracARIMA.R:328:simARMA0 <- function(n,H) {

### `R/WhittleEst.R` :

CetaFGN(eta, ...)
 \--> specFGN(eta,m, ...)

CetaARIMA(eta, p,q, ...)
 \--> specARIMA(eta,p,q,m, ...)

Qeta(eta, model = c("fGn","ARIMA"), .....)
 |--> specFGN(eta,m, ...)
 \--> specARIMA(eta,p,q,m, ...)

Qmin(etatry)
 \--> Qeta(eta, model = c("fGn","ARIMA"), .....)

per(z)

Unfinished main function {was "main program"}:

   mainARIMA(x, model = c("fGn","ARIMA"),


### `R/polyFEXP.R` :

polyFEXP.R:35:lxplot <- function(yper,spec) {
polyFEXP.R:35:llplot <- function(yper,spec) {
polyFEXP.R:...   _new main function __

