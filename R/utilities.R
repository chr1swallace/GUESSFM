
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}
makestr <- function(x) {paste(sort(unique(x)),collapse="%")}
##' internal function: xscale
##'
##' Linear transformation to map x to a new range, so min(x)==xrange[1], and max(x)==xrange[2]
##' @title xscale
##' @param x numeric vector to be transformed
##' @param torange new range 
##' @param xrange optional. full range of x vector if you are scaling only part of it. calculated from x if not supplied.
##' @return transformed numeric vector
##' @author Chris Wallace
##' @examples
##' x<-1:10
##' cbind(x,snpmod:::xscale(x,c(0.1,1)),snpmod:::xscale(x,c(1319,20578)))
##' ## Now use only part of x
##' x<-x[2:7]
##' ## WRONG
##' cbind(x,snpmod:::xscale(x,c(0.1,1)),snpmod:::xscale(x,c(1319,20578)))
##' ## RIGHT
##' cbind(x,snpmod:::xscale(x,c(0.1,1),xrange=c(1,10)),snpmod:::xscale(x,c(1319,20578),xrange=c(1,10)))
xscale <- function(x,torange,xrange=c(min(x),max(x))) {
  x0 <- (x - xrange[1]) / (xrange[2] - xrange[1]) # from 0 -> 1
  x0 * (torange[2] - torange[1]) + torange[1]
  
##   LD$A <- (n-2) * (LD$A - mn) / (mx-mn) + 1.5

}

cummean <- function(x) {
  cumsum(x) / 1:length(x)
}

