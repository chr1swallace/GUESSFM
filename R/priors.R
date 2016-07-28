##' (Beta) Binomial prior for number of SNPs in a model
##'
##' A binomial prior for the number of SNPs in a model may be
##' considered too peaked if there is relatively little prior
##' information about the number of causal SNPs, and, particularly if
##' the posterior model choice is sensitive to this prior, it can be
##' useful to consider a prior with greater spread.  One such choice
##' is the beta binomial model, implemented here, under which the
##' number of SNPs follows a binomial distribution with parameters n,
##' p while p follows a beta distribution with parameters chosen so
##' that the mean and the overdispersion (relative to a binomial
##' distribution) of the number of SNPs is as specified by
##' \code{expected} and \code{overdispersion}, respectively.  An
##' overdispersion of 1 corresponds to a binomial prior.
##' 
##' @title snpprior
##' @param x number of SNPs in a model (defaults to 1:length(groups), ie returns a vector)
##' @param n total number of SNPs or SNP groups available
##' @param expected expected number of SNPs in a model
##' @param overdispersion overdispersion parameter.  Setting this to 1
##' gives a binomial prior.  Values < 1 are nonsensical: if you really
##' believe the prior should be underdispersed relative to a binomial
##' distribution, consider using a hypergeometric prior.
##' @param pi0 prior probability that no SNP is associated
##' @param truncate optional, if supplied priors will be adjusted so models with x>truncate have prior 0
##' @param overdispersion.warning by default, prior distributions should be binomial or beta-binomial (overdispersed).  If you give an overdispersion <1, snpprior will stop with an error.  Set overdispersion.warning=FALSE to override this.
##' @param ## 
##' @param value 
##' @return prior probabilities as a numeric vector
##' @export
##' @author Chris Wallace
##' @examples
##' n<-100 # 100 SNPs in region
##' x <- 1:10 # consider prior for up to 10 causal SNPs
##' xbar <- 3 # expect around 3 causal
##'
##' ## a binomial prior
##' y <- snpprior(x, n, xbar)
##' plot(x, y, type="h")
##' 
##' ## is equivalent to
##' y1.0 <- snpprior(x, n, xbar, overdispersion=1.0)
##' points(x, y1.0, col="red")
##'
##' ##larger values of overdispersion change the distribution:
##' y1.1 <- snpprior(x, n, xbar, overdispersion=1.1)
##' y1.5 <- snpprior(x, n, xbar, overdispersion=1.5)
##' y2.0 <- snpprior(x, n, xbar, overdispersion=2.0)
##' points(x, y1.1, col="orange")
##' points(x, y1.5, col="pink")
##' points(x, y2.0, col="green")
snpprior <- function(x=0:10, n, expected, overdispersion=1, pi0=NA, truncate=NA, overdispersion.warning=TRUE                     
                     ## , value=c("prob","odds")
                           ) {
## ##' @param groups groups of SNPs, from which at most one SNP should be selected
  if(overdispersion < 1 & overdispersion.warning)
    stop("overdispersion parameter should be >= 1")
  x <- as.integer(x)
  if(any(x<0))
    stop("x should be an integer vector >= 0")
  if(!is.na(truncate))
    x <- seq(min(x), max(min(truncate,n),x))
  if(any(x>n))
    stop("max x should be <= n")
##  value <- match.arg(value)
  p <- expected/n
    rho <- (overdispersion - 1)/(n-1)
    
    prob <- if(rho==0) {
                log(dbinom(x, size=n, prob=p))
            } else {
                log(dbetabinom(x, size=n, prob=p, rho=rho))
            }
#  if(value=="prob") {
    ## otherwise value=="odds"
  if(!is.na(pi0) && 0 %in% x)
    prob[x==0] <- log(pi0)
  prob <- prob - lchoose(n,x)
  if(!is.na(truncate))
    prob <- prob - logsum(prob)
  names(prob) <- as.character(x)
  return(exp(prob))
}

change.prior <- function(d, pr) {
  M <- d@models
#  M$ABF <- 
}
