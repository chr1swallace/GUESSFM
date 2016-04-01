##' Add Bayes factors across datasets for models evaluated in each
##'
##' @title snpmod.add
##' @param ... two or more objects of class snpmod containing the same models.  Bayes factors will be summed for each model across the datasets.
##' @param priors named numeric vector of prior probabilities for models of given sizes.  Elemenets are the probabilities.  Names are the size.
##' @param BFcolumn column name that contains the (log) Bayes factor.
##' @param .marginals if TRUE, evaluate marginal SNP probabilities.  This can be time consuming, so can be disabled for repeated use of this function to add over multiple objects.  Not intended to be set by users
##' @return an object of class snpmod
##' @author Chris Wallace
snpmod.add <- function(...,priors=snpprior(0:20,n=932,expected=3),BFcolumn="logABF",.marginals=TRUE) {
## e1 <- dx[[1]]
## e2 <- dx[[2]]
  L <- list(...)
  print(length(L))
  if(length(L)==1)
    stop("require at least two items to add, you only gave one")
  if(length(L)>2) {
    tmp <- snpmod.add(L[[1]],L[[2]],priors=priors,BFcolumn=BFcolumn,.marginals=FALSE)
    for(i in 3:length(L)) {
      tmp <- snpmod.add(tmp,L[[i]],priors=priors,BFcolumn=BFcolumn,.marginals=FALSE)
    }
    return(marg.snps(tmp))
  }

  e1 <- L[[1]]
  e2 <- L[[2]]
    
## find shared and distinct models
  m1 <- e1@models$str
  m2 <- e2@models$str
  Mint <- intersect(m1,m2)
  U1 <- setdiff(m1,m2)
  U2 <- setdiff(m2,m1)
  if(length(U1) || length(U2))
    stop("both snpmod objects must contain the same set of models")
  message(length(Mint)," models found.")

## matches
mint1 <- match(Mint,e1@models$str)
mint2 <- match(Mint,e2@models$str)

## size
  size=nchar(gsub("[^%]","",Mint))+1
  size[ Mint=="" ] <- 0
  prior <- priors[as.character(size)]

## combine
newmodels <- data.frame(str=Mint,
                        jeffreys=e1@models[ mint1, BFcolumn ] + e2@models[ mint2, BFcolumn ],
                        size=size,
                        prior=prior,
                        stringsAsFactors=FALSE)
  newmodels$logPP <- log(newmodels$prior) + newmodels$jeffreys
  newmodels$logPP <- newmodels$logPP - logsum(newmodels$logPP)
  newmodels$PP <- exp(newmodels$logPP)
  newmodels <- newmodels[ order(newmodels$PP,decreasing=TRUE), ]
  newmodels$rank <- 1:nrow(newmodels)

  newmod <- e1
  newmod@models <- newmodels
  newmod@model.snps <- strsplit(newmodels$str,"%")
  if(!.marginals)
    return(newmod)
  message("calculating marginal probabilities of inclusion, this may take a while")
  return(marg.snps(newmod))
}
                   
                   
