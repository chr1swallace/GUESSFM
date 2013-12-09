##' Add Bayes factors across datasets for models evaluated in each
##'
##' @title snpmod.add
##' @param ... two or more objects of class snpmod containing the same models.  Bayes factors will be summed for each model across the datasets.
##' @param priors named numeric vector of prior probabilities for models of given sizes.  Elemenets are the probabilities.  Names are the size.
##' @param BFcolumn column name that contains the (log) Bayes factor.
##' @param marginals if TRUE, evaluate marginal SNP probabilities.  This can be time consuming, so can be disabled for repeated use of this function to add over multiple objects.  Not intended to be set by users
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
                   
snpmod.add.old <- function(...,priors=snpprior(0:20,n=932,expected=3),dataset.numbers=c(1,2)) {
 e1 <- dx[[1]]
 e2 <- dx[[2]]
    
## find shared and distinct models
 m1 <- e1@models$str
m2 <- e2@models$str
Mint <- intersect(m1,m2)
U1 <- setdiff(m1,m2)
U2 <- setdiff(m2,m1)
N <- length(Mint) + length(U1) + length(U2)
message(N," models found: ",length(Mint)," shared, ",length(U1),"/",length(U2)," unique to datasets 1/2.")

## minimum BF for datasets
BF1 <- 0 #min(e1@models$jeffreys)
BF2 <- 0 #min(e2@models$jeffreys)

## matches
mint1 <- match(Mint,e1@models$str)
mint2 <- match(Mint,e2@models$str)
munq1 <- match(U1,e1@models$str)
munq2 <- match(U2,e2@models$str)

## labels
if(!("dataset" %in% colnames(e1@models))) { e1@models$dataset <- dataset.numbers[[1]] }
if(!("dataset" %in% colnames(e2@models))) { e2@models$dataset <- dataset.numbers[[2]] }
Lint <- paste(e1@models[ mint1, "dataset"],
              e2@models[ mint2, "dataset"],
              sep="")
L1 <- e1@models[ munq1, "dataset"] 
L2 <- e2@models[ munq2, "dataset"] 

## combine
newmodels <- data.frame(str=c(Mint,U1,U2),
                        jeffreys=c(e1@models[ mint1, "jeffreys" ] + e2@models[ mint2, "jeffreys" ],
                          e1@models[ munq1, "jeffreys" ] + BF2,
                          BF1 + e2@models[ munq2, "jeffreys" ]),
                        visits=c(e1@models[ mint1, "visits" ] +
                          e2@models[ mint2, "visits" ],
                          e1@models[ munq1, "visits" ],
                          e2@models[ munq2, "visits" ]),
                        dataset=c(Lint,L1,L2),
                        stringsAsFactors=FALSE)
newmodels$size <- nchar(gsub("[^%]","",newmodels$str))+1
newmodels$size[ newmodels$str=="" ] <- 0
newmodels$prior <- priors[as.character(newmodels$size)]
newmodels$logPP <- log(newmodels$prior) + newmodels$jeffreys
newmodels$logPP <- newmodels$logPP - logsum(newmodels$logPP)
newmodels$PP <- exp(newmodels$logPP)
newmodels <- newmodels[ order(newmodels$PP,decreasing=TRUE), ]
newmodels$rank <- 1:nrow(newmodels)

newmod <- e1
newmod@models <- newmodels
newmod@model.snps <- strsplit(newmodels$str,"%")
message("calculating marginal probabilities of inclusion, this may take a while")
return(marg.snps(newmod))
}
                   
