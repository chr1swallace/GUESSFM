##' Describe and extract the union of models contained in two snpmod objects
##'
##' @title Model unions
##' @param e1 snpmod object 1 or a character vectors of model strings extracted from a snpmod
##' @param e2 snpmod object 2 or a character vectors of model strings extracted from a snpmod
##' @param dataset.labels labels used to refer to e1 and e2 in a detailed union of models. Ignored unless \code{detail==TRUE}.
##' @param detail default FALSE. If TRUE, returns a more detailed data.frame showing union of models and from which dataset they came.
##' @return If \code{detail=FALSE}, invisibly returns a character vector of the unique model strings found in e1 and e2.  If \code{detail=TRUE}, invisibly returns a data.frame containing that vector and some information on which dataset each was found in.
##' @export
##' @author Chris Wallace
model.union <- function(e1,e2,dataset.labels=c(1,2),detail=FALSE) {
## e1 <- dx[[1]]
## e2 <- dx[[2]]
## find shared and distinct models
  if(is(e1,"snpmod")) {
    m1 <- e1@models$str
  } else {
    m1 <- e1
  }
  if(is(e2,"snpmod")) {
    m2 <- e2@models$str
  } else {
    m2 <- e2
  }
  Mint <- intersect(m1,m2)
U1 <- setdiff(m1,m2)
U2 <- setdiff(m2,m1)
N <- length(Mint) + length(U1) + length(U2)
message(N," models found: ",length(Mint)," shared, ",length(U1),"/",length(U2)," unique to datasets 1/2.")
if(detail) {

## ## minimum BF for datasets
## BF1 <- 0 #min(e1@models$jeffreys)
## BF2 <- 0 #min(e2@models$jeffreys)

## matches
mint1 <- match(Mint,m1)
mint2 <- match(Mint,m2)
munq1 <- match(U1,m1)
munq2 <- match(U2,m2)

## labels
  m1 <- data.frame(model=m1,dataset=dataset.labels[[1]])
  m2 <- data.frame(model=m2,dataset=dataset.labels[[2]])
  
## if(!("dataset" %in% colnames(e1@models))) { e1@models$dataset <- dataset.numbers[[1]] }
## if(!("dataset" %in% colnames(e2@models))) { e2@models$dataset <- dataset.numbers[[2]] }
Lint <- paste(m1[ mint1, "dataset"],
              m2[ mint2, "dataset"],
              sep="")
L1 <- m1[ munq1, "dataset"] 
L2 <- m2[ munq2, "dataset"] 

## combine
newmodels <- data.frame(str=c(Mint,U1,U2),
##                         jeffreys=c(e1@models[ mint1, "jeffreys" ] + e2@models[ mint2, "jeffreys" ],
##                           e1@models[ munq1, "jeffreys" ] + BF2,
##                           BF1 + e2@models[ munq2, "jeffreys" ]),
##                         visits=c(e1@models[ mint1, "visits" ] +
##                           e2@models[ mint2, "visits" ],
##                           e1@models[ munq1, "visits" ],
##                           e2@models[ munq2, "visits" ]),
                        dataset=c(Lint,L1,L2),
                        stringsAsFactors=FALSE)
newmodels$size <- nchar(gsub("[^%]","",newmodels$str))+1
newmodels$size[ newmodels$str=="" ] <- 0

invisible(newmodels)
} else {
    invisible(unique(c(m1,m2)))
}
}                 
