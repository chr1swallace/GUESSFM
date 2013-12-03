
snpmod.add <- function(e1,e2,priors=snpprior(0:20,n=932,expected=3),dataset.numbers=c(1,2)) {
## e1 <- dx[[1]]
## e2 <- dx[[2]]
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
                   
