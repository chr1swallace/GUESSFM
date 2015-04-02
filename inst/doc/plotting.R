## ------------------------------------------------------------------------
mydir <-  system.file("extdata",package="GUESSFM")
(load(file.path(mydir,"simdata.RData")))

## ------------------------------------------------------------------------
pp.nsnp(d,plot=TRUE)
pp.nsnp(d,plot=TRUE,expected=3)

## ------------------------------------------------------------------------
groups <- as(sp,"groups")
groups

library(snpStats)
data(for.exercise, package="snpStats") # for SNP information
summx <- guess.summ(sm,groups=groups,snps=snp.support,position="position")
summx <- scalepos(summx,position="position")

library(ggplot2)
signals <- signal.plot(summx)
signals
snps <- ggsnp(summx)
snps
lds <- ggld(X, summx)
lds
chr <- ggchr(summx)
chr

