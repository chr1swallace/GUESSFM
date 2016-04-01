## ------------------------------------------------------------------------
library(snpStats)
data(for.exercise, package="snpStats")
set.seed(12346)
X <- snps.10[,101:120]
n <- nrow(X)
causal <- c("rs1555897","rs7069505")
Y <- rnorm(n,mean=as.numeric(X[,causal[1]]))*sqrt(0.2) +
  rnorm(n,mean=as.numeric(X[,causal[2]]))*sqrt(0.2) +
  rnorm(n)*sqrt(0.6)

## ------------------------------------------------------------------------
library(GUESSFM)
summary(col.summary(X))
X <- impute.missing(X)
summary(col.summary(X))

## ------------------------------------------------------------------------
ld <- show.ld(X=X)

## ------------------------------------------------------------------------
mydir <-  system.file("extdata",package="GUESSFM")

## what files were created?
list.files(mydir)

## read output with GUESSFM
d <- read.snpmod(mydir)

## examine the best models and SNPs with greatest marginal support within the tagged data.
best.models(d)
best.snps(d)

## ------------------------------------------------------------------------
sel=rownames(best.snps(d))
ld[c(sel,causal),c(sel,causal)]

## ------------------------------------------------------------------------
(load(file.path(mydir,"tags.RData")))
tags
tagsof(tags,causal)
taggedby(tags,sel)

## ------------------------------------------------------------------------
dx<-expand.tags(d,tags)

## ------------------------------------------------------------------------
best <- best.models(dx,cpp.thr=0.9)
library(speedglm)
abf <- abf.calc(y=Y,x=X,models=best$str,family="gaussian")
sm <- abf2snpmod(abf,expected=3)

## ------------------------------------------------------------------------
best.snps(d)
best.snps(dx)
best.snps(sm)

## ------------------------------------------------------------------------
sp <- snp.picker(sm,X)
summary(sp)
plot(sp)

