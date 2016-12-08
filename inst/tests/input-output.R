## load datasets from R2GUESS
library(R2GUESS)
data(data.X, package="R2GUESS")
data(data.Y.Hopx, package="R2GUESS")

## make into SnpMatrix format
library(snpStats)
X <- new("SnpMatrix",matrix(as.raw(as.matrix(data.X)+1),nrow=nrow(data.X),dimnames=dimnames(data.X)))
Y <- data.Y.Hopx[,3]

## dummy run using run.bvs
mydir <- tempfile()
run.bvs(X,Y,gdir=mydir,tag.r2=NA)
list.files(mydir)
list.files("/home/chris/RP/GUESSFM/inst/extdata")
system(paste0("cp ",mydir,"/* /home/chris/RP/GUESSFM/inst/extdata"))
system(paste0("rm /home/chris/RP/GUESSFM/inst/extdata/*sweeps*"))

## read output with GUESSFM
gfm <- read.snpmod(mydir)
best.models(gfm)
"V616" %in% rownames(best.snps(gfm))

