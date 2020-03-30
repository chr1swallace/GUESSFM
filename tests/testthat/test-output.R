context("Check runs ok")

library(snpStats)
data(testdata, package="snpStats")
z <- Autosomes
cs <- col.summary(z)
z <- z[,which(cs[,"Call.rate"]==1 & !is.na(cs[,"RAF"]) & abs(cs[,"z.HWE"])<5)]

rownames(z) <- paste0("R",1:nrow(z))
colnames(z) <- paste0("V",1:ncol(z))
y <- rnorm(nrow(z))*sqrt(0.5) + rnorm(nrow(z),mean=as(z[,50],"numeric"))*sqrt(0.5)
names(y) <- rownames(z)

tmp <- tempdir()
run.bvs(X=z,Y=y,
          gdir=tmp,family="gaussian",tag.r2=NA)
files <- list.files(tmp)

## specific files
fX <- file.path(tmp,grep("^X_",files,value=TRUE))[1]
fY <- file.path(tmp,grep("^Y_",files,value=TRUE))[1]
fdecode <- file.path(tmp,grep("^decode_[0-9]",files,value=TRUE))[1]
fsdecode <- file.path(tmp,grep("^decode_samples_[0-9]",files,value=TRUE))

hX <- scan(fX,n=2)
hY <- scan(fY,n=2)
dX <- read.table(fX,skip=2)
dY <- read.table(fY,skip=2)
ddecode <- read.table(fdecode)

context("output dimensions")

test_that("output has correct dimensions", {
  expect_equal(hX,dim(z))
  expect_equal(hY,c(length(y),1))
  expect_equal(dim(dX),hX)
  expect_equal(dim(dY),hY)
})

test_that("reader works", {
    expect_s4_class(read.snpmod(tmp),"snpmod")
})

