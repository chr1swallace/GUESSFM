
cond.best <- function(X,Y,best=NULL,p.thr=1e-6, ...) {
  if(is.null(best)) {    
    cond <- snp.rhs.tests(Y ~ 1, snp.data=X, data=data.frame(Y=Y, row.names=rownames(X)), ...)
    p.thr <- 1
  } else {
    cond <- snp.rhs.tests(as.formula(paste("Y ~", paste(best, collapse="+"))),
                          snp.data=X,
                          data=cbind(data.frame(Y=Y, row.names=rownames(X)),as(X[,best],"numeric")), ...)
  }
  p <- p.value(cond)
  pmin <- min(p, na.rm=TRUE)
  if(pmin<p.thr) {
    newbest <- colnames(X)[ which.min(p) ]
    message(newbest," ",format.pval(pmin))
    return(newbest)
  }
  return(NULL)
}

backend.guess <- function(gX, gY, gdir, nsweep, nchains, best, nsave, nexp, nexp.sd, guess.command) {
  ## print decode file
  decode.file <- paste0(gdir,"/decode_",nsweep)
  decode.file2 <- paste0(gdir,"/decode_samples_",nsweep)
  message("creating decode files under ",gdir)
  if(!file.exists(gdir))
    dir.create(gdir)
  decode <- matrix((1:ncol(gX)) - 1,ncol=1,dimnames=list(varname=colnames(gX),"varnum")) # 1-based
  write.table(decode,
              file=decode.file,
              append=FALSE,quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE)
  cat(rownames(gX),file=decode.file2,sep="\n")

  cols.guess <- which(colnames(gX) %in% best) - 1 # 0-based
  message("conditional analysis suggests a starting model with 0-based columns: ",paste(cols.guess,collapse=" "))
  message("aka: ", paste(best, collapse=" "))
  
  ##   library(glmnet)
  ##   m <- cv.glmnet(x=gX,y=gY,family=family,nfolds=3)
  ##   plot(m)
  ##   cf <- coef(m)
  ##   cols.glmnet <- which(cf[,1]!=0) # includes intercept, 1-based
  ##   cols.guess <- cols.glmnet[-1] - 2 # excludes intercept, 0-based
  ##   message("lasso suggests a starting model with 0-based columns: ",paste(cols.guess,collapse=" "))
  
  ## print intput files
  message("creating input files under ",gdir)
  if(!file.exists(gdir))
    dir.create(gdir)
  x.file <- paste0(gdir,"/X_",nsweep)
  y.file <- paste0(gdir,"/Y_",nsweep)
  init.file <- paste0(gdir,"/init_",nsweep)
  cat(length(cols.guess),cols.guess, sep="\n", file=init.file)
  cat(nrow(gX),"\n",ncol(gX),"\n",sep="",file=x.file)
  cat(nrow(gX),"\n",1,"\n",sep="",file=y.file)
  write.table(gX,
              file=x.file,
              append=TRUE,quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
  cat(gY, file=y.file,
      sep="\n",append=TRUE)
  
  ## par file
  par.file=file.path(gdir,"par.xml")
  if(!file.exists(par.file))
    file.copy(system.file("Par_file_example.xml",package="GUESSFM"), par.file)
  
  com <- sprintf("%s -history -X %s -Y %s -nsweep %s -burn_in %s -out %s/out -par %s/par.xml -top %s -init %s -Egam %s -Sgam %s -n_chain %s",
                 guess.command,x.file,y.file,nsweep,round(nsweep/11),gdir,gdir,nsave,init.file,nexp,nexp.sd,nchains)
#  if(grid) {
    message("running GUESS with command")
    message(com)
#    system(paste("/home/chrisw/local/bin/qCom.sh -N GUESS",com))
#  } else {
#    message("running GUESS locally with command")
#    message(com)
    system(sprintf("%s > %s/log",com,gdir), wait=FALSE)
#  }
}
##' Bayesian variable selection
##'
##' Wrapper for GUESS, an alternative to using R2GUESS
##' @param X snpMatrix object, holding genotypes
##' @param Y vector or matrix of phenotypes, nrow(Y)==nrow(X)
##' @param gdir directory where all the GUESS input/output files will be. If it doesn't exist, it will be created.
##' @param sub optional number < nrow(X). If supplied, only the subset of samples defined by 1:sub will be used.
##' @param strat optional covariate vector. If supplied GUESS will be run on residuals from Y ~ strat.
##' @param family family for Y ~ strat regression.  default "gaussian".
##' @param tag.r2 r squared value at which to tag to avoid numerical instability.  Default of 0.99 has worked well in our experience.
##' @param nexp expected number of causal variants in region
##' @param nsave number of models to save, see documentation for GUESS 
##' @param nsweep number of sweeps, see documentation for GUESS
##' @param nchains number of chains, see documentation for GUESS
##' @param guess.command Command to run GUESS, if GUESS is not on your PATH.
##' @param backend "GUESS" (default) to run GUESS, or "dummy" to just generate the files for checking
##' @param boot EXPERIMENTAL. Optional integer > 0 indicates that a number of bootstrap replicates should be run to understand variability in output
##' @param as.is EXPERIMENTAL if boot=TRUE, and as.is=TRUE, explicitly exclude the unchanged data structure from the set of bootstrapped datasets
##' @param dominance EXPERIMENTAL. If TRUE, include a dominance effect at each SNP
##' @export
##' @return nothing.  side effect is to set GUESS running in the background.  This takes a while (typically several hours).
run.bvs <- function(X,Y,gdir="test",sub=NA,
                    strat=NULL,family="gaussian", nsweep=55000,nchains=3,
                    nexp=3,boot=0,tag.r2=0.99, nsave=1000, as.is=TRUE,dominance=FALSE,
                    guess.command="GUESS",
                    backend=c("guess","dummy")) { 

  backend <- match.arg(backend)
  if(!file.exists(gdir))
    dir.create(gdir,recursive=TRUE)
  
  ## tag
  colnames(X) <- make.names(colnames(X))
  if(!is.na(tag.r2)) {
    tags <- tag(X,tag.threshold=tag.r2)
    X <- X[, unique(tags@tags)]
  }

  ## matrix phenotypes
  if(!is.matrix(Y))
    Y <- matrix(Y,ncol=1)
  
  ## format snp data  
  N <- as(X,"numeric")
  p <- 0
  
  use <- complete.cases(N) & complete.cases(Y)
  if(!is.null(strat))
      use <- use & complete.cases(strat)
  use <- which(use)
  if(!is.na(sub))
    use <- use[1:sub]
  gX <- N[use,]
  gY <- Y[use]
  n <- length(use)
  v <- apply(gX,2,var)
  use.cols <- which(v>0)
  gX <- gX[,use.cols]
  m <- length(use.cols)


  if(!is.null(strat)) {
      strat <- strat[use]
      m0 <- glm(gY ~ strat, family=family)
      gY <- residuals(m0)
  }
      
  message("using ",n," samples ",m," SNPs, ",p," confounders.")
  
  ## IF GUESS: what are the best regressors? (95,220?)
  if(backend=="guess") {
    best <- NULL
    for(j in 1:ncol(Y)) {
      while(length(newbest <- cond.best(X, Y[,j], best, family=family)))
        best <- c(best,newbest)
    }
  
    ## prior
    if(dominance) {
        p <- 1.5*nexp/m
    } else {
        p <- nexp/m
    }
    nexp.sd <- sqrt(m*p*(1-p)*1) # + overdispersion=0
    message("setting prior parameters nexp, sd: ",nexp, " ",nexp.sd)

  }

  ## dominance effects
  if(dominance) {
      D <- matrix(ifelse(gX>1.5,1,0),nrow=nrow(gX),dimnames=dimnames(gX))
      v <- apply(D,2,var)
      D <- D[,v>0]
      cr <- cor(D)
      wh <- which(cr>0.99,arr.ind=TRUE)
      wh <- wh[wh[,1]>wh[,2],]
      if(nrow(wh))
          D <- D[,-wh[,2]]
      cr <- cor(gX,D)
      wh <- which(cr>0.99,arr.ind=TRUE)
      if(nrow(wh))
          D <- D[,-wh[,2]]
      colnames(D) <- paste0(colnames(D),".D")
      gX <- cbind(gX,D)
  }
  
  for(i in 0:boot) {    
    if(i==0 & !as.is)
      next
    
    if(i>0) {
      use <- sample(1:nrow(gX),size=nrow(gX),replace=TRUE)
      gdir2 <- tempfile(paste0(gdir,"_boot_"), tmpdir=".")
      dir.create(gdir2,recursive=TRUE)
    } else {
      use <- 1:nrow(gX)
      gdir2 <- gdir
    }
      
    switch(backend,
           guess=backend.guess(gX=gX[use,], gY=gY[use], gdir=gdir2,
               nsweep, nchains, best, nsave, nexp, nexp.sd, guess.command=guess.command),
           dummy=NULL## ,
##            sbams=backend.sbams(X=gX[use,],Y=gY[use],gdir=gdir2)
           )
  }
    
}
 


