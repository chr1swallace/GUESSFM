##' Stepwise regression on a snpMatrix
##'
##' Calls snp.rhs.test to find next best predictor to add, possibly given a list of existing chosen predictors.
##' @title cond.best
##' @param X snpMatrix object
##' @param Y phenotype vector
##' @param best SNPs chosen so far, NULL if none
##' @param stepwise.p.thr maximum p value to continue adding predictors
##' @param stepwise.max.predictors maximum predictors to choose, NA if continue until stepwise.p.thr is met
##' @param ... arguments passed to snp.rhs.test
##' @return new predictor, NULL if none found that meets conditions
##' @export
##' @author Chris Wallace
cond.best <- function(X,Y,best=NULL,stepwise.p.thr=1e-3,stepwise.max.predictors=NA, ...) {
 if(!is.na(stepwise.max.predictors) & length(best)>=stepwise.max.predictors)
   return(NULL)
 ## cs <- col.summary(X)
 ## X <- X[,cs[,"MAF"]>0.05]
 data <- data.frame(Y=Y, row.names=rownames(X))
 if(is.null(best)) {
   Xtest <- X
   cond <- single.snp.tests(phenotype=data$Y, snp.data=Xtest)
   p <- p.value(cond,1)
   stepwise.p.thr <- 1
  } else {
    data <- cbind(data,as(X[,best],"numeric"))
    LD <- ld(X,X[,best],stats="R.squared")
    maxLD <- apply(LD,1,max,na.rm=TRUE)
    drop <- unique(c(names(maxLD)[which(is.na(maxLD) | maxLD>0.5)],best))
    Xtest <- X[,setdiff(colnames(X),drop)]
    if(!ncol(Xtest))
      return(NULL)
    cond <- snp.rhs.tests(as.formula(paste("Y ~", paste(best, collapse="+"))),
                          snp.data=Xtest,
                          data=data,...) # binomial by default
    
    p <- p.value(cond)
  }
  pmin <- min(p, na.rm=TRUE)
  if(pmin<stepwise.p.thr) {
    newbest <- colnames(Xtest)[ which.min(p) ]
    cs <- col.summary(Xtest[,newbest])
    message(newbest,"\tMAF=",signif(cs[1,"MAF"],2),"\tp=",format.pval(pmin))
    return(newbest)
  }
  return(NULL)
}

backend.guess <- function(gX, gY, gdir, nsweep, nchains, best, nsave, nexp, nexp.sd, guess.command,wait=FALSE) {
    
  opt.bak <- options(scipen=1000000000)
  ## print decode file
  decode.file <- file.path(gdir,paste0("decode_",nsweep))
  decode.file2 <- file.path(gdir,paste0("decode_samples_",nsweep))
  message("creating decode files under ",gdir)
  if(!file.exists(gdir))
    dir.create(gdir)
  decode <- matrix((1:ncol(gX)) - 1,ncol=1,dimnames=list(varname=colnames(gX),"varnum")) # 0-based
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
  message("writing to ",init.file,":")
  cat(length(cols.guess),cols.guess, sep="\n")
  cat(length(cols.guess),cols.guess, sep="\n", file=init.file)
  cat(nrow(gX),"\n",ncol(gX),"\n",sep="",file=x.file)
  cat(nrow(gX),"\n",1,"\n",sep="",file=y.file)
  ## uncertain genotypes?
  cs <- col.summary(gX)
  if(any(cs[,"Certain.calls"]<1)) {
    P <- g2post(gX)
    P <- P[,2] + 2*P[,3]
    N <- matrix(P,nrow(gX),ncol(gX))    
    write.table(N,file=x.file, append=TRUE, quote = FALSE, sep = " ", row.names=FALSE,col.names=FALSE)
  } else {  
    write.SnpMatrix(gX, file=x.file, as.alleles= FALSE, append = TRUE, quote = FALSE, sep = " ", eol = "\n", na = "NA",row.names=FALSE,col.names=FALSE)
  }
  cat(gY, file=y.file,
      sep="\n",append=TRUE)
  
  ## par file
  par.file=file.path(gdir,"par.xml")
  if(!file.exists(par.file))
    file.copy(system.file("Par_file_example.xml",package="GUESSFM"), par.file)

  com <- sprintf("%s -history -X %s -Y %s -nsweep %s -burn_in %s -out %s/out -par %s/par.xml -top %s -init %s -Egam %s -Sgam %s -n_chain %s > %s/log",
                 guess.command,x.file,y.file,nsweep,round(nsweep/11),gdir,gdir,nsave,init.file,nexp,nexp.sd,nchains,gdir)
  options(opt.bak)
    message("running GUESS with command")
    message(com)
  if(!is.null(guess.command))
    system(com, wait=wait)
  return(com)
}
##' Bayesian variable selection
##'
##' Wrapper for GUESS, an alternative to using R2GUESS
##' @param X snpMatrix object, holding genotypes
##' @param Y vector or matrix of phenotypes, nrow(Y)==nrow(X)
##' @param gdir directory where all the GUESS input/output files will
##'     be. If it doesn't exist, it will be created.
##' @param sub optional number < nrow(X). If supplied, only the subset
##'     of samples defined by 1:sub will be used.
##' @param covars optional matrix or vector of covariates. If supplied
##'     GUESS will be run on residuals from glm(Y ~ .,
##'     data=as.data.frame(covars)).
##' @param family family for Y ~ covars regression.  default
##'     "gaussian".
##' @param nsweep number of sweeps, see documentation for GUESS
##' @param nchains number of chains, see documentation for GUESS
##' @param nexp expected number of causal variants in region
##' @param tag.r2 r squared value at which to tag to avoid numerical
##'     instability.  Default of 0.99 has worked well in our
##'     experience.
##' @param nsave number of models to save, see documentation for GUESS
##' @param guess.command Command to run GUESS.  This is normally
##'     automatically set to the version of GUESS installed by
##'     R2GUESS, but you may override with a full path to a system
##'     version of GUESS if you prefer.
##' @param wait logical.  default FALSE. if TRUE, run.bvs will wait
##'     for GUESS to finish, rather than running in background
##' @param ... GUESS starts from a stepwise solution found by
##'     cond.best.  Use ... to pass arguments firectly to cond.best to
##'     influence the p value threshold or number of predictors at
##'     which the stepwise search stops.
##' @export
##' @return nothing.  side effect is to set GUESS running in the
##'     background.  This takes a while (typically several hours).
run.bvs <- function(X,Y,gdir="test",sub=NA,
                    covars=NULL,family="gaussian", nsweep=55000,nchains=3,
                    nexp=3,tag.r2=0.99, nsave=1000, 
                    guess.command=NULL,
                    wait=FALSE,
                    ...) { 

    if(is.null(guess.command)) {
        pack.root <- system.file(package = "R2GUESS")
        ESS.directory <- file.path(pack.root, "bin", .Platform$r_arch)
        guess.executable <- ifelse(.Platform$OS.type == "unix", "GUESS", 
                                   "GUESS.exe")
        guess.command <- file.path(ESS.directory, guess.executable)
    }
    
  if(!file.exists(gdir))
    dir.create(gdir,recursive=TRUE)
  
  ## tag
  colnames(X) <- make.names(colnames(X))
  if(!is.na(tag.r2)) {
    tfile <- file.path(gdir,"tags.RData")
    tags <- tag(X,tag.threshold=tag.r2)
    X <- X[, unique(tags@tags)]
    message("saving tags object to ",tfile)
    save(tags,file=tfile)
  }

  ## matrix phenotypes
  if(!is.matrix(Y))
    Y <- matrix(Y,ncol=1)
  
  ## format snp data  
  p <- 0
  rs <- row.summary(X)
  
  use <- rs[,"Call.rate"]==1 & complete.cases(Y)
  if(!is.null(covars)) {
      use <- use & complete.cases(covars)
      p <- ncol(covars)
    }
  if(any(!use) || !is.na(sub)) {
    use <- which(use)
    if(!is.na(sub))
      use <- use[1:sub]
    gX <- X[use,]
    gY <- Y[use]
    if(!is.null(covars))
      covars <- covars[use,]
  } else {
    gX <- X
    gY <- Y
  }
  n <- nrow(gX)
  cs <- col.summary(gX)
  use.cols <- which(!is.na(cs[,"z.HWE"]) & cs[,"MAF"]>0)
  m <- length(use.cols)
  if(m < ncol(gX))
    gX <- gX[,use.cols]

  if(!is.null(covars)) {
    if(!is.data.frame(covars) & !is.matrix(covars)) # deal with vectors, without messing up data.frames
      covars <- as.data.frame(as.matrix(covars))
    m0 <- glm(gY ~ ., data=covars,family=family)
    gY <- residuals(m0)
    family <- "gaussian"
  }
      
  message("using ",n," samples ",m," SNPs, ",p," covariates.")

  ## any uncertain genotypes?
  ## cs <- col.summary(X)
  ## if(length(wh <- which(cs[,"Certain.calls"] < cs[,"Calls"]))) {
  ##   N <- round(as(X,"numeric"))
  ##   X <- new("SnpMatrix",N+1)
  ## }
  
  ## IF GUESS: what are the best regressors? (95,220?)
  best <- NULL
  cs <- col.summary(X)
  wh <- use.cols[which(cs[,"MAF"]>0.05 & cs[,"Certain.calls"]>0.9)]
  while(length(newbest <- cond.best(X[use,wh], gY, best, family=family, ...))) {
    best <- c(best,newbest)
  }
  if(length(best)<3) {
    wh <- setdiff(wh,which(colnames(X) %in% best))
    best <- c(best, sample(colnames(X)[wh],3-length(best)))
  }
    
    ## prior
    ## TODO: allow set overdispersion parameter, relax requirement that overdisp >=1 in snpprior
  p.each <- nexp/m
  nexp.sd <- sqrt(m*p.each*(1-p.each)*1) # + overdispersion=0
  message("setting prior parameters nexp, sd: ",nexp, " ",nexp.sd)

  backend.guess(gX=gX, gY=gY, gdir=gdir,
                nsweep, nchains, best, nsave, nexp, nexp.sd, guess.command=guess.command,wait=wait)
}
 


