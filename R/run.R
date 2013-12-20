
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


run.bvs <- function(X,Y,gdir="test",
                    strat=NULL,sub=NA,grid=FALSE, family="gaussian", nsweep=55000,nchains=3,
                    nexp=3,boot=0,tag.r2=0.99, nsave=1000, as.is=TRUE,dominance=FALSE,
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
    library(snpStats)
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
               nsweep, nchains, best, nsave, nexp, nexp.sd, grid),
           dummy=NULL,
           sbams=backend.sbams(X=gX[use,],Y=gY[use],gdir=gdir2))
  }
    
}
 


backend.sbams <- function(X,Y,gdir="test",grid=FALSE) {

  library(snpStatsWriter)
  data.file=paste0(gdir,"/data")
  grid.file=paste0(gdir,"/grid")
  out.file=paste0(gdir,"/output")
  write.sbams(X=X, response=Y, file=data.file)
  prior.corr <- 0.9
  omega <- c(0.001,0.005,0.01,0.02,0.05)
  write.table(cbind(omega, omega*(1-prior.corr)/prior.corr),
              file=grid.file,
              col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  com <- paste("/home/chrisw/local/bin/sbams -d",data.file, "-g", grid.file, "-o", out.file)
  message(com)
  com <- sprintf("%s &> %s/log",com,gdir)
  if(grid) {
    message("running sbams on Q")
    system(paste("/home/chrisw/local/bin/qCom.sh -N sbams",com))
  } else {
    message("running sbams locally")
    system(com, wait=FALSE)
  }
}
backend.guess <- function(gX, gY, gdir, nsweep, nchains, best, nsave, nexp, nexp.sd, grid) {
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
  par.file=paste0(gdir,"/par.xml")
  if(!file.exists(par.file))
    file.copy("/home/chrisw/GUESS_v1.1/Example/Input/Par_file_example.xml", par.file)
  ##  a <- xmlNode("arg", attrs = c(default="T"), xmlNode("name", "foo"), xmlNode("defau
  ## a$children[[3]] <- xmlNode("duncan")
  
  com <- sprintf("/home/chrisw/GUESS_v1.1/Main/GUESS -history -X %s -Y %s -nsweep %s -burn_in %s -out %s/out -par %s/par.xml -top %s -init %s -Egam %s -Sgam %s -n_chain %s",
                 x.file,y.file,nsweep,round(nsweep/11),gdir,gdir,nsave,init.file,nexp,nexp.sd,nchains)
  if(grid) {
    message("running GUESS on Q with command")
    message(com)
    system(paste("/home/chrisw/local/bin/qCom.sh -N GUESS",com))
  } else {
    message("running GUESS locally with command")
    message(com)
    system(sprintf("%s > %s/log",com,gdir), wait=FALSE)
  }
}

