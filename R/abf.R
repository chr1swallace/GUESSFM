##' Uses the BIC approximation to calculate approximate Bayes factors for specified models
##'
##' @title Calculate Approximate Bayes Factors
##' @param y response vector
##' @param x explanatory variables, a SnpMatrix object, data.frame or numeric matrix
##' @param models vector of models to consider
##' @param family family argument to pass to glm.  Currently only "binomial" and "gaussian" are implemented for the \code{glm.fit} method.
##' @param backend currently ignored 
##' @param method use glm.fit or glm to fit the models.  Default is glm.fit which should be faster, but is less forgiving about the odd missing value etc, so if your code is giving potentially glm related errors, try running with \code{method="glm"} in the first instance.
##' @param R2 matrix giving pairwise r-squared measures of LD from which tag SNPs will be calculated when not directly available
##' @param snp.data if R2 is missing, it is calculated from this SnpMatrix object
##' @param return.R2 if true, return the calculated R2 matrix.  Useful if you are analysing several strata of a population and you wish to avoid repeating the calculation.
##' @param verbose print lots of progress messages if TRUE.  Default is FALSE.
##' @return a data.frame containing model name, ABF, and an indicator of whether the ABF was calculating directly or via a tag SNP
##' @author Chris Wallace
abf.calc <- function(y,x,models,family="binomial",
                     q=NULL,method=c("glm.fit","glm"),
                     R2=NULL,snp.data=NULL,return.R2=FALSE,verbose=FALSE) { # raftery, wen

  method <- match.arg(method)
  message("Calculating BICs using ",method)
  
  ## initialize return object
  models.orig <- models
  lBF <- structure(rep(as.numeric(NA),length(models)),names=models.orig)

  ## models that can be fitted directly
  snps <- strsplit(models,"%")
  cols.ok <- c("1",colnames(x))
  use <- sapply(snps, function(snp.names) all(snp.names %in% cols.ok))
  message(sum(use)," of ", length(models), " models can be evaluated directly")
  if(!all(use)) {
    snps <- snps[use]
    models <- models[use]
  }  

  ## models that need to be assessed through tagging
  if(!all(use) && (!is.null(snp.data) || !is.null(R2))) {
    message("attempting to find tag models for those with missing SNPs")    
    models.missing <- setdiff(models.orig,models)
    snps.missing <- setdiff(unlist(strsplit(models.missing,"%")),colnames(x))
    allsnps <- c(snps.missing,colnames(x))
    if(is.null(R2))
      R2 <- ld(snp.data[,snps.missing],
               snp.data[,intersect(colnames(x),colnames(snp.data))],
               stat="R.squared")
    snps.tags <- colnames(R2)[ apply(R2,1,which.max) ]
    names(snps.tags) <- snps.missing
    models.alt <- models.missing
    for(snp in snps.missing) {
      models.alt <- sub(paste0("(%|^)",snp,"(%|$)"), paste0("\\1",snps.tags[snp],"\\2"), models.alt)
    }
    names(models.alt) <- models.missing
    models <- unique(c(models.alt,models))
    message("including tags, a total of ",length(models)," need evaluation")
}

  ## do the fitting
  models <- c("1",models)
  snps <- strsplit(models,"%")
  
  ## prepare data
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(!is.null(q))
      qm <- model.matrix(~q)
  
  if(method=="glm.fit") {
    if(is(x,"data.frame"))
      x <- as.matrix(x)
    x2<-cbind(one=1,x[,intersect(unique(unlist(snps)),colnames(x))])
    comp <- complete.cases(x2) & !is.na(y)
    if(!is.null(q))
        comp <- comp & complete.cases(qm)
    if(!all(comp)) {
        message("Dropping ",sum(!comp)," samples due to incompleteness. ",sum(comp)," remain.")
      x2 <- x2[comp,]
      y2 <- y[comp]
        if(!is.null(q)) 
            qm <- qm[use,drop=FALSE]
    } else {
      y2 <- y
    }
    logn <- log(nrow(x))
    family <- switch(family,
                     "gaussian"=gaussian(link="identity"),
                     "binomial"=binomial(link="logit"))
    print(family)
    snps <- lapply(snps,setdiff,"1")
       bics <- mclapply(seq_along(models), function(i) {
      if(verbose && i %% 100 == 0)
        cat(i,"\t")
      k=length(snps[[i]])+1
      if(!is.null(q)) {
          model <- glm.fit(cbind(x2[, snps[[i]] ],qm), y2, family=family)
      } else {
          model <- glm.fit(x2[,c("one",snps[[i]])], y2, family=family)
      }
      class(model) <- c(class(model),"glm")
      BIC(model)
    })    
  }
  
  if(method=="glm") {
    df <- as.data.frame(x[,intersect(unique(unlist(snps)),colnames(x))])
    if(!is.null(q))
        df$q <- q
    df$y <- y                      
    comp <- complete.cases(df)
    if(!all(comp)) {
         message("dropping ",sum(!comp)," SNPs due to incompleteness")
     df2 <- df[comp,]
    } else {
      df2 <- df
    }
    ## is baseline already included?
    l <- sapply(snps,length)
    for(wh in which(l==0))
        snps[[wh]] <- "1"
  
    bics <- mclapply(seq_along(models), function(i) {
        if(verbose && i %% 100 == 0)
            cat(i,"\t")
        if(!is.null(q)) 
            f <- as.formula(paste("y ~ q +",paste(snps[[i]],collapse="+")))
        else
            f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
        BIC(glm(f, data=df2, family=family))
    })
  }
  
  ## calculate ABF
  bics <- unlist(bics)
  ##  use <- use[-1] # drop base model
  models <- models[-1]
  lBF.fits <- (bics[-1] - bics[1]) / (-2)
  names(lBF.fits) <- models

  ## remap to original model ids
  ## direct
  mint <- intersect(names(lBF),models)
  lBF[ mint ] <- lBF.fits[ mint ]
  ## tags
  mint <- intersect(names(lBF),models.missing)
  lBF[ mint ] <- lBF.fits[ models.alt[ mint ] ]

  lBF <- data.frame(model=models.orig,tag=!(models.orig %in% names(lBF)),lBF=lBF[models.orig],stringsAsFactors=FALSE)
  if(return.R2 && exists("R2"))
    return(list(lBF=lBF,R2=R2))
  return(lBF)
  
##   tmp <- new("snpmod")
##   tmp@models=data.frame(str=models[-1],
##     size=sapply(snps,length)[-1],
##     logABF=lBF)
##   tmp@model.snps <- snps[-1]
##   return(marg.snps(tmp))
}

abf.tag.notworking <- function(models,sm,BFcolumn="jeffreys",xR2=NULL,snp.data=NULL,return.R2=FALSE,verbose=FALSE) { # raftery, wen

  ## initialize return object
  lBF <- matrix(NA,nrow=length(models),ncol=4,dimnames=list(models,c(BFcolumn,"r2.min","r2.mean","r2.max")))

  ## models that have been fitted
  Mint <- intersect(models,sm@models$str)
  lBF[Mint,BFcolumn] <- sm@models[ match(Mint,sm@models$str), BFcolumn ]

  ## find tags for remainder
  Mdiff <- setdiff(models,sm@models$str)
  snps.missing <- strsplit(Mdiff,"%")
  snps.eval <- strsplit(Mint,"%")
  usnps.missing <- unique(unlist(snps.missing))
  usnps.eval <- unique(unlist(snps.eval))
  if(is.null(R2))
    R2 <- ld(snp.data[,usnps.missing],snp.data[,usnps.eval],stat="R.squared")

  
      snps.tags <- colnames(R2)[ apply(R2,1,which.max) ]
    names(snps.tags) <- snps.missing
    models.alt <- models.missing
    for(snp in snps.missing) {
      models.alt <- sub(paste0("(%|^)",snp,"(%|$)"), paste0("\\1",snps.tags[snp],"\\2"), models.alt)
    }

  snps.missing <- setdiff(unlist(strsplit(models.missing,"%")),colnames(x))

  ## prepare data.frame
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(is(x,"data.frame"))
    x <- as.matrix(x)

  if(method=="glm.fit") {
    x2<-cbind(one=1,x[,intersect(unique(unlist(snps)),colnames(x))])
    comp <- complete.cases(x2) & !is.na(y)
    if(!all(comp)) {
      x2 <- x2[comp,]
      y2 <- y[comp]
    } else {
      y2 <- y
    }
    logn <- log(nrow(x))
    family <- switch(family,
                     "gaussian"=gaussian(link="identity"),
                     "binomial"=binomial(link="logit"))
    print(family)
    snps <- lapply(snps,setdiff,"1")
  }
  
  if(method=="glm") {
    df <- as.data.frame(x[,intersect(unique(unlist(snps)),colnames(x))])
    df$y <- y                      
    comp <- complete.cases(df)
    if(!all(comp)) {
      df2 <- df[comp,]
    } else {
      df2 <- df
    }
  }
  
  ## generate bics
  message("calculating ABFs")
  if(method=="glm.fit") {
    bics <- mclapply(seq_along(models), function(i) {
      if(verbose && i %% 100 == 0)
        cat(i,"\t")
      k=length(snps[[i]])+1
      glm.fit(x2[,c("one",snps[[i]])], y2, family=family)$deviance + k*logn
    })
  }
  if(method=="glm") {
    ## is baseline already included?
  l <- sapply(snps,length)
  for(wh in which(l==0))
    snps[[wh]] <- "1"
  
  bics <- mclapply(seq_along(models), function(i) {
      if(verbose && i %% 100 == 0)
        cat(i,"\t")
      f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
      BIC(glm(f, data=df2, family=family))
    })
  }
  
  bics <- unlist(bics)
  use <- use[-1] # drop base model
  models <- models[-1]
  cat("use: ",sum(use),", bics: ",length(bics)-1,"\n")
  lBF[use] <- (bics[-1] - bics[1]) / (-2)

  if(!all(use) && (!is.null(snp.data) || !is.null(R2))) {
    message("attempting to find tag models for those with missing SNPs")    
    models.missing <- setdiff(models.orig,models)
    snps.missing <- setdiff(unlist(strsplit(models.missing,"%")),colnames(x))
    allsnps <- c(snps.missing,colnames(x))
    if(is.null(R2))
      R2 <- ld(snp.data[,snps.missing],
               snp.data[,intersect(colnames(x),colnames(snp.data))],
               stat="R.squared")
    snps.tags <- colnames(R2)[ apply(R2,1,which.max) ]
    names(snps.tags) <- snps.missing
    models.alt <- models.missing
    for(snp in snps.missing) {
      models.alt <- sub(paste0("(%|^)",snp,"(%|$)"), paste0("\\1",snps.tags[snp],"\\2"), models.alt)
    }

    ## have any of these been done already?
    wh <- which(models.alt %in% models)
    if(length(wh)) {
      message(length(wh)," tagged by models already evaluated")
      lBF[ models.missing[wh] ] <- lBF[ models.alt[wh] ]                                       
    }

    ## do any need to be run fresh?
    models.todo <- unique(setdiff(models.alt,models))
    message("remainder can be tagged by ",length(models.todo)," needing evaluation")
    snps <- strsplit(models.todo,"%")

    if(method=="glm.fit") {
      snps <- lapply(snps,setdiff,"1")
      x2<-cbind(one=1,x[,intersect(unique(unlist(snps)),colnames(x))])
      comp2 <- complete.cases(x2) & !is.na(y)
      if(!all(comp2) || !all(comp)) {
        keep <- sample(which(comp2),size=sum(comp),replace=sum(comp) > sum(comp2))
        x2 <- x2[comp,]
        y2 <- y[comp]
      } # otherwise existing y2 and this x2 both ok
      bics.todo <- mclapply(seq_along(models.todo), function(i) {
        k=length(snps[[i]])+1
        glm.fit(x2[,c("one",snps[[i]])], y2, family=family)$deviance + k*logn
      })
    }
    if(method=="glm") {
     ## is baseline already included?
  l <- sapply(snps,length)
  for(wh in which(l==0))
    snps[[wh]] <- "1"
  
   df <- as.data.frame(x[,intersect(unique(unlist(snps)),colnames(x))])
      df$y <- y                      
      comp2 <- complete.cases(df)
      if(!all(comp2) || !all(comp)) {
        keep <- sample(which(comp2),size=sum(comp),replace=sum(comp) > sum(comp2))
        df2 <- df[keep,]
      } # otherwise existing y2 and this x2 both ok
      bics.todo <- mclapply(seq_along(models.todo), function(i) {
        f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
        BIC(glm(f, data=df2, family=family))
      })
    }

    
    lbf.todo <- (unlist(bics.todo) - bics[1]) / (-2)
    names(lbf.todo) <- models.todo
    
    m <- match(models.alt,models.todo)
    lBF[ models.missing[ !is.na(m) ] ] <- lbf.todo[ m[!is.na(m)]  ]    
    
  }

  lBF <- data.frame(model=models.orig,tag=models.orig %in% models,lBF=lBF,stringsAsFactors=FALSE)
  if(return.R2 && exists("R2"))
    return(list(lBF=lBF,R2=R2))
  return(lBF)
  
##   tmp <- new("snpmod")
##   tmp@models=data.frame(str=models[-1],
##     size=sapply(snps,length)[-1],
##     logABF=lBF)
##   tmp@model.snps <- snps[-1]
##   return(marg.snps(tmp))
}

abf2snpmod <- function(abf,prior=snpprior(x=0:20,expected=3,n=932,truncate=20,overdispersion=2)) {
  tmp <- new("snpmod")
  msize <- nchar(gsub("[^%]","",abf$model))
  mprior <- prior[as.character(msize)]
  mpp <- log(mprior) + abf$lBF
  mpp <- mpp - logsum(mpp)
  tmp@models <- data.frame(str=abf$model,
                           logABF=abf$lBF,
                           by.tagging=abf$tag,
                           size=msize,
                           prior=mprior,
                           lPP=mpp,
                           PP=exp(mpp),
                           stringsAsFactors=FALSE)
  tmp@model.snps <- strsplit(tmp@models$str,"%")  
  marg.snps(tmp)                      
}
