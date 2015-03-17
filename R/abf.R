##' Uses the BIC approximation to calculate approximate Bayes factors for specified models
##'
##' The central idea of GUESSFM is to use GUESS to rapidly survery the model space for a tagged version of the data, and select a set of plausible models.  From their, the models tagged by the top model from GUESS should be evaluated using abf.calc.
##'
##' The use of a parallel file enables abf.calc to be run in two ways.  If you name a file that does not exist, objects will be saved to that file to enable subsets of models to be fitted in a parallel fashion.  If you name a file that exists, it is assumed to be the joined results of such model fits and will be loaded.  Without a parallel file, all models will be fitted, which may take a long time.  See vignette for more information.
##'
##' @title Calculate Approximate Bayes Factors
##' @param y response vector
##' @param x explanatory variables, a SnpMatrix object, data.frame or numeric matrix
##' @param models vector of models to consider
##' @param family family argument to pass to glm.  Currently only "binomial" and "gaussian" are implemented for the \code{glm.fit} method.
##' @param q optional vector of covariates to include in all models
##' @param method use the speedglm library (if available), glm.fit or glm to fit the models.  Default is speedglm which should be faster, but, like glm.fit, is less forgiving about the odd missing value etc, so if your code is giving potentially glm related errors, try running with \code{method="glm"} in the first instance.
##' @param R2 matrix giving pairwise r-squared measures of LD from which tag SNPs will be calculated when not directly available
##' @param snp.data if R2 is missing, it is calculated from this SnpMatrix object
##' @param return.R2 if true, return the calculated R2 matrix.  Useful if you are analysing several strata of a population and you wish to avoid repeating the calculation.
##' @param verbose print lots of progress messages if TRUE.  Default is FALSE.
##' @param parallel.file optional file name to enable parallelisation. 
##' @return a data.frame containing model name, ABF, and an indicator of whether the ABF was calculating directly or via a tag SNP
##' @author Chris Wallace
abf.calc <- function(y,x,models,family="binomial",
                     q=NULL,method=c("speedglm","glm.fit","glm"),
                     R2=NULL,snp.data=NULL,return.R2=FALSE,verbose=FALSE,
                     parallel.file=NULL) { # raftery, wen
  
  method <- match.arg(method)
  message("Calculating BICs using ",method)
  
  models[ models=="" ] <- "1"
  models.orig <- models
  
  ## models that can be fitted directly
  if(verbose)
    message("Finding models which can be evaluated directly...")
  snps <- strsplit(models,"%")
  cols.ok <- c("1",colnames(x))
  use <- sapply(snps, function(snp.names) all(snp.names %in% cols.ok))
  message(sum(use)," of ", length(models), " models can be evaluated directly")
  if(!all(use)) {
    snps <- snps[use]
    models <- models[use]
  }  
  
  ## models that need to be assessed through tagging
  models.alt <- NULL
  if(!all(use) && (!is.null(snp.data) || !is.null(R2))) {
    message("Attempting to find tag models for those with missing SNPs...")    
    models.missing <- setdiff(models.orig,models)
    snps.missing <- setdiff(unlist(strsplit(models.missing,"%")),colnames(x))
    allsnps <- c(snps.missing,colnames(x))
    if(is.null(R2))
      R2 <- ld(snp.data[,snps.missing],
               snp.data[,intersect(colnames(x),colnames(snp.data))],
               stats="R.squared")
    snps.tags <- colnames(R2)[ apply(R2,1,which.max) ]
    names(snps.tags) <- snps.missing
    models.alt <- models.missing
    for(snp in snps.missing) {
      models.alt <- sub(paste0("(%|^)",snp,"(%|$)"), paste0("\\1",snps.tags[snp],"\\2"), models.alt)
    }
    names(models.alt) <- models.missing
    models <- unique(c(models.alt,models))
    message("Including tags, a total of ",length(models)," need evaluation")
  }
  
  models <- c("1",models)
  snps <- strsplit(models,"%")
  
  ## load prepared results?
  if(!is.null(parallel.file) && file.exists(parallel.file)) {
    message("Found ",parallel.file,", loading results")
    (load(parallel.file))
  } else { ## prepare models
    results <- switch(method,
                      speedglm=abf.speedglm.fit(x,y,q,family,snps,parallel.file),
                      glm.fit=abf.glm.fit(x,y,q,family,snps,parallel.file),
                      glm=abf.glm(x,y,q,family,snps))
    if(is.null(results))
      if(return.R2 && !is.null(R2)) {
        return(R2)
      } else {
        return(NULL)
      }
  }
  
  ## summarize ABF
  ##  use <- use[-1] # drop base model
  message("calculating ABF")
  models <- models[-1]
  lBF.fits <- (results$bics[-1] - results$bics[1]) / (-2)
  names(lBF.fits) <- models
  names(results$coeff) <- c("one",models)
##   results <- results[-1]
##   names(results) <- models
  
  ## initialize return object
  lBF <- structure(rep(as.numeric(NA),length(models.orig)),names=models.orig)
  coeff <- structure(vector("list",length(models.orig)),names=models.orig)
  
  ## directly fitted
  mint <- intersect(names(lBF),models)
  message("evaluating direct models, ",length(mint)," found.")
  if(length(mint)) {
    lBF[ mint ] <- lBF.fits[ mint ]
    coeff[ mint ] <- results$coeff[ mint ]
  }
  ## tags
  if(!is.null(models.alt)) {
    models.alt.missing <- models.alt[models.missing]
    message("evaluating tag models: ",length(models.missing)," can be found through ",length(unique(models.alt.missing))," models.")
    if(length(models.alt.missing)) {
      lBF[ models.missing ] <- lBF.fits[ models.alt.missing ]
      coeff[ models.missing ] <- results$coeff[ models.alt.missing ]
    }
  }
  
  ## mint <- intersect(names(lBF.fits),models.alt[models.missing])
  ## if(length(mint)) {
  ##   m0 <- models.missing[match(mint,models.alt)]
  ##   lBF[ m0 ] <- lBF.fits[ mint  ]
  ##   coeff[ m0 ] <- results$coeff[ mint ]
  ## }
  
  lBF.df <- data.frame(model=models.orig,
                       tag=!(models.orig %in% names(lBF)),
                       lBF=lBF[models.orig],
                       stringsAsFactors=FALSE)
  if(return.R2 && exists("R2"))
    return(list(lBF=lBF.df,coeff=coeff,R2=R2))
  return(list(lBF=lBF.df,coeff=coeff))
  
  ##   tmp <- new("snpmod")
  ##   tmp@models=data.frame(str=models[-1],
  ##     size=sapply(snps,length)[-1],
  ##     logABF=lBF)
  ##   tmp@model.snps <- snps[-1]
  ##   return(marg.snps(tmp))
}



abf2snpmod <- function(abf,prior=snpprior(x=0:20,expected=3,n=932,truncate=20,overdispersion=1)) {
  tmp <- new("snpmod")
  msize <- nchar(gsub("[^%]","",abf$model)) + 1
  mprior <- prior[as.character(msize)]
  mpp <- log(mprior) + abf$lBF
  mpp <- mpp - logsum(mpp)
  message("creating data.frame")
  tmp@models <- data.frame(str=abf$model,
                           logABF=abf$lBF,
                           by.tagging=abf$tag,
                           size=msize,
                           prior=mprior,
                           lPP=mpp,
                           PP=exp(mpp),
                           stringsAsFactors=FALSE)
  tmp@model.snps <- strsplit(tmp@models$str,"%")
  message("calculating marginal SNP inclusion probabilities")
  marg.snps(tmp)                      
}

abf.manual.join <- function(parallel.file, ...) {
  (load(parallel.file))
  L <- as.list(...)
  for(results.file in L) {
    (load(results.file))    
  }
  
}

abf.manual <- function(parallel.file,targets,bic.file,coeff.file,verbose=FALSE) {
  
  if(!is.numeric(targets))
    stop("targets must be a numeric vector indexing which models to fit")
  message("Loading data from ",parallel.file)
  ## to circumvent R CMD check NOTE:
  x2 <- y2 <- NULL
  (load(parallel.file))
##  snps <- snps[targets]
##   l <- sapply(snps,length)
##   snps[[which(l==0)]] <- "one"
  message("Fitting ",length(targets)," models")
  results <- abf.glm.fit(x=x2,y=y2,q=q,family=family,snps=snps[targets])
  bics <- results$bics
  coeff <- results$coeff
  message("Saving results to ",bic.file," and ",coeff.file)
  save(targets,bics,file=bic.file)
  save(targets,coeff,file=coeff.file)
  
}

abf.speedglm.fit <- function(x,y,q,family,snps,parallel.file=NULL,verbose=FALSE) {
  
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(is(x,"data.frame"))
    x <- as.matrix(x)
  if(!is.null(q)) {
    if(is.vector(q)) {
      qm <- model.matrix(~q)
    } else {
      qm <- cbind(one=1,as.matrix(q))
    }
  } else {
    qm <- matrix(1,nrow=nrow(x),ncol=1,dimnames=list(NULL,"one"))
  }
  x2<-x[,intersect(unique(unlist(snps)),colnames(x))]
  comp <- complete.cases(x2) & !is.na(y)
  if(!is.null(q))
    comp <- comp & complete.cases(qm)
  if(!all(comp)) {
    message("Dropping ",sum(!comp)," samples due to incompleteness. ",sum(comp)," remain.")
    x2 <- x2[comp,]
    y2 <- y[comp]
    qm <- qm[comp,,drop=FALSE]
    q <- q[comp]
  } else {
    y2 <- y
  }
  logn <- log(nrow(x))
##   if(verbose)
##     print(family)
  snps <- lapply(snps,setdiff,"1")
  
  ## check
  allsnps <- unique(unlist(snps))
  if(!all(allsnps %in% colnames(x2)))
    stop("Not all SNPs found")
  
  if(is.null(parallel.file)) {
  if(is.character("family"))
    family <- switch(family,
                     "gaussian"=gaussian(link="identity"),
                     "binomial"=binomial(link="logit"))
    results <- lapply(seq_along(snps), function(i) {
      if(verbose && i %% 100 == 0)
        cat(i,"\t")
      k=length(snps[[i]])+ncol(qm)
      model <- speedglm.wfit(y2, cbind(x2[, snps[[i]] ],qm), family=family)
      model0 <- speedglm.wfit(y2, qm, family=family)
#      class(model) <- class(model0) <- c(class(model),"glm")
      A1 <- AIC(model) - 2*k + k*log(length(y2))
      A0 <- AIC(model0) - 2*ncol(qm) + ncol(qm)*log(length(y2))
      list(BIC=A1-A0,
           coeff=cbind(beta=model$coefficients,
             se=sqrt(diag(vcov(model)))))
    })    
  } else {
    if(file.exists(parallel.file)) {
      load(parallel.file)
    } else {
      message("Saving objects in ",parallel.file)
      save(snps, x2, y2, qm, q, family, file=parallel.file)
      message("Please fit the models using abf.manual and rerun with parallel.file")
      return(NULL)
    }
  }
#  print(results)
  bics <- unlist(lapply(results, "[[", "BIC"))
  coeff  <- lapply(results, "[[", "coeff")
  return(list(bics=bics,coeff=coeff))
  
}
abf.glm.fit <- function(x,y,q,family,snps,parallel.file=NULL,verbose=FALSE) {
  
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(is(x,"data.frame"))
    x <- as.matrix(x)
  if(!is.null(q)) {
    if(is.vector(q)) {
      qm <- model.matrix(~q)
    } else {
      qm <- cbind(one=1,as.matrix(q))
    }
  } else {
    qm <- matrix(1,nrow=nrow(x),ncol=1,dimnames=list(NULL,"one"))
  }
  x2<-x[,intersect(unique(unlist(snps)),colnames(x))]
  comp <- complete.cases(x2) & !is.na(y)
  if(!is.null(q))
    comp <- comp & complete.cases(qm)
  if(!all(comp)) {
    message("Dropping ",sum(!comp)," samples due to incompleteness. ",sum(comp)," remain.")
    x2 <- x2[comp,]
    y2 <- y[comp]
    qm <- qm[comp,,drop=FALSE]
    q <- q[comp]
  } else {
    y2 <- y
  }
  logn <- log(nrow(x))
##   if(verbose)
##     print(family)
  snps <- lapply(snps,setdiff,"1")
  
  ## check
  allsnps <- unique(unlist(snps))
  if(!all(allsnps %in% colnames(x2)))
    stop("Not all SNPs found")
  
  if(is.null(parallel.file)) {
  if(is.character("family"))
    family <- switch(family,
                     "gaussian"=gaussian(link="identity"),
                     "binomial"=binomial(link="logit"))
    results <- lapply(seq_along(snps), function(i) {
      if(verbose && i %% 100 == 0)
        cat(i,"\t")
      k=length(snps[[i]])+1
      model <- glm.fit(cbind(x2[, snps[[i]] ],qm), y2, family=family)
      model0 <- glm.fit(qm, y2, family=family)
      class(model) <- class(model0) <- c(class(model),"glm")
      list(BIC=BIC(model) - BIC(model0),
           coeff=cbind(beta=model$coefficients,
             se=sqrt(diag(vcov(model)))))
    })    
  } else {
    if(file.exists(parallel.file)) {
      load(parallel.file)
    } else {
      message("Saving objects in ",parallel.file)
      save(snps, x2, y2, qm, q, family, file=parallel.file)
      message("Please fit the models using abf.manual and rerun with parallel.file")
      return(NULL)
    }
  }
#  print(results)
  bics <- unlist(lapply(results, "[[", "BIC"))
  coeff  <- lapply(results, "[[", "coeff")
  return(list(bics=bics,coeff=coeff))
  
}

abf.glm <- function(x,y,q,family,snps,verbose=FALSE) {
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(!is.null(q))
    qm <- model.matrix(~q)
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
  
  results <- mclapply(seq_along(snps), function(i) {
    if(verbose && i %% 100 == 0)
      cat(i,"\t")
    if(!is.null(q)) {
      f <- as.formula(paste("y ~ ",paste(snps[[i]],collapse="+"), "+ q"))
    } else {
      f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
    }
    model <- glm(f, data=df2, family=family)
    list(BIC=BIC(model),
         coeff=cbind(beta=model$coefficients,
           se=sqrt(diag(vcov(model)))))
  })
  bics <- unlist(lapply(results, "[[", "BIC"))
  coeff  <- unlist(lapply(results, "[[", "coeff"))
  return(list(bics=bics,coeff=coeff))

}


