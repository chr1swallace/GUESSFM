utils:::globalVariables(c("is.sparse","cp","control","ll.speedglm"))

##' Uses the BIC approximation to calculate approximate Bayes factors for specified models
##'
##' The central idea of GUESSFM is to use GUESS to rapidly survery the model space for a tagged version of the data, and select a set of plausible models.  From their, the models tagged by the top model from GUESS should be evaluated using abf.calc.
##'
##' The use of a parallel file enables abf.calc to be run in two ways.  If you name a file that does not exist, objects will be saved to that file to enable subsets of models to be fitted in a parallel fashion.  If you name a file that exists, it is assumed to be the joined results of such model fits and will be loaded.  Without a parallel file, all models will be fitted, which may take a long time, particularly for glms.  See vignette for more information.
##'
##' After abf.calc, you may want to use \code{\link{abf2snpmod}} to generate a snpmod with the same structed returned by \code{\link{read.snpmod}}.
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
##' @param parallel.dir optional directory name to enable manual parallelisation. 
##' @return a data.frame containing model name, ABF, and an indicator of whether the ABF was calculating directly or via a tag SNP
##' @param approx.lm default FALSE. If TRUE and family=="binomial", fit lm on the residuals from a binomial logistic regression of y on covars.  This is the same approximation used by GUESS, but its use here is experimental and should not be used unless you are prepared to sanity check results yourself (eg by running binomial for a selection of models and comparing ABF from the two approaches)
##' @export
##' @author Chris Wallace
abf.calc <- function(y,x,models,family="binomial",
                     q=NULL,method=c("speedglm","glm.fit","glm"),
                     R2=NULL,snp.data=NULL,return.R2=FALSE,verbose=FALSE,
                     parallel.dir=NULL,
                     approx.lm=FALSE) { # raftery, wen
  
  method <- match.arg(method)
    message("Calculating BICs using ",method)

    ## speedglm should check input type for y, but doesn't, so do it here
    if(is.null(nrow(x)))
        stop(paste("x should be a matrix: samples=rows, snps=columns"))
     if(!is(y,"vector") || nrow(x)!=length(y))
        stop(paste("y should be a vector of length equal to the rows in x"))
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
#    snps.direct <- snps[use]
    models.direct <- models[use]
  }  else {
      models.direct <- models
  }
  names(models.direct) <- models.direct
  
  ## models that need to be assessed through tagging
  models.indirect <- NULL
  if(!all(use) && (!is.null(snp.data) || !is.null(R2))) {
    message("Attempting to find tag models for those with missing SNPs...")    
    models.missing <- setdiff(models.orig,models.direct)
    snps.missing <- setdiff(unlist(strsplit(models.missing,"%")),colnames(x))
    allsnps <- c(snps.missing,colnames(x))
    if(is.null(R2))
      R2 <- ld(snp.data[,snps.missing],
               snp.data[,intersect(colnames(x),colnames(snp.data))],
               stats="R.squared")
    snps.tags <- colnames(R2)[ apply(R2,1,which.max) ]
    names(snps.tags) <- snps.missing
    models.indirect <- models.missing
    for(snp in snps.missing) {
      models.indirect <- sub(paste0("(%|^)",snp,"(%|$)"), paste0("\\1",snps.tags[snp],"\\2"), models.indirect)
    }
    names(models.indirect) <- models.missing
    models <- unique(c(models.indirect,models.direct))
    message("Including tags, a total of ",length(models)," need evaluation")
  }
  
  models <- c(one="1",models)
  snps <- c(one="1",snps) #strsplit(models,"%")

    if(approx.lm && family=="binomial") {
        m <- glm(y ~ ., data=as.data.frame(q))
        y <- residuals(m)
        family <- "gaussian"
        q <- NULL
    }
    
  ## load prepared results?
  parallel.results <- presults(parallel.dir)
  if(!is.null(parallel.dir) && file.exists(parallel.results)) {
      results <- abf.fit.parallel.gather(parallel.dir)
  } else {
      ## prepare models and run
      results <- abf.speedglm.fit(x,y,q,family,snps,parallel.dir)
    ## switch(method,
    ##                   speedglm=abf.speedglm.fit(x,y,q,family,snps,parallel.dir),
    ##                   glm.fit=abf.glm.fit(x,y,q,family,snps,parallel.dir),
    ##                   glm=abf.glm(x,y,q,family,snps))
    if(!is.list(results)) {
      ret <- list(N=results)
      if(return.R2 && !is.null(R2))
        ret$R2 <- R2
      return(ret)
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
  #mint <- intersect(names(lBF),models)
  message("evaluating direct models, n=",length(models.direct))
  if(length(models.direct)) {
    lBF[ models.direct ] <- lBF.fits[ models.direct ]
    coeff[ models.direct ] <- results$coeff[ models.direct ]
  }
  ## tags
  if(!is.null(models.indirect)) {
      message("evaluating tag models: ",length(models.indirect)," tagged by ",length(unique(models.indirect))," unique models.")
      lBF[ names(models.indirect) ] <- lBF.fits[ models.indirect ]
      coeff[ names(models.indirect) ] <- results$coeff[ models.indirect ]
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
  return(list(lBF=lBF.df,coeff=coeff,nsamp=nrow(x),nsnp=ncol(x)))
  
  ##   tmp <- new("snpmod")
  ##   tmp@models=data.frame(str=models[-1],
  ##     size=sapply(snps,length)[-1],
  ##     logABF=lBF)
  ##   tmp@model.snps <- snps[-1]
  ##   return(marg.snps(tmp))
}


##' Convert an abf object to a snpmod
##'
##' @title abf2snpmod
##' @param abf object returned by abf.calc
##' @param nsnps number of SNPs in the region, optional, but required if not found in output of abf
##' @inheritParams snpprior
##' @return a snpmod
##' @export
##' @author Chris Wallace
##' @seealso \link{snpprior}
abf2snpmod <- function(abf,expected,overdispersion=1,nsnps=NULL) {
  tmp <- new("snpmod")
  msize <- nchar(gsub("[^%]","",abf$lBF$model)) + 1
  msize[abf$lBF$model=="1"] <- 0
  if(!is.null(nsnps))
      abf$nsnp <- nsnps
  if(is.null(abf$nsnp))
      stop("please supply the number of SNPs in the region")
  prior <- snpprior(x=0:max(msize),expected=expected,
                    n=abf$nsnp,
                    truncate=max(c(msize,20)),
                    overdispersion=overdispersion)
  mprior <- prior[as.character(msize)]
  mpp <- log(mprior) + abf$lBF$lBF
  mpp <- mpp - logsum(mpp)
  message("creating snpmod data.frame")
  tmp@models <- data.frame(str=abf$lBF$model,
                           logABF=abf$lBF$lBF,
                           by.tagging=abf$lBF$tag,
                           size=msize,
                           prior=mprior,
                           lPP=mpp,
                           PP=exp(mpp),
                           stringsAsFactors=FALSE)
  tmp@model.snps <- strsplit(tmp@models$str,"%")
  message("calculating marginal SNP inclusion probabilities")
  marg.snps(tmp)                      
}

fill.abf.list <- function(n,files) {
  L <- numeric(n)
  targets <- NULL # R CMD check warning
  for(f in files) {
    obj <- (load(f))
    obj <- setdiff(obj,"targets")
    L[targets] <- get(obj)
  }
  return(L)
}

abf.manual.join <- function(parallel.dir, ...) {
  if(!file.exists(parallel.dir))
    stop("parallel.dir not found: ",parallel.dir)
  bic.files <- list.files(parallel.dir,pattern="^bic-.*.RData",full.names=TRUE)
  coeff.files <- sub("bic","coeff",bic.files)
  parallel.data <- pdata(parallel.dir)
  parallel.results <- presults(parallel.dir)
  if(!file.exists(parallel.data))
    stop("parallel data file not found, did you run abf.calc?")
  if(!length(bic.files))
    stop("no BIC/coeff files found in ",parallel.dir)
  targets <- bics <- NULL # R CMD check warning
  (load(parallel.data))

  ## bics 
  message("joining bics")
  BICS <- numeric(length(snps))
  for(f in bic.files) {
    (load(f))
    BICS[targets] <- bics
  }
    
  ## coeffs
  message("joining coeffs")
  COEFFS <- vector("list",length(snps))
  for(f in coeff.files) {
    (load(f))
    coeff <- lapply(coeff, function(x) {x[!grepl("^one|^q|Intercept",rownames(x)),,drop=FALSE]})
    COEFFS[targets] <- coeff
  }
  models <- unlist(lapply(snps,paste,collapse="%"))
  names(BICS) <- names(COEFFS) <- models
  results <- list(bics=BICS,coeff=COEFFS)
  save(results, models, file=parallel.results)
}

abf.manual <- function(parallel.dir,start,end) {
  
  if(!is.numeric(start) | !is.numeric(end))
    stop("start and end must be integers indexing which models to fit")
  parallel.data <- pdata(parallel.dir)
  message("Loading data from ",parallel.data)
  ## to circumvent R CMD check NOTE:
  x2 <- y2 <- NULL
  (load(parallel.data))
  targets <- start:end
  if(any(targets<0))
    stop("targets must be positive integers")
  if(min(targets)>length(snps))
    stop("targets out of range")
  if(max(targets)>length(snps))
    targets <- targets[targets<=length(snps)]
  t1 <- min(targets)
  t2 <- max(targets)
  bic.file <- file.path(parallel.dir,paste0("bic-",t1,"-",t2,".RData"))
  coeff.file <- file.path(parallel.dir,paste0("coeff-",t1,"-",t2,".RData"))
##  snps <- snps[targets]
##   l <- sapply(snps,length)
##   snps[[which(l==0)]] <- "one"
  message("Fitting ",length(targets)," models")
  results <- abf.speedglm.fit(x=x2,y=y2,q=q,family=family,snps=snps[targets])
  bics <- results$bics
  coeff <- results$coeff
  message("Saving results to ",bic.file," and ",coeff.file)
  save(targets,bics,file=bic.file)
  save(targets,coeff,file=coeff.file)
  
}
##' Internal function for fitting abf
##'
##' @title abf.speedglm.fit
##' @param x design matrix
##' @param y response vector
##' @param q matrix of covariates
##' @param family glm family
##' @param snps list of snps models
##' @param parallel.dir optional directory used for manual parallelisation
##' @param verbose if TRUE, print progress. Default is FALSE
##' @return a list of bics and coeffs, or number of models to be fit if parallel.dir is not NULL
##' @author Chris Wallace
abf.speedglm.fit <- function(x,y,q,family,snps,parallel.dir=NULL,verbose=FALSE) {
  
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  if(is(x,"data.frame"))
    x <- as.matrix(x)
  if(is.vector(q) || is.data.frame(q))
    q <- as.data.frame(q)
  x2<-x[,intersect(unique(unlist(snps)),colnames(x))]
  comp <- complete.cases(x2) & !is.na(y)
  if(!is.null(q))
    comp <- comp & complete.cases(q)
  if(!all(comp)) {
    message("Dropping ",sum(!comp)," samples due to incompleteness. ",sum(comp)," remain.")
    x2 <- x2[comp,]
    y2 <- y[comp]
    if(!is.null(q))
      q <- q[comp,,drop=FALSE]
  } else {
    y2 <- y
  }
  if(!is.null(q)) {
    if(is.data.frame(q))
      qm <- model.matrix(~., data=q)
    if(is.vector(q))
      qm <- model.matrix(~q)
    if(is.matrix(q))
      qm <- cbind(one=1,as.matrix(q))
  } else {
    qm <- matrix(1,nrow=nrow(x2),ncol=1,dimnames=list(NULL,"one"))
  }

  logn <- log(nrow(x2))
##   if(verbose)
##     print(family)
  snps <- lapply(snps,setdiff,c("1","one"))
  
  ## check
  allsnps <- unique(unlist(snps))
  if(!all(allsnps %in% colnames(x2)))
    stop("Not all SNPs found")
  
  if(is.null(parallel.dir)) { # run it
    if(is.character(family))
      family <- switch(family,
                       "gaussian"=gaussian(link="identity"),
                       "binomial"=binomial(link="logit"))
    model0 <- speedglm.wfit(y2, qm, family=family)
    k0=ncol(qm)
    A0 <- AIC(model0) - 2*k0 + k0*log(length(y2))
    aadj <- log(length(y2)) - 2
    results <- mclapply(seq_along(snps), function(i) {
      ## if(verbose && i %% 100 == 0)
      ##   cat(i,"\t")
      k=length(snps[[i]]) + k0
      model <- speedglm.wfit(y2, cbind(x2[, snps[[i]], drop=FALSE ],qm), family=family)
      A1 <- AIC(model) + k * aadj # - 2*k + k*log(length(y2))
      list(BIC=A1-A0,
           coeff=cbind(beta=model$coefficients,
                       se=sqrt(diag(vcov(model)))))
    })    
  } else {
      abf.fit.parallel.setup(parallel.dir, snps, x2, y2, q, family)
      return(length(snps))      
  }
  abf.fit.postprocess(results)
}

abf.fit.parallel.gather <- function(parallel.dir) {
    parallel.results <- presults(parallel.dir)
    if(!file.exists(parallel.results))
        stop("results file not found: ",parallel.results)
    message("Loading results from ",parallel.results)
    results <- NULL # avoid R CMD check warning
    load(parallel.results)
    return(results)
}
##' Internal function: create parallel.dir, save necessary files
##'
##' @param parallel.dir directory to create
##' @param snps snps to model
##' @param x2 genotype data
##' @param y2 phenotype data
##' @param q covariates
##' @param family glm family
##' @param ... items to save in data.RData under parallel
##' @return No return value
##' @author chris
abf.fit.parallel.setup <- function(parallel.dir, snps, x2, y2, q, family) {
    parallel.data <- pdata(parallel.dir)
    if(!file.exists(parallel.dir)) {
        message("creating parallel dir ",parallel.dir)
        dir.create(parallel.dir)
    }
    message("Saving objects in ",parallel.data)
    save(snps, x2, y2, q, family, file=parallel.data)
    message("Please fit the models using abf.manual and rerun with parallel.dir")
} 
abf.fit.postprocess <- function(results) {
    bics <- unlist(lapply(results, "[[", "BIC"))
    coeff  <- lapply(results, "[[", "coeff")
    return(list(bics=bics,coeff=coeff))
}

fasterglm.fit <- function (y, X, intercept = TRUE, row.chunk = NULL, 
    family = gaussian(), start = NULL, etastart = NULL, mustart = NULL, 
    offset = NULL, acc = 1e-08, maxit = 25, k = 2, sparselim = 0.9, 
    camp = 0.01, eigendec = TRUE, tol.values = 1e-07, tol.vectors = 1e-07, 
    tol.solve = .Machine$double.eps, sparse = NULL, method = c("eigen", 
        "Cholesky", "qr"), trace = FALSE, ...) 
{
    nobs <- NROW(y)
    nvar <- ncol(X)
    offset <- rep.int(0, nobs)
    weights <- rep(1, nobs)
    col.names <- dimnames(X)[[2]]
    method <- match.arg(method)
    fam <- family$family
    link <- family$link
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (is.null(sparse)) 
        sparse <- is.sparse(X = X, sparselim, camp)
    if (is.null(start)) {
        if (is.null(mustart)) 
            eval(family$initialize)
        eta <- if (is.null(etastart)) 
            family$linkfun(mustart)
        else etastart
        mu <- mustart
        start <- rep(0, nvar)
    }
    else {
        eta <- offset + as.vector(if (nvar == 1) 
            X * start
        else {
            if (sparse) 
                X %*% start
            else tcrossprod(X, t(start))
        })
        mu <- linkinv(eta)
    }
    iter <- 0
    dev <- sum(dev.resids(y, mu, weights))
    tol <- 1
    if ((fam == "gaussian") & (link == "identity")) 
        maxit <- 1
    C_Cdqrls <- getNativeSymbolInfo("Cdqrls", PACKAGE = getLoadedDLLs()$stats)
    while ((tol > acc) & (iter < maxit)) {
        iter <- iter + 1
        beta <- start
        dev0 <- dev
        varmu <- variance(mu)
        mu.eta.val <- mu.eta(eta)
        z <- (eta - offset) + (y - mu)/mu.eta.val
        W <- (mu.eta.val * mu.eta.val)/varmu
        XTX <- cp(X, W, row.chunk, sparse)
        XTz <- if (sparse) 
            t(X) %*% (W * z)
        else t(crossprod((W * z), X))
        if (iter == 1 & method != "qr") {
            variable <- colnames(X)
            ris <- if (eigendec) 
                control(XTX, , tol.values, tol.vectors, , method)
            else list(rank = nvar, pivot = 1:nvar)
            ok <- ris$pivot[1:ris$rank]
            if (eigendec) {
                XTX <- ris$XTX
                X <- X[, ok]
                XTz <- XTz[ok]
                start <- start[ok]
            }
            beta <- start
        }
        if (method == "qr") {
            ris <- .Call(C_Cdqrls, XTX, XTz, tol.values, FALSE)
            start <- if (ris$rank < nvar) 
                ris$coefficients[ris$pivot]
            else ris$coefficients
        }
        else {
            start <- solve(XTX, XTz, tol = tol.solve)
        }
        eta <- if (sparse) 
            drop(X %*% start)
        else drop(tcrossprod(X, t(start)))
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        tol <- max(abs(dev0 - dev)/(abs(dev) + 0.1))
        if (trace) 
            cat("iter", iter, "tol", tol, "\n")
    }
    wtdmu <- if (intercept) 
        sum(y)/nobs
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs
    nulldf <- n.ok - as.integer(intercept)
    rank <- ris$rank
    dfr <- nobs - rank
    aic.model <- aic(y, nobs, mu, weights, dev) + k * rank
    ll.nuovo <- ll.speedglm(fam, aic.model, rank)
    res <- (y - mu)/mu.eta(eta)
    resdf <- n.ok - rank
    RSS <- sum(W * res * res)
    var_res <- RSS/dfr
    dispersion <- if (fam %in% c("poisson", "binomial")) 
        1
    else var_res
    if (method == "qr") {
        coefficients <- start
        coefficients[coefficients == 0] = NA
        ok <- ris$pivot[1:rank]
    }
    else {
        coefficients <- rep(NA, nvar)
        start <- as(start, "numeric")
        coefficients[ok] <- start
    }
    names(coefficients) <- col.names
    rval <- list(coefficients = coefficients, logLik = ll.nuovo, 
        iter = iter, tol = tol, family = family, link = link, 
        df = dfr, XTX = XTX, dispersion = dispersion, ok = ok, 
        rank = rank, RSS = RSS, method = method, aic = aic.model, 
        sparse = sparse, deviance = dev, nulldf = nulldf, nulldev = nulldev, 
        ngoodobs = n.ok, n = nobs, intercept = intercept, convergence = (!(tol > 
            acc)))
    class(rval) <- "speedglm"
    rval
}


## abf.glm.fit <- function(x,y,q,family,snps,parallel.file=NULL,verbose=FALSE) {
  
##   if(is(x,"SnpMatrix"))
##     x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
##   if(is(x,"data.frame"))
##     x <- as.matrix(x)
##   if(!is.null(q)) {
##     if(is.vector(q)) {
##       qm <- model.matrix(~q)
##     } else {
##       qm <- cbind(one=1,as.matrix(q))
##     }
##   } else {
##     qm <- matrix(1,nrow=nrow(x),ncol=1,dimnames=list(NULL,"one"))
##   }
##   x2<-x[,intersect(unique(unlist(snps)),colnames(x))]
##   comp <- complete.cases(x2) & !is.na(y)
##   if(!is.null(q))
##     comp <- comp & complete.cases(qm)
##   if(!all(comp)) {
##     message("Dropping ",sum(!comp)," samples due to incompleteness. ",sum(comp)," remain.")
##     x2 <- x2[comp,]
##     y2 <- y[comp]
##     qm <- qm[comp,,drop=FALSE]
##     q <- q[comp]
##   } else {
##     y2 <- y
##   }
##   logn <- log(nrow(x))
## ##   if(verbose)
## ##     print(family)
##   snps <- lapply(snps,setdiff,"1")
  
##   ## check
##   allsnps <- unique(unlist(snps))
##   if(!all(allsnps %in% colnames(x2)))
##     stop("Not all SNPs found")
  
##   if(is.null(parallel.file)) {
##   if(is.character("family"))
##     family <- switch(family,
##                      "gaussian"=gaussian(link="identity"),
##                      "binomial"=binomial(link="logit"))
##     results <- mclapply(seq_along(snps), function(i) {
##       if(verbose && i %% 100 == 0)
##         cat(i,"\t")
##       k=length(snps[[i]])+1
##       model <- glm.fit(cbind(x2[, snps[[i]] ],qm), y2, family=family)
##       model0 <- glm.fit(qm, y2, family=family)
##       class(model) <- class(model0) <- c(class(model),"glm")
##       list(BIC=BIC(model) - BIC(model0),
##            coeff=cbind(beta=model$coefficients,
##              se=sqrt(diag(vcov(model)))))
##     })    
##   } else {
##     if(file.exists(parallel.file)) {
##       load(parallel.file)
##     } else {
##       message("Saving objects in ",parallel.file)
##       save(snps, x2, y2, qm, q, family, file=parallel.file)
##       message("Please fit the models using abf.manual and rerun with parallel.file")
##       return(NULL)
##     }
##   }
## #  print(results)
##   bics <- unlist(lapply(results, "[[", "BIC"))
##   coeff  <- lapply(results, "[[", "coeff")
##   return(list(bics=bics,coeff=coeff))
  
## }

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


pdata <- function(d) file.path(d,"data.RData")
presults <- function(d) file.path(d,"results.RData")
