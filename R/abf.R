abf.calc <- function(y,x,models,family="binomial",backend=c("BIC"),R2=NULL,snp.data=NULL,return.R2=FALSE) { # raftery, wen
  ## initialize return object
  lBF <- structure(rep(as.numeric(NA),length(models)),names=models)

  ## add a base model
  models.orig <- models
  models <- c("1",models)
  snps <- strsplit(models,"%")

  ## is baseline already included?
  l <- sapply(snps,length)
  for(wh in which(l==0))
    snps[[wh]] <- "1"
  
  ## restrict to models that can be fitted
  cols.ok <- c("1",colnames(x))
  use <- sapply(snps, function(snp.names) all(snp.names %in% cols.ok))
  message(sum(use)-1," of ", length(models)-1, " models can be attempted")
  if(!all(use)) {
    snps <- snps[use]
    models <- models[use]
  }  
  
  ## prepare data.frame
  if(is(x,"SnpMatrix"))
    x <- matrix(as(x,"numeric"),nrow=nrow(x),dimnames=dimnames(x))
  df <- as.data.frame(x[,intersect(unique(unlist(snps)),colnames(x))])
  df$y <- y                      

  ## generate bics
  message("calculating ABFs")
  bics <- mclapply(seq_along(models), function(i) {
    f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
    BIC(glm(f, data=df, family=family))
  })

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
    df <- as.data.frame(x[,intersect(unique(unlist(snps)),colnames(x))])
    df$y <- y                      
    bics.todo <- mclapply(seq_along(models.todo), function(i) {
      f <- as.formula(paste("y ~",paste(snps[[i]],collapse="+")))
      BIC(glm(f, data=df, family=family))
    })
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
