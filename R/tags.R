add.neighbours <- function(d) {

  
  
}

expand.tags <- function(d, tags) {
    best <- d@model.snps
    B <- length(best)
    bsnps <- unique(unlist(best))
    wh <- which(make.names(tags@tags) %in% bsnps)
    proxies <- split(make.names(tags@snps[wh]),make.names(tags@tags[wh]))
    table(sapply(proxies,length))
    message("expanding tags for ",B," models over ",length(proxies)," tag SNPs, tagging a total of ",length(unlist(proxies)), " SNPs.")

    ## prepare expanded models
    pm <- mclapply(as.list(1:B), function(i) {
        ##  r2 <- ld(XX[,setdiff(colnames(X),best[[i]])],XX[,best[[i]]],stat="R.squared")
        ## all proxy models
        if(d@models[i,"size"]==0) {
            pm.str <- ""
        } else {
            pm <- do.call(expand.grid,proxies[best[[i]]])
            pm.str <- apply(pm,1,makestr)
        }
        return(pm.str)
    })

    ## tie these to their index models
    npm <- sapply(pm,length)
    tmp <- d@models[rep(i,each=npm),]
    colnames(tmp) <- sub("str","index.str",colnames(tmp))
    neighb <- cbind(data.frame(str=unlist(pm.str),
                            stringsAsFactors=FALSE),
                 tmp)
        
    neighb$index <- neighb$str==neighb$index.str
    rownames(neighb) <- NULL

#    neighb <- do.call("rbind",neighb)
##     vars <- colnames(neighb[[1]])
##     neighb2 <- structure(vector("list",length(vars)),names=vars)
##     for(i in seq_along(vars))
##         neighb2[[i]] <- unlist(lapply(neighb,"[[",i))
##     neighb <- as.data.frame(neighb2, stringsAsFactors=FALSE)

    ## normalise, PP should sum to 1
    neighb$logPP <- neighb$logPP - logsum(neighb$logPP)
    neighb$PP <- exp(neighb$logPP)
    neighb <- neighb[ order(neighb$PP,decreasing=TRUE),]
    neighb$rank <- 1:nrow(neighb)
    ##    return(neighb)

    d@models <- neighb
    d@model.snps <- strsplit(neighb$str,"%")
    return(marg.snps(d))

}


##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters,
##' then cuts the tree at 1-tag.threshold
##' @title tag
##' @param X SnpMatrix object
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param snps colnames of the SnpMatrix object to be used
##' @param samples optional, subset of samples to use
##' @param strata optional, return a list of tag vectors, one for each stratum defined by as.factor(strata)
##' @param quiet if FALSE (default), show progress messages
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL, strata=NULL,quiet=FALSE,method="complete") {
  if(!is(X,"SnpMatrix"))
    X <- as(X,"SnpMatrix")
  
  if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  if(!is.null(strata)) {
    strata <- factor(strata)
    tags <- lapply(levels(strata), function(l) {
      wh <- which(strata==l)
      if(!length(wh))
        return(NULL)
      return(tag(X[wh,], tag.threshold=tag.threshold, method=method))
    })
    return(tags)
  }
  
  r2 <- myr2(X)
  D <- as.dist(1-r2)
  hc <- hclust(D, method=method)
  clusters <- cutree(hc, h=1-tag.threshold)
  
  snps.use <- names(clusters)[!duplicated(clusters)]
  groups <- split(names(clusters),clusters)
  
  ## now process each group, picking best tag
  n <- sapply(groups,length)
  names(groups)[n==1] <- unlist(groups[n==1])
  for(i in which(n>1)) {
    g <- groups[[i]]
##     cat(i,g,"\n")
##     print(r2[g,g])
    a <- apply(r2[g,g],1,mean)
    names(groups)[i] <- g[ which.max(a) ]
  }
  groups <- new("groups",groups=groups,tags=names(groups))
  
  ## check
  r2 <- myr2(X[,groups@tags])
  diag(r2) <- 0
##   if(max(r2)==1) 
##     stop("max r2 still 1!")
  if(!quiet)
    message("max r2 is now",max(r2),"\n")
  return(as(groups,"tags"))
}
  
group.tags <- function(tags, keep) {
    groups <- tags[ names(tags) %in% keep ]
    groups <- split(names(groups), groups)
}
myr2 <- function(X) {
  r2 <- ld(X,
           depth=ncol(X)-1,
           symmetric=TRUE,
           stat="R.squared")
  if(any(is.na(r2))) {
    r2.na <- as(is.na(r2),"matrix")
    use <- rowSums(r2.na)>0
    ## work around for r2=NA bug.  
    r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
    r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
  }
  diag(r2) <- 1
  return(r2)
}

