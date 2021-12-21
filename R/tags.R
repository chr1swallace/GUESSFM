##' We can use tagging to reduce the number of models.  This function expands each tagged model to every possible model that could be created from its component (tag) SNPs.
##'
##' @title Expand tags for a snpmod object
##' @param d snpmod object
##' @param tags tags object
##' @import data.table
##' @return snpmod object with tags expanded
##' @export
expand.tags <- function(d, tags) {
    best <- d@model.snps
    B <- length(best)
    bsnps <- unique(unlist(best))
    wh <- which(make.names(tags(tags)) %in% bsnps)
    if(!length(wh))
      stop("none of the supplied tags are amongst the best SNPs in d")
    proxies <- split(make.names(snps(tags)[wh]),make.names(tags(tags)[wh]))    

    ## check - all best SNPs should be in names of proxies
    if(!all(bsnps %in% names(proxies)))
      stop("not all model SNPs found in tags object")
##    table(sapply(proxies,length))
    message("expanding tags for ",B," models over ",length(proxies)," tag SNPs, tagging a total of ",length(unlist(proxies)), " SNPs.")

    ## map snp -> num to save memory, hopefully
    allsnps <- unique(c(names(proxies),unlist(proxies)))
    allsnps <- structure(seq_along(allsnps), names=allsnps)
    ## check - names of proxies should correspond to 1..length(proxies)
    if(!all(allsnps[names(proxies)]==seq_along(proxies)))
        stop("cannot map proxies to integers in expand.tags - this shouldn't happen")
    ## Bnum <- lapply(best, function(b) allsnps[b])
    ## Pnum <- lapply(proxies, function(p) allsnps[p])

    ## modified from 
    ## https://stat.ethz.ch/pipermail/r-help/2006-February/087972.html
    ## expand.grid.str  <- function(id, vars) {
    ##     nv <- length(vars)
    ##     lims <- sapply(vars,length)
    ##     stopifnot(length(lims) > 0, id <= prod(lims), length(names(vars)) == nv)
    ##     res <- structure(vector("list",nv), .Names = names(vars))
    ##     if (nv > 1) for(i in nv:2) {
    ##                     f <- prod(lims[1:(i-1)])
    ##                     res[[i]] <- vars[[i]][(id - 1)%/%f + 1]
    ##                     id <- (id - 1)%%f + 1
    ##                 }
    ##     res[[1]] <- vars[[1]][id]
    ##     res <- do.call("cbind",res)
    ##     apply(res,1,makestr)
    ## }

  
    ## prepare expanded models
    pm <- mclapply(as.list(1:B), function(i) {
        ##  r2 <- ld(XX[,setdiff(colnames(X),best[[i]])],XX[,best[[i]]],stat="R.squared")
        ## all proxy models
        if(d@models[i,"size"]==0) {
            pm.str <- ""
        } else {
            #pm.num <- apply(do.call(expand.grid,Pnum[ best[[i]] ]), 1, function(x) paste(allsnps[sort(x)],collapse="%"))
           # microbenchmark({
                tmp <- t(do.call(expand.grid,
                                 c(proxies[best[[i]]], stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)))
                pm.str <- character(ncol(tmp))
                for(j in seq_along(pm.str))
                    pm.str[j] <- paste(sort(as.vector(tmp[,j])),collapse="%")
            ## }, {
            ##     tmp <- t(do.call(expand.grid,
            ##                      c(proxies[best[[i]]], stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)))
            ##     pm.str3 <- apply(tmp,2,function(x) makestr(x))
            ## }, {
            ##     pm.str2 <- apply(do.call(expand.grid,
            ##                              c(proxies[best[[i]]],
            ##                                stringsAsFactors=FALSE,
            ##                                KEEP.OUT.ATTRS=FALSE)),
            ##                      1, makestr)
            ## })
            ## all.equal(pm.str3,pm.str)
            ## system.time({
            ##     toexpand <- proxies[best[[i]]]
            ##     l <- sapply(toexpand,length)
            ##     pm.str2 <- expand.grid.str(1:prod(l), toexpand)
            ## })
        }
        #message(i)
        #gc()
        return(pm.str)
    })

    ## tie these to their index models
    npm <- sapply(pm,length)
    neighb <- as.data.table(d@models)[rep(1:length(pm),times=npm),]
    ##colnames(neighb) <- sub("str","index.str",colnames(neighb))
    setnames(neighb,sub("str","index.str",names(neighb)))
    neighb[,str:=unlist(pm)]
    neighb[,index:=str==index.str]

#    neighb <- do.call("rbind",neighb)
##     vars <- colnames(neighb[[1]])
##     neighb2 <- structure(vector("list",length(vars)),names=vars)
##     for(i in seq_along(vars))
##         neighb2[[i]] <- unlist(lapply(neighb,"[[",i))
##     neighb <- as.data.frame(neighb2, stringsAsFactors=FALSE)

    ## normalise, PP should sum to 1
    neighb[,logPP:=logPP - logsum(logPP)]
    neighb[,PP:=exp(neighb$logPP)]
    neighb <- neighb[ order(neighb$PP,decreasing=TRUE),]
    neighb[,rank:=1:nrow(neighb)]
    d@models <- as.data.frame(neighb)
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
##' @param method method used for heirarchical clustering.  See hclust for options.
##' @param split.at for large numbers of SNPs, tag first splits them into subsets of size split.at with similar MAF, then groups tag groups in high LD between subsets.  If you wish to avoid this, at additional computational cost, set split.at=NULL
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL, strata=NULL,quiet=FALSE,method="single",split.at=500) {
  if(!is(X,"SnpMatrix"))
    X <- as(X,"SnpMatrix")
  
  if(!is.null(snps) && !is.null(samples)) {
    X <- X[samples,snps]
  } else {
    if(!is.null(snps))
      X <- X[,snps]
    if(!is.null(samples))
      X <- X[samples,]
  }
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

  cs <- col.summary(X)
  if(any(is.na(cs[,"z.HWE"])))
    stop("some SNPs have zero standard deviation (missing HWE stat). Please fix and rerun")
  ## X <- X[,order(cs$MAF)]
  

  ## if too large, split first, then join
  if(!is.null(split.at) && ncol(X)>split.at) {
    cs <- col.summary(X)
    Q <- quantile(cs$MAF,seq(0,1,length=ceiling(ncol(X)/split.at)+1))
    maf <- cut(cs$MAF,breaks=unique(Q),include.lowest=TRUE)
    tags <- mclapply(levels(maf), function(l) {
      wh <- which(maf==l)
      if(!length(wh))
        return(NULL)
      return(tag(X, snps=wh, tag.threshold=tag.threshold, method=method, split.at=NULL))
    })

    ## merge any high-LD groups between tag sets
    TAGS <- tags[[1]]
    for(i in 2:length(tags)) {
      tg.0 <- unique(tags(tags[[i-1]]))
      tg.1 <- unique(tags(tags[[i]]))
      r2 <-   ld(X[,tg.0],
                 X[,tg.1],
                 symmetric=TRUE,
                 stats="R.squared")
      wh <- which(r2 >= tag.threshold, arr.ind=TRUE)
      if(nrow(wh)) {
        wh <- wh[!duplicated(wh[,2]),,drop=FALSE] # just in case
        for(j in 1:nrow(wh)) {
          tags[[i]]@tags[ tags[[i]]@tags==tg.1[wh[j,2]] ] <- tg.0[wh[j,1]]
        }
      }
      TAGS <- new("tags", .Data=c(TAGS@.Data,tags[[i]]@.Data),
                  tags=c(TAGS@tags,tags[[i]]@tags))
    }
    return(TAGS)
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
  groups <- new("groups",groups,tags=names(groups))
  
  ## check
  r2 <- r2[tags(groups),tags(groups), drop=FALSE]
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
           stats="R.squared")
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

