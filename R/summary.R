##' guess.summ
##'
##' @title Summarise a snpmod object 
##' @param results an object of class snpmod, or a named list of such objects
##' @param groups object of class groups.  If supplied, all SNPs here will be summarised, and grouped according to this structure
##' @param snps data.frame giving details of the SNPs
##' @param snp.data if groups is missing, tag groups are determined from this object of class SnpMatrix
##' @param tag.thr if groups is missing, threshold at which to tag
##' @param pp.thr if groups is missing, threshold above which SNPs are selected for summary.  
##' @param method if groups is missing, method to determine tag groups using heirarchical clustering, default is "complete"
##' @return data.frame
##' @author Chris Wallace
guess.summ <- function(results, groups=NULL, snps=NULL, snp.data=NULL, tag.thr=0.8, pp.thr=0.01, method="complete") {
  if(is.null(groups) && is.null(snp.data))
    stop("must supply either SNP groups or snp.data for fixed LD threshold groups")
  if(!is.list(results))
    results <- list(trait=results)
  if(!is.null(groups)) {
    all.snps <- make.names(unlist(snps(groups)))
    guess.snps <- lapply(results, function(x) subset(x@snps, var %in% all.snps))
  } else {
    guess.snps <- lapply(results, function(x) subset(x@snps,x@snps$Marg_Prob_Incl>pp.thr))
    all.snps <- make.names(unique(unlist(lapply(guess.snps, "[[", "var"))))
  }
  if(!is.null(groups)) {
    tags <- as(groups,"tags")
  } else {
    if(length(drop <- setdiff(all.snps, colnames(snp.data)))) {
      warning("dropping ",length(drop)," SNPs: ", paste(drop, collapse=" "))
      kk <- lapply(guess.snps, function(x) subset(x, make.names(var) %in% drop))
      kk <- kk[ sapply(kk,nrow)>0 ]
      print(kk)
      all.snps <- intersect(all.snps, colnames(snp.data)) # dangerous
      ##     if(!is.null(groups)) {
      ##       g <- lapply(groups@groups,setdiff,drop)
      ##       g <- g[ sapply(g,length)>0 ]
      ##       groups <- new("groups", groups=g, tags=sapply(g, "[[", 1))
      ##     }
    }
    colnames(snp.data) <- make.names(colnames(snp.data))
    tags <- tag(snp.data[,all.snps], tag.threshold=tag.thr, method=method)
  }
  
  ## first order tags
  cl <- unique(tags(tags))
  cl <- cl[ order(snps[cl, "position.hg19"]) ]
  ## then order snps within tag groups
  tsplit <- split(snps(tags), tags(tags))[ cl ]
  tsplit <- lapply(tsplit, function(x) x[ order(snps[x, "position.hg19"]) ])
  ## put it back together
  all.snps <- unlist(tsplit)
  df.snps <- data.frame(snp=all.snps,tag=rep(names(tsplit),sapply(tsplit,length)),
                        row.names=all.snps,snpnum=1:length(all.snps),
                        stringsAsFactors=FALSE)
  
  ## add SNP info if provided
  if(!is.null(snps)) {
    df.snps <- cbind(df.snps, snps[ as.character(df.snps$snp), ])
                                        #    df.snps$position.scale <- with(df.snps, xscale(position, c(min(snpnum),max(snpnum))))
  }
  
  ## set xmin, xmax for easy plotting
  tag.names <- unique(as.character(df.snps$tag))
  tmp <- tapply(df.snps$snpnum, df.snps$tag, min)
  df.snps[names(tmp),"x.min"] <- tmp
  tmp <- tapply(df.snps$snpnum, df.snps$tag, max)
  df.snps[names(tmp),"x.max"] <- tmp
  
  ## add pp results
  df.list <- vector("list",length(results))
  for(i in 1:length(results)) {
    df <- data.frame(x=1:length(all.snps), pp=results[[i]]@snps[all.snps,"Marg_Prob_Incl"], snp=df.snps$snp, tag=df.snps$tag,
                     row.names=df.snps$snp)
    ppsum <- tapply(df$pp,df$tag,sum, na.rm=TRUE)
    df[names(ppsum),"ppsum"] <- ppsum
    df.list[[i]] <- cbind(df.snps, trait=names(results)[i])
    df.list[[i]]$pp[rownames(df)] <- df$pp
    df.list[[i]]$ppsum[rownames(df)] <- df$ppsum  
  }
  df.snps <- do.call("rbind",df.list)
  df.snps$pp[ is.na(df.snps$pp) ] <- 0
  df.snps$ppsum[ is.na(df.snps$ppsum) ] <- 0
                                        # df.snps[ is.na(df.snps) ] <- 0
  df.snps$snp <- make.names(df.snps$snp)
  return(df.snps)
}

guess.bootsum <- function(summ) {
  ## reshape wide
  marg <- do.call("rbind",with(summ, tapply(pp, snp, cummean)))
  mmarg <- melt(marg)
  ss <- unique(summ[,c("snp","tag","snpnum","position","position.scale","x.min","x.max")])
  ss$pp <- marg[rownames(ss),ncol(marg)]
  ppsum <- tapply(ss$pp,ss$tag,sum, na.rm=TRUE)
  ss[names(ppsum),"ppsum"] <- ppsum
  ss$trait <- "bootsumm"
  mplot <- ggplot(mmarg,aes(x=X2,col=X1,y=value)) + geom_point() + geom_path() + theme(legend.position="none")
  return(list(bootsumm=ss,
              plot=mplot))
}

makemod <- function(snps) {
  snps <- strsplit(snps,"%")
  all.snps <- unique(unlist(snps))
  mod <- Matrix(0,length(snps),length(all.snps),dimnames=list(NULL,all.snps))
  snum <- lapply(snps, function(s) which(all.snps %in% s))
  I <- rep(1:length(snum),times=sapply(snum,length))
  J <- unlist(snum)
  mod[cbind(I,J)] <- 1
  mod
}

marg.snps <- function(d) {
    mod <- makemod(d@models$str)
    marg.pp <- (d@models$PP %*% mod)[1,,drop=TRUE]
    d@snps <- marg.snps.vecs(d@models$str, d@models$PP)
    return(d)
}

marg.snps.vecs <- function(str,pp) {
    mod <- makemod(str)
    marg.pp <- (pp %*% mod)[1,,drop=TRUE]
    data.frame(Marg_Prob_Incl=marg.pp, var=names(marg.pp), rownames=names(marg.pp), stringsAsFactors=FALSE)
}

best.snps <- function(d,pp.thr=0.1) {
  tmp <- subset(d@snps, d@snps$Marg_Prob_Incl>pp.thr)
  return(tmp[order(tmp$Marg_Prob_Incl,decreasing=TRUE),])
}
best.models <- function(d,pp.thr=0.01,cpp.thr=NA) {
  if(!is.na(cpp.thr)) {
    d@models <- d@models[ order(d@models$PP,decreasing=TRUE),] # just in case
    cpp <- cumsum(d@models$PP)
    wh <- which(cpp<=cpp.thr)
  } else {
    wh <- which(d@models$PP>pp.thr)
    wh <- wh[ order(d@models$PP[wh],decreasing=TRUE) ]
  }
  return(cbind(d@models[wh,], snps=unlist(lapply(d@model.snps[wh],makestr))))
}
