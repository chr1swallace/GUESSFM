##' Group SNPs by r-squared and their influence on the posterior
##'
##' In the presence of LD, the posterior is diluted over correlated
##' SNPs.  We deal with this by picking groups of SNPs according to r
##' squared and their influence on the posterior, as measured by the
##' marginal probability of inclusion (MPI).  The aim is to identify
##' groups of SNPs, at most one of which should be included in any
##' model.
##'
##' The algorithm is rather simplistic, and works according to a
##' series of tuning parameters detailed below.  It appears to work
##' for the datasets I have tried, but I would be very interested to
##' hear about its performance on alternative datasets, whether
##' positive or negative!  It may be possible to set sensible defaults
##' for specific parameters according to the type of dataset.
##' 
##' By default, the method checks explicitly that SNPs in the same
##' group do not occur together in models.  This requires creation of
##' a model matrix, which is a time consuming step. It can be avoided
##' by setting skip.shared.models=TRUE, but this is not advised as SNP
##' groups should be sets of SNPs, either none or one of which are
##' needed in the model.
##' 
##' @title snp.picker
##' @param d object of class \code{snpmod} containing SNPs to be grouped
##' @param data genotype data as a \code{SnpMatrix} object which will be used to determine r-squared
##' @param start.thr MPI threshold to identify new groups.  As long as at least one ungrouped SNP remains with MPI>start.thr the algorithm will attempt to form new groups
##' @param nochange.thr the threshold below which MPI will be regarded as unchanged
##' @param nochange.run the number of SNPs with MPI < nochange.thr for the algorithm to decide a group is complete
##' @param r2.gap if MPI falls below nochange.thr, then at the next SNP rises above start.thr and if that SNP is > r2.gap away, the group is terminated
##' @param shared.models the maximum proportion of models for the group's index SNP in which another group SNP may belong
##' @param skip.shared.models if TRUE ignore the limit on shared.models
##' @return object of class \code{snppicker}
##' @author Chris Wallace
##' @export
snp.picker <- function(d,data,start.thr=0.01,nochange.thr=0.001,nochange.run=3,r2.gap=0.1,shared.models=0.03,skip.shared.models=FALSE) {
  groups <- plotsdata <- list()
  i <- 0
  a <- d@snps
  a$var <- as.character(a$var)
  a <- subset(a,var!="1") # ignore null model
  max.r2 <- 0.5
  ## check all snps in a are in data
  if(!all(a$var %in% colnames(data)))
    stop("not all SNPs in d found in data")
  
  ## create model matrix
  if(!skip.shared.models) {
    M <- modelMatrix(d)
  } else {
    cs <- structure(numeric(length(a$var)),names=a$var)
  }
  
  while(any(a$Marg_Prob_Incl>start.thr)) {
    
    (i <- i+1)
    ## best snp
    wh <- which.max(a$Marg_Prob_Incl)
    (snp <- a$var[wh])

    ## r2 with best snp
    ## system.time({
    ##   r2 <- ld(data[,setdiff(colnames(data),snp)], data[,snp], stats="R.squared")[,1];
    ##   r2[snp] <- 1})
    r2 <- ld(data, data[,snp], stats="R.squared")[,1]
    if(all(is.na(r2))) {
      r2[is.na(r2)] <- 0
      r2[snp] <- 1
    }
    
    ## models that contain best snp
    if(!skip.shared.models) {
      Mbest <- M[ M[,snp]==1,,drop=FALSE]
      ## count of snps in models - proportion of models in which each SNP is present in addition to best snp
      cs <- structure(colSums(Mbest),names=colnames(Mbest))
      cs[snp] <- 0
    }
    ## summary df
    df <- data.frame(r2=1-r2[a$var], mpi=a$Marg_Prob_Incl,model.count=cs[a$var]/nrow(Mbest))
    df <- df[order(df$r2),]
    dfp <- df
    df <- df[df$r2<max.r2 & !is.na(df$r2),,drop=FALSE]
    if(length(wh <- which(df$model.count>shared.models)))
      df <- df[-wh,]
    print(df)

    ## process
    r2lim <- NULL
    if(nrow(df)==1) {     ## singletons
      r2lim <- 0
    } else {
      ## group into r2 bins
      df$r2bin <- round(df$r2,2)
      df2 <- data.frame(r2bin=unique(df$r2bin),
                        mpibin=tapply(df$mpi,df$r2bin,sum),
                        mcountbin=tapply(df$model.count,df$r2bin,max))
      wh <- which(df2$mpibin < nochange.thr)
      if(!length(wh)) { # use all
        r2lim <- max(df2$r2bin)
      } else {
        for(k in wh) { ## find break point
          ## k is candidate for first position beyond snp group
          if(all(df2$mpibin[k+seq(0,min(nochange.run-1, nrow(df2)-k))] < nochange.thr)) {
            ## have hit our run of low change, stop
            r2lim <- df2$r2bin[ k-1 ]
            break
          }
          if((df2$r2bin[k+1] - df2$r2bin[k])>r2.gap && df2$mpi[k+1]>start.thr) {
            ## k is low, then big gap to next SNP, which has a high enough mpi to be a new hit in its own right
            r2lim <- df2$r2bin[ k-1 ]
            break
          }                          
          if((df2$r2bin[k+1] - df2$r2bin[k])>r2.gap && df2$mpi[k+1]<nochange.thr*3) {
           ## k is low, then big gap to next SNP, which is nearly low enough to be ignored
            r2lim <- df2$r2bin[ k-1 ]
            break
          }
        }
      }
      if(is.null(r2lim)) {
          message("outside given rules! Stopping on row ",wh[1]-1)
          print(df2)
          r2lim <- df2$r2bin[ wh[1]-1 ]    
      }
      if(length(r2lim)==0)
          r2lim <- 0
    }    
    ## store
    wh <- which(round(df$r2,2)<=r2lim)
    snp.group <- rownames(df)[wh]
    groups[[i]] <- cbind(a[snp.group,,drop=FALSE],r2=df$r2[wh],cmpi=cumsum(df$mpi[wh]))
    cat(i,snp,length(snp.group),max(groups[[i]]$cmpi),"\n")
    a <- subset(a, !(a$var %in% snp.group))
    dfp$changepoint <- 1:nrow(dfp)==max(wh)
    plotsdata[[i]] <- dfp    

  }
  
  return(new("snppicker",groups=groups,plotsdata=plotsdata))
}
##' Check if SNP groups can be combined.
##'
##' snp.picker tries to automatically group SNPs according the r2 and
##' the model posterior probs, but may sometimes make two separate
##' groups when one combined group might be considered more
##' appropriate.  To check whether two groups can be combined, we need
##' to check that SNPs from both groups are rarely included in the
##' same models.
##'
##' If the output shows that pp["all"] >> pp["any"], and the SNPs have
##' been shown to be in some LD, then it is likely that the groups can
##' be appropriately merged, using \code{groups.merge}.
##' @title Check if snp groups can be combined 
##' @param d object of class snpmod
##' @param test.groups list of character vectors, each vector giving the SNPs that comprise one group
##' @return shows the posterior probability of models containing exactly one, any, or all of the SNPs.  
##' @author Chris Wallace
##' @export
check.merge <- function(d,test.groups) {
  inmod <- lapply(test.groups, function(query) {
    unlist(lapply(d@model.snps,function(target) as.numeric(any(query %in% target))))
  })
  inmod <- do.call("rbind",inmod)
  #print(dim(inmod)
  rownames(inmod) <- test.groups@tags
  cs <- colSums(inmod)
  inmod <- rbind(inmod,
                 any=ifelse(cs>0,1,0),
                 all=ifelse(cs==length(test.groups),1,0))
  return(inmod %*% d@models$PP)
}
##' Merge groups in an object of class groups
##'
##' Groups identified by the index tag SNPs given by \code{tags} will
##' be merged, and the number of groups in the object reduced.
##' @title merge groups
##' @param groups object of class groups
##' @param tags index tags indicating the groups which should be merged
##' @return object of class groups
##' @export
groups.merge <- function(groups,tags) {
  if(length(tags)<2)
    stop("makes no sense to merge <2 groups")
  if(!all(tags %in% groups@tags))
    stop("can't find all index tags in groups@tags")
  wh <- which(groups@tags %in% tags)
  drop <- wh[-1]
  keep <- wh[1]
  groups@.Data[[keep]] <- unlist(snps(groups)[wh])
  groups@.Data <- snps(groups)[-drop]
  groups@tags <- tags(groups)[-drop]
  return(groups)
}

##   df$cmpi <- cumsum(df$mpi)
##   ## group at points of inflection, where cmpi<=1
##   split.1 <- which(df$cmpi>1)[1] - 1

##   wh.max.r2 <- which(df$r2>max.r2)[1]-1
##   ij <- expand.grid(i=1:wh.max.r2,j=1:wh.max.r2)
##   ij <- ij[ ij$i<=ij$j , ]
##   ij$p <- sapply(1:nrow(ij), function(k) {
##     t.test(df$mpi[ij[k,"i"] : ij[k,"j"] ], mu=thr, alternative="greater")$p.value
##   })
  
                
  
##   runs <- unlist(lapply(1:wh.max.r2, function(i) {
##     lapply(i:wh.max.r2, function(j) {
##       i:j
##     })
##   }))
                 
  
##   library(bcp)
##   cp <- bcp(log(df$mpi/(1-df$mpi)),return.mcmc=TRUE)
##   plot(cp)
##   wh <- which(cp$posterior.prob<0.2)[1] - 1
##   wh <- 1:wh
## ##       x <- zoo(df$cmpi,order.by=seq_along(df$cmpi))
## ##       df$ma <- c(as.numeric(rollmean(x, k=3, align="left")),NA,NA)
## ##       ggplot(df,aes(x=r2)) + geom_point(aes(y=cmpi)) + geom_path(aes(y=ma))
## ##       df$dma <- c(diff(df$ma),NA)
## ##       library(cpm)
## ##       max.r2 <- 0.5
## ##       eps <- 1e-3
## ##       cp <- detectChangePoint(df$cmpi,cpmType="Mann-Whitney")
## ##       if(cp$changeDetected && df$r2[cp$changePoint - 1] < max.r2) {
## ##         wh <- 1:(cp$changePoint-1)
## ## ##         lt <- (df$mpi[wh]<eps)
## ## ##         keep <- which(rev(lt))
## ##         ## ? step backwards whilst diff < eps?        
## ##       } else {
## ##         wh <- 1
## ##       }
  
                                                  
## }

  

      
      
##       findipiterplot(df$r2,df$cmpi,index=0)

      
      
##       df$mpid <- c(NA,diff(df$mpi))
##         wh <- which(df$mpid>0)[1] - 1
##         df$group <- ifelse(rownames(df) %in% snp.group,snp,NA)
        
##         plots[[i]] <- ggplot(df,aes(x=r2,y=pp,col=group)) + geom_point() + ggtitle(snp)
## #        dev.off()
        
## ##         library(cluster)
## ##         gap <- clusGap(data.frame(r2=r2[a$snp],pp=a$Marg_Prob_Incl), kmeans, 10, B = 100, verbose = FALSE)
## ##         k <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
## ##         km <- kmeans        (r2[a$snp],centers=k)
## ##         snp.group <- names(km$cluster)[ km$cluster == km$cluster[snp] ]
##         qplot(r2[a$snp],a$Marg_Prob_Incl,col=as.factor(km$cluster)) + ggtitle(snp)
##         groups[[i]] <- subset(a, a$snp %in% snp.group)
##     }
    
## do.call("grid.arrange",plots)
    
## ## NEED TO BASE THIS ON WHETHER SNPS ARE ALTERNATIVE OR CUMULATIVE IN MODELS
    
## }

snp.picker.old <- function(d,data,start.thr=0.01,nochange.thr=0.001,nochange.run=3,r2.gap=0.1) {
  groups <- plotsdata <- list()
  i <- 0
  a <- d@snps
  a$var <- as.character(a$var)
  max.r2 <- 0.5
  ## check all snps in a are in data
  if(!all(a$var %in% colnames(data)))
    stop("not all SNPs in a found in data")
  while(any(a$Marg_Prob_Incl>start.thr)) {
    i <- i+1
    ## best snp
    wh <- which.max(a$Marg_Prob_Incl)
    (snp <- a$var[wh])
    r2 <- ld(data[,setdiff(colnames(data),snp)], data[,snp], stats="R.squared")[,1]
    r2[snp] <- 1
    ## if(mean(is.na(r2))>0.5) {
##       o <- order(a$Marg_Prob_Incl,decreasing=TRUE)
##       wh <- o[2]
##       snp <- a$var[wh]
##       print(a[wh,])
##       r2 <- ld(data[,setdiff(colnames(data),snp)], data[,snp], stats="R.squared")[,1]
##       r2[snp] <- 1
##     }  
    
    df <- data.frame(r2=1-r2[a$var], mpi=a$Marg_Prob_Incl, row.names=a$var)
    df <- df[order(df$r2),]
    dfp <- df
    df <- df[df$r2<max.r2 & !is.na(df$r2),,drop=FALSE]
    print(df)

    if(nrow(df)==1) {     ## singletons
      r2lim <- 0
    } else {

      ## group into r2 bins
      df$r2bin <- round(df$r2,2)
      df2 <- data.frame(r2bin=unique(df$r2bin),
                        mpibin=tapply(df$mpi,df$r2bin,sum))
      wh <- which(df2$mpibin < nochange.thr)
      if(!length(wh)) { # use all
        r2lim <- max(df2$r2bin)

      } else {

        for(k in wh) { ## find break point
          ## k is candidate for first position beyond snp group
          if(all(df2$mpibin[k+seq(0,min(nochange.run-1, nrow(df2)-k))] < nochange.thr)) {
            ## have hit our run of low change, stop
            r2lim <- df2$r2bin[ k-1 ]
            break
          }
          if((df2$r2bin[k+1] - df2$r2bin[k])>r2.gap && df2$mpi[k+1]>start.thr) {
            ## k is low, then big gap to next SNP, which has a high enough mpi to be a new hit in its own right
            r2lim <- df2$r2bin[ k-1 ]
            break
          }                          
        }
      }
    }
    
    wh <- which(round(df$r2,2)<=r2lim)
    snp.group <- rownames(df)[wh]
    groups[[i]] <- cbind(a[snp.group,,drop=FALSE],r2=df$r2[wh],cmpi=cumsum(df$mpi[wh]))
    cat(i,snp,length(snp.group),max(groups[[i]]$cmpi),"\n")
    a <- a[!(a$var %in% snp.group),]
    dfp$changepoint <- 1:nrow(dfp)==max(wh)
    plotsdata[[i]] <- dfp
  }
  return(new("snppicker",groups=groups,plotsdata=plotsdata))
}



##   df$cmpi <- cumsum(df$mpi)
##   ## group at points of inflection, where cmpi<=1
##   split.1 <- which(df$cmpi>1)[1] - 1

##   wh.max.r2 <- which(df$r2>max.r2)[1]-1
##   ij <- expand.grid(i=1:wh.max.r2,j=1:wh.max.r2)
##   ij <- ij[ ij$i<=ij$j , ]
##   ij$p <- sapply(1:nrow(ij), function(k) {
##     t.test(df$mpi[ij[k,"i"] : ij[k,"j"] ], mu=thr, alternative="greater")$p.value
##   })
  
                
  
##   runs <- unlist(lapply(1:wh.max.r2, function(i) {
##     lapply(i:wh.max.r2, function(j) {
##       i:j
##     })
##   }))
                 
  
##   library(bcp)
##   cp <- bcp(log(df$mpi/(1-df$mpi)),return.mcmc=TRUE)
##   plot(cp)
##   wh <- which(cp$posterior.prob<0.2)[1] - 1
##   wh <- 1:wh
## ##       x <- zoo(df$cmpi,order.by=seq_along(df$cmpi))
## ##       df$ma <- c(as.numeric(rollmean(x, k=3, align="left")),NA,NA)
## ##       ggplot(df,aes(x=r2)) + geom_point(aes(y=cmpi)) + geom_path(aes(y=ma))
## ##       df$dma <- c(diff(df$ma),NA)
## ##       library(cpm)
## ##       max.r2 <- 0.5
## ##       eps <- 1e-3
## ##       cp <- detectChangePoint(df$cmpi,cpmType="Mann-Whitney")
## ##       if(cp$changeDetected && df$r2[cp$changePoint - 1] < max.r2) {
## ##         wh <- 1:(cp$changePoint-1)
## ## ##         lt <- (df$mpi[wh]<eps)
## ## ##         keep <- which(rev(lt))
## ##         ## ? step backwards whilst diff < eps?        
## ##       } else {
## ##         wh <- 1
## ##       }
  
                                                  
## }

  

      
      
##       findipiterplot(df$r2,df$cmpi,index=0)

      
      
##       df$mpid <- c(NA,diff(df$mpi))
##         wh <- which(df$mpid>0)[1] - 1
##         df$group <- ifelse(rownames(df) %in% snp.group,snp,NA)
        
##         plots[[i]] <- ggplot(df,aes(x=r2,y=pp,col=group)) + geom_point() + ggtitle(snp)
## #        dev.off()
        
## ##         library(cluster)
## ##         gap <- clusGap(data.frame(r2=r2[a$snp],pp=a$Marg_Prob_Incl), kmeans, 10, B = 100, verbose = FALSE)
## ##         k <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"])
## ##         km <- kmeans        (r2[a$snp],centers=k)
## ##         snp.group <- names(km$cluster)[ km$cluster == km$cluster[snp] ]
##         qplot(r2[a$snp],a$Marg_Prob_Incl,col=as.factor(km$cluster)) + ggtitle(snp)
##         groups[[i]] <- subset(a, a$snp %in% snp.group)
##     }
    
## do.call("grid.arrange",plots)
    
## ## NEED TO BASE THIS ON WHETHER SNPS ARE ALTERNATIVE OR CUMULATIVE IN MODELS
    
## }
