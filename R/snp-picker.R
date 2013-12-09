
library(gridExtra)
snp.picker <- function(d,data,pp.thr=0.01,nochange.thr=0.001,nochange.run=3,r2.gap=0.1) {
  groups <- plotsdata <- list()
  i <- 0
  a <- d@snps
  a$var <- as.character(a$var)
  max.r2 <- 0.5
  ## check all snps in a are in data
  if(!all(a$var %in% colnames(data)))
    stop("not all SNPs in a found in data")
  while(any(a$Marg_Prob_Incl>pp.thr)) {
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
    
    df <- data.frame(r2=1-r2[a$var], mpi=a$Marg_Prob_Incl)
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
          if((df2$r2bin[k+1] - df2$r2bin[k])>r2.gap && df2$mpi[k+1]>pp.thr) {
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
