######################################################################

##                          display                                 ##

######################################################################

setMethod("show", signature="snpmod",
          function(object) {
            nmod <- nrow(object@models)
            nsnp <- nrow(object@snps)
            maxpp <- max(object@models$PP)
            minpp <- min(object@models$PP)
            spp <- sum(object@models$PP)
            message("snpmod object, containing information on ",nmod," models / ",nsnp," SNPs.")
            message(sprintf("PP ranges from %4.3f-%4.3f (sum: %4.3f).",minpp,maxpp,spp))
          })

setMethod("show", signature="snppicker",
          function(object) {
            ngroup <- length(object@groups)
            nsnps <- if(ngroup==0) { 0 } else { sapply(object@groups,nrow) }
            message("snppicker object, containing ",sum(nsnps)," SNPs grouped into ",ngroup," groups.")
          })

setMethod("show", signature="tags",
          function(object) {
            ntags <- length(unique(object@tags))
            nsnps <- length(object@snps)
            message("tags object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

setMethod("show", signature="groups",
          function(object) {
            ntags <- length(object@tags)
            nsnps <- length(unlist(object@groups))
            message("groups object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

######################################################################

##                          plotting                                ##

######################################################################

setMethod("plot", signature(x="snppicker",y="missing"),
          function(x, do.plot=TRUE) {
            index <- unlist(lapply(x@groups, function(z) rownames(z)[1]))
            plots <- lapply(seq_along(x@plotsdata), function(i) {
              pdata <- x@plotsdata[[i]]
              nsel <- which(pdata$changepoint==TRUE)
              tmpi <- sum(pdata[1:nsel,"mpi"])
              pdata$cmpi <- cumsum(pdata$mpi)
              ggplot(pdata, aes(x=r2,y=cmpi)) + 
#                geom_hline(yintercept=thr,col="lightblue",linetype="dashed") +
                  geom_vline(mapping=aes(xintercept=r2), data=subset(pdata,changepoint==TRUE),col="red",linetype="dashed") +
                  geom_hline(mapping=aes(yintercept=cmpi), data=subset(pdata,changepoint==TRUE),col="red",width=0.2) +
                    geom_point() + geom_path() +
                      geom_point(aes(y=mpi),col="blue") + geom_path(aes(y=mpi),col="blue") +
                        ylim(0,1.2) + xlab("rsq with index SNP") + ylab("MPI/cum.MPI") +
                          ggtitle(paste(nsel,"SNPs; cum.MPI =",signif(tmpi,2))) + theme_bw()
            })
            names(plots) <- index
            if(do.plot)
              print(do.call("tracks",c(plots,label.text.cex=0.5)))
            invisible(plots)            
          })

################################################################################

## convert between groups and tags

################################################################################

setAs("groups", "tags",      
      def=function(from) {
        if(!length(from@groups))
          return(new("tags")) # return empty object
        new("tags",tags=rep(from@tags,times=sapply(from@groups,length)),snps=unlist(from@groups))
      })
setAs("tags", "groups",
      def=function(from) {
        if(!length(tags@tags))
          return(new("groups"))
        new("groups",tags=unique(from@tags),groups=split(from@snps,from@tags))
      })
setAs("snppicker","groups",
      def=function(from) {
        if(!length(from@groups))
          return(new("groups")) # return empty object
        new("groups",
            tags=unlist(lapply(from@groups,function(x) x[1,"var"])),
            groups=lapply(from@groups, "[[", "var"))
      })
setAs("snppicker","tags",
      def=function(from) {
        if(!length(from@groups))
          return(new("tags")) # return empty object
        g <- new("groups",
            tags=unlist(lapply(from@groups,function(x) x[1,"var"])),
            groups=lapply(from@groups, "[[", "var"))
        as(g,"tags")
      })
          
######################################################################

##                        concatenate                               ##

######################################################################

setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {
  union(as(x,"groups"),as(y,"groups")) })
setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
  as(union(as(x,"groups"),as(y,"groups")),"tags") })
setMethod("union",signature(x="groups",y="groups"),definition=function(x,y) {
  if(!length(x@groups))
    return(y)
  if(!length(y@groups))
    return(x)
  ## find intersecting groups, and make unions
  int <- matrix(FALSE,length(x@groups),length(y@groups))
  for(i in seq_along(x@groups)) {
    for(j in seq_along(y@groups)) {      
      int[i,j] <- any(x@groups[[i]] %in% y@groups[[j]])
    }
  }
  wh <- which(int,arr.ind=TRUE)
  if(nrow(wh)) {
    ## merge y into x
    for(i in 1:nrow(wh)) {
      xi <- wh[i,1]
      yi <- wh[i,2]
      x@groups[[xi]] <- unique(c(x@groups[[xi]],y@groups[[yi]]))
    }
    ## check - duplicate y's mean x groups must be merged
    dups <- wh[,2][ duplicated(wh[,2]) ]
    if(length(dups)) {
      ## do the merge
      for(ydup in dups) {
        xi <- wh[ wh[,2]==ydup, 1]
        x@groups[[ xi[1] ]] <- unique(unlist(x@groups[xi]))
      }
      ## drop NULLs now we don't need wh any more
      to.drop <- vector("list",length(dups))
      for(i in seq_along(dups)) {
        yi <- dups[[i]]
        to.drop[[i]] <- (wh[ wh[,2]==ydup, 1])[-1]
      }
      to.drop <- unique(unlist(to.drop))
      x@groups <- x@groups[ -to.drop ]
      x@tags <- x@tags[ -to.drop ]
    }
    
    add <- setdiff(seq_along(y@groups),unique(wh[,2]))
    if(!length(add))
      return(x)
    return(new("groups",groups=c(x@groups,y@groups[add]),tags=c(x@tags,y@tags[add])))
  }
  ## no overlap: simply concatenate
  new("groups",groups=c(x@groups,y=y@groups),tags=c(x@tags,y@tags))
})

######################################################################

##                      add two snpmods                             ##

######################################################################

## setMethod("+",signature(e1="snpmod",e2="snpmod"), function(e1,e2) snpmod.add(e1,e2))
