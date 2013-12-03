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
            nsnps <- sapply(object@groups,nrow)
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
          function(x) {            
            plots <- lapply(seq_along(x@plotsdata), function(i) {
              pdata <- x@plotsdata[[i]]
              nsel <- which(pdata$changepoint==TRUE)
              cmpi <- sum(pdata[1:nsel,"mpi"])
              ggplot(pdata, aes(x=r2,y=cumsum(mpi))) + 
#                geom_hline(yintercept=thr,col="lightblue",linetype="dashed") +
                  geom_vline(mapping=aes(xintercept=r2), data=subset(pdata,changepoint==TRUE),col="red",linetype="dashed") +
                    geom_point() + geom_path() +
                      geom_point(aes(y=mpi),col="blue") + geom_path(aes(y=mpi),col="blue") +
                        ylim(0,1.2) + xlab("rsq with index SNP") + ylab("MPI/cum.MPI") +
                          ggtitle(paste(nsel,"SNPs; cum.MPI =",signif(cmpi,2)))
            })
            print(do.call("tracks",plots))
            invisible(plots)            
          })

################################################################################

## convert between groups and tags

################################################################################

setAs("groups", "tags",
      def=function(from) {
        new("tags",tags=rep(from@tags,times=sapply(from@groups,length)),snps=unlist(from@groups))
      })
setAs("tags", "groups",
      def=function(from) {
        new("groups",tags=unique(from@tags),groups=split(from@snps,from@tags))
      })
setAs("snppicker","groups",
      def=function(from) {
        new("groups",
            tags=unlist(lapply(from@groups,function(x) x[1,"var"])),
            groups=lapply(from@groups, "[[", "var"))
      })
setAs("snppicker","tags",
      def=function(from) {
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
  ## find intersecting groups, and make unions
  int <- matrix(FALSE,length(x@groups),length(y@groups))
  for(i in seq_along(x@groups)) {
    for(j in seq_along(y@groups)) {      
      int[i,j] <- any(x@groups[[i]] %in% y@groups[[j]])
    }
  }
  wh <- which(int,arr.ind=TRUE)
  if(nrow(wh)) {
    for(i in 1:nrow(wh)) {
      xi <- wh[i,1]
      yi <- wh[i,2]
      x@groups[[xi]] <- unique(c(x@groups[[xi]],y@groups[[yi]]))
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
