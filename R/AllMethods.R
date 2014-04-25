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
            ngroup <- length(object@.Data)
            nsnps <- if(ngroup==0) { 0 } else { sapply(object@.Data,length) }
            message("snppicker object, containing ",sum(nsnps)," SNPs grouped into ",ngroup," groups.")
          })

setMethod("show", signature="tags",
          function(object) {
            ntags <- length(unique(object@tags))
            nsnps <- length(object@.Data)
            message("tags object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

setMethod("show", signature="groups",
          function(object) {
            ntags <- length(object@tags)
            nsnps <- length(unlist(object@.Data))
            message("groups object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

######################################################################

##          detailed summary                                        ##

######################################################################

setMethod("summary",signature="snppicker",
          function(object){
            ngroups <- length(object@.Data)
            cmpi <- unlist(lapply(object@groups, function(x) max(x$cmpi)))
            nsnp <- unlist(lapply(object@groups, nrow))
            index <- unlist(lapply(object@groups, function(x) x[1,"var"]))
            maxr2 <- unlist(lapply(object@groups, function(x) x[nrow(x),"r2"]))
            data.frame(SNP.index=index,
                       SNP.count=nsnp,
                       min.R2=1-maxr2,
                       sum.MPI=cmpi)
          })
setMethod("summary",signature="groups",
          function(object) {
            data.frame(tag=object@tags,
                       nsnp=sapply(object@.Data,length))
          })

######################################################################

##                          Plotting                                ##

######################################################################

setMethod("plot", signature(x="snppicker",y="missing"),
          function(x, do.plot=TRUE) {
            index <- unlist(lapply(x@groups, function(z) rownames(z)[1]))
            df <- lapply(x@plotsdata, function(xx) {
              xx$cmpi <- cumsum(xx$mpi)
              xx <- cbind(subset(xx, cmpi<1.25),index=rownames(xx)[1])
            })
            df <- do.call("rbind",df)
            df <- subset(df,!is.na(r2))
            cols <- c("#2171B5", "#08306B")
            cols <- c("#238B45", "#00441B")
            red <-  "#FF7F00"
            red <- "#666666"
            plots <-

              ggplot(df,aes(x=r2,y=cmpi)) + 
                geom_vline(mapping=aes(xintercept=r2), data=subset(df,changepoint==TRUE),col=red) +
                  geom_hline(mapping=aes(yintercept=cmpi), data=subset(df,changepoint==TRUE),col=red,width=0.2) +
                    geom_point(col=cols[2],pch=3) + geom_path(col=cols[2],linetype="longdash") +
                      geom_point(aes(y=mpi),col=cols[1],pch=3) + geom_path(aes(y=mpi),col=cols[1],linetype="longdash") +
                        ylim(0,1.2) + xlab("rsq with index SNP") + ylab("cumulative MPI/MPI") +
                          facet_grid(index ~ .) +  theme(strip.text.y = element_text(angle = 0, size = 6,
                                                                        vjust = 0.5)) 

            if(do.plot)
              print(plots)
            invisible(plots)            
          })

################################################################################

## convert between snppicker, groups and tags

################################################################################

setAs("groups", "tags",      
      def=function(from) {
        if(!length(from))
          return(new("tags")) # return empty object
        new("tags",unlist(from@.Data),tags=rep(from@tags,times=sapply(from@.Data,length)))
      })
setAs("tags", "groups",
      def=function(from) {
        if(!length(from@tags))
          return(new("groups"))
        new("groups",split(from@.Data,from@tags),tags=unique(from@tags))
      })
setAs("snppicker","groups",
      def=function(from) {
        if(!length(from@groups))
          return(new("groups")) # return empty object
        new("groups",
            lapply(from@groups, "[[", "var"),
            tags=unlist(lapply(from@groups,function(x) x[1,"var"])))
      })
setAs("snppicker","tags",
      def=function(from) {
        if(!length(from@groups))
          return(new("tags")) # return empty object
        g <- new("groups",
                 lapply(from@groups, "[[", "var"),
                 tags=unlist(lapply(from@groups,function(x) x[1,"var"])))
        as(g,"tags")
      })
setMethod("[",signature=c(x="groups",i="character",j="missing",drop="missing"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            new("groups",x@.Data[wh],tags=x@tags[wh])
          })
setMethod("[",signature=c(x="groups",i="numeric",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
setMethod("[",signature=c(x="groups",i="logical",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
setMethod("[[",signature=c(x="groups",i="numeric"),
          function(x,i) {
            x@.Data[[i]]
          })
setMethod("[[",signature=c(x="groups",i="logical"),
          function(x,i) {
            x@.Data[[i]]
          })
setMethod("[[",signature=c(x="groups",i="character"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            x@.Data[[wh]]
          })

######################################################################

##                         conversion                               ##

######################################################################

setMethod("convert",signature=c(object="groups"), function(object) {
  object@.Data=object@groups
  return(object)
})
setMethod("convert",signature=c(object="tags"), function(object) {
  object@.Data=object@snps
  return(object)
})

######################################################################

##                           manipulate snppicker                                 ##

######################################################################

setMethod("[",signature=c(x="snppicker",i="ANY",j="missing",drop="missing"),
          function(x,i) {
            new("snppicker",
                groups=x@groups[i],
                plotsdata=x@plotsdata[i])
          })

setMethod("[[",signature=c(x="snppicker",i="ANY"),
          function(x,i) {
            groups=x@groups[[i]]
          })

######################################################################

##           snpin: is a snp in a  group?                  ##

######################################################################

#' @rdname snpin-methods
#' @aliases snpin,character,snppicker-method
setMethod("snpin",signature(x="character",y="snppicker"),definition=function(x,y) {
  snpin(x,as(y,"groups"))
})
#' @rdname snpin-methods
#' @aliases snpin,character,tags-method
setMethod("snpin",signature(x="character",y="tags"),definition=function(x,y) {
  snpin(x,as(y,"tags"))
})
#' @rdname snpin-methods
#' @aliases snpin,character,groups-method
setMethod("snpin",signature(x="character",y="groups"),definition=function(x,y) {
  if(!length(y@.Data))
    return(NULL)
  names(y@.Data) <- paste0("group",seq_along(y@.Data))
  ret <- sapply(y@.Data,function(yg) x %in% yg)
  if(nrow(ret))
    rownames(ret) <- x
  return(ret)
})

######################################################################

##                        concatenate                               ##

######################################################################

#' @rdname union-methods
#' @aliases union,snppicker,snppicker-method
setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {
  union(as(x,"groups"),as(y,"groups")) })
#' @rdname union-methods
#' @aliases union,tags,tags-method
setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
  as(union(as(x,"groups"),as(y,"groups")),"tags") })
#' @rdname union-methods
#' @aliases union,groups,groups-method
setMethod("union",signature(x="groups",y="groups"),definition=function(x,y) {
  if(!length(x))
    return(y)
  if(!length(y))
    return(x)
  ## find intersecting groups, and make unions
  int <- matrix(FALSE,length(x),length(y))
  for(i in seq_along(x)) {
    for(j in seq_along(y)) {      
      int[i,j] <- any(x[[i]] %in% y[[j]])
    }
  }
  wh <- which(int,arr.ind=TRUE)
  if(nrow(wh)) {
    ## merge y into x
    for(i in 1:nrow(wh)) {
      xi <- wh[i,1]
      yi <- wh[i,2]
      x[[xi]] <- unique(c(x[[xi]],y[[yi]]))
    }
    ## check - duplicate y's mean x groups must be merged
    dups <- wh[,2][ duplicated(wh[,2]) ]
    if(length(dups)) {
      ## do the merge
      for(ydup in dups) {
        xi <- wh[ wh[,2]==ydup, 1]
        x[[ xi[1] ]] <- unique(unlist(x[xi]))
      }
      ## drop NULLs now we don't need wh any more
      to.drop <- vector("list",length(dups))
      for(i in seq_along(dups)) {
        yi <- dups[[i]]
        to.drop[[i]] <- (wh[ wh[,2]==ydup, 1])[-1]
      }
      to.drop <- unique(unlist(to.drop))
      x@.Data <- x@.Data[ -to.drop ]
      x@tags <- x@tags[ -to.drop ]
    }
    
    add <- setdiff(seq_along(y),unique(wh[,2]))
    if(!length(add))
      return(x)
    return(new("groups",c(x,y[add]),tags=c(x@tags,y@tags[add])))
  }
  ## no overlap: simply concatenate
  new("groups",c(x,y=y),tags=c(x@tags,y@tags))
})

######################################################################

##                      manipulate snpmods                             ##

######################################################################

## setMethod("snpdrop",signature(x="snpmod",y="character"),
##           function(x,y) {
##             wh <- which(sapply(x@model.snps,function(x) any(x %in% y)))
##             message(length(wh), " / ", length(x@model.snps), " models (",
##                     format.pval(100*length(wh)/length(x@model.snps)), "%) will be dropped.")
##             if(!length(wh))
##               return(d)
##             x@models <- x@models[-wh,]
##             x@model.snps <- x@model.snps[-wh]  
##             return(marg.snps(x))  
##           })

## setMethod("+",signature(e1="snpmod",e2="snpmod"), function(e1,e2) snpmod.add(e1,e2))
