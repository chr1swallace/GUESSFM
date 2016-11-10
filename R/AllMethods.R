######################################################################

##                          display                                 ##

######################################################################
##' Show
##'
##' Methods to briefly summarise objects to screen. \code{show} is
##' called implicitly when you type an object's name at the command
##' line.  Use \code{summary}, where available, to get more details.
##' @param object the thing to show
##' @export
##' @rdname show-methods 
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

##' @rdname show-methods
setMethod("show", signature="snppicker",
          function(object) {
            ngroup <- length(object@groups)
            nsnps <- if(ngroup==0) { 0 } else { sapply(object@groups,nrow) }
            message("snppicker object, containing ",sum(nsnps)," SNPs grouped into ",ngroup," groups.")
          })
##' @rdname show-methods
setMethod("show", signature="ppnsnp",
          function(object) {
            L <- object@.Data
            names(L) <- object@traits
            show(L)
          })

##' @rdname show-methods
setMethod("show", signature="tags",
          function(object) {
            ntags <- length(unique(object@tags))
            nsnps <- length(object@.Data)
            message("tags object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

##' @rdname show-methods
setMethod("show", signature="groups",
          function(object) {
            ntags <- length(object@tags)
            nsnps <- length(unlist(object@.Data))
            message("groups object, containing ",nsnps," SNPs in ",ntags," groups.")
          })

######################################################################

##          detailed summary                                        ##

######################################################################


##' Summaries
##' 
##' Print summary of an object
##'
##' @param object the thing to summarise
##' @export
##' @rdname summary
setMethod("summary",signature="snppicker",
          function(object){
            ngroups <- length(object@groups)
            cmpi <- unlist(lapply(object@groups, function(x) max(x$cmpi)))
            nsnp <- unlist(lapply(object@groups, nrow))
            index <- unlist(lapply(object@groups, function(x) x[1,"var"]))
            maxr2 <- unlist(lapply(object@groups, function(x) x[nrow(x),"r2"]))
            data.frame(SNP.index=index,
                       SNP.count=nsnp,
                       min.R2=1-maxr2,
                       gMMPI=cmpi)
          })
##' @rdname summary
setMethod("summary",signature="groups",
          function(object) {
            data.frame(tag=object@tags,
                       nsnp=sapply(object@.Data,length))
          })

######################################################################

##                          Plotting                                ##

######################################################################

##' Plots
##'
##' Descriptive plots showing how the sets in snppicker objects were
##' generated.  Uses \code{ggplot2}, so you can customize using
##' \code{theme} etc.
##'
##' @rdname plot-methods
##' @param x ppsnp or snppicker object
##' @param do.plot if TRUE (default) print the plot on the current
##' device, otherwise invisibly return the plot object.
##' @export
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
              ggplot(df,aes(x=1-r2,y=cmpi)) + 
                geom_vline(mapping=aes(xintercept=1-r2), data=subset(df,changepoint==TRUE),col=red) +
                  geom_hline(mapping=aes(yintercept=cmpi), data=subset(df,changepoint==TRUE),col=red) +
                    geom_point(col=cols[2],pch=3) + geom_path(col=cols[2],linetype="longdash") +
                      geom_point(aes(y=mpi),col=cols[1],pch=3) + geom_path(aes(y=mpi),col=cols[1],linetype="longdash") +
                        ylim(0,1.2) + scale_x_reverse("rsq with index SNP") + ylab("cumulative MPI/MPI") +
                          facet_grid(index ~ .) +  theme(strip.text.y = element_text(angle = 0, size = 6,
                                                                        vjust = 0.5)) 

            if(do.plot)
              print(plots)
            invisible(plots)            
          })

##' @rdname plot-methods
##' @export
setMethod("plot", signature(x="ppnsnp",y="missing"),
          function(x) print(x@plot))
          

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

################################################################################

## subsets

################################################################################

##' @rdname groups-subset
setMethod("taggedby",signature=c(object="tags",i="character"),
          function(object,i) {
            wh <- which(object@tags %in% i)
            if(!length(wh))
              stop("tags not found")
            ret <- data.frame(tag=object@tags[wh],snps=object@.Data[wh])
            return(ret[order(ret$tag),])  
          })
##' @rdname groups-subset
setMethod("tagsof",signature=c(object="tags",i="character"),
          function(object,i) {
            wh <- which(object@.Data %in% i)
             if(!length(wh))
              stop("SNPs not found")
           ret <- data.frame(tag=object@tags[wh],snps=object@.Data[wh])
            return(ret[order(ret$tag),])  
          })
##' @rdname groups-subset
setMethod("taggedby",signature=c(object="groups",i="character"),
          function(object,i) {
            taggedby(as(object,"tags"),i)
          })
##' @rdname groups-subset
setMethod("tagsof",signature=c(object="groups",i="character"),
          function(object,i) {
            tagsof(as(object,"tags"),i)
          })

##' Subset groups or tags objects
##'
##' '[' will extract another object of the same class.  '[[' will extract a single element.
##' @param x groups or tags object
##' @param i numeric, logical or character vector to index SNPs or tags
##' @return subsetted groups or tags object
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="character",j="missing",drop="missing"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            new("groups",x@.Data[wh],tags=x@tags[wh])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="numeric",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="groups",i="logical",j="missing",drop="missing"),
          function(x,i) {
            new("groups",x@.Data[i],tags=x@tags[i])
          })
##' @rdname groups-subset
setMethod("[",signature=c(x="tags",i="character",j="missing",drop="missing"),
          function(x,i) {
            wh <- sapply(i,function(ii) which(snps(x)==ii))
            tags(tags)[wh]
          })
##' @rdname groups-subset
##' @export
setMethod("[[",signature=c(x="groups",i="numeric"),
          function(x,i) {
            x@.Data[[i]]
          })
##' @rdname groups-subset
setMethod("[[",signature=c(x="groups",i="logical"),
          function(x,i) {
            x@.Data[[i]]
          })
##' @rdname groups-subset
setMethod("[[",signature=c(x="groups",i="character"),
          function(x,i) {
            wh <- which(x@tags %in% i)
            x@.Data[[wh]]
          })

######################################################################

##                         conversion                               ##

######################################################################

#' @rdname conversion
setMethod("convert",signature=c(object="groups"), function(object) {
  object@.Data=object@groups
  return(object)
})
#' @rdname conversion
setMethod("convert",signature=c(object="tags"), function(object) {
  object@.Data=object@snps
  return(object)
})

######################################################################

##                           manipulate snppicker                                 ##

######################################################################

##' Subset snppicker object
##'
##' This works a bit like subsetting a list.  If you use the '['
##' construct you will get another snppicker object with a subset of
##' the elements.  The ability to plot that subset should remain.  If
##' you use the '[[' construct it will extract just that element, and
##' return an object of class groups.  No plotting ability is
##' retained, but a single groups object is easy to inspect and
##' manipulate.
##' @param x snppicker object
##' @param i numeric, logical or character vector to index SNPs
##' @return subsetted snppicker object
##' @rdname snppicker-subset
setMethod("[",signature=c(x="snppicker",i="ANY",j="missing",drop="missing"),
          function(x,i) {
            new("snppicker",
                groups=x@groups[i],
                plotsdata=x@plotsdata[i])
          })

##' @rdname snppicker-subset
setMethod("[[",signature=c(x="snppicker",i="ANY"),
          function(x,i) {
            groups=x@groups[[i]]
          })

######################################################################

##           snpin: is a snp in a  group?                  ##

######################################################################

#' @rdname snpin
setMethod("snpin",signature(x="character",y="snppicker"),definition=function(x,y) {
  snpin(x,as(y,"groups"))
})

#' @rdname snpin
setMethod("snpin",signature(x="character",y="tags"),definition=function(x,y) {
  snpin(x,as(y,"groups"))
})

#' @rdname snpin
setMethod("snpin",signature(x="character",y="groups"),definition=function(x,y) {
  if(!length(y@.Data) || !length(x))
    return(NULL)
  names(y@.Data) <- paste0("group",seq_along(y@.Data))
  ret <- sapply(y@.Data,function(yg) x %in% yg)
  if(is.null(dim(ret)))
    ret <- matrix(ret,nrow=1)
  rownames(ret) <- x
  return(ret)
})

######################################################################

##                        concatenate                               ##

######################################################################
## two groups in x (8, 11) match to one group in y (1) with dups
## x.7 also matches to y.1

#' @rdname union
setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {
  xgr <- as(x,"groups")
  ygr <- as(y,"groups")
  int <- .group.intersection(xgr,ygr)
  wh <- which(int,arr.ind=TRUE)
  M <- new("snppicker",
               plotsdata=c(x@plotsdata,y@plotsdata),
               groups=c(x@groups,y@groups))
  if(!length(wh)) { # no overlap, concatenate
    return(M)
  }
  Mkeep <- rep(TRUE,length(M@groups))
  for(i in 1:nrow(wh)) {
    ix <- wh[i,1]
    iy <- wh[i,2]
    iy.M <- iy + length(x@groups)
    gx <- M@groups[[ ix ]]
    gy <- M@groups[[ iy ]]
    gint <- intersect(rownames(gx),rownames(gy))
    pintx <- sum(gx[gint,"Marg_Prob_Incl"])
    pinty <- sum(gy[gint,"Marg_Prob_Incl"])
    pnintx <- sum(gx[setdiff(rownames(gx),gint),"Marg_Prob_Incl"])
    pninty <- sum(gy[setdiff(rownames(gy),gint),"Marg_Prob_Incl"])
    if(pnintx > pintx || pninty > pinty) { # unmerge
      if(pintx/(pintx + pnintx) > pinty/(pinty + pninty)) { # keep int in x
        M@groups[[iy.M]] <- M@groups[[iy.M]][ !(rownames(M@groups[[iy.M]]) %in% gint), ]
      } else { # keep int in y
        M@groups[[ix]] <- M@groups[[ix]][ !(rownames(M@groups[[ix]]) %in% gint), ]
      }
    } else { # merge
      tmp <- rbind(M@groups[[ix]],
                   M@groups[[iy.M]][ !(rownames(M@groups[[iy.M]]) %in% gint), ])
      M@groups[[ix]] <- tmp[ order(tmp$Marg_Prob_Incl, decreasing=TRUE), ]
      Mkeep[[ iy.M ]] <- FALSE
    }
    cat(length(M@groups[[ix]]$var), length(unique(M@groups[[ix]]$var)),
        length(M@groups[[iy.M]]$var), length(unique(M@groups[[iy.M]]$var)))
        
  }
  new("snppicker",
      plotsdata=M@plotsdata[which(Mkeep)],
      groups=M@groups[which(Mkeep)])
})

#' @rdname union
setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
  ugr <- union(as(x,"groups"),as(y,"groups"))
  as(ugr,"tags") })

.group.intersection <- function(x,y) {
  int <- matrix(FALSE,length(x),length(y))
  for(i in seq_along(x)) {
    for(j in seq_along(y)) {      
      int[i,j] <- any(x[[i]] %in% y[[j]])
    }
  }
  return(int)
}

#' @rdname union
setMethod("union",signature(x="groups",y="groups"),definition=function(x,y) {
  if(!length(x))
    return(y)
  if(!length(y))
    return(x)
  ## find intersecting groups, and make unions
  int <- .group.intersection(x,y)
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

## split

######################################################################

## # @rdname union-methods
## # @aliases union,snppicker,snppicker-method
## setMethod("union",signature(x="snppicker",y="snppicker"),definition=function(x,y) {
##   union(as(x,"groups"),as(y,"groups")) })
## # @rdname union-methods
## # @aliases union,tags,tags-method
## setMethod("union",signature(x="tags",y="tags"),definition=function(x,y) {
##   as(union(as(x,"groups"),as(y,"groups")),"tags") })
## # @rdname split-methods
## # @aliases split,groups,groups-method
## setMethod("split",signature(x="groups",y="groups"),definition=function(x,y) {


######################################################################

##                      accessors                                   ##

######################################################################

##' @rdname accessors
setMethod("snps",signature(object="groups"), function(object) { object@.Data })
##' @rdname accessors
setMethod("tags",signature(object="groups"), function(object) { object@tags })
##' @rdname accessors
setMethod("snps",signature(object="tags"), function(object) { object@.Data })
##' @rdname accessors
setMethod("tags",signature(object="tags"), function(object) { object@tags })

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

## 
################################################################################

## qc
##' @rdname qc
##' @export
setMethod("qc",signature(object="ppnsnp",data="missing"),
          function(object) {
            qc.skew <- skewness(object)
            qc.nmodes <- nmodes(object)
            ret <- data.frame(trait=object@traits, skewness=unlist(qc.skew),
                              nmodes=unlist(qc.nmodes))
            ret$flag <- ret$skewness<0 | ret$nmodes>1
            return(ret)
          })
##' @rdname qc
##' @export
setMethod("qc",signature(object="snpmod",data="SnpMatrix"),
          function(object,data) {
            data <- data[,rownames(object@snps)]
            LD <- ld(data, stats="R.squared", symmetric=TRUE, depth=ncol(data)-1)
            diag(LD) <- 0
            M <- object@models
            ss <- strsplit(M$str,"%")
            maxld <- numeric(length(ss))
            n <- sapply(ss,length)
            wh <- which(n > 1)
            maxld[ wh ] <- unlist(lapply(ss[wh],function(s) max(LD[s,s],na.rm=TRUE)))
            ret <- data.frame(model=M$str, nsnps=n, PP=M$PP, maxr2=maxld, stringsAsFactors=FALSE)
            ret.s <- snps.from.correlated.models(ret)
            invisible(list(summary=ret.s$summary,SNPs=ret.s$SNPs,models=ret))
          })
##' @rdname qc
##' @export
##' @author Chris Wallace
setMethod("qc",signature(object="list",data="ANY"),
          function(object) {
            lapply(object,qc,data=data)
          })
