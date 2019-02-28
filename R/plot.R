utils:::globalVariables(c("A", "B", "cPP", "height", "jeffreys",
"logPP", "n", "phenotype", "position.plot", "pp", "ppsum", "ppthr",
"size", "trait", "value", "variable", "x", "X1", "X2", "xend", "xmax",
"xmin", "x.min", "x.scale", "y", "ymin","ymax","yend"))

##' Plot two summaries of the models, a diffusion plot and a summary of the prior and posterior number of SNPs in the models
##'
##' @title Summary plots
##' @param results object of class snpmod, or a named list of such objects
##' @return a list of two ggplot2 objects, which may be print()ed to the current graphics device
##' @export
plotsummary <- function(results) {
  if(!is.list(results))
    results <- list(trait=results)
  plots <- list(diffusion=plot_diffuse(results),
                pp.nsnp=pp.nsnp(results,plot=TRUE))                
}
##' create a track from a guess.summ object showing snp ids
##'
##' @title ggsnp
##' @param summx summary of a snpmod + groups object
##' @param snp.name column name for snp labels, default \code{snp}
##' @return a ggplot object
##' @author Chris Wallace
##' @export
ggsnp <- function(summx,snp.name="snp") {
  summx <- unique(summx[,c(snp.name,"tag","x.scale")])
  ggplot(summx,aes_string(x="x.scale",y=1,label=snp.name,xintercept="x.scale")) +
    geom_vline(aes(xintercept=x.scale,col=tag),size=0.2,alpha=0.5) +
      geom_text(angle=90,vjust=0.5,hjust=0,size=4,fontface="bold") +
        theme(legend.position="none",axis.text=element_blank(),
              axis.ticks=element_blank(),axis.title=element_blank()) + ylim(1,1.4)
}

plot.bf.pp <- function(results, nm="unknown",maxrank=1000) {
  if(!is.list(results))
    results <- structure(list(results),names=nm)
  plots <- structure(vector("list",length(results)), names=names(results))
  for(nm in names(results)) {
    df <- results[[nm]]$models 
    df$phenotype <- nm
    df <- subset(df, rank<=maxrank)
    df$jeffreys[ is.infinite(df$jeffreys) ] <- 1000
    plots[[nm]] <- ggplot(df, aes(x=jeffreys,y=logPP,col=as.factor(size))) + geom_point()
  }
  do.call("grid.arrange", plots)
}

plot.fdr <- function(summ,causal) {
  ppgrid <- seq(0,1,0.01)
  called <- lapply(ppgrid,function(i) {
    df <- subset(summ,pp>i, select=c("trait","snp"))
    df$true <- df$snp %in% causal
    df$ppthr <- rep(i,nrow(df))
    return(df)
  })
  called <- do.call(rbind,called)
  true <- melt(with(called, tapply(true,list(trait,ppthr),sum)), varnames=c("trait","ppthr"))
  false <- melt(with(called, tapply(!true,list(trait,ppthr),sum)), varnames=c("trait","ppthr"))
  colnames(true) <- sub("value","true",colnames(true))
  colnames(false) <- sub("value","false",colnames(false))
  result <- merge(true,false)
#  result[ is.na(result) ] <- 0

  result$prop.called.true <- with(result, true/(true + false))
  result$prop.true.called <- result$true/length(causal)

  mresult <- melt(result[,c("trait","ppthr","prop.called.true","prop.true.called")],id.vars=c("trait","ppthr"))
  p <- ggplot(mresult, aes(x=ppthr,y=value,col=variable)) + geom_point() + geom_path() + facet_wrap(~trait) 
  
  return(list(data=mresult, plot=p))
}
##' A summary plot to indicate how diffuse the posterior is.
##'
##' This is measured by how rapidly the cumulative posterior
##' probability increases with number of models.  For one dataset,
##' this may be informative in an absolute sense, but comparing
##' different datasets, by supplying a list of results, can be useful
##' too.
##'
##' If a list of snpmod objects are given, multiple diffusion plots
##' are overlaid, with colours indicating labels.
##' @title Diffusion plot
##' @param results an object of class snpmod, or a named list of such objects
##' @param maxrank truncate the x axis at this value
##' @param thin thin data to display only \code{thin} points per line
##' @return the ggplot2 object, which by default is print()ed to the current graphics device
##' @author Chris Wallace
##' @export
plot_diffuse <- function(results, maxrank=1000, thin=500) {
  if(!is.list(results))
    results <- list(trait=results)
  df <- lapply(results, function(x) x@models)
  for(i in 1:length(results)) {
    if(!("rank" %in% colnames(df[[i]]))) {
      df[[i]]$rank <- 1:nrow(df[[i]])
    }
    df[[i]]$phenotype <- names(results)[i]
    df[[i]]$cPP <- cumsum(df[[i]]$PP)
    df[[i]]$cPP <- df[[i]]$cPP/max(df[[i]]$cPP)
  }
  df.max <- lapply(df,function(x) x[nrow(x),])
  df.max <- do.call("rbind",df.max)
  df <- do.call("rbind",df)
  df <- subset(df, rank<=maxrank)
  np <- levels(df$phenotype <- factor(df$phenotype))
  use <- unlist(tapply(1:nrow(df), df$phenotype, function(x) x[ seq(1,length(x), length.out=thin) ],simplify=FALSE)) ## WORKS
  
  ggplot(df[use,], aes(x=rank,y=cPP, col=phenotype)) + geom_path() + ylim(0,1) +
    geom_hline(aes(yintercept=cPP,col=phenotype),data=df.max,linetype="dotted")
}
##' Rotated LD plot
##'
##' Generates a plot of r2, rotating through 45 degrees so it can be
##' aligned with other results, similar to Haploview.
##'
##' The positions in the summ object are used to align the snps to
##' the correct coordinates.
##' @title ggld
##' @param data object of class SnpMatrix from which LD will be calculated 
##' @inheritParams signal.plot
##' @return ggplot output
##' @author Chris Wallace
##' @export
##' @family plotting GUESSFM results
ggld <- function(data, summ) {
use <- !duplicated(summ$snpnum)
snps.num <- structure(summ$snpnum[use],names=rownames(summ)[use])
all.snps <- names(snps.num)
  LD <- melt(as(ld(data[,all.snps], stats="R.squared", depth=length(all.snps)-1, symmetric=TRUE),"matrix"))
  LD$X1 <- snps.num[as.character(LD$X1)] - 1
  LD$X2 <- snps.num[as.character(LD$X2)] - 1
n <- length(all.snps)
if(n<2)
    return(NULL)
  offset <- n/sqrt(2)
  LD$A <- with(LD, -(n-X1-X2)/sqrt(2))
  LD$B <- with(LD, -(n-X1-X2)/sqrt(2) - X1*sqrt(2))
  LD <- LD[ LD$X1>LD$X2, ]

  ## align A to 1:n
  LD$A <- xscale(LD$A, torange=c(1.5, n-0.5))
  ## scale by pos
LD$A <- xscale(LD$A,
               xrange=c(min(summ$snpnum),max(summ$snpnum)),
               torange=c(min(summ$x.scale),max(summ$x.scale)))
  
  ## length of ticks
  tlength <- (max(LD$B) - min(LD$B))/50
  
  hm.colours <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                  "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  hmap <- ggplot(LD, aes(x=A,y=B)) + geom_point(aes(fill=value), pch=23, colour="white", size=2) +
    scale_fill_gradientn(colours = rev(grey.colors(20))) +
      labs(x = NULL, y = NULL) +
#scale_x_continuous(expand=c(0,0),breaks=NULL) +
        scale_y_continuous(expand=c(0,0),breaks = NULL) +
          theme(plot.margin = unit(rep(0, 4), "lines"),
                legend.position="none",
                panel.grid=element_blank()) +
                  geom_segment(data     = data.frame(x=sort(unique(summ$x.scale)),
                                 y=max(LD$B) + tlength, yend=max(LD$B) + 2*tlength), 
                               aes(x    = x, 
                                   y    = y,
                                   xend = x, 
                                   yend = yend  ))      #top-ticks
}  
##' Generate main results plot: the sets of SNPs, their groups and marginal posterior probabilities of inclusion
##'
##' @title Signal Plot
##' @param summ data.frame generated by \code{\link{guess.summ}()}
##' @param w half the width of a block, default 0.2
##' @param highlight list of SNP ids which, if found in summ$snp, are highlighted
##' @return ggplot2 object, which by default is print()ed to the current graphics device
##' @family plotting GUESSFM results
##' @export
##' @family plotting GUESSFM results
signal.plot <- function(summ,w=0.2,highlight=NULL) {
  if(!("x.scale" %in% colnames(summ)) || any(is.na(summ$x.scale))) {
    stop("Missing x co-ordinates.  Do you need to run scalepos or fix some missing values?")
  }
  summ$xmin <- summ$xmin.scale - w
  summ$xmax <- summ$xmax.scale + w
  summ$x <- summ$x.scale
  summ$tag <- as.factor(summ$tag)
  p <- ggplot(summ, aes(x=x, col=tag)) +
    geom_vline(aes(xintercept=x,col=tag),alpha=0.5,size=0.2) +
      geom_hline(yintercept=1,col="grey",size=0.2,linetype="dashed") +
        geom_point(aes(y=pp),size=1) + 
      ylab("PP") +
        theme(legend.position="none",
              strip.text.y = element_text(size=9,angle=0),
              panel.grid=element_blank()) +
#        geom_vline(data=subset(summ, !is.na(x.min)), mapping=aes(xintercept=snpnum-0.5),col="grey",linetype="dashed") +
          geom_segment(data=subset(summ, !is.na(xmin)),
                       mapping=aes(x=xmin,xend=xmax,y=ppsum,yend=ppsum, col=tag),size=0.5) +
                         geom_rect(data=subset(summ, !is.na(x.min)), mapping=aes(xmin=xmin,xmax=xmax,ymin=0,ymax=ppsum,fill=tag), alpha=0.1,size=0.2)
  if(length(highlight)) {
    summ$highlight <- ifelse(summ$snp %in% highlight, summ$x, NA)    
    if(any(!is.na(summ$highlight)))
      p <- p + geom_vline(data=subset(summ,!is.na(highlight)),
                          mapping=aes(xintercept=highlight),col="black",linetype="dotted",size=0.5)
  }
  if(length(unique(summ$trait))>1)
    p <- p + facet_grid(trait ~ .)
  return(p)
}

##' Plot lines connecting each SNP's physical position to their representative position in the summ object
##'
##' @title ggchr
##' @inheritParams signal.plot
##' @return ggplot2 object
##' @author chris
##' @export
##' @family plotting GUESSFM results
ggchr <- function(summ) {
  ggplot(summ, aes(x=x.scale,xend=position.plot,y=0,yend=1, colour=tag)) +
    ggplot2::geom_segment() +
  theme(legend.position="none",
          strip.text.y = element_text(size=9,angle=0),
          axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
              panel.grid=element_blank())
}
##' Plot a bed file, overlaying lines showing the positions of SNPs in a summ object
##'
##' @title Plot a bed file
##' @param bed data.frame representing data read in from a bed file
##' @inheritParams signal.plot
##' @return ggplot object
##' @author chris
##' @export
##' @family plotting GUESSFM results
ggbed <- function(bed,summ) {
  ggplot(bed, aes(xmin=start,xmax=end,ymin=0,ymax=1,alpha=height)) + ggplot2::geom_rect() + facet_grid(name ~ ., margins=TRUE) +
    geom_vline(data=summ, mapping=aes(xintercept=position.plot,col=tag)) +
    theme(legend.position="none",
          strip.text.y = element_text(size=9,angle=0),
          axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
              panel.grid=element_blank())
}



##' Scale SNP positions in a summ object in order to be plotted by signal.plot and friends
##' 
##' @title Scale SNP positions
##' @inheritParams signal.plot
##' @param position character string: name of the column used for SNP position, default is "position"
##' @return summ data.frame with additional columns, x.scale, xmin.scale, xmax.scale
##' @export
##' @family plotting GUESSFM results
scalepos <- function(summ,position="position",prange=NULL) {
    if(!(position %in% colnames(summ)))
        stop("Position column not found: ",position)       
    summ$position.plot <- summ[,position]
    pr <- if(is.null(prange)) {
              c(min(summ$position.plot),max(summ$position.plot))
          } else {
              prange
          }
    xr <- c(min(summ$snpnum),max(summ$snpnum))
    summ$x.scale <- xscale(summ$snpnum,torange=pr,xrange=xr)
    summ$xmin.scale <- xscale(summ$x.min,torange=pr,xrange=xr)
    summ$xmax.scale <- xscale(summ$x.max,torange=pr,xrange=xr)
    return(summ)
}
##' Summarize the posterior model support by the number of SNPs contained in a model
##'
##' @param results object of class snpmod
##' @param plot if TRUE, print a plot to current plotting device
##' @param expected optional, if specified, this is the prior expectation of the number of causal SNPs underlying the trait in the studied region.  Specifying will cause a line for the prior to be added to the plot
##' @param overdispersion optional, default=1.  Ignored unless expected is given.  The overdispersion of the beta binomial distribution used for the prior relative to a binomial.
##' @return a named list containing
##' \itemize{
##'  \item{"pp"}{a named vector summarizing the posterior for each count of SNPs}
##'  \item{"plot"}{a ggplot object}
##' } 
##' @family plotting GUESSFM results
##' @export
##' @seealso \link{snpprior}
pp.nsnp <- function(results,plot=FALSE,expected=NULL,overdispersion=1) {
  if(!is.list(results)) ## single trait    
    results <- list(trait=results)
  if(is.null(names(results))) ## unnamed list
    names(results) <- paste0("trait",1:length(results))
  df <- pp.nsnps <- vector("list",length(results))
  names(df) <- names(pp.nsnps) <- names(results)
  for(i in seq_along(df)) {
    d <- results[[i]]
    pp.nsnps[[i]] <- tapply(d@models$PP,d@models$size,sum)
    df[[i]] <- data.frame(n=as.numeric(names(pp.nsnps[[i]])),pp=pp.nsnps[[i]],trait=names(results)[i])
  }
  df <- do.call("rbind",df)
  if(!is.null(expected)) {
    nsnp <- nrow(results[[1]]@snps)
    rho <- (overdispersion - 1)/(nsnp-1)
    p <- expected/nsnp
    nind <- 0:max(df$n)
    pp <- if(rho==0) {
              dbinom(nind, size=nsnp, prob=p)
          } else {
              dbetabinom(nind, size=nsnp, prob=p, rho=rho)
          }
    df <- rbind(df,data.frame(n=nind,
                              pp=pp,
                              trait="Prior"))
  }
    p <- ggplot(df,aes(x=n,y=pp,col=trait)) + geom_point() + geom_path() + 
    labs(x="Number of SNPs in model",
           y=if(is.null(expected)) {"Posterior probability"} else {"Prior/posterior probability"}) +
             scale_x_continuous(breaks=seq(0,max(df$n),by=2))
  if(plot)
    print(p)
##  print(pp.nsnps)
  new("ppnsnp",.Data=pp.nsnps,plot=p,traits=names(pp.nsnps))
#  invisible(list(pp=pp.nsnps,plot=p))
}
##' Add lines to a ggplot, indicating position of snps
##'
##' @title addlines
##' @inheritParams signal.plot
##' @param alpha transparency, default=0.5
##' @param size line width, default=0.2
##' @param xvar name of variable in summ which contains the x coordinates at which vertical lines should be plotted
##' @return object that can be added to a ggplot
##' @author Chris Wallace
##' @export
##' @family plotting GUESSFM results
addlines <- function(summ,alpha=0.5,size=0.2,xvar="position.plot") {
  geom_vline(data=summ, mapping=aes_string(xintercept=xvar,col="tag"),alpha=alpha,size=size) 
}

##' Calculate D' and R^2 for a SnpMatrix object and display as a heatmap
##'
##' Side effect: displays plot on current graphics device
##' @title show.ld
##' @return invisibly returns matrix with D' in upper.tri() entries and R^2 in lower.tri() entries
##' @author Chris Wallace
##' @export
##' @param X a SnpMatrix object
##' @param snps optional character vector of column names of X (SNPs)
##' for which LD should be calculated
##' @param samples optional character vector of row names of X
##' (Samples) for which LD should be calculated
##' @param lines.limit draw faint lines to make separation between SNPs clearer if the number of SNPs is < lines.limit
show.ld <- function(X, snps=colnames(X), samples=rownames(X),
                    lines.limit=20) {
  my.colours <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
                  "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  if(!is(X,"SnpMatrix"))
    X <- as(X,"SnpMatrix")
  breaks <- NULL
  groups <- NULL
  if(is.list(snps)) {
    breaks <- sapply(snps,length)
    groups <- rep(1:length(snps), times=breaks)
    snps <- unlist(snps)
  }
  LD <- ld(X[samples,snps],
           depth=length(snps)-1,
           stat=c("D.prime","R.squared"),
           symmetric=TRUE)
  ld <- as.matrix(LD$R.squared)
  wh <- which(upper.tri(ld))
  ld[wh] <- as.matrix(LD$D.prime)[ wh ]
 
 diag(ld) <- 1
 if(!is.null(groups))
   colnames(ld) <- rownames(ld) <- paste(groups,colnames(ld),sep="/")
 
  df <- melt(as.matrix(ld))
  df$X1 <- factor(df$X1, levels=colnames(ld))
  df$X2 <- factor(df$X2, levels=colnames(ld))
 n <- max(as.numeric(df$X1))
 p <- ggplot(df) + geom_tile(mapping=aes(x=X1,y=X2,fill=value), linetype=1) +
    scale_fill_gradientn(colours = my.colours) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text = element_text(colour="black"),
            axis.title = element_blank())
 if(n<lines.limit) {
   p <- p + geom_hline(yintercept=(0:n)+0.5, col="grey60") +
     geom_vline(xintercept=(0:n)+0.5, col="grey60")
 }
 
if(!is.null(breaks)) {
  x <- c(0,cumsum(breaks),max(as.numeric(df$X1))) + 0.5
  n <- length(x)
  xf <- rbind(data.frame(x=x[-n], y=x[-n], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-n], y=x[-n], xend=x[-n], yend=x[-1]),
              data.frame(x=x[-1], y=x[-1], xend=x[-1], yend=x[-n]),
              data.frame(x=x[-1], y=x[-1], xend=x[-n], yend=x[-1]))
              
  show(p + geom_segment(aes(x=x,xend=xend,y=y,yend=yend), data=xf, size=1))
 } else {
   show(p)
 }
  ## pheatmap(ld,
  ## cluster_rows=FALSE, cluster_cols=FALSE)
  invisible(ld)
}
## myr2 <- function(X) {
##   r2 <- ld(X,
##            depth=ncol(X)-1,
##            symmetric=TRUE,
##            stat="R.squared")
##   if(any(is.na(r2))) {
##     r2.na <- as(is.na(r2),"matrix")
##     use <- rowSums(r2.na)>0
##     ## work around for r2=NA bug.
##     r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
##     r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
##   }
##   diag(r2) <- 1
##   return(r2)
## }

mod2group <- function(str,groups) {  
  x.snps <- strsplit(str,"%")
  x.groups <- lapply(x.snps,snpin,groups)
  G <- sapply(x.groups,function(g) {
    if(is.null(g))
      return("N")
    match <- apply(g,1,any)
    ret <- numeric(nrow(g))
    if(any(match))
      ret[match] <- apply(g[match,,drop=FALSE],1,which)
    return(paste(sort(ret),collapse="-"))
  })
  G[str=="1"] <- "-" # distinguish null model from ungrouped snps
  return(G)
}
numorneg <- function(x) {
  suppressWarnings(n <- as(x,"numeric"))
  if(any(is.na(n)))
    n[is.na(n)] <- -1
  return(n)
}
gsumm <- function(x,groups) {
  n <- max(numorneg(unlist(strsplit(x$group,"-"))))
  counts <- tapply(x$PP,x$group,sum)
  df <- data.frame(pattern=names(counts), PP=counts)
  df <- df[order(df$PP,decreasing=TRUE), ]
  df$ymax <- cumsum(df$PP)
  df$ymin <- c(0,df$ymax[-nrow(df)])
  cn <- lapply(strsplit(rownames(df),"-"),numorneg)  
  df2 <- df[rep(1:nrow(df), times=sapply(cn,length)),]
  df2$xmin <- unlist(cn)
  df2$xmax <- df2$xmin+1
  return(df2)
}

##' @importFrom data.table as.data.table melt setnames
##' @importFrom cowplot plot_grid
NULL

##' Plot pattern of SNP group inclusion
##'
##' @title pattern.plot
##' @param SM snpmod or list of snpmods
##' @param groups groups object
##' @return a ggplot object, by default printed to current graphics device
##' @author Chris Wallace
##' @export
pattern.plot <- function(SM,groups,r2=NULL) {
  if(!is.list(SM))
    SM <- list(trait=SM)
  BM <- lapply(SM,function(x) best.models(x,cpp.thr=0.99)[,c("str","PP")])
  for(i in seq_along(BM)) 
    BM[[i]]$group <- mod2group(BM[[i]]$str,groups)
  G <- lapply(BM,gsumm)
  if(is.null(names(G)))
      names(G) <- paste0("trait",seq_along(G))
  for(i in names(G)) 
    G[[i]]$trait <- i
  G <- do.call("rbind",G)
  G$xmin <- G$xmin - 1
  G$xmax <- G$xmax - 1
  
  p <- ggplot(G,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)) +
    geom_rect() +
      geom_vline(aes(xintercept=xmin),col="gray60", size=0.2) +
        geom_hline(aes(yintercept=ymax),col="gray60",size=0.2) +
          scale_x_continuous(breaks=sort(unique(G$xmin))+0.5,labels=sort(unique(G$xmin))+1,limits=c(min(G$xmin),max(G$xmin)+1),expand=c(0,0)) +
            scale_y_continuous(breaks=c(0,1), expand=c(0,0), limits=c(0,1)) +
              theme_bw() +
                theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.position="none", panel.spacing = unit(1, "lines")) +
                  facet_grid(trait~.) + xlab("SNP group index") + ylab("Cumulative Model Posterior Probility") + 
      scale_fill_manual(values=c("FALSE"="grey","TRUE"="grey20"))
  if(is.null(r2))
      return(p)

  maxr2 <- calc.maxmin(r2,groups,fun=max)
  
  LD <- as.data.table(melt(maxr2))
  n <- ncol(maxr2)
  offset <- n/sqrt(2)
  setnames(LD,c("X1","X2","R2"))
  LD$A <- with(LD, -(n-X1-X2)/sqrt(2))
  LD$B <- with(LD, -(n-X1-X2)/sqrt(2) - X1*sqrt(2))
  LD <- LD[ LD$X1>LD$X2, ]
  if(!nrow(LD))
      return(p)
  ## align A to 1:n
  LD$A <- xscale(LD$A, torange=c(0.5,n-0.5))
  tlength <- (max(LD$B) - min(LD$B))/50
  hmap <- ggplot(LD, aes(x=A,y=B)) +
    geom_point(aes(fill=R2,col=R2), pch=23,  size=4) +
    scale_colour_gradient(low="grey90",high="lightblue",breaks=seq(0,1,by=0.1)) +
    scale_fill_gradient(## colours = rev(grey.colors(10)),
        low="grey90",high="lightblue",breaks=seq(0,1,by=0.1)) +
    labs(x = NULL, y = NULL) +
                                        #scale_x_continuous(expand=c(0,0),breaks=NULL) +
    geom_text(aes(label=sub("0.",".",signif(R2,1))),col="black",data=LD[R2>0.2,]) +
    scale_y_continuous(expand=c(0,0),breaks = NULL) +
    theme(plot.margin = unit(rep(0, 4), "lines"),
                                        #            legend.position="none",
          panel.grid=element_blank()) +
    geom_segment(data     = data.frame(x=1:n-0.5,
                                       y=max(LD$B) + tlength, yend=max(LD$B) + 2*tlength), 
                 aes(x    = x, 
                     y    = y,
                     xend = x, 
                     yend = yend  )) +
    scale_x_continuous(breaks=NULL,limits=c(-1,n),expand=c(0,0))
plot_grid(p,hmap,ncol=1,align="v",axis="blr",rel_heights=c(2,1))

}
