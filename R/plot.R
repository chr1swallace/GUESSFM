summary.plots <- function(results) {
  if(!is.list(results))
    results <- list(trait=results)
  plots <- list(diffusion=plot.diffuse(results),
                pp.nsnp=pp.nsnp(results))                
}


plot.bf.pp <- function(results, nm="unknown",maxrank=1000) {
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
##' This is measured by how rapidly the cumulative posterior probability increases with number of models.  For one dataset, this may be informative in an absolute sense, but comparing different datasets, by supplying a list of results, can be useful too.
##' @title Diffusion plot
##' @param results an object of class snpmod, or a list of such objects
##' @param maxrank truncate the x axis at this value
##' @return no return value
##' @author Chris Wallace
plot.diffuse <- function(results, maxrank=1000) {
  if(!is.list(results))
    results <- list(trait=results)
  df <- lapply(results, function(x) x@models)
  for(i in 1:length(results)) {
    df[[i]]$phenotype <- names(results)[i]
    df[[i]]$cPP <- cumsum(df[[i]]$PP)
  }
  df <- do.call("rbind",df)
  df <- subset(df, rank<=maxrank)
  cat(levels(df$phenotype <- factor(df$phenotype)))
  ggplot(df, aes(x=rank,y=cPP, col=phenotype)) + geom_path()
}

ggld <- function(data, df.snps) {
  all.snps <- unique(as.character(df.snps$snp))
  snps.num <- structure(1:length(all.snps),names=all.snps)
  LD <- melt(as(ld(data[,all.snps], stats="R.squared", depth=length(all.snps)-1, symmetric=TRUE),"matrix"))
  LD$X1 <- snps.num[as.character(LD$X1)] - 1
  LD$X2 <- snps.num[as.character(LD$X2)] - 1
  n <- length(all.snps)
  offset <- n/sqrt(2)
  LD$A <- with(LD, -(n-X1-X2)/sqrt(2))
  LD$B <- with(LD, -(n-X1-X2)/sqrt(2) - X1*sqrt(2))
  LD <- LD[ LD$X1>LD$X2, ]

  ## align A to 1:n
  LD$A <- xscale(LD$A, torange=c(1.5, n-0.5))
  ## scale by pos
  LD$A <- xscale(LD$A, xrange=c(min(df.snps$snpnum),max(df.snps$snpnum)), torange=c(min(df.snps$position.plot),max(df.snps$position.plot)))
  
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
                  geom_segment(data     = data.frame(x=sort(unique(df.snps$x.scale)),
                                 y=max(LD$B) + tlength, yend=max(LD$B) + 2*tlength), 
                               aes(x    = x, 
                                   y    = y,
                                   xend = x, 
                                   yend = yend  ))      #top-ticks
}  

signal.plot <- function(summ,w=0.2,highlight=NULL) {
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

ggbed <- function(bed,summ) {
  mn <- min(summ$position.plot, na.rm=TRUE)
  mx <- max(summ$position.plot, na.rm=TRUE)
  ## silly, can use + xlim(mn,mx) if needed
  ##  bed <- subset(bed, end>mn | start < mx)

  ggplot(bed, aes(xmin=start,xmax=end,ymin=0,ymax=1,alpha=height)) + ggplot2::geom_rect() + facet_grid(name ~ ., margins=TRUE) +
    geom_vline(data=summ, mapping=aes(xintercept=position.plot,col=tag)) +
    theme(legend.position="none",
          strip.text.y = element_text(size=9,angle=0),
          axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
              panel.grid=element_blank())

}


getgenes <- function (summ) {
  ## these lines based on plotGenes function in LDheatmap
    minRange <- min(summ$position)
    maxRange <- max(summ$position)
    require(rtracklayer)
    require(GenomicRanges)
    session <- browserSession()
    genome(session) <- "hg19"
    query1 <- ucscTableQuery(session, "knownGene", GRangesForUCSCGenome("hg19", 
        "chr10", IRanges(minRange, maxRange)))
    t <- getTable(query1)
    if (!dim(t)[1]) {
        print("The genetic region of the data does not correspond to any genes in the UCSC genome browser")
        return()
    }
   t[, "name"] <- as.character(t[, "name"])
     t[, "gene_name"] <- ""
    tbl <- "kgXref"
    query2 <- ucscTableQuery(session, "knownGene", GRangesForUCSCGenome("hg19", 
        "chr10", IRanges(minRange, maxRange)), table = tbl, 
        names = t[, "name"])
    t1 <- getTable(query2)
    t1[, "kgID"] <- as.character(t1[, "kgID"])
    t1[, "geneSymbol"] <- as.character(t1[, "geneSymbol"])
    for (i in 1:dim(t)[1]) {
        gene_name <- t1[t1[, "kgID"] == t[i, "name"], "geneSymbol"]
        if (length(gene_name) != 0) 
            t[i, "gene_name"] <- gene_name
    }
  ##   if ("kgColor" %in% tableNames(query1)) {
##         query3 <- ucscTableQuery(session, "knownGene", GRangesForUCSCGenome(genome, 
##             chromosome, IRanges(minRange, maxRange)), table = "kgColor", 
##             names = t[, "name"])
##         color_tbl <- getTable(query3)
##     }
    
    return(t)
  }

data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start, end = df$end))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                     "DataFrame")
        elementMetadata(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

scalepos <- function(summ,pos="position.hg19") {
  summ$position.plot <- summ[,pos]
  pr <- c(min(summ$position.plot),max(summ$position.plot))
  xr <- c(min(summ$snpnum),max(summ$snpnum))
  summ$x.scale <- xscale(summ$snpnum,torange=pr,xrange=xr)
  summ$xmin.scale <- xscale(summ$x.min,torange=pr,xrange=xr)
  summ$xmax.scale <- xscale(summ$x.max,torange=pr,xrange=xr)
  return(summ)
}

pp.nsnp <- function(results,plot=FALSE) {
  if(!is.list(results))
    results <- list(trait=results)
  df <- vector("list",length(results))
  for(i in seq_along(df)) {
    d <- results[[i]]
    pp.nsnps <- tapply(d@models$PP,d@models$size,sum)
    df[[i]] <- data.frame(n=as.numeric(names(pp.nsnps)),pp=pp.nsnps,trait=names(results)[i])
  }
  df <- do.call("rbind",df)
  p <- ggplot(df,aes(x=n,y=pp,col=trait)) + geom_point() + geom_path() + 
    xlab("Number of SNPs in model") + ylab("Posterior Probability") + scale_x_continuous(breaks=seq(0,max(df$n),by=2))
  print(p)
  print(pp.nsnps)
  invisible(list(pp=pp.nsnps,plot=p))
}

addlines <- function(summ) {
  geom_vline(data=summ, mapping=aes(xintercept=position.plot,col=tag),alpha=0.5,size=0.2) 
}
