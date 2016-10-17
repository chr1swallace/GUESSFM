##' Create an object of class ESS given a base file name
##'
##' A sensible base file name is the file ending in
##' _sweeps_features.txt.  If you used run.bvs() to run GUESS, then
##' this function should find all the correct files.  If you didn't,
##' then this function probably won't work, and you should use
##' as.ESS.object from the R2GUESS package directly and use ess2snpmod
##' to allow R2GUESS functions to work.
##' @title ess.read
##' @param f base file name, a character string
##' @param ... optional arguments passed to as.ESS.object
##' @return object of class ESS
##' @export
##' @author Chris Wallace
read.ess <- function(f,...) {
    f <- basefile(f)
    message("Reading from base file ",f)
    DIR <- paste0(dirname(f),"/")
    n <- gsub(".*/out_|_sweeps_features.txt","",f)
    suff <- function(str) {
        paste(str,n,sep="_")
    }
    decode <- read.decode(paste0(DIR,"/",suff("decode")))
    as.ESS.object(dataY=suff("Y"),dataX=suff("X"),file.par="par.xml",command=FALSE,
                  path.input=DIR,path.output=DIR,path.par=DIR,
                  root.file.output=sprintf("out_%s_sweeps",n),
                  label.X=decode,...)
}
##' Convert an ESS object to a snpmod object
##'
##' R2GUESS and GUESSFM are aimed at slightly different things.  That
##' means different aspects of the output are important.  R2GUESS's
##' major object class is ESS.  For GUESSFM it is snpmod.  This
##' function takes an object of class ESS and creates the
##' corresponding snpmod object.  The two contain slightly different
##' information, so this is not a one to one transformation.  You
##' cannot do the reverse, and create an ESS object from a snpmod one,
##' instead use the read.ess() function.
##' @title ess2snpmod
##' @param ess an object of class ESS
##' @return an object of class snpmod
##' @export
##' @author Chris Wallace
ess2snpmod <- function(ess## ,decode=NULL
                       ) {
    ## ESS is just a list.  Extract the paths we need
#    if(is.null(decode)) {
        read.snpmod(ess$path.input)
#    } else {
#        reader(f,decode)
#    }   
}

##' Read a decode file, internal function
##'
##' NB decode file is written 0-based, but ESS input is read 1-based, so do that switch here.
##' @title read.decode
##' @param dfile filename
##' @return named character vector, names=snp numbers (1-based), elements=snp names
##' @author Chris Wallace
read.decode <- function(dfile) {

    ## decode snp numbers to names
    if(!file.exists(dfile)) 
        stop("decode file ",dfile," not found")
    decode <- read.table(dfile,header=FALSE,as.is=TRUE,col.names=c("varname","varnum"))
    structure(decode$varname, names=decode$varnum+1) ## num -> name

}

basefile <- function(f,patt="^out_.*.txt") {
    if(file.exists(f) && file.info(f)$isdir) {
        ofiles <- list.files(f, pattern=patt,full.names=TRUE)
        if(!length(ofiles))
            return(NULL)
        ofiles <- ofiles[ order( file.info(ofiles)$mtime, decreasing=TRUE ) ]
        f <- gsub("_features.*|_output.*|_sweeps.*", "", ofiles[[1]])
        message(f," is a directory - reading most recent filestub: ",f)
    } else {
        message("reading supplied filestub: ",f)
    }
    return(f)
}

reader <- function(f,decode,offset) {
  files <- paste0(f,c("_output_marg_prob_incl.txt","_output_best_visited_models.txt"))
  if(any(!file.exists(files))) {
    message("not both marg_prob_incl.txt and _best_visited_models.txt exist")
    return(NULL)
  }

  pp <- try(read.table(files[[1]], as.is=TRUE, header=TRUE))
  if(inherits(pp, "try-error"))
      return(NULL)
  rownames(pp) <- pp$var <- make.names(decode[ pp$Predictor ])

  models <- try(read.table(files[[2]], as.is=TRUE, header=FALSE, skip=1, comment.char="", fill=TRUE)) # sometimes file is empty
  if(inherits(models, "try-error"))
      return(NULL)

  cnames <- c("rank","visits","size","logPP","PP","jeffreys")
  maxsize <- ncol(models) - length(cnames)
  colnames(models) <- c(cnames,paste0("s",1:maxsize))  
  models <- subset(models, !is.na(size))
  if(max(models$size)>maxsize) { # run again
      maxsize <- max(models$size)
      models <- read.table(files[[2]], as.is=TRUE, header=FALSE, skip=1, comment.char="", fill=TRUE,
                           col.names=c("rank","visits","size","logPP","PP","jeffreys",paste0("s",1:maxsize)))
      models <- subset(models, !is.na(size))
  }
  ## force numerics
  models$logPP <- as.numeric(models$logPP)
  models$PP <- as.numeric(models$PP)
  models$jeffreys <- as.numeric(models$jeffreys)
  if(any(is.nan(models$PP)))
    models$PP <- exp(models$logPP - logsum(models$logPP))

  maxsize <- max(models$size)
  models <- models[, 1:(6+max(models$size))]
  message(nrow(models), " models, with pp > ",signif(models[nrow(models), "PP"],3),
          ", and, cumulatively, to ",signif(sum(models$PP),3))

  model.snps <- apply(models[,-c(1:6)], 1, paste, collapse=",")
  model.snps <- strsplit(model.snps, ",")
  model.snps <- lapply(lapply(model.snps, setdiff, "NA"), as.numeric)
  model.snpnames <- lapply(model.snps, function(x) make.names(decode[ x ]))
  models <- models[,1:6]
  models$str <- unlist(lapply(model.snpnames,makestr))
  if(offset!=0) {
    pp$Predictor <- pp$Predictor - offset
    model.snps <- lapply(model.snps, function(x) x-offset)
    ##   for(j in 1:max(models$size)) {
    ##     v <- paste0("s",j)
    ##     models[,v] <- models[,v] - offset
    ##   }
  }
  return(new("snpmod",
             snps=pp, models=models, model.snps=model.snpnames))

}

##' Read in the output files and create an object of class snpmod
##'
##' @title read.snpmod
##' @param f filename
##' @param offset number of covariates
##' @return object of class snpmod
##' @author Chris Wallace
##' @export
read.snpmod <- function(f,offset=0) {
    f.decode <- basefile(f,patt="decode_[0-9]")
    decode <- read.decode(f.decode)

    f <- basefile(f)
    reader(f=f,decode=decode,offset=offset)
}

##' guess.read, for backwards compatability
##'
##' Just a wrapper for read.snpmod.
##' @title guess.read
##' @param ... arguments to be passed to read.snpmod
##' @return see read.snpmod
##' @author Chris Wallace
##' @export
guess.read <- function(...) { read.snpmod(...) }
