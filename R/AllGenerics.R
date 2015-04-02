## generics from S3
setGeneric("length")
setGeneric("plot")
setGeneric("summary")
setOldClass("ggplot")

#' Create a union of groups, snppicker or tags objects
#'
#' First, objects are converted to class groups.  Then any groups
#' which overlap are identified, and replaced by their union.  Groups
#' which do not overlap are unchanged.  The combined set of groups is
#' returned.
#' @rdname union
#' @author Chris Wallace
#' @export
#' @param x object of class \code{groups}, \code{snppicker} or \code{tags}
#' @param y object of same class as x
#' @return object of class groups
setGeneric("union", function(x,y) standardGeneric("union"))

#' Check whether a snp is in a snppicker, groups or tags object
#'
#' @rdname snpin
#' @author Chris Wallace
#' @export
#' @param x character vector of SNPs to be checked
#' @param y object of class snppicker, groups or tags
#' @return logical matrix of nrow equal to length(x) and ncol equal to number of groups in y
setGeneric("snpin", function(x,y) standardGeneric("snpin"))
#setGeneric("snpdrop", function(x,y) standardGeneric("drop"))

##' Convert from old definitions of groups, tags classes to new
##' 
##' DON'T USE THIS FUNCTION UNLESS YOU HAVE OBJECTS STORED FROM A
##' PREVIOUS PRE-DEVELOPMENT VERSION!
##' 
#' @rdname conversion
##' @param object GUESSFM object from pre-development version
##' @return new S4 structure
setGeneric("convert",function(object) standardGeneric("convert"))

##' Accessor functions
##'
##' \code{snps} shows list of snps in an object of class groups and
##' returns a list of character vectors
##'
##' @title Accessors for groups objects
##' @param object object from which items should be extracted
##' @return a list of character vectors
##' @author Chris Wallace
##' @rdname accessors
##' @export
setGeneric("snps",function(object) standardGeneric("snps"))

##' @details \code{tags} shows tags from an object of class groups and
##' returns a character vector of tag SNPs
##' @rdname accessors
##' @export
setGeneric("tags",function(object) standardGeneric("tags"))

##' \code{tagsof} shows tags for a named character vector of SNPs
##' @rdname groups-subset
##' @export
##' @param object tags or groups object
##' @return data.frame of tags and their tagged SNPs
setGeneric("tagsof", function(object,i) standardGeneric("tagsof"))

##' \code{taggedby} shows SNPs tagged by a named character vector of tag SNPs
##' @rdname groups-subset
##' @export
setGeneric("taggedby", function(object,i) standardGeneric("taggedby"))

##' qc a GUESSFM run
##'
##' With all genetic data, we use some QC measures to determine "bad"
##' SNPs.  The qc() functions in GUESSFM attempt to flag features that
##' experience suggests is related to spurious differential SNP calls
##' between cases and controls.
##'
##' The function \code{\link{pp.nsnp}} generates a posterior
##' distribution for the number of SNPs in a model.  We expect this
##' posterior distribution to have some right skew (as does the
##' binomial or beta binomial prior) and be unimodal.  Experience
##' suggests that a posterior that does not have these properties may
##' have favoured models with "bad" SNPs.  Running \code{qc} on the
##' object returned by \code{pp.nsnp} will flag these issues.
##'
##' You can also call \code{qc} directly on a \code{snpmod} object.
##' This may take a little longer, and attempts to estimate the
##' maximum r squared between SNPs in any model.  GUESS has a prior
##' which should enforce that highly correlated SNPs are not both
##' placed in a model.  Sometimes it may be that two correlated SNPs
##' are indeed required to model a trait, but experience with imputed
##' data suggests that when a majority of models above a given size
##' contain highly correlated SNPs, there is a problem with
##' differential genotype calling which requires further
##' investigation.
##' 
##' @title Quality control
##' @return data.frame of traits in pp.nsnp together with qc measures
##' or data.frame of models and associated size and max r squared.
##' @export
##' @author Chris Wallace
##' @rdname qc
##' @param object snpmod object or object returned by \code{pp.nsnp}
##' @param data SnpMatrix data for LD calculation if object is a snpmod
setGeneric("qc", function(object,data) standardGeneric("qc"))
